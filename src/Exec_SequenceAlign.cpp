#include <cstring>
#include <cctype> // isdigit, isalpha
#include "Exec_SequenceAlign.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "Trajout_Single.h"
#include "StringRoutines.h"

// EXPERIMENTAL ALPHA CODE

/** CONSTRUCTOR */
Exec_SequenceAlign::Exec_SequenceAlign() :
  Exec(GENERAL),
  queryOffset_(0),
  sbjctOffset_(0)
{
  SetHidden(true);
}

/// \return true if character is a valid sequence align char
static inline bool is_valid_seq(char c) {
  return (isalpha( c ) || c == '-');
}

/** Advance past any digit to the sequence alignment. */
int Exec_SequenceAlign::advance_past_rnum(std::string::const_iterator& it,
                                          std::string const& line)
{
  int rnum = -1;
  std::string rnumStr;
  bool in_number = false;
  while (it != line.end()) {
    if (isdigit( *it )) { // TODO are negative res numbers allowed?
      if (!in_number && rnum == -1)
        in_number = true;
      if (in_number)
        rnumStr += *it;
    } else if (isspace( *it )) {
      if (in_number) {
        in_number = false;
        rnum = convertToInteger( rnumStr );
      }
    } else if (is_valid_seq( *it )) {
      // At first alpha. Break.
      break;
    }
    ++it;
  }
  mprintf("DEBUG: alpha starts at %li, rnum= %i\n", it-line.begin(), rnum);
  return rnum;
}

/** Read sequence alignment. */
int Exec_SequenceAlign::read_blast(std::string& queryMaskStr,
                                   std::string& sbjctMaskStr,
                                   ResNumCharMap& queryToSbjct,
                                   std::string const& blastfile)
const
{
  mprintf("\tReading BLAST alignment from '%s'\n", blastfile.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( blastfile )) return 1;
  // Seek down to first Query (but not Query:) line.
  const char* ptr = infile.Line();
  bool atFirstQuery = false;
  while (ptr != 0) {
    if (*ptr == 'Q') {
      if ( strncmp(ptr, "Query ", 6) == 0 ) {
        atFirstQuery = true;
        break;
      }
    }
    ptr = infile.Line();
  }
  if (!atFirstQuery) {
    mprinterr("Error: 'Query' not found.\n");
    return 1;
  }

  // Read alignment.
  //typedef std::vector<char> Carray;
  typedef std::vector<int> Iarray;
  Iarray queryResidues;
  Iarray sbjctResidues;
  std::string Query; // Query residues
  std::string Align; // Alignment line
  std::string Sbjct; // Sbjct residues
  int qres = -1; // Current query residue #
  int sres = -1; // Current subject residue #

  while (ptr != 0) {
    Query.clear();
    Sbjct.clear();
    Align.clear();
    if (ptr[0] == 'Q') {
      // We are at a Query line
      Query.assign( ptr );
      // Get the alignment line
      const char* aline = infile.Line();
      if (aline == 0) {
        mprinterr("Error: Missing alignment line after Query, line %i\n", infile.LineNumber());
        return 1;
      }
      Align.assign( aline );
    } else {
      mprintf("Warning: No Query/Align lines, line %i.\n", infile.LineNumber());
    }

    if (ptr[0] == 'S') {
      // Already at Sbjct line. This can happen when there are no more Query lines.
      Sbjct.assign( ptr);
    } else {
      // Get subject line
      const char* sline = infile.Line();
      if (sline == 0 || sline[0] != 'S') {
        mprintf("Warning: Missing Subject line after alignment, line %i\n", infile.LineNumber());
      } else {
        Sbjct.assign( sline );
      }
    }

    // DEBUG
    mprintf("DEBUG: Query: %s\n", Query.c_str());
    mprintf("DEBUG: Align: %s\n", Align.c_str());
    mprintf("DEBUG: Sbjct: %s\n", Sbjct.c_str());

    if (!Query.empty() && !Align.empty() && !Sbjct.empty()) {
      if (Query.size() < 6 || Sbjct.size() < 6) {
        mprinterr("Error: Query and/or Sbjct line sizes are short (line %i)\n", infile.LineNumber());
        return 1;
      }
      // Process Query/Sbjct alignment
      std::string::const_iterator qit = Query.begin() + 5;
      std::string::const_iterator sit = Sbjct.begin() + 5;
      // Get and advance past any residue numbers
      int rn = advance_past_rnum(qit, Query);
      if (qres == -1) {
        qres = rn;
        if (rn == -1) {
          mprintf("Warning: No initial residue number for query. Setting to 1.\n");
          qres = 1;
        }
        mprintf("Query starting from residue %i\n", qres);
      } else if (rn != -1) {
        if (qres != rn)
          mprintf("Warning: Residue # mismatch for Query; got %i, expected %i\n", rn, qres);
      }
      rn = advance_past_rnum(sit, Sbjct);
      if (sres == -1) {
        sres = rn;
        if (rn == -1) {
          mprintf("Warning: No initial residue number for subject. Setting to 1.\n");
          sres = 1;
        }
        mprintf("Subject starting from residue %i\n", sres);
      } else if (rn != -1) {
        if (sres != rn)
          mprintf("Warning: Residue # mismatch for Sbjct; got %i, expected %i\n", rn, sres);
      }
      // Ensure query and subject are starting from the same column
      long int coloffset = qit - Query.begin();
      long int sbjoffset = sit - Sbjct.begin();
      if (sbjoffset != coloffset) {
        mprinterr("Error: Query column (%li) != Sbjct solumn (%li), line %i\n",
                   coloffset+1, sbjoffset+1, infile.LineNumber());
        return 1;
      }
      std::string::const_iterator ait = Align.begin() + coloffset;
      // Loop over query/subject
      while (qit != Query.end() || sit != Sbjct.end()) {
        if (is_valid_seq( *qit )) {
          if (!is_valid_seq( *sit )) {
            mprinterr("Error: Invalid Sbjct char %c corresponding to Query %c, line %i\n",
                      *sit, *qit, infile.LineNumber());
            return 1;
          }
          if (isalpha(*qit) && isalpha(*sit)) {
            // 1 to 1 correspondence
            int qResnum = qres+queryOffset_;
            int sResnum = sres+sbjctOffset_;
            mprintf("DEBUG: qres %i %c to sres %i %c (%c)\n", qResnum, *qit, sResnum, *sit, *ait);
            queryResidues.push_back( qResnum );
            sbjctResidues.push_back( sResnum );
            queryToSbjct.insert( ResNumCharPair(qResnum, ResNumChar(*qit, sResnum, *sit)) );
            qres++;
            sres++;
          } else if (isalpha(*qit) && *sit == '-') {
            // Query, no sbjct
            mprintf("DEBUG: qres %i %c\n", qres, *qit);
            qres++;
          } else if (isalpha(*sit) && *qit == '-') {
            // Sbjct, no query
            mprintf("DEBUG: sres %i %c\n", sres, *sit);
            sres++;
          }
        }
        if (qit != Query.end()) ++qit;
        if (sit != Sbjct.end()) ++sit;
        if (ait != Align.end()) ++ait;
        // TODO check end numbers if present?
        if (!is_valid_seq(*qit)) {
          if (is_valid_seq(*sit)) {
            mprinterr("Error: Query %c has ended before Sbjct %c, line %i\n", *qit, *sit, infile.LineNumber());
            return 1;
          }
        }
        if (!is_valid_seq(*sit)) {
          if (is_valid_seq(*qit)) {
            mprinterr("Error: Sbjct %c has ended before Query %c, line %i\n", *sit, *qit, infile.LineNumber());
            return 1;
          }
        }
      } // END loop over lines
    } // END process alignment

    // Scan to next Query or Sbjct
    ptr = infile.Line();
    while (ptr != 0) {
      if (*ptr == 'Q') {
        if ( strncmp(ptr, "Query", 5) == 0 ) break;
      } else if (*ptr == 'S') {
        if ( strncmp(ptr, "Sbjct", 5) == 0 ) break;
      }
      ptr = infile.Line();
    }
  } // END loop over Alignment file
  infile.CloseFile();

  std::string queryResStr = ArrayToRangeExpression( queryResidues, 0 );
  std::string sbjctResStr = ArrayToRangeExpression( sbjctResidues, 0 );
  mprintf("\tQuery Residues (%zu): %s\n", queryResidues.size(), queryResStr.c_str());
  mprintf("\tSbjct Residues (%zu): %s\n", sbjctResidues.size(), sbjctResStr.c_str());
  mprintf("DEBUG: Query to Sbjct map:\n");
  for (ResNumCharMap::const_iterator it = queryToSbjct.begin();
                                     it != queryToSbjct.end(); ++it)
    mprintf("\t%8i %c %8i %c\n", it->first, it->second.ResChar(), it->second.MappedNum(), it->second.MappedChar());

  if (queryResidues.size() != sbjctResidues.size()) {
    mprinterr("Error: There is not a 1 to 1 correspondence between Query and Sbjct residues.\n");
    return 1;
  }
  queryMaskStr = ":" + queryResStr;
  sbjctMaskStr = ":" + sbjctResStr;
  return 0;
}

/** Print help text. */
void Exec_SequenceAlign::Help() const {
  mprintf("\t%s blastfile <file> out <file>\n"
          "\t[{pdb | mol2}] [<trajout args>] [smaskoffset <#>] [qmaskoffset <#>]\n"
          "  blastfile: File containing sequence alignment.\n"
          "  out: File to write resulting structure to.\n"
          "  [{pdb | mol2}] Format of structure (default pdb).\n"
          "  [smaskoffset]|[qmaskoffset]: Resdiue offset for output.\n"
          "Given a reference structure and BLAST-like sequence alignment of format:\n"
          "    Query  1    MIV...\n"
          "                MI+...\n"
          "    Sbjct  1    MIM...\n"
          "where the Query sequence corresponds to the reference, output a structure\n"
          "with the Sbjct sequence. All hydrogens will be removed. Non-matching residues\n"
          "will include only backbone and CB (if applicable) atoms.\n",
          DataSetList::RefArgs);
}

/** Execute. */
Exec::RetType Exec_SequenceAlign::Execute(CpptrajState& State, ArgList& argIn) {
  mprintf("Warning: THIS COMMAND IS NOT FULLY IMPLEMENTED.\n");
  std::string blastfile = argIn.GetStringKey("blastfile");
  if (blastfile.empty()) {
    mprinterr("Error: 'blastfile' must be specified.\n");
    return CpptrajState::ERR;
  }

  bool map_residues = true;
  ReferenceFrame qref = State.DSL().GetReferenceFrame(argIn);
  TrajectoryFile::TrajFormatType fmt = TrajectoryFile::UNKNOWN_TRAJ;
  std::string outfilename;
  
  if (qref.error()) {
    mprinterr("Error: Must specify reference structure for query.\n");
    return CpptrajState::ERR;
  } else if (qref.empty()) {
    map_residues = false;
  } else {
    // Get residue mapping keywords.
    outfilename = argIn.GetStringKey("out");
    if (outfilename.empty()) {
      mprinterr("Error: Must specify output file if reference is specified.\n");
      return CpptrajState::ERR;
    }
    // Default to PDB. TODO only allow PDB/Mol2?
    fmt = TrajectoryFile::WriteFormatFromArg(argIn, TrajectoryFile::PDBFILE);
  }
  sbjctOffset_ = argIn.getKeyInt("smaskoffset", 0);
  queryOffset_ = argIn.getKeyInt("qmaskoffset", 0);

  // Load blast file
  ResNumCharMap queryToSbjct;
  std::string queryMaskStr, sbjctMaskStr;
  if (read_blast(queryMaskStr, sbjctMaskStr, queryToSbjct, blastfile))
    return CpptrajState::ERR;

  mprintf("\tReading BLAST alignment from '%s'\n", blastfile.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( blastfile )) return CpptrajState::ERR;
  // Seek down to first Query line.
  const char* ptr = infile.Line();
  bool atFirstQuery = false;
  while (ptr != 0) {
    if (*ptr == 'Q') {
      if ( strncmp(ptr, "Query", 5) == 0 ) {
        atFirstQuery = true;
        break;
      }
    }
    ptr = infile.Line();
  }
  if (!atFirstQuery) {
    mprinterr("Error: 'Query' not found.\n");
    return CpptrajState::ERR;
  }

  // Read alignment. Replacing query with subject.
  typedef std::vector<char> Carray;
  typedef std::vector<int> Iarray;
  Carray Query; // Query residues
  Carray Sbjct; // Sbjct residues
  Iarray Smap;  // Smap[Sbjct index] = Query index
  while (ptr != 0) {
    const char* qline = ptr;           // query line
    const char* aline = infile.Line(); // alignment line
    const char* sline = infile.Line(); // subject line
    if (aline == 0 || sline == 0) {
      mprinterr("Error: Missing alignment line or subject line after Query:\n");
      mprinterr("Error:  %s", qline);
      return CpptrajState::ERR;
    }
    for (int idx = 12; qline[idx] != ' '; idx++) {
      if (qline[idx] == '-') {
        // Sbjct does not have corresponding res in Query
        Smap.push_back(-1);
        Sbjct.push_back( sline[idx] );
      } else if (sline[idx] == '-') {
        // Query does not have a corresponding res in Sbjct
        Query.push_back( qline[idx] );
      } else {
        // Direct Query to Sbjct map
        Smap.push_back( Query.size() );
        Sbjct.push_back( sline[idx] );
        Query.push_back( qline[idx] );
      }
    }
    // Scan to next Query 
    ptr = infile.Line();
    while (ptr != 0) {
      if (*ptr == 'Q') {
        if ( strncmp(ptr, "Query", 5) == 0 ) break;
      }
      ptr = infile.Line();
    }
  }
  // DEBUG
  int smaskoffset=1; // FIXME
  int qmaskoffset=1; // FIXME
  std::string SmaskExp, QmaskExp;
  if (State.Debug() > 0) mprintf("  Map of Sbjct to Query:\n");
  for (int sres = 0; sres != (int)Sbjct.size(); sres++) {
    if (State.Debug() > 0)
      mprintf("%-i %3s %i", sres+smaskoffset, Residue::ConvertResName(Sbjct[sres]),
              Smap[sres]+qmaskoffset);
    const char* qres = "";
    if (Smap[sres] != -1) {
      qres = Residue::ConvertResName(Query[Smap[sres]]);
      if (SmaskExp.empty())
        SmaskExp.assign( integerToString(sres+smaskoffset) );
      else
        SmaskExp.append( "," + integerToString(sres+smaskoffset) );
      if (QmaskExp.empty())
        QmaskExp.assign( integerToString(Smap[sres]+qmaskoffset) );
      else
        QmaskExp.append( "," + integerToString(Smap[sres]+qmaskoffset) );

    }
    if (State.Debug() > 0) mprintf(" %3s\n", qres);
  }
  mprintf("Smask: %s\n", SmaskExp.c_str());
  mprintf("Qmask: %s\n", QmaskExp.c_str());
  // Check that query residues match reference.
  for (unsigned int sres = 0; sres != Sbjct.size(); sres++) {
    int qres = Smap[sres];
    if (qres != -1) {
      if (Query[qres] != qref.Parm().Res(qres).SingleCharName()) {
        mprintf("Warning: Potential residue mismatch: Query %s reference %s\n",
                Residue::ConvertResName(Query[qres]), qref.Parm().Res(qres).c_str());
      }
    }
  }
  // Build subject using coordinate from reference.
  //AtomMask sMask; // Contain atoms that should be in sTop
  Topology sTop;
  Frame sFrame;
  Iarray placeHolder; // Atom indices of placeholder residues.
  for (unsigned int sres = 0; sres != Sbjct.size(); sres++) {
    int qres = Smap[sres];
    NameType SresName( Residue::ConvertResName(Sbjct[sres]) );
    if (qres != -1) {
      Residue const& QR = qref.Parm().Res(qres);
      Residue SR(SresName, sres+1, ' ', QR.ChainID());
      if (Query[qres] == Sbjct[sres]) { // Exact match. All non-H atoms.
        for (int qat = QR.FirstAtom(); qat != QR.LastAtom(); qat++)
        {
          if (qref.Parm()[qat].Element() != Atom::HYDROGEN) {
            sTop.AddTopAtom( qref.Parm()[qat], SR );
            sFrame.AddXYZ( qref.Coord().XYZ(qat) );
            //sMask.AddAtom(qat);
          }
        }
      } else { // Partial match. Copy only backbone and CB.
        for (int qat = QR.FirstAtom(); qat != QR.LastAtom(); qat++)
        {
          if ( qref.Parm()[qat].Name().Match("N" ) ||
               qref.Parm()[qat].Name().Match("CA") ||
               qref.Parm()[qat].Name().Match("CB") ||
               qref.Parm()[qat].Name().Match("C" ) ||
               qref.Parm()[qat].Name().Match("O" ) )
          {
            sTop.AddTopAtom( qref.Parm()[qat], SR );
            sFrame.AddXYZ( qref.Coord().XYZ(qat) );
          }
        }
      }
    } else {
      // Residue in query does not exist for subject. Just put placeholder CA for now.
      Vec3 Zero(0.0);
      placeHolder.push_back( sTop.Natom() );
      sTop.AddTopAtom( Atom("CA", "C "), Residue(SresName, sres+1, ' ', "") );
      sFrame.AddXYZ( Zero.Dptr() );
    }
  }
  //sTop.PrintAtomInfo("*");
  mprintf("\tPlaceholder residue indices:");
  for (Iarray::const_iterator p = placeHolder.begin(); p != placeHolder.end(); ++p)
    mprintf(" %i", *p + 1);
  mprintf("\n");
  // Try to give placeholders more reasonable coordinates.
  if (!placeHolder.empty()) {
    Iarray current_indices;
    unsigned int pidx = 0;
    while (pidx < placeHolder.size()) {
      if (current_indices.empty()) {
        current_indices.push_back( placeHolder[pidx++] );
        // Search for the end of this segment
        for (; pidx != placeHolder.size(); pidx++) {
          if (placeHolder[pidx] - current_indices.back() > 1) break;
          current_indices.push_back( placeHolder[pidx] );
        }
        // DEBUG
        mprintf("\tSegment:");
        for (Iarray::const_iterator it = current_indices.begin();
                                    it != current_indices.end(); ++it)
          mprintf(" %i", *it + 1);
        // Get coordinates of residues bordering segment.
        int prev_res = sTop[current_indices.front()].ResNum() - 1;
        int next_res = sTop[current_indices.back() ].ResNum() + 1;
        mprintf(" (prev_res=%i, next_res=%i)\n", prev_res+1, next_res+1);
        Vec3 prev_crd(sFrame.XYZ(current_indices.front() - 1));
        Vec3 next_crd(sFrame.XYZ(current_indices.back()  + 1));
        prev_crd.Print("prev_crd");
        next_crd.Print("next_crd");
        Vec3 crd_step = (next_crd - prev_crd) / (double)(current_indices.size()+1);
        crd_step.Print("crd_step");
        double* xyz = sFrame.xAddress() + (current_indices.front() * 3);
        for (unsigned int i = 0; i != current_indices.size(); i++, xyz += 3) {
          prev_crd += crd_step;
          xyz[0] = prev_crd[0];
          xyz[1] = prev_crd[1];
          xyz[2] = prev_crd[2];
        }
        current_indices.clear();
      }
    }
  }
  //Topology* sTop = qref.Parm().partialModifyStateByMask( sMask );
  //if (sTop == 0) return CpptrajState::ERR;
  //Frame sFrame(qref.Coord(), sMask);
  // Write output traj
  Trajout_Single trajout;
  if (trajout.PrepareTrajWrite(outfilename, argIn, State.DSL(), &sTop, CoordinateInfo(), 1, fmt))
    return CpptrajState::ERR;
  if (trajout.WriteSingle(0, sFrame)) return CpptrajState::ERR;
  trajout.EndTraj();
  return CpptrajState::OK;
}
