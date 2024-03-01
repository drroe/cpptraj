#include "AmberParamFile.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "ParameterHolders.h"
#include "StringRoutines.h"
#include "ParameterSet.h"
#include "Constants.h"
#include <cstdio> // sscanf

const int AmberParamFile::MAXSYMLEN = 16;

/// CONSTRUCTOR
AmberParamFile::AmberParamFile() :
  debug_(0)
{}

/** Set debug level */
void AmberParamFile::SetDebug(int d) {
  debug_ = d;
}

/** Read symbols delimited by - and space. */
int AmberParamFile::read_symbols(const char* ptrIn, std::vector<std::string>& symbols, int nsymbols)
{
  int isymbol = 0;
  bool char_has_been_read = false;
  for (const char* ptr = ptrIn; *ptr != '\0'; ++ptr)
  {
    if (*ptr == '-') {
      isymbol++;
      char_has_been_read = false;
    } else if (*ptr == ' ' && isymbol + 1 == nsymbols && char_has_been_read) {
      return (ptr - ptrIn);
    } else {
      symbols[isymbol] += *ptr;
      if (*ptr != ' ') char_has_been_read = true;
    }
  }
  return -1;
}

/// Hold a set of nonbonded parameters
class AmberParamFile::NonbondSet {
  public:
    NonbondSet(std::string const& n) : name_(n) {}

    std::string name_;          ///< Name of set parameters
    ParmHolder<LJparmType> LJ_; ///< Hold LJ 6-12 parameters
};

/// Hold an off-diagonal NB modification
class AmberParamFile::OffdiagNB {
  public:
    OffdiagNB(NameType const& AT1, NameType const& AT2, double sig1, double eps1, double sig2, double eps2) :
      types_(2), LJ1_(sig1, eps1), LJ2_(sig2, eps2)
    {
      types_.AddName(AT1);
      types_.AddName(AT2);
    }

    TypeNameHolder types_;
    LJparmType LJ1_;
    LJparmType LJ2_;
};

/** Read input for atom symbols and masses. */
int AmberParamFile::read_atype(ParameterSet& prm, const char* ptr)
const
{
  // Format (A2,2X,F10.2x,f10.2)
  if (debug_ > 1) mprintf("DEBUG: Atype: %s\n", ptr);
  char kndsym[MAXSYMLEN];
  double amass = 0;
  double atpol = 0;
  int nscan = sscanf(ptr, "%s %lf %lf", kndsym, &amass, &atpol);
  ParameterHolders::RetType ret;
  if (nscan == 3) {
    ret = prm.AT().AddParm( TypeNameHolder(kndsym),
                            AtomType(amass, atpol),
                            true );
  } else if (nscan == 2) {
    // Only mass
    ret = prm.AT().AddParm( TypeNameHolder(kndsym),
                            AtomType(amass),
                            true );
  } else {
    mprinterr("Error: Expected atom type, mass, polarizability, got only %i columns.\n", nscan);
    return 1;
  }
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining atom type %s\n", kndsym);
  return 0;
}

/** Read input for bond. */
int AmberParamFile::read_bond(ParameterSet& prm, const char* ptr)
const
{
  // Bond parameters
  // IBT , JBT , RK , REQ
  // FORMAT(A2,1X,A2,2F10.2)
  if (debug_ > 1) mprintf("DEBUG: Bond: %s\n", ptr);
  std::vector<std::string> symbols(2);
  int pos = read_symbols(ptr, symbols, 2);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for bond from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), ptr+pos);
  double RK, REQ;
  int nscan = sscanf(ptr+pos, "%lf %lf", &RK, &REQ);
  if (nscan != 2) {
    mprinterr("Error: Expected RK, REQ, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(2);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  ParameterHolders::RetType ret = prm.BP().AddParm(types, BondParmType(RK, REQ), true);
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining bond type %s - %s\n", *(types[0]), *(types[1]));

  return 0;
}

/** Read input for angle. */
int AmberParamFile::read_angle(ParameterSet& prm, const char* ptr)
const
{
  // Angle parameters
  // ITT , JTT , KTT , TK , TEQ
  // FORMAT(A2,1X,A2,1X,A2,2F10.2)
  if (debug_ > 1) mprintf("DEBUG: Angle: %s\n", ptr);
  std::vector<std::string> symbols(3);
  int pos = read_symbols(ptr, symbols, 3);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for angle from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), ptr+pos);
  double TK, TEQ;
  int nscan = sscanf(ptr+pos, "%lf %lf", &TK, &TEQ);
  if (nscan != 2) {
    mprinterr("Error: Expected TK, TEQ, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(3);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  ParameterHolders::RetType ret = prm.AP().AddParm(types, AngleParmType(TK, TEQ*Constants::DEGRAD), true);
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining angle type %s - %s - %s\n", *(types[0]), *(types[1]), *(types[2]));

  return 0;
}

/** Read input for dihedral.
  * IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
  * FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
  * If IPT .eq. 'X ' .and. LPT .eq. 'X ' then any dihedrals in the
  * system involving the atoms "JPT" and and "KPT" are assigned 
  * the same parameters.  This is called the general dihedral type
  * and is of the form "X "-"JPT"-"KPT"-"X ".
  * IDIVF is the factor by which the torsional barrier is divided.
  * Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
  * details. Basically, the actual torsional potential is
  *   (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
  * If PN .lt. 0.0 then the torsional potential is assumed to have more
  * than one term, and the values of the rest of the terms are read from
  * the next cards until a positive PN is encountered.  The negative value
  * of pn is used only for identifying the existence of the next term and 
  * only the absolute value of PN is kept.
  */
int AmberParamFile::read_dihedral(ParameterSet& prm, const char* ptr,
                                  std::vector<std::string>& last_symbols,
                                  bool first_char_is_space)
const
{
  // Dihedral parameters
  if (debug_ > 1)
    mprintf("DEBUG: Dihedral: %s\n", ptr);
  std::vector<std::string> symbols(4);
  int pos = 0;
  if (first_char_is_space) {
    // Assume no symbols on this line. Use previous symbols.
    if (last_symbols.empty()) {
      mprinterr("Error: No symbols in dihedral line: %s\n", ptr);
      return 1;
    }
    symbols = last_symbols;
    // Advance past the whitespace
    while (ptr[pos] == ' ') pos++;
  } else {
    pos = read_symbols(ptr, symbols, 4);
    if (pos < 0) {
      mprinterr("Error: Could not read symbols for dihedral from %s\n", ptr);
      return 1;
    }
  }
  //mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
  int IDIVF;
  double PK, PHASE, PN;
  char sSCEE[128];
  char sSCNB[128];
  // TODO note when PN is negative and expect more terms?
  int nscan = sscanf(ptr+pos, "%i %lf %lf %lf %s %s", &IDIVF, &PK, &PHASE, &PN, sSCEE, sSCNB);
  if (nscan < 4) {
    mprinterr("Error: Expected IDIVF, PK, PHASE, PN, got only %i elements\n", nscan);
    return 1;
  }
  if (PN < 0.0) PN = -PN;
  double scee = 1.2; // AMBER DEFAULT
  double scnb = 2.0; // AMBER DEFAULT
  if (nscan == 6) {
    // Check for SCEE/SCNB (GLYCAM)
    if (sSCEE[0] == 'S' && sSCEE[1] == 'C' && sSCEE[2] == 'E' && sSCEE[3] == 'E')
     sscanf( sSCEE, "SCEE=%lf", &scee);
    if (sSCNB[0] == 'S' && sSCNB[1] == 'C' && sSCNB[2] == 'N' && sSCNB[3] == 'B')
     sscanf( sSCNB, "SCNB=%lf", &scnb);
  }
  TypeNameHolder types(4);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  types.AddName( symbols[3] );
  ParameterHolders::RetType ret =
    prm.DP().AddParm(types, DihedralParmType(PK / (double)IDIVF, PN, PHASE*Constants::DEGRAD, scee, scnb), true);
  if (ret == ParameterHolders::UPDATED) {
    mprintf("Warning: Redefining dihedral type %s - %s - %s - %s (PN=%g)\n",
            *(types[0]), *(types[1]), *(types[2]), *(types[3]), PN);
    //mprintf("DEBUG: %s\n", ptr);
  }
  last_symbols = symbols;
  return 0;
}

/** Read input for improper. */
int AmberParamFile::read_improper(ParameterSet& prm, const char* ptr)
const
{
  // Improper parameters
  // IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
  // FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
  if (debug_ > 1) mprintf("DEBUG: Improper: %s\n", ptr);
  std::vector<std::string> symbols(4);
  int pos = read_symbols(ptr, symbols, 4);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for improper from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
  double PK, PHASE, PN;
  int nscan = sscanf(ptr+pos, "%lf %lf %lf", &PK, &PHASE, &PN);
  if (nscan != 3) {
    mprinterr("Error: Expected PK, PHASE, PN, got only %i elements\n", nscan);
    return 1;
  }
  if (PN < 0.0) {
    mprintf("Warning: Improper for %s-%s-%s-%s has negative phase (%g)\n",
            symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), PN);
    PN = -PN;
  }
  TypeNameHolder types(4);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  types.AddName( symbols[3] );
  //types.SortImproperByAlpha("X"); // FIXME wildcard should be a static var
  ParameterHolders::RetType ret =
    prm.IP().AddParm(types, DihedralParmType(PK, PN, PHASE*Constants::DEGRAD), true);
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining improper type %s - %s - %s - %s\n",
            *(types[0]), *(types[1]), *(types[2]), *(types[3]));

  return 0;
}

/** Read input for LJ 10-12 hbond. */
int AmberParamFile::read_lj1012(ParameterSet& prm, const char* ptr)
const
{
  // Lennard-Jones 10-12 hydrogen bonding term.
  // According to the docs, the format should be:
  // KT1 , KT2 , A , B , ASOLN , BSOLN , HCUT , IC
  // FORMAT(2X,A2,2X,A2,2x,5F10.2,I2)
  // In practive (in e.g. parm91.dat), the actual format appears to be:
  // KT1, KT2, A, B, HCUT
  char KT1[MAXSYMLEN], KT2[MAXSYMLEN];
  double A, B, HCUT;
  //double ASOLN, BSOLN, HCUT;
  //int IC;
  if (debug_ > 1) mprintf("DEBUG: LJ 10-12: %s\n", ptr);
  //int nscan = sscanf(ptr, "%s %s %lf %lf %lf %lf %lf %i\n",
  int nscan = sscanf(ptr, "%s %s %lf %lf %lf", KT1, KT2, &A, &B, &HCUT);
  if (nscan < 5) {
    mprinterr("Error: Expected at least KT1, KT2, A, B, HCUT, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(2);
  types.AddName( KT1 );
  types.AddName( KT2 );
  ParameterHolders::RetType ret = prm.HB().AddParm(types, HB_ParmType(A, B, HCUT), true);
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining LJ 10-12 hbond type %s %s\n", *(types[0]), *(types[1]));
  return 0;
}

/** Read radius/depth input for LJ 6-12 nonbond. */
int AmberParamFile::read_nb_RE(NonbondSet& nbset, const char* ptr)
const
{
  // ***** ONLY IF KINDNB .EQ. 'RE' *****
  // LTYNB , R , EDEP
  if (debug_ > 1) mprintf("DEBUG: Nonbond: %s\n", ptr);
  //char LTYNB[MAXSYMLEN];
  //double R, EDEP;
  //double R14, E14;
  // This section is a little tricky. Apparently CHARMM-style Amber FF
  // files can have 14 LJ params here. Try to detect this.
  ArgList nbargs( ptr, " " );
  if (nbargs.Nargs() < 3) {
    mprinterr("Error: Expected at least TYPE, R, DEPTH, got %i elements.\n", nbargs.Nargs());
    return 1;
  }
  bool has_14 = false;
  if (nbargs.Nargs() >= 5) {
    if (validDouble( nbargs[3] ) && validDouble( nbargs[4] )) {
      has_14 = true;
    }
  }
  if (has_14) mprintf("DEBUG: NB HAS 1-4 LJ PARAMS.\n"); // FIXME save these
  double R = convertToDouble( nbargs[1] );
  double EDEP = convertToDouble( nbargs[2] );
  TypeNameHolder types( nbargs[0] );
  ParameterHolders::RetType ret = nbset.LJ_.AddParm( types, LJparmType(R, EDEP), true );
  if (ret == ParameterHolders::UPDATED)
    mprintf("Warning: Redefining LJ 6-12 type %s\n", *(types[0]));
  return 0;
}

/** Read LJ off-diagonal modifications */
int AmberParamFile::read_ljedit(Oarray& Offdiag, const char* ptr)
const
{
  if (debug_ > 1) mprintf("DEBUG: LJedit: %s\n", ptr);
  // Lennard-Jones sigma and epsilon of the first atom type when it
  // interacts with anything under the normal rules, then the sigma
  // and epsilon of the second atom type when it interacts with the first.
  char AT1[MAXSYMLEN], AT2[MAXSYMLEN];
  double sig1, eps1, sig2, eps2;
  int nscan = sscanf(ptr, "%s %s %lf %lf %lf %lf", AT1, AT2, &sig1, &eps1, &sig2, &eps2);
  if (nscan != 6) {
    mprinterr("Error: Expected AT1, AT2, SIG1, EPS1, SIG2, EPS2, got %i elements.\n", nscan);
    return 1;
  }
  Offdiag.push_back( OffdiagNB(AT1, AT2, sig1, eps1, sig2, eps2) );
  return 0;
}

/** Read CMAP section */
int AmberParamFile::read_cmap(ParameterSet& prm, CmapType& currentCmapFlag, std::string const& line) const {
  //if (ptr[0] == '%' && ptr[1] == 'F' && ptr[2] == 'L' &&
  //    ptr[3] == 'A' && ptr[4] == 'G')
  int currentCmapIdx;
  if (prm.CMAP().empty())
    currentCmapIdx = -1;
  else
    currentCmapIdx = prm.CMAP().size() - 1;
  ArgList argline(line);
  if (argline.Nargs() > 1)
  {
    if (argline[0] == "%FLAG") {
      if (argline[1] == "CMAP_COUNT") {
        // New CMAP term. Ignore the index for now. If a previous CMAP
        // was read make sure its OK.
        if (currentCmapIdx > -1) {
          if (!prm.CMAP()[currentCmapIdx].CmapIsValid()) {
            mprinterr("Error: Previous CMAP term is not valid.\n");
            return 1;
          }
          if (!prm.CMAP()[currentCmapIdx].CmapNresIsValid()) {
            mprintf("Warning: # expected CMAP %i residue names %i not equal to actual number %zu\n",
                    currentCmapIdx+1,
                    prm.CMAP()[currentCmapIdx].NcmapResNames(),
                    prm.CMAP()[currentCmapIdx].ResNames().size());
          }
        }
        int cmapcount = argline.getKeyInt("CMAP_COUNT", -1);
        mprintf("DEBUG: Cmap count: %i\n", cmapcount);
        prm.CMAP().push_back( CmapGridType() );
        if ( cmapcount != (int)prm.CMAP().size() )
          mprintf("Warning: CMAP term is not in numerical order. CMAP_COUNT %i, expected %zu\n",
                  cmapcount, prm.CMAP().size());
        return 0;
      } else if (argline[1] == "CMAP_TITLE") {
        currentCmapFlag = CMAP_TITLE;
        return 0;
      } else if (argline[1] == "CMAP_RESLIST") {
        currentCmapFlag = CMAP_RESLIST;
        int nres = argline.getKeyInt("CMAP_RESLIST", -1);
        if (nres < 1) {
          mprinterr("Error: Bad CMAP # residues: %s\n", line.c_str());
          return 1;
        }
        prm.CMAP()[currentCmapIdx].SetNumCmapRes( nres );
        return 0;
      } else if (argline[1] == "CMAP_RESOLUTION") {
        int cmapres = argline.getKeyInt("CMAP_RESOLUTION", -1);
        mprintf("DEBUG: Cmap res: %i\n", cmapres);
        if (cmapres < 1) {
          mprinterr("Error: Bad CMAP resolution: %s\n", line.c_str());
          return 1;
        }
        prm.CMAP()[currentCmapIdx].SetResolution( cmapres );
        return 0;
      } else if (argline[1] == "CMAP_PARAMETER") {
        currentCmapFlag = CMAP_PARAMETER;
        return 0;
      }
    }
  }

  if (currentCmapFlag == CMAP_PARAMETER) {
    double terms[8];
    int nterms = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf",
                        terms, terms+1, terms+2, terms+3,
                        terms+4, terms+5, terms+6, terms+7);
    for (int i = 0; i != nterms; i++) {
      //mprintf("DEBUG: cmap term %f\n", terms[i] );
      prm.CMAP()[currentCmapIdx].AddToGrid( terms[i] );
    }
  } else if (currentCmapFlag == CMAP_RESLIST) {
    ArgList resnames( line );
    std::string rn = resnames.GetStringNext();
    while (!rn.empty()) {
      prm.CMAP()[currentCmapIdx].AddResName( rn );
      rn = resnames.GetStringNext();
    }
  } else if (currentCmapFlag == CMAP_TITLE) {
    prm.CMAP()[currentCmapIdx].SetTitle( line );
  }
  return 0;
}

/** Assign nonbond parameters from a NonbondSet to ParameterSet */
int AmberParamFile::assign_nb(ParameterSet& prm, NonbondSet const& nbset) const {
  // Nonbonds
  for (ParmHolder<LJparmType>::const_iterator it = nbset.LJ_.begin();
                                              it != nbset.LJ_.end(); ++it)
  {
    ParmHolder<AtomType>::iterator at = prm.AT().GetParam( it->first );
    if (at == prm.AT().end()) {
      mprinterr("Error: Nonbond parameters defined for previously undefined type '%s'.\n",
                *(it->first[0]));
      return 1;
    } 
    at->second.SetLJ().SetRadius( it->second.Radius() );
    at->second.SetLJ().SetDepth( it->second.Depth() );
  }
  return 0;
}

/** Assign off-diagonal nonbond parameters from an OffdiagNB array to ParameterSet */
int AmberParamFile::assign_offdiag(ParameterSet& prm, Oarray const& Offdiag) const {
  // Do off diagonal NB mods
  if (!Offdiag.empty()) {
    for (Oarray::const_iterator od = Offdiag.begin(); od != Offdiag.end(); ++od)
    {
      if (debug_ > 1) mprintf("DEBUG: Off diag %s %s\n", *(od->types_[0]), *(od->types_[1]));
      // Set type 1
      ParmHolder<AtomType>::iterator it = prm.AT().GetParam( od->types_[0] );
      if (it == prm.AT().end()) {
        mprinterr("Error: Off-diagonal nonbond parameters defined for previously undefined type '%s'.\n",
                  *(od->types_[0]));
        return 1;
      }
      it->second.SetLJ().SetRadius( od->LJ1_.Radius() );
      it->second.SetLJ().SetDepth( od->LJ1_.Depth() );
      // Set off-diagonal for type1-type2
      // FIXME different combine rules?
      prm.NB().AddParm(od->types_, od->LJ1_.Combine_LB(od->LJ2_), false);
    }
  } // END off-diagonal NB mods
  return 0;
}

/** Read parametrers from Amber frcmod file. */
int AmberParamFile::ReadFrcmod(ParameterSet& prm, FileName const& fname, int debugIn) const
{
  // Set wildcard character for dihedrals and impropers
  prm.DP().SetWildcard('X');
  prm.IP().SetWildcard('X');
  // Read title
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read anything from Amber FF file %s\n", fname.full());
    return 1;
  }
  std::vector<std::string> last_symbols;
  last_symbols.reserve(4);
  std::string title(ptr);
  mprintf("\tTitle: %s\n", title.c_str());
  prm.SetParamSetName( title );
  NonbondSet nbset(title);
  Oarray Offdiag;
  // Read file
  SectionType section = UNKNOWN;
  CmapType currentCmapFlag = CMAP_INITIAL;
  ptr = infile.Line();
  while (ptr != 0) {
    bool first_char_is_space = (*ptr == ' ');
    // Advance to first non-space char
    while (*ptr == ' ' && *ptr != '\0') ++ptr;
    // Is this a recognized section keyword?
    if (*ptr != '\0') {
      std::string line(ptr);
      if      (line.compare(0, 4, "MASS") == 0) section = ATYPE;
      else if (line.compare(0, 4, "BOND") == 0) section = BOND;
      else if (line.compare(0, 4, "ANGL") == 0) section = ANGLE;
      else if (line.compare(0, 4, "DIHE") == 0) section = DIHEDRAL;
      else if (line.compare(0, 4, "IMPR") == 0) section = IMPROPER;
      else if (line.compare(0, 4, "HBON") == 0) section = LJ1012;
      else if (line.compare(0, 4, "CMAP") == 0) section = CMAP;
      else if (line.compare(0, 6, "LJEDIT") == 0) {
        section = LJEDIT;
        prm.SetHasLJparams( true );
      } else if (line.compare(0, 4, "NONB") == 0) {
        section = NONBOND;
        prm.SetHasLJparams( true );
        // TODO check RE
      } else {
        //mprintf("DEBUG: Section %i: %s\n", (int)section, ptr);
        int err = 0;
        if (section == ATYPE)
          err = read_atype(prm, ptr);
        else if (section == BOND)
          err = read_bond(prm, ptr);
        else if (section == ANGLE)
          err = read_angle(prm, ptr);
        else if (section == DIHEDRAL)
          err = read_dihedral(prm, ptr, last_symbols, first_char_is_space);
        else if (section == IMPROPER)
          err = read_improper(prm, ptr);
        else if (section == LJ1012)
          err = read_lj1012(prm, ptr);
        else if (section == NONBOND)
          err = read_nb_RE(nbset, ptr);
        else if (section == LJEDIT)
          err = read_ljedit(Offdiag, ptr);
        else if (section == CMAP)
          err = read_cmap(prm, currentCmapFlag, line);
        if (err != 0) {
          mprinterr("Error: Reading line: %s\n", ptr);
          return 1;
        }
      }
    }
    ptr = infile.Line();
  }
  // Nonbonds
  if (assign_nb(prm, nbset)) return 1;
  // Off-diagonal NB modifications
  if (assign_offdiag(prm, Offdiag)) return 1;

  if (debug_ > 0) prm.Debug();
  infile.CloseFile();

  return 0;
}

/** Read parameters from Amber main FF parameter file. */
int AmberParamFile::ReadParams(ParameterSet& prm, FileName const& fname,
                               std::string const& nbsetnameIn, int debugIn) const
{
  // Set wildcard character for dihedrals and impropers
  prm.DP().SetWildcard('X');
  prm.IP().SetWildcard('X');
  // For files with > 1 set of NB params
  typedef std::vector<NonbondSet> NbSetArrayType;
  NbSetArrayType NBsets;
  // For holding equivalent NB type names
  typedef std::vector<NameType> Narray;
  typedef std::vector<Narray> XNarray;
  XNarray EquivalentNames;
  // For holding off-diagonal mods
  Oarray Offdiag;

  // Read title
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read anything from Amber FF file %s\n", fname.full());
    return 1;
  }
  std::string title(ptr);
  mprintf("\tTitle: %s\n", title.c_str());
  prm.SetParamSetName( title );
  std::vector<std::string> last_symbols;
  last_symbols.reserve(4);
  // Read file
  SectionType section = ATYPE;
  ptr = infile.Line();
  while (ptr != 0) {
    bool first_char_is_space = (*ptr == ' ');
    // Advance to first non-space char
    while (*ptr == ' ' && *ptr != '\0') ++ptr;
    //mprintf("DEBUG: First char: %c (%i)\n", *ptr, (int)*ptr);
    int read_err = 0;
    if (*ptr == '\0') {
      // Section Change
      if (section != UNKNOWN) {
        if (section == NONBOND) {
          // Do a lookahead to see if there are multiple NB sets.
          // It will either be another set, END, or LJEDIT
          ptr = infile.Line();
          while (*ptr == ' ' && *ptr != '\0') ++ptr;
          std::string nbline(ptr);
          if (nbline == "END") {
            mprintf("END\n");
            section = UNKNOWN;
          } else if (nbline == "LJEDIT") {
            section = LJEDIT;
            prm.SetHasLJparams( true );
          } // Otherwise assume another nonbond section
        } else {
          if (debug_ > 0) mprintf("SECTION %i change to %i\n", (int)section, (int)section + 1);
          section = (SectionType)((int)section + 1);
        }
      }
      // Special cases
      if (section == HYDROPHILIC) {
        // Special case: hydrophilic atom types. One line.
        // FORMAT(20(A2,2X))
        ptr = infile.Line();
        if (debug_ > 1) mprintf("DEBUG: Hydrophilic: %s\n", ptr);
        // Take advantage of the fact that we expect whitespace-delimiters
        ArgList hsymbols( ptr, " " );
        for (int iarg = 0; iarg != hsymbols.Nargs(); ++iarg)
          prm.AddHydrophilicAtomType( NameType( hsymbols[iarg] ) );
        section = (SectionType)((int)section + 1);
        // Look ahead. Older parm files have the hydrophilic section delimited
        // with a newline.
        ptr = infile.Line();
        while (*ptr == ' ' && *ptr != '\0') ++ptr;
        if (*ptr != '\0') continue;
      } else if (section == NONBOND) {
        prm.SetHasLJparams( true );
        // Special case: first read the line.
        // LABEL , KINDNB
        // FORMAT(A4,6X,A2)
        if (NBsets.empty())
          ptr = infile.Line();
        char nb_label[MAXSYMLEN], nb_kind[MAXSYMLEN];
        int nscan = sscanf(ptr, "%s %s", nb_label, nb_kind);
        if (nscan != 2) {
          mprinterr("Error: Expected nonbond label, nonbond kind, got %i elements.\n", nscan);
          return 1;
        }
        mprintf("DEBUG: NB label= %s  NB kind = %s\n", nb_label, nb_kind);
        if (nb_kind[0] != 'R' || nb_kind[1] != 'E') {
          mprinterr("Error: Nonbond parameters are not of type 'RE' (Rmin, Epsilon).\n");
          return 1;
        }
        NBsets.push_back( NonbondSet( std::string(nb_label) ) );
      }
    } else if (section == ATYPE) {
      read_err = read_atype(prm, ptr);
    } else if (section == BOND) {
      read_err = read_bond(prm, ptr);
    } else if (section == ANGLE) {
      read_err = read_angle(prm, ptr);
    } else if (section == DIHEDRAL) {
      read_err = read_dihedral(prm, ptr, last_symbols, first_char_is_space);
    } else if (section == IMPROPER) {
      read_err = read_improper(prm, ptr);
    } else if (section == LJ1012) {
      read_err = read_lj1012(prm, ptr);
    } else if (section == NB_EQUIV) {
      // EQUIVALENCING ATOM SYMBOLS FOR THE NON-BONDED 6-12 POTENTIAL PARAMETERS
      // IORG , IEQV(I) , I = 1 , 19
      // FORMAT(20(A2,2X))
      if (debug_ > 1) mprintf("DEBUG: Nonbond equiv: %s\n", ptr);
      EquivalentNames.push_back( Narray() );
      ArgList equiv_line( ptr, " " );
      for (int iarg = 0; iarg != equiv_line.Nargs(); iarg++)
        EquivalentNames.back().push_back( equiv_line[iarg] );
    } else if (section == NONBOND) {
      // ***** ONLY IF KINDNB .EQ. 'RE' *****
      read_err = read_nb_RE(NBsets.back(), ptr);
    } else if (section == LJEDIT) {
      read_err = read_ljedit(Offdiag, ptr);
    }
    if (read_err != 0) {
      mprinterr("Error: Reading line: %s\n", ptr);
      return 1;
    } 
    ptr = infile.Line();
  } // END loop over file.

  // Deal with nonbond and equivalence
  if (!NBsets.empty()) {
    int nbsetidx = 0;
    if (NBsets.size() > 1 || !nbsetnameIn.empty()) {
      // Choose from multiple nbsets
      if (nbsetnameIn.empty()) {
        mprinterr("Error: Parm set in file '%s' contains %zu nonbonded parameter sets\n"
                  "Error:  but no set has been specified.\n", fname.full(), NBsets.size());
        mprinterr("Error: Need to specify one of");
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
          mprinterr(" %s", it->name_.c_str());
        mprinterr("\n");
      } else {
        nbsetidx = -1;
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it) {
          if (it->name_ == nbsetnameIn) {
            nbsetidx = it - NBsets.begin();
            break;
          }
        }
        if (nbsetidx < 0) {
          mprinterr("Error: Nonbonded set '%s' not found.\n", nbsetnameIn.c_str());
          mprinterr("Error: Need to specify one of");
          for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
            mprinterr(" %s", it->name_.c_str());
          mprinterr("\n");
          return 1;
        }
      }
    }
    mprintf("\tUsing nonbonded parm set: %s\n", NBsets[nbsetidx].name_.c_str());
    prm.SetNbParamName( NBsets[nbsetidx].name_ );
    if (assign_nb(prm, NBsets[nbsetidx])) return 1;
    
    // Do equivalent atoms.
    for (XNarray::const_iterator equivAts = EquivalentNames.begin();
                                 equivAts != EquivalentNames.end(); ++equivAts)
    {
      // First name is the type to copy TODO check size?
      Narray::const_iterator typeName = equivAts->begin();
      ParmHolder<AtomType>::const_iterator at0 = prm.AT().GetParam( *typeName );
      if (at0 == prm.AT().end()) {
        mprinterr("Error: Equivalent atom type '%s' not found.\n", *(*typeName) );
        return 1;
      }
      ++typeName;
      for (; typeName != equivAts->end(); ++typeName) {
        ParmHolder<AtomType>::iterator at1 = prm.AT().GetParam( *typeName );
        if (at1 == prm.AT().end()) {
          mprinterr("Error: Equivalent atom type '%s' (base '%s') not found.\n", *(*typeName), *(at0->first[0]));
          return 1;
        }
        if (debug_ > 1) mprintf("DEBUG: Equiv '%s' => '%s'\n", *(at0->first[0]), *(*typeName));
        at1->second.SetLJ().SetRadius( at0->second.LJ().Radius() );
        at1->second.SetLJ().SetDepth( at0->second.LJ().Depth() );
      }
    } // END loop over EquivalentNames
  } // END nonbond parameters

  // Do off diagonal NB mods
  if (assign_offdiag(prm, Offdiag)) return 1;

  if (debug_ > 0) prm.Debug();
  infile.CloseFile();
  return 0;
}

// DataIO_AmberFF::WriteData()
int AmberParamFile::WriteParams(ParameterSet& prm, FileName const& fname, int debugIn) const
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) return 1;
  // 1 - Title
  std::string title;
  if (prm.ParamSetName().empty())
    title.assign("Parameters written from CPPTRAJ");
  else
    title = prm.ParamSetName();
  outfile.Printf("%s\n", title.c_str());
  // 2 - Atom symbols and masses
  for (ParmHolder<AtomType>::const_iterator at = prm.AT().begin(); at != prm.AT().end(); ++at)
  {
    std::string asym = at->first[0].Truncated();
    if (asym.size() > 2)
      mprintf("Warning: Atom symbol %s is larger than 2 characters, which breaks Amber FF format.\n");
    outfile.Printf("%-2s", asym.c_str());
    if (at->second.Mass() < 10.0)
      outfile.Printf(" %-10.3f", at->second.Mass());
    else if (at->second.Mass() < 100.0)
      outfile.Printf(" %-10.2f", at->second.Mass());
    else
      outfile.Printf(" %-10.1f", at->second.Mass());
    outfile.Printf(" %10.3f\n", at->second.Polarizability());

    //outfile.Printf("%-2s  %-10.2f %-10.2f\n", asym.c_str(), at->second.Mass(), at->second.Polarizability());
  }

  return 0;
}
