#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToDouble
#include "DataSet_double.h"
#include "DataSet_Mesh.h"

/** Must correspond to FieldType enum:
  * ETOT = 0, EPTOT, GMAX,        BOND,
  * ANGLE,    DIHED, VDWAALS,     EEL,        EGB,    EPB, ECAVITY, EDISPER,
  * VDW14,    EEL14, RESTRAINT,   EAMBER,     DENSITY,
  * RMS,      EKTOT, ESURF,       EAMD_BOOST, VOLUME, TEMP,
  * PRESS,    DVDL,  N_FIELDTYPES
  */
const char* DataIO_Mdout::fieldTypeStr_[] = {
  "Etot",   "EPtot",  "GMAX",  "BOND",
  "ANGLE",  "DIHED",  "VDW",   "EELEC",      "EGB",     "EPB", "ECAVITY", "EDISPER",
  "VDW1-4", "EEL1-4", "RST",   "EAMBER",     "Density",
  "RMS",    "EKtot",  "ESURF", "EAMD_BOOST", "VOLUME",  "TEMP",
  "PRESS",  "DVDL",   0
};

/** CONSTRUCTOR - populate initial map and set aspects with known fields. */
DataIO_Mdout::DataIO_Mdout() :
  imin_(-1)
{
  // Load "known" terms
  for (unsigned int idx = 0; idx != (unsigned int)N_FIELDTYPES; idx++)
  {
    SetAspects_.push_back( std::string(fieldTypeStr_[idx]) );
    SetOffsets_.push_back(0);
    EneSets_.push_back(0);
  }
  // Populate the term name to index map. In some cases, multiple term names
  // map to the same index.
  termIdxMap_.insert(NameIdxPair("Etot", ETOT));
  termIdxMap_.insert(NameIdxPair("EPtot", EPTOT));
  termIdxMap_.insert(NameIdxPair("GMAX", GMAX)); // Not necessary?
  termIdxMap_.insert(NameIdxPair("BOND", BOND));
  termIdxMap_.insert(NameIdxPair("ANGLE", ANGLE));
  termIdxMap_.insert(NameIdxPair("DIHED", DIHED));
  termIdxMap_.insert(NameIdxPair("VDWAALS", VDWAALS));
  termIdxMap_.insert(NameIdxPair("EEL", EEL));
  termIdxMap_.insert(NameIdxPair("EELEC", EEL));
  termIdxMap_.insert(NameIdxPair("EGB", EGB));
  termIdxMap_.insert(NameIdxPair("EPB", EPB));
  termIdxMap_.insert(NameIdxPair("ECAVITY", ECAVITY));
  termIdxMap_.insert(NameIdxPair("EDISPER", EDISPER));
  termIdxMap_.insert(NameIdxPair("1-4 VDW", VDW14));
  termIdxMap_.insert(NameIdxPair("1-4 NB", VDW14));
  termIdxMap_.insert(NameIdxPair("1-4 EEL", EEL14));
  termIdxMap_.insert(NameIdxPair("RESTRAINT", RESTRAINT));
  termIdxMap_.insert(NameIdxPair("EAMBER", EAMBER));
  termIdxMap_.insert(NameIdxPair("Density", DENSITY));
  termIdxMap_.insert(NameIdxPair("RMS", RMS)); // Not necessary?
  termIdxMap_.insert(NameIdxPair("EKtot", EKTOT));
  termIdxMap_.insert(NameIdxPair("ESURF", ESURF));
  termIdxMap_.insert(NameIdxPair("EAMD_BOOST", EAMD_BOOST));
  termIdxMap_.insert(NameIdxPair("VOLUME", VOLUME));
  termIdxMap_.insert(NameIdxPair("TEMP(K)", TEMP));
  termIdxMap_.insert(NameIdxPair("PRESS", PRESS));
  termIdxMap_.insert(NameIdxPair("DV/DL", DVDL));
}

// DataIO_Mdout::ID_DataFormat()
bool DataIO_Mdout::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  bool isMdout = false;
  std::string line = infile.GetLine();
  if (line[0] == '\n') {
    line = infile.GetLine();
    if (line.compare(0, 15, "          -----") == 0) {
      line = infile.GetLine();
      if (line.compare(0, 15, "          Amber") == 0)
        isMdout = true;
    }
  }
  infile.CloseFile();
  return isMdout;
}

static inline int EOF_ERROR() {
  mprinterr("Error: Unexpected EOF in MDOUT file.\n");
  return 1;
}

/** \return Known idx of given term if known, next unknown idx otherwise.
  */ 
unsigned int DataIO_Mdout::getTermIdx(std::string const& name) {
  NameIdxMap::const_iterator it = termIdxMap_.find( name );
  if (it == termIdxMap_.end()) {
    // Create new term.
    unsigned int newIdx = SetAspects_.size();
    SetAspects_.push_back( name ); // TODO replace spaces?
    SetOffsets_.push_back( 0 );
    EneSets_.push_back( 0 );
    termIdxMap_.insert( NameIdxPair(name, newIdx) );
    return newIdx;
  } else {
    // Existing term. just return the index.
    return it->second;
  }
}
/*
/// Hold a single amber energy term
class AmberEterm {
  public:
    AmberEterm(std::string const& nameIn, double valueIn) :
      name_(nameIn), value_(valueIn) {}
    std::string const& Name() const { return name_; }
    double Value() const { return value_; }
  private:
    std::string name_; ///< ASCII text before '=' describing the term.
    double value_;     ///< The value of the energy term.
};
*/
/// Array of amber energy terms
//typedef std::vector<AmberEterm> EtermArray;

/** Add given value to specified set. */
int DataIO_Mdout::AddData(unsigned int idx, double value,
                           std::string const& dsname, DataSetList& dsl)
{
  // Get the set
  DataSet* ds = 0;
  if (EneSets_[idx] == 0) {
    // Check if the set exists
    MetaData md(dsname, SetAspects_[idx]);
    ds = dsl.CheckForSet( md );
    if (ds == 0) {
      // Need to create the data set
      ds = dsl.AddSet(DataSet::XYMESH, MetaData(dsname, SetAspects_[idx]));
      if (ds == 0) {
        mprinterr("Error: Could not allocate set %s[%s]\n", dsname.c_str(),
                  SetAspects_[idx].c_str());
        return 1;
      }
      mprintf("\tCreated new set: %s\n", ds->legend());
      SetOffsets_[idx] = t0_;
      Dimension Xdim;
      if      (imin_ == 5) Xdim = Dimension(1, 1, "Set");
      else if (imin_ == 1) Xdim = Dimension(1, 1, "Nstep");
      else                 Xdim = Dimension(t0_, ntpr_*dt_, "Time"); // imin == 0
      ds->SetDim(Dimension::X, Xdim);
    } else {
      // Needs to be 1D scalar
      if (ds->Group() != DataSet::SCALAR_1D) {
        mprinterr("Error: Set %s is wrong type for MD data append.\n", ds->legend());
        return 1;
      }
      DataSet_1D const& ds1 = static_cast<DataSet_1D const&>( *ds );
      if (ds1.Size() > 0) {
        SetOffsets_[idx] = ds1.Xcrd( ds1.Size() - 1 );
      }
      mprintf("\tAppending to set %s starting from X= %g\n", ds->legend(), SetOffsets_[idx]);
    }
    EneSets_[idx] = ds;
  }

  if (EneSets_[idx]->Type() == DataSet::XYMESH) {
    double Xval;
    if (imin_ == 5)
      Xval = (double)nstep_ + SetOffsets_[idx];
    else if (imin_ == 1)
      Xval = (double)minStep_ + SetOffsets_[idx];
    else
      Xval = ((double)nstep_ * dt_) + SetOffsets_[idx];
    ((DataSet_Mesh*)EneSets_[idx])->AddXY(Xval, value);
  }

  return 0;
}
    

/** Parse the incoming line into individual energy terms. In the Amber
  * MDOUT file this takes the form of 1 or more '<name>=<value>'. There 
  * may or may not be a space after <name>/before <value>. There should
  * be spaces between different energy terms though.
  */
int DataIO_Mdout::GetAmberEterms(const char* ptr, std::string const& dsname,
                                 DataSetList& datasetlist)
{
  //mprintf("DBG: [%s]\n", ptr);
  if (ptr == 0 || ptr[0] == '|') return 0;
  const char* beg = ptr;
  //          111111111122222222223
  //0123456789012345678901234567890
  // NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =   435.99  PRESS =-10207.6
  bool eol = false;
  while (!eol) {
    // Skip leading whitespace
    while (*beg == ' ' && *beg != '\0') ++beg;
    if (*beg == '\0') {
      // Line is blank or no more terms. Bail out.
      break;
    }
    //mprintf("DBG: beg= %c\n", *beg);
    // Search for next '='
    const char* eq = beg + 1;
    while (*eq != '=' && *eq != '\0') ++eq;
    if (*eq == '\0')
      eol = true;
    else {
      // Search for end token. Start just after '='.
      const char* val = eq + 1;
      // Skip leading whitespace
      while (*val == ' ' && *val != '\0') ++val;
      if (*val == '\0') {
        eol = true;
        mprintf("Warning: EOL encountered before energy term could be read.\n");
        return 1;
      } else {
        //mprintf("DBG: val= %c\n", *val);
        // Search for next whitespace or line end.
        const char* end = val + 1;
        while (*end != ' ' && *end != '\0' && *end != '\n' && *end != '\r') ++end;
        // Term is now complete. Convert.
        std::string valstr(val, end);
        //mprintf("DBG: valstr= '%s'\n", valstr.c_str());
        std::string termName = NoTrailingWhitespace(std::string(beg,eq));
        if (!reachedNstep_ && termName != "NSTEP") return 0;
        reachedNstep_ = true;
        if (!validDouble(valstr)) {
          mprintf("Warning: Invalid number detected: %s = %s\n", termName.c_str(), valstr.c_str());
        } else {
          mprintf("DBG: %s = %s\n", termName.c_str(), valstr.c_str());
          // Special cases
          if (termName == "NSTEP")
            nstep_ = convertToInteger(valstr);
          else if (termName != "TIME(PS)") {
            unsigned int idx = getTermIdx(termName);
            AddData(idx, convertToDouble(valstr), dsname, datasetlist);
          }
        }
        beg = end;
      }
    }
  } // END loop over line
  
  return 0;
}

// DataIO_Mdout::ReadData()
int DataIO_Mdout::ReadData(FileName const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  mprintf("\tReading from mdout file: %s\n", fname.full());
  EneSets_.assign(EneSets_.size(), 0);
  reachedNstep_ = false;
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  const char* ptr = buffer.Line();
  if (ptr == 0) {
    mprinterr("Error: Nothing in MDOUT file: %s\n", fname.full());
    return 1;
  }
  // ----- PARSE THE INPUT SECTION ----- 
  imin_ = -1;              // imin for this file
  const char* Trigger = 0; // Trigger for storing energies, must be 8 chars long.
  int frame = 0;           // Frame counter for this file
  dt_ = 1.0;         // Timestep for this file (MD)
  t0_ = 0.0;         // Initial time for this file (MD)
  ntpr_ = 1;            // Value of ntpr
  int irest = 0;           // Value of irest
  while ( ptr != 0 && strncmp(ptr, "   2.  CONTROL  DATA", 20) != 0 )
    ptr = buffer.Line();
  if (ptr == 0) return EOF_ERROR();
  // Determine whether this is dynamics or minimization, get dt
  ptr = buffer.Line(); // Dashes 
  ptr = buffer.Line(); // Blank 
  ptr = buffer.Line(); // title line
  while ( strncmp(ptr, "   3.  ATOMIC", 13) != 0 ) 
  {
    ArgList mdin_args( ptr, " ,=" ); // Remove commas, equal signs
    // Scan for stuff we want
    //mprintf("DEBUG:\tInput[%i] %s", mdin_args.Nargs(), mdin_args.ArgLine());
    for (int col=0; col < mdin_args.Nargs(); col += 2) {
      int col1 = col + 1;
      if (mdin_args[col] == "imin") {
        imin_ = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: imin is %i\n", imin_);
        // Set a trigger for printing. For imin=5 this is the word minimization.
        // For imin=0 or imin=1 this is NSTEP.
        if      (imin_==0) Trigger = " NSTEP =";
        else if (imin_==1) Trigger = "   NSTEP";
        else if (imin_==5) Trigger = "minimiza";
        // Since imin0 and imin1 first trigger has no data, set frame 1 lower.
        if (imin_==1 || imin_==0) frame = -1;
      } else if (mdin_args[col] == "dt") {
        dt_ = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: dt is %f\n", dt_);
      } else if (mdin_args[col] == "t") {
        t0_ = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: t is %f\n", t0_);
      } else if (mdin_args[col] == "ntpr") {
        ntpr_ = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: ntpr is %i\n", ntpr_);
      } else if (mdin_args[col] == "irest") {
        irest = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: irest is %i\n", irest);
      }
    }
    ptr = buffer.Line();
    if (ptr == 0) return EOF_ERROR();
  }
  if (Trigger == 0) {
    mprinterr("Error: Could not determine whether MDOUT is md, min, or post-process.\n");
    return 1;
  }
  // ----- PARSE THE ATOMIC ... SECTION -----
  while ( ptr != 0 && strncmp(ptr, "   4.  RESULTS", 14) != 0 )
  {
    ptr = buffer.Line();
    // If run is a restart, set the initial time value.
    if (irest == 1) {
      if (strncmp(ptr, " begin time", 11) == 0) {
        sscanf(ptr, " begin time read from input coords = %lf", &t0_);
        if (debug_ > 0) mprintf("\t\tMD restart initial time= %f\n", t0_);
      }
    }
  }
  if (ptr == 0) return EOF_ERROR();
  // ----- PARSE THE RESULTS SECTION -----
  bool finalE = false;
  minStep_ = -1; // For imin=1 only
  if (irest == 0)
    nstep_ = 0;
  else
    nstep_ = ntpr_;
  //double time = 0.0;
  // Loop over remaining lines of the file.
  while (ptr != 0) {
    // Check for end of imin 0 or 1 run; do not record Average and Stdevs
    if ( (imin_ == 1 && (strncmp(ptr, "                    FINAL", 25) == 0 ||
                         strncmp(ptr, "   5.  TIMINGS",            14) == 0   )) ||
         (imin_ == 0 && strncmp(ptr, "      A V", 9) == 0))
      finalE = true;
    // Check for '| TI region  2' to prevent reading duplicate energies
    if ( strncmp(ptr, "| TI region  2", 14) == 0 ) {
      while (ptr != 0 && !(ptr[0] == ' ' && ptr[1] == '-'))
        ptr = buffer.Line();
      if (ptr == 0) return EOF_ERROR();
    }
    // Record set for energy post-processing
    if (imin_ == 5 && strncmp(ptr, "minimizing", 10) == 0) {
      reachedNstep_ = true;
      nstep_ = atoi( ptr + 22 );
    }
    // MAIN OUTPUT ROUTINE
    // If the trigger has been reached print output.
    // For imin0 and imin1 the first trigger will have no data.
    // If the end of the file has been reached print then exit.
    if ( strncmp(ptr, Trigger, 8) == 0 || finalE ) {
      /*if (frame > -1) {
        // Store all energies present.
        for (int i = 0; i < (int)N_FIELDTYPES; i++) {
          if (EnergyExists[i]) {
            if (inputSets[i] == 0) {
              MetaData md( dsname, Enames[i] );
              md.SetLegend( dsname + "_" + Enames[i] );
              inputSets[i] = new DataSet_double();
              inputSets[i]->SetMeta( md );
            }
            // Since energy terms can appear and vanish over the course of the
            // mdout file, resize if necessary.
            if (frame > (int)inputSets[i]->Size())
              ((DataSet_double*)inputSets[i])->Resize( frame );
            ((DataSet_double*)inputSets[i])->AddElement( Energy[i] );
          }
        }
        TimeVals.push_back( time );*/
        nstep_ += ntpr_;
      //}
      frame++;
      if (finalE) break;
    }
    // Check for NSTEP in minimization or post-processing. Values will be
    // on the next line. NOTE: NSTEP means something different for imin=5.
    if ((imin_ == 1 || imin_ == 5) && strncmp(ptr, "   NSTEP", 8) == 0) {
      reachedNstep_ = true;
      ptr = buffer.Line(); // Get next line
      //sscanf(ptr, " %6lf    %13lE  %13lE  %13lE", Energy+NSTEP, Energy+EPtot, Energy+RMS, Energy+GMAX);
      double tmpbuf[3];
      sscanf(ptr, " %i %lE %lE %lE", &minStep_, tmpbuf, tmpbuf+1, tmpbuf+2);
      AddData((unsigned int)EPTOT, tmpbuf[0], dsname, datasetlist);
      AddData((unsigned int)RMS,   tmpbuf[1], dsname, datasetlist);
      AddData((unsigned int)GMAX,  tmpbuf[2], dsname, datasetlist);
      ptr = buffer.Line();
    }
    // Tokenize line, scan through until '=' is reached; value after is target.
    GetAmberEterms(ptr, dsname, datasetlist);
    /*int ntokens = buffer.TokenizeLine(" ");
    if (ntokens > 0) {
      int nidx = 0;
      Name[0].clear();
      Name[1].clear();
      for (int tidx = 0; tidx < ntokens; tidx++) {
        const char* tkn = buffer.NextToken();
        if (tkn[0] == '=') {
          FieldType Eindex = getEindex(Name);
          tkn = buffer.NextToken();
          ++tidx;
          if (tkn == 0)
            mprintf("Warning: No numerical value, line %i column %i. Skipping.\n",
                    buffer.LineNumber(), tidx+1);
          else if (tkn[0] == '*' || tkn[0] == 'N') // Assume if number begins with N it is NaN
            mprintf("Warning: Numerical overflow detected, line %i column %i. Skipping.\n",
                     buffer.LineNumber(), tidx+1);
          else {
            if (Eindex != N_FIELDTYPES) {
              Energy[Eindex] = atof( tkn );
              EnergyExists[Eindex] = true;
            }
          }
          nidx = 0;
          Name[0].clear();
          Name[1].clear();
        } else {
          if (nidx > 1) break; // Two tokens, no '=' found. Not an E line.
          Name[nidx++].assign( tkn );
          //mprintf("DEBUG: %s\n", tkn);
        }
      }
    }*/
    // Set time

    // Read in next line
    ptr = buffer.Line();
  }
  mprintf("\t%i frames\n", frame);
  buffer.CloseFile();
  //if (datasetlist.AddOrAppendSets( Xlabel, TimeVals, inputSets )) return 1;
  return 0;
}
