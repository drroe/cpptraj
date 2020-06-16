#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToDouble
#include "DataSet_double.h"

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
DataIO_Mdout::DataIO_Mdout() {
  for (unsigned int idx = 0; idx != (unsigned int)N_FIELDTYPES; idx++)
  {
    SetAspects_.push_back( std::string(fieldTypeStr_[idx]) );
    SetOffsets_.push_back(0);
    EneSets_.push_back(0);
  }
  // Populate the term name to index map. In some cases, multiple term names
  // map to the same index.
  termIdxMap_.insert(NameIdxPair("Etot", ETOT));
  termIdxMap_.insert(NameIdxPair("EPtot") return EPtot;
  termIdxMap_.insert(NameIdxPair("GMAX") return GMAX; // Not necessary?
  termIdxMap_.insert(NameIdxPair("BOND") return BOND;
  termIdxMap_.insert(NameIdxPair("ANGLE") return ANGLE;
  termIdxMap_.insert(NameIdxPair("DIHED") return DIHED;
  termIdxMap_.insert(NameIdxPair("VDWAALS") return VDWAALS;
  termIdxMap_.insert(NameIdxPair("EEL", EEL));
  termIdxMap_.insert(NameIdxPair("EELEC", EEL));
  termIdxMap_.insert(NameIdxPair("EGB") return EGB;
  termIdxMap_.insert(NameIdxPair("EPB") return EPB;
  termIdxMap_.insert(NameIdxPair("ECAVITY") return ECAVITY;
  termIdxMap_.insert(NameIdxPair("EDISPER") return EDISPER;
  if ((Name[0]=="1-4" && Name[1]=="VDW") || (Name[0]=="1-4" && Name[1]=="NB")) return VDW14;
  if  (Name[0]=="1-4" && Name[1]=="EEL") return EEL14;
  termIdxMap_.insert(NameIdxPair("RESTRAINT") return RESTRAINT;
  termIdxMap_.insert(NameIdxPair("EAMBER") return EAMBER;
  termIdxMap_.insert(NameIdxPair("Density") return Density;
  termIdxMap_.insert(NameIdxPair("RMS") return RMS; // Not necessary?
  termIdxMap_.insert(NameIdxPair("EKtot") return EKtot;
  termIdxMap_.insert(NameIdxPair("ESURF") return ESURF;
  termIdxMap_.insert(NameIdxPair("EAMD_BOOST") return EAMD_BOOST;
  termIdxMap_.insert(NameIdxPair("VOLUME") return VOLUME;
  termIdxMap_.insert(NameIdxPair("TEMP(K)") return TEMP;
  termIdxMap_.insert(NameIdxPair("PRESS") return PRESS;
  termIdxMap_.insert(NameIdxPair("DV/DL") return DVDL;

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
    // New term. Create new term.
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

DataIO_Mdout::FieldType DataIO_Mdout::getEindex(Sarray const& Name) {
  //mprintf("DEBUG:\tgetEindex(%s,%s)\n", Name[0].c_str(), Name[1].c_str());
  return N_FIELDTYPES;
}

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

/// Array of amber energy terms
typedef std::vector<AmberEterm> EtermArray;

/** Parse the incoming line into individual energy terms. In the Amber
  * MDOUT file this takes the form of 1 or more '<name>=<value>'. There 
  * may or may not be a space after <name>/before <value>. There should
  * be spaces between different energy terms though.
  */
EtermArray GetAmberEterms(const char* ptr) {
  //mprintf("DBG: [%s]\n", ptr);
  EtermArray Terms;
  if (ptr == 0) return Terms;
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
        return Terms;
      } else {
        //mprintf("DBG: val= %c\n", *val);
        // Search for next whitespace or line end.
        const char* end = val + 1;
        while (*end != ' ' && *end != '\0' && *end != '\n' && *end != '\r') ++end;
        // Term is now complete. Convert.
        std::string valstr(val, end);
        //mprintf("DBG: valstr= '%s'\n", valstr.c_str());
        if (!validDouble(valstr)) {
          mprintf("Warning: Invalid number detected: %s\n", valstr.c_str());
        } else {
          Terms.push_back( AmberEterm(
                             NoTrailingWhitespace(std::string(beg,eq)),
                             convertToDouble(valstr)) );
        }
        beg = end;
      }
    }
  } // END loop over line
  // DEBUG
  mprintf("DBG: %s\n", ptr);
  mprintf("DBG: %zu terms:", Terms.size());
  for (EtermArray::const_iterator it = Terms.begin(); it != Terms.end(); ++it)
    mprintf(" %s = %g", it->Name().c_str(), it->Value());
  mprintf("\n");
  return Terms;
}

// DataIO_Mdout::ReadData()
int DataIO_Mdout::ReadData(FileName const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  mprintf("\tReading from mdout file: %s\n", fname.full());
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  const char* ptr = buffer.Line();
  if (ptr == 0) {
    mprinterr("Error: Nothing in MDOUT file: %s\n", fname.full());
    return 1;
  }
  // ----- PARSE THE INPUT SECTION ----- 
  int imin = -1;           // imin for this file
  const char* Trigger = 0; // Trigger for storing energies, must be 8 chars long.
  int frame = 0;           // Frame counter for this file
  double dt = 1.0;         // Timestep for this file (MD)
  double t0 = 0.0;         // Initial time for this file (MD)
  int ntpr = 1;            // Value of ntpr
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
        imin = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: imin is %i\n", imin);
        // Set a trigger for printing. For imin5 this is the word minimization.
        // For imin0 or imin1 this is NSTEP.
        if      (imin==0) Trigger = " NSTEP =";
        else if (imin==1) Trigger = "   NSTEP";
        else if (imin==5) Trigger = "minimiza";
        // Since imin0 and imin1 first trigger has no data, set frame 1 lower.
        if (imin==1 || imin==0) frame = -1;
      } else if (mdin_args[col] == "dt") {
        dt = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: dt is %f\n", dt);
      } else if (mdin_args[col] == "t") {
        t0 = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: t is %f\n", t0);
      } else if (mdin_args[col] == "ntpr") {
        ntpr = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: ntpr is %i\n", ntpr);
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
        sscanf(ptr, " begin time read from input coords = %lf", &t0);
        if (debug_ > 0) mprintf("\t\tMD restart initial time= %f\n", t0);
      }
    }
  }
  if (ptr == 0) return EOF_ERROR();
  // ----- PARSE THE RESULTS SECTION -----
  bool finalE = false;
  int nstep;
  int minStep = 0; // For imin=1 only
  if (irest == 0)
    nstep = 0;
  else
    nstep = ntpr;
  double Energy[N_FIELDTYPES];
  std::fill( Energy, Energy+N_FIELDTYPES, 0.0 );
  std::vector<bool> EnergyExists(N_FIELDTYPES, false);
  DataSetList::Darray TimeVals;
  DataSetList::DataListType inputSets(N_FIELDTYPES, 0);
  Sarray Name(2);
  double time = 0.0;
  while (ptr != 0) {
    // Check for end of imin 0 or 1 run; do not record Average and Stdevs
    if ( (imin == 1 && (strncmp(ptr, "                    FINAL", 25) == 0 ||
                        strncmp(ptr, "   5.  TIMINGS",            14) == 0   )) ||
         (imin == 0 && strncmp(ptr, "      A V", 9) == 0))
      finalE = true;
    // Check for '| TI region  2' to prevent reading duplicate energies
    if ( strncmp(ptr, "| TI region  2", 14) == 0 ) {
      while (ptr != 0 && !(ptr[0] == ' ' && ptr[1] == '-'))
        ptr = buffer.Line();
      if (ptr == 0) return EOF_ERROR();
    }
    // Record set for energy post-processing
    if (imin == 5 && strncmp(ptr, "minimizing", 10) == 0)
      nstep = atoi( ptr + 22 );
    // MAIN OUTPUT ROUTINE
    // If the trigger has been reached print output.
    // For imin0 and imin1 the first trigger will have no data.
    // If the end of the file has been reached print then exit.
    if ( strncmp(ptr, Trigger, 8) == 0 || finalE ) {
      if (frame > -1) {
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
        TimeVals.push_back( time );
        nstep += ntpr;
      }
      frame++;
      if (finalE) break;
    }
    // Check for NSTEP in minimization or post-processing. Values will be
    // on the next line. NOTE: NSTEP means something different for imin=5.
    if ((imin == 1 || imin == 5) && strncmp(ptr, "   NSTEP", 8) == 0) {
      ptr = buffer.Line(); // Get next line
      //sscanf(ptr, " %6lf    %13lE  %13lE  %13lE", Energy+NSTEP, Energy+EPtot, Energy+RMS, Energy+GMAX);
      sscanf(ptr, " %i %lE %lE %lE", &minStep, Energy+EPtot, Energy+RMS, Energy+GMAX);
      EnergyExists[EPtot] = true;
      EnergyExists[RMS] = true;
      EnergyExists[GMAX] = true;
      ptr = buffer.Line();
    }
    // Tokenize line, scan through until '=' is reached; value after is target.
    GetAmberEterms(ptr);
    int ntokens = buffer.TokenizeLine(" ");
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
    }
    // Set time
    switch (imin) {
      case 5: time = (double)nstep + t0; break;
      case 1: time = (double)minStep + t0; break;
      case 0: time = ((double)nstep * dt) + t0; break;
    }
    // Read in next line
    ptr = buffer.Line();
  }
  mprintf("\t%i frames\n", frame);
  buffer.CloseFile();
  std::string Xlabel;
  if      (imin == 5) Xlabel.assign("Set");
  else if (imin == 1) Xlabel.assign("Nstep");
  else                Xlabel.assign("Time"); // imin == 0
  if (datasetlist.AddOrAppendSets( Xlabel, TimeVals, inputSets )) return 1;
  return 0;
}
