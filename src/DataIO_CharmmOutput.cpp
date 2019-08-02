#include <cstdio> // sscanf
#include <cctype> // isspace
#include "DataIO_CharmmOutput.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_double.h"

/// CONSTRUCTOR
DataIO_CharmmOutput::DataIO_CharmmOutput() { }

// DataIO_CharmmOutput::ID_DataFormat()
bool DataIO_CharmmOutput::ID_DataFormat(CpptrajFile& infile)
{
  // Assume file set up for read.
  if (infile.OpenFile()) return false;
  bool isCharmmOut = false;
  // Scan first 3 lines. That should be plenty.
  for (int num = 0; num != 3; num++)
  {
    ArgList line( infile.GetLine() );
    if (line.Nargs() >= 3 && line[0] == "Chemistry" && line[2] == "HARvard") {
      isCharmmOut = true;
      break;
    } else if (line.Nargs() >=1 && line[0] == "CHARMM>") {
      isCharmmOut = true;
      break;
    }
  }
  infile.CloseFile();

  return isCharmmOut;
}

// DataIO_CharmmOutput::ReadHelp()
void DataIO_CharmmOutput::ReadHelp() { }

// DataIO_CharmmOutput::processReadArgs()
int DataIO_CharmmOutput::processReadArgs(ArgList& argIn)
{
  return 0;
}

// DataIO_CharmmOutput::ReadData()
int DataIO_CharmmOutput::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  const char* ptr = buffer.Line();
  // Figure out if this contains DYNA/ENER/MINI
  bool DYNA = false;
  bool ENER = false;
  bool MINI = false;
  unsigned int headerOffset = 0;
  // Out of all present, prefer DYNA
  while (ptr != 0) {
    if ( ptr[0] == 'D' && ptr[1] == 'Y' && ptr[2] == 'N' && ptr[3] == 'A' ) {
      // Reached DYNA which is preferred. Break now.
      DYNA = true;
      break;
    }
    if ( ptr[0] == 'E' && ptr[1] == 'N' && ptr[2] == 'E' && ptr[3] == 'R' ) {
      // Reached ENER but DYNA/MINI is preferred, so keep looking.
      ENER = true;
    }
    if ( ptr[0] == 'M' && ptr[1] == 'I' && ptr[2] == 'N' && ptr[3] == 'I' ) {
      // Reached MINI but DYNA is preferred, so keep looking.
      MINI = true;
    }
    ptr = buffer.Line();
  }
  if (!DYNA && !ENER && !MINI) {
    mprinterr("Error: 'DYNA'/'ENER'/'MINI' not found in output file '%s'.\n", fname.full());
    return 1;
  }
  if (DYNA) {
    if (ENER) mprintf("Warning: Both DYNA and ENER found. Only reading output from DYNA.\n");
    if (MINI) mprintf("Warning: Both DYNA and MINI found. Only reading output from DYNA.\n");
    ENER = false;
    MINI = false;
    headerOffset = 14;
  } else if (MINI) {
    if (ENER) mprintf("Warning: Both MINI and ENER found. Only reading output from MINI.\n");
    DYNA = false;
    ENER = false;
    headerOffset = 15;
    // Rewind and go to MINI 
    buffer.CloseFile();
    if (buffer.OpenFileRead( fname )) return 1;
    ptr = buffer.Line();
    while (ptr != 0) {
      if ( ptr[0] == 'M' && ptr[1] == 'I' && ptr[2] == 'N' && ptr[3] == 'I' ) break;
      ptr = buffer.Line();
    }
  } else if (ENER) {
    DYNA = false;
    MINI = false;
    headerOffset = 16;
    // Rewind and go to ENER
    buffer.CloseFile();
    if (buffer.OpenFileRead( fname )) return 1;
    ptr = buffer.Line();
    while (ptr != 0) {
      if ( ptr[0] == 'E' && ptr[1] == 'N' && ptr[2] == 'E' && ptr[3] == 'R' ) break;
      ptr = buffer.Line();
    }
  }

  // Figure out what terms we have. Format is 'DYNA <Type>: E0 E1 ...'
  // Terms start at character headerOffset.
  typedef std::vector<std::string> Sarray;
  Sarray Terms;
  Sarray LineHeaders;
  int timeIdx = -1;
  typedef std::vector<int> Iarray;
  Iarray nTermsInLine;
  while (ptr != 0 && ptr[1] != '-') {
    // Determine line header.
    // NOTE DYN gets no header in DYNA> section.
    std::string headerLine(ptr+5, 8);
    ArgList header(headerLine, " :");
    LineHeaders.push_back( header[0] + ">" );
    ArgList line( ptr+headerOffset );
    for (int col = 0; col < line.Nargs(); col++) {
      if (line[col] == "Time") timeIdx = (int)Terms.size();
      Terms.push_back( line[col] );
    }
    // Next line
    ptr = buffer.Line();
    nTermsInLine.push_back( line.Nargs() );
  }
  mprintf("\t%zu terms:", Terms.size());
  DataSetList::DataListType inputSets;
  inputSets.reserve( Terms.size() );
  bool isRestart = false;
  for (Sarray::const_iterator it = Terms.begin(); it != Terms.end(); ++it) {
    mprintf(" %s", it->c_str());
    MetaData md(dsname, *it);
    DataSet* ds = dsl.CheckForSet( md );
    if (DYNA) {
      if (ds != 0 && !isRestart)
        isRestart = true;
      inputSets.push_back( new DataSet_double() );
      inputSets.back()->SetMeta( md );
    } else {
      if (ds == 0) {
        // New ENER set
        ds = dsl.AddSet( DataSet::DOUBLE, md );
        if (ds == 0) return 1;
      } else {
        mprintf("\tAppending to existing ENER set '%s'\n", ds->legend());
        if (ds->Type() != DataSet::DOUBLE) {
          mprinterr("Error: Set '%s' is not of type DOUBLE.\n", ds->legend());
          return 1;
        }
      }
      inputSets.push_back( ds );
    }
    // TODO make default width.prec 12.5?
  }
  mprintf("\n");
  if (isRestart)
    mprintf("\tWill append to existing set; skipping step 0.\n");
  if (debug_ > 0) {
    mprintf("DEBUG: Time index: %i\n", timeIdx);
    mprintf("DEBUG: [%s]\n", ptr);
    mprintf("DEBUG: %zu lines:\n", nTermsInLine.size());
    for (unsigned int i = 0; i != nTermsInLine.size(); i++)
      mprintf("DEBUG:\t\t[%u] '%s' %i terms.\n",
              i+1, LineHeaders[i].c_str(), nTermsInLine[i]);
  }
  if (DYNA && timeIdx == -1) {
    mprinterr("Error: Reading energy from DYNAmics but no Time field detected.\n");
    return 1;
  }
  // Read data.
  int set = 0;
  int lastStep = -1;
  bool readFile = true;
  bool isRepD = false;
  while (readFile) {
    // Scan to next DYNA> section
    while (ptr != 0) {
      if (DYNA) {
        if (ptr[0] == 'D' && ptr[1] == 'Y' && ptr[2] == 'N' && ptr[3] == 'A' && ptr[4] == '>' )
          break;
        if (ptr[0] == 'R' && ptr[1] == 'E' && ptr[2] == 'P' && ptr[3] == 'D' && ptr[4] == '>' )
          isRepD = true;
      } else if (MINI) {
        if (ptr[0] == 'M' && ptr[1] == 'I' && ptr[2] == 'N' && ptr[3] == 'I' && ptr[4] == '>' )
          break;
      } else { // ENER
        if (ptr[0] == 'E' && ptr[1] == 'N' && ptr[2] == 'E' && ptr[3] == 'R' && ptr[4] == '>' )
          break;
      }
      ptr = buffer.Line();
    }
    if (ptr == 0)
      readFile = false;
    else {
      int step;
      double dvals[5];
      bool ignoreStep = false;
      bool repeatedStep = false;
      std::fill(dvals, dvals+5, 0.0); // DEBUG
      //mprintf("%i [%s]\n", buffer.LineNumber(), ptr); // DEBUG
      int idx = 0; // Index into inputSets
      for (unsigned int i = 0; i != nTermsInLine.size(); i++) {
        // Determine if this line is present. First line should always
        // be present, DYN
        bool lineIsPresent = true;
        if (i > 0) {
          // Advance to first non-whitespace character
          const char* key = ptr + 5;
          while (*key != '\0' && isspace(*key)) ++key;
          lineIsPresent = (LineHeaders[i].compare(0, LineHeaders[i].size(),
                                                  key, LineHeaders[i].size()) == 0);
          //mprintf("\t%s lineIsPresent= %i {%s}\n", LineHeaders[i].c_str(), (int)lineIsPresent, key); // DEBUG
          if (!lineIsPresent)
            mprintf("Warning: Missing expected term %s at line %i\n", LineHeaders[i].c_str(), buffer.LineNumber());
        } else {
          // Line 0 should have the step
          sscanf(ptr+5, "%9i", &step);
          //mprintf("DEBUG: Step %i\n", step);
          if (step == lastStep) {
            // If REPD, dynamics are restarted after exchange, so the output
            // from the last step is repeated and should be ignored. Steps
            // can also be repeated after e.g. velocity scaling, so those
            // values should overwrite previous ones.
            if (isRepD)
              ignoreStep = true;
            else {
              repeatedStep = true;
              set--;
            }
          } else if (step == 0 && isRestart) {
            // First step of restarted run is pretty much same as last step of previous.
            ignoreStep = true;
          }
        }
        if (ignoreStep) {
          // Only advance if the line was actually present.
          if (lineIsPresent) ptr = buffer.Line();
        } else if (lineIsPresent) {
          int nvals = sscanf(ptr+headerOffset,"%13lf%13lf%13lf%13lf%13lf",
                             dvals, dvals+1, dvals+2, dvals+3, dvals+4);
          // SANITY CHECK
          if (nvals != nTermsInLine[i]) {
            mprinterr("Error: Number of terms in line %i (%i) != expected terms (%i)\n",
                      buffer.LineNumber(), nvals, nTermsInLine[i]);
            return 1; // TODO carry on?
          }
          //mprintf("DEBUG: %8s", LineHeaders[i].c_str());
          if (repeatedStep) {
            for (int val = 0; val < nvals; val++) {
              DataSet_double& dset = static_cast<DataSet_double&>( *(inputSets[idx+val]) );
              dset[set] = dvals[val];
            }
          } else {
            for (int val = 0; val < nvals; val++) {
              //mprintf(" %13.5f", dvals[val]); // DEBUG
              inputSets[idx+val]->Add(set, dvals + val);
            }
            //mprintf("\n"); // DEBUG
          }
          ptr = buffer.Line();
        }
        idx += nTermsInLine[i];
      } // END loop over DYNA> lines
      if (!ignoreStep) set++;
      lastStep = step;
    } // END start DYNA> section
  }
  if (isRepD)
    mprintf("\tREPD run detected. Repeated steps were ignored.\n");
  buffer.CloseFile();

  if (DYNA) {
    // Separate out time values.
    DataSetList::DataListType dataSets;
    dataSets.reserve( inputSets.size() - 1 );
    for (int idx = 0; idx != (int)inputSets.size(); idx++) {
      if (idx != timeIdx)
        dataSets.push_back( inputSets[idx] );
    }
    DataSetList::Darray const& timeVals = ((DataSet_double*)inputSets[timeIdx])->Data();
    if (dsl.AddOrAppendSets( "Time", timeVals, dataSets )) return 1;
    // Delete the time value data set
    delete inputSets[timeIdx];
  }
  
  return 0;
}

// DataIO_CharmmOutput::WriteHelp()
void DataIO_CharmmOutput::WriteHelp() { }

// DataIO_CharmmOutput::processWriteArgs()
int DataIO_CharmmOutput::processWriteArgs(ArgList& argIn)
{
  return 0;
}

// DataIO_CharmmOutput::WriteData()
int DataIO_CharmmOutput::WriteData(FileName const& fname, DataSetList const& dsl)
{
  mprinterr("Error: Charmm output write not supported.\n");
  return 1;
}
