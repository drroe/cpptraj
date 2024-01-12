#include "DataIO_LeapRC.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataIO_AmberFF.h"
#include "DataIO_AmberFrcmod.h"
#include "DataIO_AmberLib.h"
#include "DataIO_AmberPrep.h"
#include <cstdlib> //getenv

/// CONSTRUCTOR
DataIO_LeapRC::DataIO_LeapRC()
{

}

// DataIO_LeapRC::ID_DataFormat()
bool DataIO_LeapRC::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  bool isLeaprc = false;
  // Scan the first 5 non-blank non-comment lines
  int nLinesScanned = 0;
  while (nLinesScanned < 5 && !isLeaprc) {
    const char* ptr = infile.NextLine();
    if (ptr == 0) break;
    if (ptr[0] == '\0' || ptr[0] == '#') continue;
    nLinesScanned++;
    ArgList line(ptr, " \t");
    if (line.Nargs() > 0) {
      if (line[0] == "logFile" || line[0] == "logfile")
        isLeaprc = true;
      else if (line[0] == "source")
        isLeaprc = true;
      else if (line[0] == "addAtomTypes" || line[0] == "addatomtypes")
        isLeaprc = true;
    }
  }
  infile.CloseFile();
  return isLeaprc;
}

// DataIO_LeapRC::ReadHelp()
void DataIO_LeapRC::ReadHelp()
{

}

// DataIO_LeapRC::processReadArgs()
int DataIO_LeapRC::processReadArgs(ArgList& argIn)
{

  return 0;
}

/** LEaP loadAmberParams command. */
int DataIO_LeapRC::LoadAmberParams(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  // TODO detect this better
  ArgList fargs( filename, "." );
  if (fargs.hasKey("frcmod")) {
    mprintf("\tLoading force field modifications from '%s'\n", filename.c_str());
    DataIO_AmberFrcmod infile;
    if (infile.ReadData(amberhome_ + "parm/" + filename, dsl, dsname)) {
      mprinterr("Error: Could not load force field modifications from '%s'\n", filename.c_str());
      return 1;
    }
  } else {
    mprintf("\tLoading force field from '%s'\n", filename.c_str());
    DataIO_AmberFF infile;
    if (infile.ReadData(amberhome_ + "parm/" + filename, dsl, dsname)) {
      mprinterr("Error: Could not load force field from '%s'\n", filename.c_str());
      return 1;
    }
  }
  return 0;
}

/** LEaP loadOff command. */
int DataIO_LeapRC::LoadOFF(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  DataIO_AmberLib infile;
  if (infile.ReadData(amberhome_ + "lib/" + filename, dsl, dsname)) {
    mprinterr("Error: Could not load library file '%s'\n", filename.c_str());
    return 1;
  }
  return 0;
}

/** LEaP loadAmberPrep command. */
int DataIO_LeapRC::LoadAmberPrep(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  DataIO_AmberPrep infile;
  if (infile.ReadData(amberhome_ + "lib/" + filename, dsl, dsname)) {
    mprinterr("Error: Could not load prep file '%s'\n", filename.c_str());
    return 1;
  }
  return 0;
}

/** LEaP addAtomTypes command. */
int DataIO_LeapRC::AddAtomTypes(NHarrayType& atomHybridizations, BufferedLine& infile)
const
{
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // Some entries (like LP and EP) are not required to have elements.
          // Set the hybridization index to 1 or 2.
          int hidx;
          if (aline.Nargs() == 3)
            hidx = 2;
          else if (aline.Nargs() == 2)
            hidx = 1;
          else {
            mprinterr("Error: Malformed addAtomTypes entry %s\n", tmp.c_str());
            return 1;
          }
          AtomType::HybridizationType ht;
          if (aline[hidx] == "sp3")
            ht = AtomType::SP3;
          else if (aline[hidx] == "sp2")
            ht = AtomType::SP2;
          else if (aline[hidx] == "sp")
            ht = AtomType::SP;
          else {
            mprintf("Warning: Unknown hybridization in addAtomTypes entry %s\n", tmp.c_str());
            ht = AtomType::UNKNOWN_HYBRIDIZATION;
          }
          atomHybridizations.push_back( NHpairType(aline[0], ht) );
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addAtomTypes command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
      //mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addAtomTypes command.\n");
    return 1;
  }
  mprintf("\tRead %zu atom hybridizations.\n", atomHybridizations.size());
  return 0;
}

// DataIO_LeapRC::ReadData()
int DataIO_LeapRC::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  // First, need to determine where the Amber FF files are
  const char* env = getenv("AMBERHOME");
  if (env != 0)
    amberhome_ = std::string(env) + "/dat/leap/";
  else {
    mprintf("Warning: AMBERHOME is not set. Determining FF file location based on leaprc file.\n");
    // Try to guess based on where the leaprc file is
    FileName leapcmddir( fname.DirPrefix_NoSlash() );
    if (leapcmddir.Base() == "oldff") {
      FileName leapcmddir2( leapcmddir.DirPrefix_NoSlash() );
      amberhome_ = leapcmddir2.DirPrefix();
    } else
      amberhome_ = leapcmddir.DirPrefix();
  }
  mprintf("\tForce field files located in '%s'\n", amberhome_.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) {
    mprinterr("Error: Could not open leaprc file '%s'\n", fname.full());
    return 1;
  }
  NHarrayType atomHybridizations;
  int err = 0;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (ptr[0] != '\0' && ptr[0] != '#') {
      ArgList line( ptr, " \t" );
      std::string param_fname;
      if (line.Contains("loadAmberParams"))
        err = LoadAmberParams( line.GetStringKey("loadAmberParams"), dsl, dsname );
      else if (line.Contains("loadamberparams"))
        err = LoadAmberParams( line.GetStringKey("loadamberparams"), dsl, dsname );
      else if (line.Contains("loadOff"))
        err = LoadOFF( line.GetStringKey("loadOff"), dsl, dsname );
      else if (line.Contains("loadoff"))
        err = LoadOFF( line.GetStringKey("loadoff"), dsl, dsname );
      else if (line.Contains("loadAmberPrep"))
        err = LoadAmberPrep( line.GetStringKey("loadAmberPrep"), dsl, dsname );
      else if (line.Contains("loadamberprep"))
        err = LoadAmberPrep( line.GetStringKey("loadamberprep"), dsl, dsname );
      else if (line.Contains("addAtomTypes") || line.Contains("addatomtypes"))
        err = AddAtomTypes(atomHybridizations, infile);
    }
    if (err != 0) break;
    ptr = infile.Line();
  }
  infile.CloseFile();
  return err;
}

// DataIO_LeapRC::WriteHelp()
void DataIO_LeapRC::WriteHelp()
{

}

// DataIO_LeapRC::processWriteArgs()
int DataIO_LeapRC::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_LeapRC::WriteData()
int DataIO_LeapRC::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
