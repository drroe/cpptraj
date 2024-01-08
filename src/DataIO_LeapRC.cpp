#include "DataIO_LeapRC.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataIO_AmberFF.h"
#include "DataIO_AmberFrcmod.h"
#include "DataIO_AmberLib.h"
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
    amberhome_ = leapcmddir.DirPrefix();
  }
  mprintf("\tForce field files located in '%s'\n", amberhome_.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) {
    mprinterr("Error: Could not open leaprc file '%s'\n", fname.full());
    return 1;
  }
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
