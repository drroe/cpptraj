#include "DataIO_AmberFF.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_Parameters.h"

/// CONSTRUCTOR
DataIO_AmberFF::DataIO_AmberFF()
{

}

// DataIO_AmberFF::ID_DataFormat()
bool DataIO_AmberFF::ID_DataFormat(CpptrajFile& infile)
{
  //if (infile.OpenFile()) return false;
  //std::string line = infile.GetLine(); // Title
  //infile.CloseFile();
  //bool isLib = (line == "!!index array str");
  return false;
}

// DataIO_AmberFF::ReadHelp()
void DataIO_AmberFF::ReadHelp()
{

}

// DataIO_AmberFF::processReadArgs()
int DataIO_AmberFF::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberFF::ReadData()
int DataIO_AmberFF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  // Allocate data set
  MetaData md( dsname );
  DataSet* ds = dsl.CheckForSet( md );
  if (ds != 0) {
    if (ds->Type() != DataSet::PARAMETERS) {
      mprinterr("Error: Set '%s' does not have parameters, cannot append.\n", ds->legend());
      return 1;
    }
    mprintf("\tAdding to existing set %s\n", ds->legend());
  } else {
    ds = dsl.AddSet( DataSet::PARAMETERS, md );
    if (ds == 0) return 1;
  }
  DataSet_Parameters& prm = static_cast<DataSet_Parameters&>( *ds ); 

  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  std::string title = infile.GetLine();
  mprintf("\tTitle: %s\n", title.c_str());

  infile.CloseFile();
  return 0;
}

// DataIO_AmberFF::WriteHelp()
void DataIO_AmberFF::WriteHelp()
{

}

// DataIO_AmberFF::processWriteArgs()
int DataIO_AmberFF::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberFF::WriteData()
int DataIO_AmberFF::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
