#include "DataIO_AmberFF.h"
#include "CpptrajStdio.h"
#include "AmberParamFile.h"
#include "DataSet_Parameters.h"

/// CONSTRUCTOR
DataIO_AmberFF::DataIO_AmberFF()
{
  SetValid( DataSet::PARAMETERS );
  // TODO Topology too?
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
  mprintf("\tnbset <nonbond set name> : Nonbonded set name to use when multiple nonbond parameter sets.\n");
}

// DataIO_AmberFF::processReadArgs()
int DataIO_AmberFF::processReadArgs(ArgList& argIn)
{
  nbsetname_ = argIn.GetStringKey("nbset");
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

  AmberParamFile infile;
  int err = infile.ReadParams(prm, fname, nbsetname_, debug_);
  if (err != 0) {
    mprinterr("Error: Could not read '%s'\n", fname.full());
  }

  return err;
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
  AmberParamFile outfile;
  std::vector<DataSet_Parameters*> toWrite;
  for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
  {
    if ( (*it)->Type() == DataSet::PARAMETERS )
      toWrite.push_back( (DataSet_Parameters*)(*it) );
  }
  if (toWrite.empty()) {
    mprinterr("Error: No parameter sets to write.\n");
    return 1;
  } else if (toWrite.size() == 1) {
    return outfile.WriteParams( *(toWrite.front()), fname, debug_ );
  } else {
    // Create a combined parameter set
    ParameterSet prm;
    for (std::vector<DataSet_Parameters*>::const_iterator it = toWrite.begin();
                                                          it != toWrite.end(); ++it)
    {
      ParameterSet::UpdateCount UC;
      prm.UpdateParamSet( *(*it), UC, debug_, debug_ ); // FIXME verbose
    }
    return outfile.WriteParams( prm, fname, debug_ );
  }

  return 1;
}
