#include "DataIO_Coords.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"
#include "TrajectoryFile.h"

/// CONSTRUCTOR
DataIO_Coords::DataIO_Coords()
{

}

// DataIO_Coords::ID_DataFormat()
bool DataIO_Coords::ID_DataFormat(CpptrajFile& infile)
{
  // Needs to be either a topology format or a coords format
  ParmFile::ParmFormatType parm_format = ParmFile::DetectFormat(infile.Filename());
  TrajectoryFile::TrajFormatType traj_format = TrajectoryFile::DetectFormat(infile.Filename());
  if (parm_format != ParmFile::UNKNOWN_PARM ||
      traj_format != TrajectoryFile::UNKNOWN_TRAJ)
    return true;
  return false;
}

// DataIO_Coords::ReadHelp()
void DataIO_Coords::ReadHelp()
{

}

// DataIO_Coords::processReadArgs()
int DataIO_Coords::processReadArgs(ArgList& argIn)
{

  return 0;
}

/// \return True if this is a COORDS set we can append to
static inline bool can_append(DataSet::DataType typeIn) {
  return (typeIn == DataSet::COORDS ||
          typeIn == DataSet::FRAMES);
}

// DataIO_Coords::ReadData()
int DataIO_Coords::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  DataSet* dset = 0;
  if (!dsname.empty()) {
    // Is this set already present?
    dset = dsl.GetDataSet( dsname );
    if (dset != 0) {
      if (!can_append(dset->Type())) {
        mprinterr("Error: Cannot append coordinates to existing set '%s'\n", dset->legend());
        return 1;
      } else
        mprintf("\tAppending to set '%s'\n", dset->legend());
    }
  }

  return 0;
}

// DataIO_Coords::WriteHelp()
void DataIO_Coords::WriteHelp()
{

}

// DataIO_Coords::processWriteArgs()
int DataIO_Coords::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Coords::WriteData()
int DataIO_Coords::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
