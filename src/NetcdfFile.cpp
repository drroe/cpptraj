#include "NetcdfFile.h"
#ifdef BINTRAJ
#  include <netcdf.h>
#  include "NC_Routines.h"
#endif
#include "CpptrajStdio.h"
#include "Constants.h"
#include "Version.h"

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(File::Name const& fname) {
  NCTYPE nctype = NC_UNKNOWN;
# ifdef BINTRAJ
  // NOTE: Do not use checkNCerr so this fails silently. Allows routine to
  //       be used in file autodetection.
  int tmp_ncid;
  if ( nc_open( fname.full(), NC_NOWRITE, &tmp_ncid ) != NC_NOERR )
    return NC_UNKNOWN;
  nctype = GetNetcdfConventions( tmp_ncid );
  nc_close( tmp_ncid );
# else
  mprintf("Error: Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
# endif
  return nctype;
}

#ifdef BINTRAJ
// DEFINES
#define NCENSEMBLE "ensemble"
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_LENGTHS "cell_lengths"
#define NCCELL_ANGULAR "cell_angular"
#define NCCELL_ANGLES "cell_angles"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCFRC "forces"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5
#define NCREMD_DIMENSION "remd_dimension"
#define NCREMD_DIMTYPE "remd_dimtype"
#define NCREMD_INDICES "remd_indices"
#define NCREMD_REPIDX "remd_repidx"
#define NCREMD_CRDIDX "remd_crdidx"
#define NCEPTOT "eptot"
#define NCBINS "bins"

// CONSTRUCTOR
NetcdfFile::NetcdfFile() :
  ncid_(-1),
  ncatom_(-1),
  ncatom3_(-1),
  ncframe_(-1),
  remd_dimension_(0),
  TempVID_(-1),
  coordVID_(-1),
  velocityVID_(-1),
  frcVID_(-1),
  spatialVID_(-1),
  cell_spatialVID_(-1),
  cell_angularVID_(-1),
  cellAngleVID_(-1),
  cellLengthVID_(-1),
  timeVID_(-1),
  indicesVID_(-1),
  repidxVID_(-1),
  crdidxVID_(-1),
  ensembleDID_(-1),
  frameDID_(-1),
  atomDID_(-1),
  spatialDID_(-1),
  labelDID_(-1),
  cell_spatialDID_(-1),
  cell_angularDID_(-1)
{
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  start_[3] = 0;
  count_[0] = 0;
  count_[1] = 0;
  count_[2] = 0;
  count_[3] = 0;
}

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(int ncidIn) {
  NCTYPE nctype = NC_UNKNOWN;
  std::string attrText = NC::GetAttrText(ncidIn, "Conventions");
  if (attrText == "AMBERENSEMBLE")
    nctype = NC_AMBERENSEMBLE;
  else if (attrText == "AMBER")
    nctype = NC_AMBERTRAJ;
  else if (attrText == "AMBERRESTART")
    nctype = NC_AMBERRESTART;
  else if (attrText.empty()) 
    mprinterr("Error: Could not get conventions from NetCDF file.\n");
  else {
    mprinterr("Error: NetCDF file has unrecognized conventions \"%s\".\n",
              attrText.c_str());
    mprinterr("Error:   Expected \"AMBER\", \"AMBERRESTART\", or \"AMBERENSEMBLE\".\n");
  }
  return nctype;
}
#endif /* BINTRAJ */
