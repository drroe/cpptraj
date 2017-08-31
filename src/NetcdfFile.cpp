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
  type_(NC_UNKNOWN),
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

// NetcdfFile::NC_setupRead()
int NetcdfFile::NC_setupRead(File::Name const& fnameIn) {
  if (IsOpen()) Close();
  type_ = GetNetcdfConventions(fnameIn);
  return (Setup(fnameIn, File::READ));
}

// NetcdfFile::NC_setupWrite()
int NetcdfFile::NC_setupWrite(File::Name const& fnameIn, NCTYPE type, int natomIn,
                              CoordinateInfo const& coordInfo, std::string const& title)
{
  type_ = type;
  ncatom_ = natomIn;
  cInfo_ = coordInfo;
  nc_title_ = title;
  return (Setup(fnameIn, File::WRITE));
}

int NetcdfFile::NC_open() { return Open(); }

// -----------------------------------------------------------------------------
// NetcdfFile::InternalOpen()
int NetcdfFile::InternalOpen() {
  if ( Access() == File::READ ) {
    if (NC::CheckErr( nc_open( Filename().full(), NC_NOWRITE, &ncid_ ) )) return 1;
  } else {
    if (NC::CheckErr( nc_open( Filename().full(), NC_WRITE,   &ncid_ ) )) return 1;
  }
  return 0;
}

// NetcdfFile::InternalClose()
void NetcdfFile::InternalClose() {
  if (ncid_ == -1) return;
  bool err = NC::CheckErr( nc_close(ncid_) );
  if (Debug() > 0 && !err)
    mprintf("Successfully closed ncid %i\n", ncid_);
  ncid_ = -1;
}

// NetcdfFile::InternalSetup()
int NetcdfFile::InternalSetup() {
  if (type_ == NC_UNKNOWN) {
    mprinterr("Internal Error: NetcdfFile::InternalSetup(): Type is unknown.\n");
    return 1;
  }
  int err;
  if (Access() == File::WRITE)
    err = CreateNewFile();
  else
    err = SetupExistingFile();
  if (Debug() > 0) {
    WriteVIDs();
    mprintf("DEBUG: %s\n", cInfo_.InfoString().c_str());
  }
  return err;
}

// -----------------------------------------------------------------------------
int NetcdfFile::CreateNewFile() { return 1; } // TODO

// -----------------------------------------------------------------------------

// NetcdfFile::CheckConventionsVersion()
void NetcdfFile::CheckConventionsVersion() {
  std::string attrText = NC::GetAttrText(ncid_, "ConventionVersion");
  if ( attrText != "1.0")
    mprintf("Warning: NetCDF file has ConventionVersion that is not 1.0 (%s)\n", attrText.c_str());
}

// NetcdfFile::SetupExistingFile()
int NetcdfFile::SetupExistingFile() {
  Open();
  // This will warn if conventions are not 1.0 
  CheckConventionsVersion();
  // Title
  nc_title_ = NC::GetAttrText(ncid_, "title");
  // Get frame dimension ID if present.
  frameDID_ = NC::GetDimInfo( ncid_, NCFRAME, ncframe_ );
  // Get atoms info
  atomDID_ = NC::GetDimInfo( ncid_, NCATOM, ncatom_ );
  if (atomDID_ != -1)
    ncatom3_ = ncatom_ * 3;
  // Get coord info
  coordVID_ = -1;
  if ( nc_inq_varid(ncid_, NCCOORDS, &coordVID_) == NC_NOERR ) {
    if (Debug() > 0) mprintf("\tNetCDF file has coordinates.\n");
    std::string attrText = NC::GetAttrText(ncid_, coordVID_, "units");
    if (attrText != "angstrom")
      mprintf("Warning: NetCDF file has length units of %s - expected angstrom.\n",
              attrText.c_str());
  }
  // Get spatial info
  int spatial;
  spatialDID_ = NC::GetDimInfo( ncid_, NCSPATIAL, spatial );
  if (spatialDID_ == -1) return 1;
  if (spatial != 3) {
     mprinterr("Error: Expected 3 spatial dimensions, got %i\n",spatial);
    return 1;
  }
  if ( NC::CheckErr(nc_inq_varid(ncid_, NCSPATIAL, &spatialVID_)) ) {
    mprintf("Warning: Could not get spatial VID. File may not be Amber NetCDF compliant.\n");
    mprintf("Warning: Assuming spatial variables are 'x', 'y', 'z'\n");
  } else {
    start_[0] = 0;
    count_[0] = 3;
    char xyz[3];
    if (NC::CheckErr(nc_get_vara_text(ncid_, spatialVID_, start_, count_, xyz))) {
      mprinterr("Error: Getting spatial variables.\n");
      return 1;
    }
    if (xyz[0] != 'x' || xyz[1] != 'y' || xyz[2] != 'z') {
      mprinterr("Error: NetCDF spatial variables are '%c', '%c', '%c', not 'x', 'y', 'z'\n",
                xyz[0], xyz[1], xyz[2]);
      return 1;
    }
  }
  // Get velocity info
  velocityVID_ = -1;
  if ( nc_inq_varid(ncid_, NCVELO, &velocityVID_) == NC_NOERR ) {
    if (Debug() > 0) mprintf("\tNetCDF file has velocities.\n");
  }
  // Get force info
  frcVID_ = -1;
  if ( nc_inq_varid(ncid_, NCFRC, &frcVID_) == NC_NOERR ) {
    if (Debug() > 0) mprintf("\tNetCDF file has forces.\n");
  }
  // Get overall replica and coordinate indices
  crdidxVID_ = -1;
  if ( nc_inq_varid(ncid_, NCREMD_REPIDX, &repidxVID_) == NC_NOERR ) {
      if (Debug() > 0) mprintf("\tNetCDF file has overall replica indices.\n");
    if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_CRDIDX, &crdidxVID_)) ) {
      mprinterr("Error: Getting overall coordinate index variable ID.\n");
      return 1;
    }
  } else
    repidxVID_ = -1;
  // Determine if Netcdf file contains time; set up timeVID and check units.
  timeVID_=-1;
  if ( nc_inq_varid(ncid_, NCTIME, &timeVID_) == NC_NOERR ) {
    std::string attrText = NC::GetAttrText(ncid_, timeVID_, "units");
    if (attrText!="picosecond")
      mprintf("Warning: NetCDF file has time units of %s - expected picosecond.\n",
              attrText.c_str());
    // Check for time values which have NOT been filled, which was possible
    // with netcdf trajectories created by older versions of ptraj/cpptraj.
    if (ncframe_ > 0 && type_ == NC_AMBERTRAJ) {
      float time;
      start_[0] = 0; count_[0] = 1;
      if (NC::CheckErr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time))) {
        mprinterr("Error: Getting time value for NetCDF file.\n");
        return 1;
      }
      if (time == NC_FILL_FLOAT) {
        mprintf("Warning: NetCDF file time variable defined but empty. Disabling.\n");
        timeVID_ = -1;
      }
    }
  }
  // Setup box variable IDs
  Box nc_box;
  if ( nc_inq_varid(ncid_, NCCELL_LENGTHS, &cellLengthVID_) == NC_NOERR ) {
    if (NC::CheckErr( nc_inq_varid(ncid_, NCCELL_ANGLES, &cellAngleVID_) )) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (Debug() > 0) mprintf("\tNetCDF Box information found.\n");
    // If present, get box lengths/angles.
    start_[0]=0; 
    start_[1]=0; 
    start_[2]=0;
    start_[3]=0;
    switch (type_) {
      case NC_AMBERRESTART:
        count_[0]=3;
        count_[1]=0;
        count_[2]=0;
        break;
      case NC_AMBERTRAJ:
        count_[0]=1; 
        count_[1]=3;
        count_[2]=0;
        break;
      case NC_AMBERENSEMBLE:
        count_[0]=1; // NOTE: All ensemble members must have same box type
        count_[1]=1; // TODO: Check all members?
        count_[2]=3;
        break;
      case NC_UNKNOWN: return 1; // Sanity check
    }
    count_[3]=0;
    double boxCrd[6]; /// XYZ ABG
    if ( NC::CheckErr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, boxCrd )) )
    {
      mprinterr("Error: Getting cell lengths.\n");
      return 1;
    }
    if ( NC::CheckErr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, boxCrd+3)) )
    {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (Debug() > 0) mprintf("\tNetCDF Box: XYZ={%f %f %f} ABG={%f %f %f}\n",
                              boxCrd[0], boxCrd[1], boxCrd[2], boxCrd[3], boxCrd[4], boxCrd[5]);
    nc_box.SetBox( boxCrd );
  }
  // Determine if Netcdf file contains temperature; set up temperature VID.
  TempVID_=-1;
  if ( nc_inq_varid(ncid_,NCTEMPERATURE,&TempVID_) == NC_NOERR ) {
    if (Debug()>0) mprintf("\tNetCDF file has replica temperatures.\n");
  } 
  // Determine if Netcdf file contains multi-D REMD info.
  int dimensionDID = -1;
  ReplicaDimArray remdDim;
  if ( nc_inq_dimid(ncid_, NCREMD_DIMENSION, &dimensionDID) == NC_NOERR) {
    dimensionDID = NC::GetDimInfo(ncid_, NCREMD_DIMENSION, remd_dimension_);
    if (dimensionDID = -1) return 1; // SANITY CHECK
    if (remd_dimension_ < 1) {
      mprinterr("Error: Number of REMD dimensions is less than 1!\n");
      return 1;
    }
    // Start and count for groupnum and dimtype, allocate mem
    start_[0]=0; 
    start_[1]=0; 
    start_[2]=0;
    count_[0]=remd_dimension_; 
    count_[1]=0; 
    count_[2]=0;
    int* remd_dimtype = new int[ remd_dimension_ ];
    // Get dimension types
    int dimtypeVID;
    if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_DIMTYPE, &dimtypeVID)) ) {
      mprinterr("Error: Getting dimension type variable ID for each dimension.\n");
      return 1;
    }
    if ( NC::CheckErr(nc_get_vara_int(ncid_, dimtypeVID, start_, count_, remd_dimtype)) ) {
      mprinterr("Error: Getting dimension type in each dimension.\n");
      return -1;
    }
    // Get VID for replica indices
    if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_INDICES, &indicesVID_)) ) {
      mprinterr("Error: Getting replica indices variable ID.\n");
      return 1;
    }
    // Print info for each dimension
    for (int dim = 0; dim < remd_dimension_; ++dim)
      remdDim.AddRemdDimension( remd_dimtype[dim] );
    delete[] remd_dimtype;
  }
  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;
  if (Debug() > 1) NC::Debug(ncid_);

  cInfo_ = CoordinateInfo(remdDim, nc_box, velocityVID_ != -1,
                          TempVID_ != -1, timeVID_ != -1, frcVID_ != -1);
  cInfo_.SetCrd( coordVID_ != -1 );

  Close();
  return 0;
}
// -----------------------------------------------------------------------------
// NetcdfFile::WriteIndices()
void NetcdfFile::WriteIndices() const {
  mprintf("DBG: Start={%zu, %zu, %zu, %zu} Count={%zu, %zu, %zu, %zu}\n",
         start_[0], start_[1], start_[2], start_[3],
         count_[0], count_[1], count_[2], count_[3]);
}

// NetcdfFile::WriteVIDs()
void NetcdfFile::WriteVIDs() const {
  rprintf("TempVID_=%i  coordVID_=%i  velocityVID_=%i frcVID_=%i  cellAngleVID_=%i"
          "  cellLengthVID_=%i  indicesVID_=%i\n",
          TempVID_, coordVID_, velocityVID_, frcVID_, cellAngleVID_, cellLengthVID_, indicesVID_);
}


 
#endif /* BINTRAJ */
