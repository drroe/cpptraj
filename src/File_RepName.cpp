#include "File_RepName.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

using namespace File;

// RepName CONSTRUCTOR
File::RepName::RepName(Name const& fname, int debugIn) :
  extChar_('.')
{
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: FileName=[%s]\n", fname.full());
  if ( fname.Ext().empty() ) {
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n"
              "Error:   detection of replica trajectories. Expected filename format is\n"
              "Error:   <Prefix>.<#> (with optional compression extension), examples:\n"
              "Error:   Rep.traj.nc.000,  remd.x.01.gz etc.\n", fname.base());
    return;
  }
  // Split off everything before replica extension
  size_t found = fname.Full().rfind( fname.Ext() );
  Prefix_.assign( fname.Full().substr(0, found) );
  ReplicaExt_.assign( fname.Ext() ); // This should be the numeric extension
  // Remove leading '.'
  if (ReplicaExt_[0] == '.') ReplicaExt_.erase(0,1);
  CompressExt_.assign( fname.Compress() );
  if (debugIn > 1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix_.c_str(), ReplicaExt_.c_str(), CompressExt_.c_str());
  }
  // CHARMM replica numbers are format <name>_<num>
  if ( !validInteger(ReplicaExt_) ) {
    size_t uscore = fname.Full().rfind('_');
    if (uscore != std::string::npos) {
      Prefix_.assign( fname.Full().substr(0, uscore) );
      ReplicaExt_.assign( fname.Full().substr(uscore+1) );
      extChar_ = '_';
      if (debugIn > 0)
        mprintf("\tREMDTRAJ: CHARMM style replica names detected, prefix='%s' ext='%s'\n",
                Prefix_.c_str(), ReplicaExt_.c_str());
    }
  }
  // Check that the numerical extension is valid.
  if ( !validInteger(ReplicaExt_) ) {
    mprinterr("Error: Replica extension [%s] is not an integer.\n", ReplicaExt_.c_str());
    Prefix_.clear(); // Empty Prefix_ indicates error.
    return;
  }
  ExtWidth_ = (int)ReplicaExt_.size();
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: Numerical Extension width=%i\n", ExtWidth_);
  // Store lowest replica number
  lowestRepnum_ = convertToInteger( ReplicaExt_ );
  // TODO: Do not allow negative replica numbers?
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n", lowestRepnum_);
}

/** \return Replica file name for given offset from lowest replica number. */
File::Name File::RepName::RepFilename(int offset) const {
  Name trajFilename;
  trajFilename.SetName_NoExpansion( Prefix_ + extChar_ +
                                        integerToString(lowestRepnum_ + offset, ExtWidth_) +
                                        CompressExt_ );
  return trajFilename;
}

