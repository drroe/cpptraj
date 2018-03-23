#include <cstdio>     // FILE, fopen
#include <cerrno>     // fileErrMsg, errno
#include <cstring>    // fileErrMsg, strerror
#include "File.h"
#include "File_RepName.h"
#include "File_WordExp.h"
#include "CpptrajStdio.h"

File::NameArray File::ExpandToFilenames(std::string const& fnameArg) {
  NameArray fnames;
  Sarray names = WordExp( fnameArg );
  for (Sarray::const_iterator it = names.begin(); it != names.end(); ++it)
  {
    Name fn;
    if (fn.SetName_NoExpansion( *it ))
      mprinterr("Internal Error: Could not set file name '%s'\n", it->c_str());
    else
      fnames.push_back( fn );
  }
  return fnames;
}

static std::string fileErrMsg_ = std::string("");

void File::ErrorMsg(const char* fname) {
  mprinterr("Error: '%s': %s\n", fname, fileErrMsg_.c_str());
}

const char* File::StrError() {
  return strerror(errno);
}

bool File::Exists(Name const& fn) {
  if (!fn.empty()) {
    FILE* infile = fopen(fn.full(), "rb");
    if (infile==0) {
      fileErrMsg_.assign( strerror(errno) );
      return false;
    }
    fclose(infile);
    return true;
  }
  return false;
}

//bool File::Exists(std::string const& fname) {
//  return File::Exists( Name(fname) );
//}

/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
File::NameArray File::SearchForReplicas(Name const& fname, int debug) {
  NameArray replica_filenames;
  if (!File::Exists(fname)) {
    mprinterr("Error: '%s' does not correspond to a file.\n", fname.full());
    return replica_filenames;
  }
  RepName repName(fname, debug);
  if (repName.Error()) return replica_filenames;
  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  if (File::Exists( repName.RepFilename( -1 ) )) {
    mprintf("Warning: Replica# found lower than file specified with trajin.\n"
            "Warning:   Found \"%s\"; 'trajin remdtraj' requires lowest # replica.\n",
            repName.RepFilename( -1 ).full());
  }
  // Add lowest replica filename, search for and add all replicas higher than it.
  replica_filenames.push_back( fname );
  int rep_offset = 0;
  bool search_for_files = true;
  Name trajFilename;
  while (search_for_files) {
    ++rep_offset;
    trajFilename = repName.RepFilename( rep_offset );
    //mprintf("\t\tChecking for %s\n", trajFilename.full());
    if (File::Exists( trajFilename ))
      replica_filenames.push_back( trajFilename );
    else
      search_for_files = false;
  }
  return replica_filenames;
}

#ifdef MPI
/** Each rank searches for replica based on lowest replica number. */
File::NameArray File::SearchForReplicas(Name const& fname, bool trajCommMaster,
                                        int ensRank, int ensSize, int debug)
{
  NameArray replica_filenames;
  RepName repName(fname, debug);
  if (repName.Error()) return replica_filenames;
  // TODO check for lower replica number?
  Name replicaFilename = repName.RepFilename( ensRank );
  // Only traj comm masters actually check for files.
  if (trajCommMaster) {
    if (!File::Exists( replicaFilename )) {
      File::ErrorMsg( replicaFilename.full() );
      rprinterr("Error: File '%s' not accessible.\n", replicaFilename.full());
      return replica_filenames;
    }
  }
  // At this point each rank has found its replica. Populate filename array.
  for (int offset = 0; offset < ensSize; ++offset)
    replica_filenames.push_back( repName.RepFilename( offset ) );
  return replica_filenames;
}
#endif
