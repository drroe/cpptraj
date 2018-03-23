#ifndef _WIN32
#   include <wordexp.h>
#endif
#include <sys/stat.h> // stat
#include <cstdio>     // FILE, fopen
#include <cerrno>     // fileErrMsg, errno
#include <cstring>    // fileErrMsg, strerror
#include "File.h"
#include "File_RepName.h"
#include "CpptrajStdio.h"

#ifndef _WIN32
/** Print error messages from the wordexp() function. */ // TODO do not duplicate
static void WexpErr(int err) {
  switch ( err ) {
    case WRDE_BADCHAR :
      mprinterr("Error: Illegal occurrence of newline or one of |, &, ;, <, >, (, ), {, }.\n");
      break;
    //case WRDE_BADVAL
    case WRDE_CMDSUB :
      mprinterr("Error: Command substitution is not allowed in file names.\n");
      break;
    case WRDE_NOSPACE :
      mprinterr("Error: Out of memory.\n");
      break;
    case WRDE_SYNTAX :
      mprinterr("Error: Bad syntax (unbalanced parentheses, unmatched quotes.\n");
      break;
  }
}
#endif /* _WIN32 */

// ===== File::Base ============================================================
File::Base::Base() :
  file_size_(0),
  debug_(0),
  access_(READ),
  compressType_(NO_COMPRESSION),
  isOpen_(false),
  isStream_(false),
  isPresent_(false)
{}

File::Base::Base(int d) :
  file_size_(0),
  debug_(d),
  access_(READ),
  compressType_(NO_COMPRESSION),
  isOpen_(false),
  isStream_(false),
  isPresent_(false)
{}

/** Copy constructor. Always copy closed. */
File::Base::Base(Base const& rhs) :
  file_size_(rhs.file_size_),
  debug_(rhs.debug_),
  access_(rhs.access_),
  compressType_(rhs.compressType_),
  isOpen_(false),
  isStream_(rhs.isStream_),
  isPresent_(rhs.isPresent_)
{}

/** Assignment. Always assign closed. */
File::Base& File::Base::operator=(Base const& rhs) {
  if (this != &rhs) {
    file_size_ = rhs.file_size_;
    debug_ = rhs.debug_;
    access_ = rhs.access_;
    compressType_ = rhs.compressType_;
    isOpen_ = false;
    isStream_ = rhs.isStream_;
    isPresent_ = rhs.isPresent_;
  }
  return *this;
}

const char* File::Base::AccessTypeName_[] = {
    "read", "write", "append", "update"
};

static inline File::CompressType SetCompressTypeFromName(File::Name const& fname) {
  if (fname.Compress() == ".gz")
    return File::GZIP;
  else if (fname.Compress() == ".bz2")
    return File::BZIP2;
  else
    return File::NO_COMPRESSION;
}

// File::Base::Setup()
/** Set the file name, determine if file is or should be compressed, then
  * call InternalSetup().
  */
int File::Base::Setup(Name const& fnameIn, AccessType accessIn)
{
  if (isOpen_) Close();
  access_ = accessIn;
  isOpen_ = false;
  file_size_ = 0;
  compressType_ = NO_COMPRESSION;
  isPresent_ = false;
  if (fnameIn.empty()) {
    // Empty file name is stream
    isStream_ = true;
    if (access_ == READ)
      fname_.SetName_NoExpansion("STDIN");
    else
      fname_.SetName_NoExpansion("STDOUT");
  } else {
    isStream_ = false;
    fname_ = fnameIn;
    isPresent_ = Exists(fname_);
    if (isPresent_) {
      // File exists. Get stats.
      struct stat frame_stat;
      if (stat(fname_.full(), &frame_stat) == -1) {
        mprinterr("Error: Could not find file status for %s\n", fname_.full());
        if (debug_>0) 
          perror("     Error from stat: ");
        return 1;
      }
      file_size_ = (unsigned int)frame_stat.st_size;
      // Access-specific setup for existing file.
      if (access_ == WRITE) {
        // TODO check for overwrite?
        compressType_ = SetCompressTypeFromName(fname_);
      } else {
        // ID compression by magic number (first 3 bytes).
        FILE* fIn = fopen(fname_.full(), "rb");
        if ( fIn == 0 ) { 
          mprinterr("Error: Could not open %s for hex signature read.\n", fname_.full());
          return 1;
        }
        unsigned char magic[3];
        magic[0] = 0; 
        magic[1] = 0; 
        magic[2] = 0;
        size_t numread = fread(magic, 1, 3, fIn);
        fclose(fIn);
        if (numread == 0) {
          if (access_ == READ)
            mprinterr("Error: File %s is empty\n", fname_.full());
        } else if (numread < 3 ) {
          mprintf("Warning: Could only read first %zu bytes of file %s.\n", numread, fname_.full());
        } else {
          if (debug_>0) mprintf("\tHex sig: %x %x %x\n", magic[0],magic[1],magic[2]);
          // Check compression
          if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8))
            compressType_ = GZIP;
          else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68))
            compressType_ = BZIP2;
          else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3))
            compressType_ = ZIP;
        }
      }
    } else {
      // File does not exist.
      if (access_ == READ) {
        mprinterr("Error: File '%s' does not exist.\n", fname_.full());
        ErrorMsg(fname_.full());
      } else
        compressType_ = SetCompressTypeFromName(fname_);
    }
  } // END file is not stream
  if (debug_ > 0) {
    mprintf("\tFILE INFO:");
    if (isStream_)
      mprintf(" STREAM\n");
    else
      mprintf(" %s\n", fname_.full());
    mprintf("\t  Size= %u\n", file_size_);
    const char* compTypeStr[4] = {"None", "Gzip", "Bzip", "Zip "};
    mprintf("\t  Compression= %s\n", compTypeStr[compressType_]);
    mprintf("\t  IsPresent= %i\n", (int)isPresent_);
  }
  return InternalSetup();
}

// File::Base::Open()
int File::Base::Open() {
  if (isOpen_) Close();
  if (InternalOpen()) return 1;
  isOpen_ = true;
  return 0;
}

// File::Base::Open()
int File::Base::Open(Name const& fnameIn, AccessType accessIn) {
  if (Setup(fnameIn, accessIn)) return 1;
  if (InternalOpen()) return 1;
  isOpen_ = true;
  return 0;
}

// File::Base::Close()
void File::Base::Close() {
  InternalClose();
  isOpen_ = false;
}

// =============================================================================
File::NameArray File::ExpandToFilenames(std::string const& fnameArg) {
  NameArray fnames;
#ifdef _WIN32
  fnames.push_back( fnameArg );
#else
  if (fnameArg.empty()) return fnames;
  wordexp_t expanded;
  int err = wordexp( fnameArg.c_str(), &expanded, WRDE_NOCMD );
  WexpErr( err );
  if ( err == 0 ) {
    for (unsigned int i = 0; i != expanded.we_wordc; i++) {
      if (expanded.we_wordv[i] == 0)
        mprinterr("Internal Error: Bad expansion at %i\n", i);
      else {
        Name fn;
        fn.SetName_NoExpansion( expanded.we_wordv[i] );
        fnames.push_back( fn );
      }
    }
    wordfree( &expanded );
  }
#endif /* _WIN32 */
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
