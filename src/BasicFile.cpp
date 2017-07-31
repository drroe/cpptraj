#include "BasicFile.h"
#include "CpptrajStdio.h"
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "FileIO_Mpi.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

using namespace File;

BasicFile::BasicFile() :
  IO_(0),
  isDos_(0),
  uncompressed_size_(0),
  fileType_(UNKNOWN_TYPE)
{}

const char* BasicFile::FileTypeName_[] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};


int BasicFile::InternalSetup() {
  Reset();
  if (Debug() > 0)
    mprintf("BasicFile: Setting up %s for %s.\n", Filename().full(), accessStr());
  int err = 0;
  switch ( Access() ) {
    case UPDATE:
    case READ:   err = SetupRead(); break;
    case WRITE:  err = SetupWrite(); break;
    case APPEND: err = SetupAppend(); break;
  }
  return err;
}

int BasicFile::SetupRead() {
  // If nameIn is empty assume reading from STDIN desired. 
  if (IsStream()) {
    // file type must be STANDARD for streams
    fileType_ = STANDARD;
    if (SetupFileIO( STANDARD )) return 1;
  } else {
    // Check if file exists. If not, fail silently
    // FIXME this check also happens in Base::Setup() - consolidate
    if (!File::Exists( Filename() )) return 1;
    if (SetupFileIO( UNKNOWN_TYPE )) return 1;
  }
  // Check for DOS line endings 

  if (Debug() > 0)
    rprintf("\t[%s] is type %s with access READ\n", Filename().full(), FileTypeName_[fileType_]);
  return 0;
}

/** Set up the IO based on given file type. */
int BasicFile::SetupFileIO( FileType typeIn ) {
  if (IO_ != 0) delete IO_;
  IO_ = 0;
  fileType_ = typeIn;
  if (fileType_ == UNKNOWN_TYPE) {
    if ( Compression() == GZIP )
      fileType_ = GZIPFILE;
    else if ( Compression() == BZIP2 )
      fileType_ = BZIP2FILE;
    else
      fileType_ = STANDARD;
  }
  switch (fileType_) {
    case STANDARD  :
      IO_ = new FileIO_Std();
      break;
    case GZIPFILE  : 
#     ifdef HASGZ
      IO_ = new FileIO_Gzip(); 
      break;
#     else
      mprinterr("Error: Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#     endif
    case BZIP2FILE :
#     ifdef HASBZ2 
      IO_ = new FileIO_Bzip2();
      break;
#     else
      mprinterr("Error: Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#     endif
    case MPIFILE   : 
#     ifdef MPI
      IO_ = new FileIO_Mpi();
      break;
#     else
      mprinterr("Error: Compiled without MPI support. Recompile with -DMPI\n");
      return 0;
#     endif
    //case ZIPFILE   : return (new ZipFile()); break;
    default        : 
      mprinterr("Error: Unrecognized file type.\n");
      return 1;
  }
  return 0;
}

/** Close file if open, reset all file information. */
void BasicFile::Reset() {
  Close();
  if (IO_!=0) delete IO_;
  IO_ = 0;
  isDos_ = 0;
  uncompressed_size_ = 0;
  fileType_ = UNKNOWN_TYPE;
}

