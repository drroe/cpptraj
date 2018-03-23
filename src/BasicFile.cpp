#include "BasicFile.h"
#include "CpptrajStdio.h"
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "FileIO_Mpi.h"
#  include "FileIO_MpiShared.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

using namespace File;

const char* BasicFile::FileTypeName_[] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE", "MPISHARED"
};

BasicFile::BasicFile() :
  IO_(0),
  isDos_(0),
  uncompressed_size_(0U),
  fileType_(UNKNOWN_TYPE)
{}

BasicFile::BasicFile(int d) : Base(d),
  IO_(0),
  isDos_(0),
  uncompressed_size_(0U),
  fileType_(UNKNOWN_TYPE)
{}

BasicFile::~BasicFile() {
  Reset();
}

/** Copy constructor. Always copy file closed. Allocate IO if necessary. */
BasicFile::BasicFile(BasicFile const& rhs) :
  Base(rhs),
  IO_(0),
  isDos_(rhs.isDos_),
  uncompressed_size_(rhs.uncompressed_size_)
{
  if (rhs.IO_ != 0)
    SetupFileIO( rhs.fileType_ );
}

/** Assignment. Always assign closed. */
BasicFile& BasicFile::operator=(BasicFile const& rhs) {
  if (this != &rhs) {
    Reset();
    Base::operator=(rhs);
    isDos_ = rhs.isDos_;
    uncompressed_size_ = rhs.uncompressed_size_;
    if (rhs.IO_ != 0)
      SetupFileIO( rhs.fileType_ );
  }
  return *this;
}
    
// BasicFile::UncompressedSize()
unsigned int BasicFile::UncompressedSize() const {
  if (Compression() == NO_COMPRESSION)
    return Size();
  else
    return uncompressed_size_;
}

/** Set up IO object, determine DOS line endings, etc.
  * \return first line size, or -1 if setup fails.
  */
int BasicFile::BasicSetup() {
  Reset();
  if (Debug() > 0)
    mprintf("BasicFile: Setting up %s for %s.\n", Filename().full(), accessStr());
  if (Compression() != NO_COMPRESSION && Access() == APPEND) {
    mprinterr("Error: Appending to compressed files is not supported.\n");
    return -1;
  }

  uncompressed_size_ = Size();
  // FIXME : determine line endings in streams?
  unsigned int lineSize = 0;
  if (IsStream()) {
    // file type must be STANDARD for streams
    if (SetupFileIO( STANDARD )) return -1;
  } else {
    if (SetupFileIO( UNKNOWN_TYPE )) return -1;
    // Check if file exists.
    if (IsPresent()) {
      // File exists.
      uncompressed_size_ = IO_->Size( Filename().full() );
      // Check for binary file, DOS line endings, first line size.
      if (IO_->Open( Filename().full(), "rb" ) != 0) return -1;
      char bufchar;
      while ( IO_->Read(&bufchar, 1) == 1 ) {
        if ( bufchar < 7 ) {
          // ASCII code less than 7 means file probably is binary. Do not
          // try to figure anything else out.
          mprintf("\t'%s' appears to be data.\n", Filename().full());
          break;
        }
        ++lineSize;
        if ( bufchar == '\n' ) break;
        if ( bufchar == '\r' ) {
          isDos_ = 1;
          if ( IO_->Read(&bufchar, 1) == 1 && bufchar == '\n' )
            ++lineSize;
          break;
        }
      }
      IO_->Close();
    } else {
      // File does not exist. If READ, fail silently.
      if (Access() == READ) return -1;
    }
  }
  
  if (Debug() > 0) {
    rprintf("\t[%s] is type %s with access %s\n", Filename().full(), FileTypeName_[fileType_],
            accessStr());
    rprintf("\t  isDos= %i\n\t  uncompressed_size_ = %u\n\t  firstLineSize = %u\n",
            isDos_, UncompressedSize(), lineSize);
  }
  return (int)lineSize;
}

#ifdef MPI
/** Open the file using MPI file routines. */
int BasicFile::ParallelOpenFile(AccessType accessIn, Parallel::Comm const& commIn, bool sharedWrite)
{
  if (IO_ == 0) {
    mprinterr("Internal Error: CpptrajFile has not been set up.\n");
    return 1;
  }
  // This will currently only work for fileType_ STANDARD
  if (fileType_ != STANDARD) {
    mprinterr("Error: Parallel file access not supported for file type '%s'\n",
              FileTypeName_[fileType_]);
    return 1;
  }
  // This will NOT work for streams.
  if (isStream_) {
    mprinterr("Error: Parallel file access not supported for streams.\n");
    return 1;
  }
  if (isOpen_) CloseFile();
  // TODO Save serial IO object?
  if (sharedWrite)
    fileType_ = MPISHARED;
  else
    fileType_ = MPIFILE;
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;
  ((FileIO_Mpi*)IO_)->SetComm( commIn );
  return OpenFile( accessIn );
}
#endif

// BasicFile::OpenIO()
int BasicFile::OpenIO() {
  if (IO_ == 0) {
    mprinterr("Internal Error: BasicFile has not been set up.\n");
    return 1;
  }
  int err = 0;
  if (IsOpen()) Close();
  if (IsStream()) {
    switch (Access()) {
      case READ : err = IO_->OpenStream( FileIO::STDIN ); break;
      case WRITE: err = IO_->OpenStream( FileIO::STDOUT); break;
      default:
        mprinterr("Internal Error: %s access not supported for file streams.\n",
                  accessStr());
        err = 1;
    }
    if (Debug() > 0 && err == 0)
      rprintf("Opened stream %s\n", Filename().full());
  } else {
    if (Filename().empty()) {
      mprinterr("Internal Error: CpptrajFile file name is empty.\n");
      err = 1;
    } else {
      switch (Access()) {
        case READ:   err = IO_->Open(Filename().full(), "rb"); break;
        case WRITE:  err = IO_->Open(Filename().full(), "wb"); break;
        case APPEND: err = IO_->Open(Filename().full(), "ab"); break;
        case UPDATE: err = IO_->Open(Filename().full(), "r+b"); break;
      }
      if (Debug() > 0 && err == 0)
        rprintf("Opened file %s with access %s\n", Filename().full(), accessStr());
    }
  }
  if (err != 0) {
    if (Debug() > 0)
      rprinterr("Could not open %s with access %s\n", Filename().full(), accessStr());
    mprinterr("Error: File '%s': %s\n", Filename().full(), StrError());
  }
  return err;
}

// BasicFile::InternalClose()
void BasicFile::InternalClose() {
  if (IsOpen()) {
    IO_->Close();
    #   ifdef MPI
    // Restore standard IO object.
    if (IsMPI()) {
      delete IO_;
      fileType_ = STANDARD;
      IO_ = SetupFileIO( fileType_ );
      if (IO_ == 0)
        mprinterr("Internal Error: Could not reset file '%s' from parallel to serial.\n",
                  fname_.full());
    }
#   endif
  }
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
