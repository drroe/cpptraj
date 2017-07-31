#include <cstring>    // strlen
#include <cstdio>     // vsprintf
#include <cstdarg>    // va_X functions
#include <algorithm>  // std::max
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
  uncompressed_size_(0U),
  BUF_SIZE_(0U),
  fileType_(UNKNOWN_TYPE),
  linebuffer_(0)
{}

BasicFile::BasicFile(int d) : Base(d),
  IO_(0),
  isDos_(0),
  uncompressed_size_(0U),
  BUF_SIZE_(0U),
  fileType_(UNKNOWN_TYPE),
  linebuffer_(0)
{}

BasicFile::~BasicFile() {
  Reset();
}

const char* BasicFile::FileTypeName_[] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};

// -----------------------------------------------------------------------------
// BasicFile::Printf()
/** Take the formatted string and write it to file using Write.
  */
void BasicFile::Printf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(linebuffer_,format,args);
  IO_->Write(linebuffer_, strlen(linebuffer_));
  va_end(args);
}

std::string BasicFile::GetLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE_) != 0) {
    //mprinterr("Error: Getting line from %s\n", Filename().full());
    return std::string();
  }
  return std::string(linebuffer_);
}

const char* BasicFile::NextLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE_) != 0) {
    //mprinterr("Error: Reading line from %s\n", Filename().full());
    return 0;
  }
  return linebuffer_;
}

unsigned int BasicFile::UncompressedSize() const {
  if (Compression() == NO_COMPRESSION)
    return Size();
  else
    return uncompressed_size_;
}

// -----------------------------------------------------------------------------
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

int BasicFile::InternalOpen() {
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

void BasicFile::InternalClose() { if (IsOpen()) IO_->Close(); }

int BasicFile::SetupRead() {
  // FIXME : determine line endings in streams?
  unsigned int lineSize = 0;
  if (IsStream()) {
    // file type must be STANDARD for streams
    fileType_ = STANDARD;
    if (SetupFileIO( STANDARD )) return 1;
    uncompressed_size_ = Size();
  } else {
    // Check if file exists. If not, fail silently
    // FIXME this check also happens in Base::Setup() - consolidate
    if (!File::Exists( Filename() )) return 1;
    if (SetupFileIO( UNKNOWN_TYPE )) return 1;
    uncompressed_size_ = IO_->Size( Filename().full() );
    // Additional file characteristics
    if (IO_->Open( Filename().full(), "rb" ) != 0) return 1;
    char bufchar;
    while ( IO_->Read(&bufchar, 1) == 1 ) {
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
  }
  BUF_SIZE_ = std::max(1024U, lineSize + 1); // +1 for null char
  linebuffer_ = new char[ BUF_SIZE_ + 1 ];
  linebuffer_[BUF_SIZE_] = '\0';
  
  if (Debug() > 0) {
    rprintf("\t[%s] is type %s with access READ\n", Filename().full(), FileTypeName_[fileType_]);
    rprintf("\t  isDos= %i  BUF_SIZE_ = %u\n", isDos_, BUF_SIZE_);
  }
  return 0;
}

int BasicFile::SetupWrite() {
  return 1;
}

int BasicFile::SetupAppend() {
  return 1;
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
  if (linebuffer_ != 0) delete[] linebuffer_;
  linebuffer_ = 0;
  isDos_ = 0;
  uncompressed_size_ = 0;
  fileType_ = UNKNOWN_TYPE;
}

