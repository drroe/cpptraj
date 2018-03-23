#include <sys/stat.h> // stat
#include <cstdio>     // FILE, fopen
#include "File_Base.h"
#include "CpptrajStdio.h"

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
