#ifndef _WIN32
#   include <wordexp.h>
#endif
#include <sys/stat.h> // stat
#include <cstdio>     // FILE, fopen
#include <cerrno>     // fileErrMsg, errno
#include <cstring>    // fileErrMsg, strerror
#include "File.h"
#include "CpptrajStdio.h"

#ifndef _WIN32
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

// ===== File::Name ============================================================
// COPY CONSTRUCTOR
File::Name::Name( const Name& rhs ) : fullPathName_(rhs.fullPathName_),
  baseName_(rhs.baseName_), extension_(rhs.extension_),
  compressExt_(rhs.compressExt_), dirPrefix_(rhs.dirPrefix_) {}

// ASSIGNMENT
File::Name& File::Name::operator=(const Name& rhs) {
  if (this != &rhs) {
    fullPathName_ = rhs.fullPathName_;
    baseName_ = rhs.baseName_;
    extension_ = rhs.extension_;
    compressExt_ = rhs.compressExt_;
    dirPrefix_ = rhs.dirPrefix_;
  }
  return *this;
}

// Name::clear()
void File::Name::clear() {
  fullPathName_.clear();
  baseName_.clear();
  extension_.clear();
  compressExt_.clear();
  dirPrefix_.clear();
}

bool File::Name::MatchFullOrBase(std::string const& rhs) const {
  if (!fullPathName_.empty()) {
    // Prefer full filename match.
    if (fullPathName_ == rhs) return true;
    if (baseName_     == rhs) return true;
  }
  return false;
}

int File::Name::SetName(std::string const& nameIn) {
  // null filename allowed (indicates STDIN/STDOUT)
  if (nameIn.empty()) {
    clear();
    return 0;
  }
#ifndef _WIN32
  wordexp_t expanded;
  int err = wordexp( nameIn.c_str(), &expanded, WRDE_NOCMD );
  WexpErr( err );
  if (err == 0) {
    if (expanded.we_wordc < 1) { // Sanity check
      mprinterr("Internal Error: Word expansion failed.\n");
      err = 1;
    } else
      err = SetName_NoExpansion( expanded.we_wordv[0] );
    wordfree( &expanded );
  }
  return err;
#else
  SetName_NoExpansion(nameIn);
  return 0;
#endif
}

int File::Name::SetName_NoExpansion(std::string const& nameIn) {
  // null filename allowed (indicates STDIN/STDOUT)
  if (nameIn.empty()) {
    clear();
    return 0;
  }
  // Assign filename with full path
  fullPathName_.assign( nameIn );
  // Get position of last occurence of '/' to determine base filename
  size_t found = fullPathName_.find_last_of("/");
  if (found == std::string::npos) {
    baseName_ = fullPathName_;
    dirPrefix_.clear();
  } else {
    baseName_ = fullPathName_.substr(found+1);
    dirPrefix_ = fullPathName_.substr(0, found+1);
  }
  // Get the filename extension
  found = baseName_.find_last_of(".");
  if (found == std::string::npos) {
    extension_.clear();
  } else {
    extension_ = baseName_.substr(found);
  }
  // See if the extension is one of the 2 recognized compression extensions.
  // If file has a compression format extension, look for another extension.
  if ( extension_ == ".gz" || extension_ == ".bz2" ) {
    compressExt_ = extension_;
    // Get everything before the compressed extension
    std::string compressPrefix = baseName_.substr(0,found);
    // See if there is another extension
    found = compressPrefix.find_last_of(".");
    if (found == std::string::npos)
      // No other extension
      extension_.clear();
    else
      extension_ = compressPrefix.substr(found);
  } else
    compressExt_.clear();
  return 0;
}

int File::Name::Append( std::string const& suffix ) {
  if (fullPathName_.empty()) return 1;
  fullPathName_.append( suffix );
  baseName_.append( suffix );
  return 0;
}

File::Name File::Name::AppendName( std::string const& suffix ) const {
  Name out( *this );
  out.Append( suffix );
  return out;
}

//TODO make this more efficient by just modifying full and base names
File::Name File::Name::PrependName( std::string const& prefix ) const {
  Name out;
  out.SetName_NoExpansion(dirPrefix_ + prefix + baseName_);
  return out;
}

File::Name File::Name::PrependExt( std::string const& extPrefix ) const {
  Name out( *this );
  // Find location of extension.
  size_t found = out.baseName_.rfind( extension_ );
  // Remove extension.
  out.baseName_.resize( found );
  // Insert extPrefix to just before extension and re-add extension.
  out.baseName_.append( extPrefix + extension_ + compressExt_ );
  // Update full path name.
  out.fullPathName_ = dirPrefix_ + out.baseName_;
  //mprintf("DEBUG: fullPathName= '%s'\n"
  //        "       baseName=     '%s'\n"
  //        "       extension=    '%s'\n"
  //        "       compressExt=  '%s'\n"
  //        "       dirPrefix=    '%s'\n",
  //        out.fullPathName_.c_str(), out.baseName_.c_str(), out.extension_.c_str(),
  //        out.compressExt_.c_str(), out.dirPrefix_.c_str());
  return out;
}

// ===== File::Base ============================================================
File::Base::Base() :
  file_size_(0UL),
  debug_(0),
  access_(READ),
  compressType_(NO_COMPRESSION),
  isOpen_(false),
  isStream_(false)
{}

const char* File::Base::AccessTypeName_[] = {
    "read", "write", "append", "update"
};

// File::Base::Setup()
int File::Base::Setup(Name const& fnameIn, AccessType accessIn)
{
  if (isOpen_) Close();
  access_ = accessIn;
  isOpen_ = false;
  file_size_ = 0;
  compressType_ = NO_COMPRESSION;
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
    // Get basic file information
    if (Exists(fname_)) {
      // File exists. Get size and compression info.
      struct stat frame_stat;
      if (stat(fname_.full(), &frame_stat) == -1) {
        mprinterr("Error: Could not find file status for %s\n", fname_.full());
        if (debug_>0) 
          perror("     Error from stat: ");
        return 1;
      }
      file_size_ = frame_stat.st_size;
      // ID compression by magic number - open for binary read access
      FILE* fIn = fopen(fname_.full(), "rb");
      if ( fIn == 0 ) { 
        mprinterr("Error: Could not open %s for hex signature read.\n", fname_.full());
        return 1;
      }
      // Read first 3 bytes
      unsigned char magic[3];
      magic[0] = 0; 
      magic[1] = 0; 
      magic[2] = 0;
      size_t numread = fread(magic, 1, 3, fIn);
      fclose(fIn);
      if (numread == 0)
        mprintf("Warning: File %s is empty\n", fname_.full());
      else if (numread < 3 ) {
        mprinterr("Warning: Could only read first %zu bytes of file %s.\n", numread, fname_.full());
      } else {
        if (debug_>0) mprintf("\t    Hex sig: %x %x %x", magic[0],magic[1],magic[2]);
        // Check compression
        if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8)) {
          if (debug_>0) mprintf(", Gzip file.\n");
          compressType_ = GZIP;
        } else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68)) {
          if (debug_>0) mprintf(", Bzip2 file.\n");
          compressType_ = BZIP2;
        } else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3)) {
          if (debug_>0) mprintf(", Zip file.\n");
          compressType_ = ZIP;
        } else {
          if (debug_>0) mprintf(", No compression.\n");
        }
      }
    } else {
      // File does not exist. Determine compression via extension.
      if (fname_.Compress() == ".gz")
        compressType_ = GZIP;
      else if (fname_.Compress() == ".bz2")
        compressType_ = BZIP2;
    }
  } // END file is not stream
  if (debug_ >= 0) { // FIXME
    mprintf("\tFILE INFO:");
    if (isStream_)
      mprintf(" STREAM\n");
    else
      mprintf(" %s\n", fname_.full());
    mprintf("\t  Size= %li\n", file_size_);
    const char* compTypeStr[4] = {"None", "Gzip", "Bzip", "Zip "};
    mprintf("\t  Compression= %s\n", compTypeStr[compressType_]);
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

bool File::Exists(std::string const& fname) {
  return File::Exists( Name(fname) );
}
