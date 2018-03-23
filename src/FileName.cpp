#include "FileName.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // validInteger, convertToInteger, integerToString
#include "File_WordExp.h"

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
  Sarray names = WordExp( nameIn );
  int err = 1;
  if (!names.empty()) {
    if (names.size() > 1)
      mprintf("Warning: '%s' matches multiple files, only using '%s'\n",
              nameIn.c_str(), names.front().c_str());
    err = SetName_NoExpansion( names.front() );
  }
  return err;
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
