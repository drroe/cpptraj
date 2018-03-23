#ifndef INC_FILE_H
#define INC_FILE_H
#include <string>
#include <vector>
#include "FileName.h"
/// This namespace contains useful file-related routines.
namespace File {
  /// File name, path, extension etc.
  class Name;
  /// Array of file names
  typedef std::vector<Name> NameArray;
  /// Expand given expression to array of file names
  NameArray ExpandToFilenames(std::string const&);
  /// File access types
  enum AccessType { READ=0, WRITE, APPEND, UPDATE };
  /// File compression types
  enum CompressType { NO_COMPRESSION=0, GZIP, BZIP2, ZIP };
  /// Basic file; no IO routines.
  class Base;
  /// Print error message with given file name corresponding to 'false' value from 'Exists()'
  void ErrorMsg(const char*);
  /// return last file error message
  const char* StrError();
  /// \return true if file exists and is accessible.
  bool Exists(Name const&);
  //bool Exists(std::string const&); // TODO remove?
  /// Given lowest replica name search for other replica names
  NameArray SearchForReplicas(Name const&, int);
# ifdef MPI
  /// Each rank searches for replica based on lowest replica number.
  NameArray SearchForReplicas(Name const&, bool, int, int, int);
# endif
}

/** Base class that all files will inherit. */ 
class File::Base {
  public:
    Base();
    Base(int);
    virtual ~Base() {} // Virtual since class is inherited
    Base(Base const&);
    Base& operator=(Base const&);
    Name const& Filename()     const { return fname_;        }
    unsigned int Size()        const { return file_size_;    }
    int Debug()                const { return debug_;        }
    AccessType Access()        const { return access_;       }
    CompressType Compression() const { return compressType_; }
    bool IsOpen()              const { return isOpen_;       }
    bool IsStream()            const { return isStream_;     }
    bool IsPresent()           const { return isPresent_;    }
    /// \return string based on current access
    const char* accessStr()    const { return AccessTypeName_[access_]; }

    void SetDebug(int d) { debug_ = d; }
    /// Set up file for given access but do not open.
    int Setup(Name const&, AccessType);
    /// Open the file with current access.
    int Open();
    /// Setup the file with given access and open.
    int Open(Name const&, AccessType);
    /// Close the file.
    void Close();
  protected:
    /// File-specific setup, called by Setup()
    virtual int InternalSetup() = 0;
    /// File-specific open, called by Open()
    virtual int InternalOpen() = 0;
    /// File-specific close, called by Close()
    virtual void InternalClose() = 0;
  private:
    static const char* AccessTypeName_[];

    Name fname_;
    unsigned int file_size_;           ///< Actual file size
    int debug_;
    AccessType access_;         ///< Current file access
    CompressType compressType_; ///< Type of compression present
    bool isOpen_;               ///< True if file is open and ready for IO
    bool isStream_;             ///< True if file is to/from a stream.
    bool isPresent_;            ///< True if file exists.
};
#endif
