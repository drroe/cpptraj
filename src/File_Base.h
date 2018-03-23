#ifndef INC_FILE_BASE_H
#define INC_FILE_BASE_H
#include "File.h"
namespace File {
/** Base class that all files will inherit. Contains basic information such
  * as the file name, what access the file is set up for, etc. Files that
  * inherit this class must implement the following functions:
  *   InternalSetup() - Called by Setup(), prepares file with specific access.
  *   InteralOpen()   - Called by Open(), open file for IO.
  *   InternalClose() - Called by Close(), close file.
  */
class Base {
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

} /* END namespace File */
#endif
