#ifndef INC_FILE_H
#define INC_FILE_H
#include <string>
#include <vector>
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
  /// Print error message corresponding to 'false' value from 'Exists()'
  void ErrorMsg(const char*);
  /// \return true if file exists and is accessible.
  bool Exists(std::string const&); // TODO remove?
  bool Exists(Name const&);
}

/** Class to hold file name, extension, etc. */
class File::Name {
  public:
    Name() {}
    Name(std::string const& s) { SetName(s); }
    Name(const char* s) { SetName( std::string(s) ); }
    Name(const Name&);
    Name& operator=(const Name&);
    /// Set file name and extensions, perform expansion as necessary.
    int SetName(std::string const&);
    /// Set file name, no expansions.
    int SetName_NoExpansion(std::string const&);
    /// Append given string to file name but do not change extension info.
    int Append(std::string const&); //TODO this can probably replace AppendNumber some places
    /// \return File name with given string appended - do not change extension info.
    Name AppendName(std::string const&) const;
    /// \return File name with given string prepended to base file name.
    Name PrependName(std::string const&) const;
    /// \return File name with given string prepended to extension
    Name PrependExt(std::string const&) const;
    /// Clear File name
    void clear();
    /// \return true if string matches full or base file name.
    bool MatchFullOrBase(std::string const&) const;

    const std::string& Full()      const { return fullPathName_;         }
    const std::string& Base()      const { return baseName_;             }
    const std::string& Ext()       const { return extension_;            }
    const char* full()             const { return fullPathName_.c_str(); }
    const char* base()             const { return baseName_.c_str();     }
    const char* ext()              const { return extension_.c_str();    }
    const std::string& Compress()  const { return compressExt_;          }
    const std::string& DirPrefix() const { return dirPrefix_;            }
    bool empty()                   const { return fullPathName_.empty(); }
  private:
    std::string fullPathName_;
    std::string baseName_;
    std::string extension_;
    std::string compressExt_;
    std::string dirPrefix_;
};

/** Base class that all files will inherit. */ 
class File::Base {
  public:
    Base();
    Name const& Filename() const { return fname_;  }
    AccessType Access()    const { return access_; }
    int Setup(const char*, AccessType);
  protected:
    virtual int Open() = 0;
    virtual int Close() = 0;
  private:
    Name fname_;
    off_t file_size_;           ///< Actual file size
    int debug_;
    AccessType access_;         ///< Current file access
    CompressType compressType_; ///< Type of compression present
    bool isOpen_;               ///< True if file is open and ready for IO
    bool isStream_;             ///< True if file is to/from a stream.
};
#endif
