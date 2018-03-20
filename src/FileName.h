#ifndef INC_FILENAME_H
#define INC_FILENAME_H
#include <string>
namespace File {
/// File name, path, extension etc.
class Name {
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
} /* namespace File */
#endif
