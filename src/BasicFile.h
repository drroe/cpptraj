#ifndef INC_BASICFILE_H
#define INC_BASICFILE_H
#include "File_Base.h"
#include "FileIO.h"
/// Abstract base class for files that will do simple file IO
class BasicFile : public File::Base {
  public:
    /// CONSTRUCTOR
    BasicFile();
    /// CONSTRUCTOR - debug level
    BasicFile(int);
    /// DESTRUCTOR
    virtual ~BasicFile(); // Virtual since this class is inherited.
    /// COPY CONSTRUCTOR
    BasicFile(BasicFile const&);
    /// Assignment
    BasicFile& operator=(BasicFile const&);

    /// \return 1 if the file contains carriage returns in addition to newlines
    int IsDos()                 const { return isDos_;                }
    /// \return uncompressed file size (just size if file is not compressed).
    unsigned int UncompressedSize() const;
#   ifdef MPI
    bool IsMPI() const { return (fileType_ == MPIFILE || fileType_ == MPISHARED); }
    // FIXME should this be a Base routine?
    int ParallelOpenFile(AccessType, Parallel::Comm const&, bool);
    int ParallelOpenFile(AccessType t, Parallel::Comm const& c) {
      return ParallelOpenFile(t, c, false); //  Default to no shared write
    }
    int ParallelOpenFile(Parallel::Comm const& c) { return ParallelOpenFile(access_, c); }
    int ParallelOpenFile(Parallel::Comm const& c, bool b) {
      return ParallelOpenFile(access_, c, b);
    }
#   endif
  protected:
    /// Set up basic file IO and determine file characteristics.
    int BasicSetup();
    /// \return pointer to basic file IO interface
    FileIO* IO() { return IO_; }
  private:
    enum FileType { UNKNOWN_TYPE=0, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };
    static const char* FileTypeName_[];
    // -------------------------------------------
    /// Close IO_ interface
    void InternalClose();
    /// Open IO_ interface
    int InternalOpen() { return OpenIO(); }
    /// Open IO interface
    int OpenIO();

    // -------------------------------------------
    /// Set up IO_ interface
    int SetupFileIO( FileType );
    /// Reset all file variables
    void Reset();

    FileIO* IO_;                     ///< The interface to basic IO operations.
    int isDos_;                      ///< 1 if CR present, need to count them as newlines
    unsigned int uncompressed_size_; ///< If compressed, uncompressed file size
    FileType fileType_;              ///< File type (determines IO)
};
#endif
