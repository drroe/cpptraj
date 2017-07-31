#ifndef INC_BASICFILE_H
#define INC_BASICFILE_H
#include "File.h"
#include "FileIO.h"
/// Basic file
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

  private:
    enum FileType { UNKNOWN_TYPE=0, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };
    static const char* FileTypeName_[];
    // -------------------------------------------
    /// Set up file
    int InternalSetup();
    /// Open file
    int InternalOpen();
    /// Close file
    void InternalClose();
    // -------------------------------------------
    int SetupRead();
    int SetupWrite();
    int SetupAppend();
    int SetupFileIO( FileType );
    void Reset();

    FileIO* IO_;                     ///< The interface to basic IO operations.
    int isDos_;                      ///< 1 if CR present, need to count them as newlines
    unsigned int uncompressed_size_; ///< If compressed, uncompressed file size
    unsigned int BUF_SIZE_;          ///< Buffer size
    FileType fileType_;              ///< File type (determines IO)
};
#endif
