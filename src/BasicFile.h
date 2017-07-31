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

    /// \return 1 if the file contains carriage returns in addition to newlines
    int IsDos()                 const { return isDos_;                }
    /// \return uncompressed file size (just size if file is not compressed).
    unsigned int UncompressedSize() const;

    // ----- IO Routines -------------------------
    int Gets(char* buf, int num)           { return IO_->Gets(buf, num);  }
    int Write(const void* buf, size_t num) { return IO_->Write(buf, num); }
    int Read(void* buf, size_t num)        { return IO_->Read(buf, num);  }
    int Seek(off_t offset)                 { return IO_->Seek(offset);    }
    int Rewind()                           { return IO_->Rewind();        }
    int Flush()                            { return IO_->Flush();         }
    off_t Tell()                           { return IO_->Tell();          }
    /// Printf using the Write routine.
    void Printf(const char*, ...);
    /// Get next line as a string
    std::string GetLine();
    /// Get next line and return pointer to raw buffer
    const char* NextLine();
  private:
    enum FileType { UNKNOWN_TYPE=0, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };
    static const char* FileTypeName_[];
    void SetupBuffer(unsigned int);
    // -------------------------------------------
    /// Set up file
    int InternalSetup();
    /// Open file
    int InternalOpen();
    /// Close file
    void InternalClose();
    // -------------------------------------------
    int SetupFileIO( FileType );
    void Reset();

    FileIO* IO_;                     ///< The interface to basic IO operations.
    int isDos_;                      ///< 1 if CR present, need to count them as newlines
    unsigned int uncompressed_size_; ///< If compressed, uncompressed file size
    unsigned int BUF_SIZE_;          ///< Buffer size
    FileType fileType_;              ///< File type (determines IO)
    char* linebuffer_;               ///< Used for IO functions.
};
#endif
