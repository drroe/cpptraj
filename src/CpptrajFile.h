#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include "BasicFile.h"
/// Perform simple file IO with a static internal buffer.
class CpptrajFile : public BasicFile {
  public:
    CpptrajFile();
    CpptrajFile(int);
    virtual ~CpptrajFile();
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);
    // ----- IO Routines -------------------------
    int Gets(char* buf, int num)           { return IO()->Gets(buf, num);  }
    int Write(const void* buf, size_t num) { return IO()->Write(buf, num); }
    int Read(void* buf, size_t num)        { return IO()->Read(buf, num);  }
    int Seek(off_t offset)                 { return IO()->Seek(offset);    }
    int Rewind()                           { return IO()->Rewind();        }
    int Flush()                            { return IO()->Flush();         }
    off_t Tell()                           { return IO()->Tell();          }
    /// Printf using the Write routine.
    void Printf(const char*, ...);
    /// Get next line as a string
    std::string GetLine();
    /// Get next line and return pointer to raw buffer
    const char* NextLine();
  private:
    // -------------------------------------------
    int InternalSetup();
    // -------------------------------------------
    void SetupBuffer(unsigned int);

    unsigned int BUF_SIZE_;
  protected:
    char* linebuffer_;
};
#endif
