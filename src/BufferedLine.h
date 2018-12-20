#ifndef INC_BUFFEREDLINE_H
#define INC_BUFFEREDLINE_H
#include <vector>
#include "BasicFile.h"
/// Used to buffer text files that will be read line-by-line
class BufferedLine : private BasicFile {
  public:
    BufferedLine();
    ~BufferedLine();
    /// \return pointer to next line in the buffer.
    const char* Line();
    /// Convert current line into tokens
    int TokenizeLine(const char*);
    /// \return next token, null-delimited.
    const char* NextToken();
    /// \return specified token, not null-delimited.
    inline const char* Token(int);
    /// Open file for reading, set up buffer.
    int OpenFileRead( File::Name const& );
    /// Open the file (must be set up), set up buffer.
    int OpenFile();
    /// Open the file for writing, no advanced buffering TODO implement buffering?
    int OpenWrite(File::Name const&);
    /// Printf using the Write routine.
    void Printf(const char*, ...);
    /// \return current line number
    int LineNumber()          const { return nline_;          }
    /// \return pointer to buffer
    const char* Buffer()      const { return buffer_;         }
    /// \return Pointer to current buffer position.
    const char* CurrentLine() const { return bufferPosition_; }
    /// \return Next line as string.
    inline std::string GetLine();
    // Members of Base that should be public
    using Base::Filename;
    using Base::Close;
    using Base::SetDebug;
  private:
    /// Basic file IO setup. Clean and reallocate buffer for currentBufSize_.
    int InternalSetup();
    /// Default initial buffer size.
    static const size_t DEFAULT_BUFFERSIZE = 16384;

    size_t currentBufSize_; ///< Current size of buffer.
    char* buffer_;         ///< Beginning of character buffer.
    char* bufferPosition_; ///< Position in buffer/start of current line.
    /// Array of pointers to beginning and ends of tokens in current line. 
    std::vector<char*> tokens_;
    size_t tokenidx_;      ///< Current position in tokens array
    char saveChar_;        ///< Saved last char of current token
    char* lineEnd_;        ///< End of current line in buffer
    char* endBuffer_;      ///< End of character buffer
    size_t nline_;         ///< Current line number.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
const char* BufferedLine::Token(int idx) {
  if (idx < 0 || idx >= (int)tokens_.size()) return 0;
  return tokens_[idx];
}

std::string BufferedLine::GetLine() {
  const char* ptr = Line();
  if (ptr == 0) return std::string();
  return std::string(ptr);
}
#endif
