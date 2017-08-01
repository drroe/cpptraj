#include <cstdarg>    // va_X functions
#include <cstdio>     // vsprintf
#include <cstring>    // strlen
#include <algorithm>  // std::max
#include "CpptrajFile.h"
#include "CpptrajStdio.h" // DEBUG

CpptrajFile::CpptrajFile() : BUF_SIZE_(0), linebuffer_(0) {}

CpptrajFile::CpptrajFile(int d) : BasicFile(d), BUF_SIZE_(0), linebuffer_(0) {}

CpptrajFile::~CpptrajFile() {
  if (linebuffer_ != 0) delete[] linebuffer_;
}

/** Copy constructor. Copies closed. */
CpptrajFile::CpptrajFile(CpptrajFile const& rhs) :
  BasicFile(rhs),
  BUF_SIZE_(rhs.BUF_SIZE_),
  linebuffer_(0)
{
  if (BUF_SIZE_ > 0)
    linebuffer_ = new char[ BUF_SIZE_ ];
}

/** Assignment. Assigns closed. */
CpptrajFile& CpptrajFile::operator=(CpptrajFile const& rhs) {
  if (this != &rhs) {
    if (linebuffer_ != 0) delete linebuffer_;
    BasicFile::operator=(rhs);
    BUF_SIZE_ = rhs.BUF_SIZE_;
    if (BUF_SIZE_ > 0)
      linebuffer_ = new char[ BUF_SIZE_ ];
  }
  return *this;
}

// -----------------------------------------------------------------------------
// TODO check if file is open?
// CpptrajFile::Printf()
/** Take the formatted string and write it to file using Write.
  */
void CpptrajFile::Printf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(linebuffer_,format,args);
  IO()->Write(linebuffer_, strlen(linebuffer_));
  va_end(args);
}

std::string CpptrajFile::GetLine() {
  if (IO()->Gets(linebuffer_, BUF_SIZE_) != 0) {
    //mprinterr("Error: Getting line from %s\n", Filename().full());
    return std::string();
  }
  return std::string(linebuffer_);
}

const char* CpptrajFile::NextLine() {
  if (IO()->Gets(linebuffer_, BUF_SIZE_) != 0) {
    //mprinterr("Error: Reading line from %s\n", Filename().full());
    return 0;
  }
  return linebuffer_;
}

void CpptrajFile::SetupBuffer(unsigned int sizeIn) {
  BUF_SIZE_ = sizeIn;
  if (linebuffer_ != 0) delete[] linebuffer_;
  linebuffer_ = new char[ BUF_SIZE_ + 1 ];
  linebuffer_[BUF_SIZE_] = '\0';
}

// -----------------------------------------------------------------------------
int CpptrajFile::InternalSetup() {
  int firstLineSize = BasicSetup();
  if (firstLineSize < 0) return 1;
  SetupBuffer( std::max(1024, firstLineSize+1) );
  mprintf("\t  BUF_SIZE_ = %u\n", BUF_SIZE_);
  return 0;
}
