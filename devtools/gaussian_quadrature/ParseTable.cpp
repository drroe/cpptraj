#include "BufferedLine.h"
#include <string>
#include <cstdio>

int main(int argc, char** argv) {
  BufferedLine infile;
  std::string fname;
  for (int iarg = 1; iarg < argc; iarg++)
  {
    std::string arg( argv[iarg] );
    if (!arg.empty())
      fname = arg;
  }
  if (fname.empty()) {
    fprintf(stderr,"No file specified.\n");
    return 1;
  }
  fprintf(stderr,"Parsing %s\n", fname.c_str());
  return 0;
}
