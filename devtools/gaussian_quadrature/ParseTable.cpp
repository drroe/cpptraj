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
  fprintf(stderr,"Parsing '%s'\n", fname.c_str());

  if (infile.OpenFileRead( fname )) {
    fprintf(stderr,"Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  const char* ptr = infile.Line();
  std::string currentExpression;
  bool in_array = false;
  while (ptr != 0) {
    std::string currentLine( ptr );
    if (!currentLine.empty()) {
      for (std::string::const_iterator it = currentLine.begin(); it != currentLine.end(); ++it) {
        if (!isspace( *it )) {
          if (*it == '$') {
            if (in_array) {
              fprintf(stderr,"Error: Encountered new array while in an array.\n");
              fprintf(stderr,"Error: '%s'\n", currentLine.c_str());
              return 1;
            }
            in_array = true;
          }
          if (in_array) {
            currentExpression += *it;
          }
          if (*it == ';') {
            printf("Expression: '%s'\n", currentExpression.c_str());
            currentExpression.clear();
            in_array = false;
          }
        } // END not space
      } // END loop over current line
    } // END current line not empty
    ptr = infile.Line();
  }
  infile.CloseFile();
  return 0;
}
