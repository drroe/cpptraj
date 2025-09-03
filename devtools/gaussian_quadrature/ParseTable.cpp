#include "ArgList.h"
#include "BufferedLine.h"
#include <string>
#include <vector>
#include <cstdio>

//double convert_

int string_to_array(std::vector<double>& array, std::string const& str)
{
  // String should have no whitespace
  ArgList line;
  line.SetList( str, "$[]=(,);" );
  if (line.Nargs() < 3) {
    // Probably blank
    return 0;
  }
  if (line[1] == "2") // DEBUG
    line.PrintDebug();
  line.MarkArg(0);
  int order = line.getNextInteger(-1);
  if (order < 1) {
    fprintf(stderr, "Error: invalid order: %s\n", str.c_str());
    return -1;
  }
  line.MarkArg(2);
  array.clear();
  for (int n = 0; n < order; n++)
    array.push_back( line.getNextDouble(-999) );
  if (order == 2) {
    printf("DEBUG: order %i\n", order);
    for (std::vector<double>::const_iterator it = array.begin(); it != array.end(); ++it)
      printf("\t%16.8f\n", *it);
  }
  
  return order;
}

int read_file(std::string const& fname)
{
  BufferedLine infile;
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
            std::vector<double> array;
            string_to_array( array, currentExpression );
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

int main(int argc, char** argv) {
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
  read_file(fname);

  return 0;
}
