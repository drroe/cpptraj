#ifndef INC_AMBERPARAMFILE_H
#define INC_AMBERPARAMFILE_H
#include <vector>
#include <string>
class ParameterSet;
//class BufferedLine;
class FileName;
/// Used to read in Amber parameters from Amber FF/FRCMOD file.
class AmberParamFile {
  public:
    AmberParamFile();
    int ReadParams(ParameterSet&, FileName const&, std::string const&, int) const;
    int WriteParams(ParameterSet&, FileName const&, int) const;
  private:
    static int read_symbols(const char*, std::vector<std::string>&, int);
    //int ReadInput(std::string&, BufferedLine&) const;
};
#endif
