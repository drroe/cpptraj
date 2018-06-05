#ifndef INC_CHARMMPARAMFILE_H
#define INC_CHARMMPARAMFILE_H
#include "ParameterSet.h"
#include "FileName.h"
/// Used to read in CHARMM parameters from CHARMM parameter file.
class CharmmParamFile {
  public:
    CharmmParamFile() {}
    int ReadParams(ParameterSet&, FileName const&, int);
  private:
};
#endif