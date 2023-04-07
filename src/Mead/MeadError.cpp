#include "MeadError.h"
#include "../CpptrajStdio.h"
#include "../../mead/MEADexcept.h"

/** Print MEAD error message. */
int Cpptraj::Mead::ERR(const char* fxn, MEADexcept& e) {
  mprinterr("Error: MEAD error in '%s': '%s' '%s' '%s'\n", fxn,
              e.get_error1().c_str(),
              e.get_error2().c_str(),
              e.get_error3().c_str());
  return 1;
}
