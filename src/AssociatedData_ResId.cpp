#include "AssociatedData_ResId.h"
#include "CpptrajStdio.h"

/** Help text */
const char* AssociatedData_ResId::HelpText = "\0";

/** Process Args */
int AssociatedData_ResId::ProcessAdataArgs(ArgList& argIn) {
  return 0;
}

/** Print info */
void AssociatedData_ResId::Ainfo() const {
  mprintf(" (Res name %s termType %s)", *resName_, Cpptraj::Structure::terminalStr(termType_));
}
