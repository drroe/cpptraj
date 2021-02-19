#include "Exec_Output.h"
#include "CpptrajStdio.h"

// Exec_Output::Help()
void Exec_Output::Help() const
{
  mprintf("\t{to <file>|reset} [err]\n");
}

// Exec_Output::Execute()
Exec::RetType Exec_Output::Execute(CpptrajState& State, ArgList& argIn)
{
  bool for_stderr = argIn.hasKey("err");
  std::string fname = argIn.GetStringKey("to");
  bool do_reset = argIn.hasKey("reset");
  if (do_reset && !fname.empty()) {
    mprinterr("Error: Specify either 'to <file>' or 'reset', not both.\n");
    return CpptrajState::ERR;
  }
  if (!do_reset && fname.empty()) {
    mprinterr("Error: Must specify either 'to <file>' or 'reset'.\n");
    return CpptrajState::ERR;
  }
  int err = 0;
  if (do_reset) {
    if (for_stderr)
      err = ErrToFile(0);
    else
      err = OutputToFile(0);
  } else {
    if (for_stderr)
      err = ErrToFile(fname.c_str());
    else
      err = OutputToFile(fname.c_str());
  }
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}
