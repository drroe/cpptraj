#include "Exec_Source.h"
#include "CpptrajStdio.h"

// Exec_Source::Help()
void Exec_Source::Help() const
{

}

// Exec_Source::Execute()
Exec::RetType Exec_Source::Execute(CpptrajState& State, ArgList& argIn)
{
  // source <file>
  std::string fname = argIn.GetStringNext();
  if (fname.empty()) {
    mprinterr("Error: No filename given for 'source'.\n");
    return CpptrajState::ERR;
  }
  // If file is in this dir, read it.
  // Otherwise, prepend AMBERHOME/dat/leap/cmd
  std::string fileName;
  if ( File::Exists( fname ) ) {
    fileName = fname;
  } else {
    const char* env = getenv("AMBERHOME");
    if (env == 0) {
      mprinterr("Error: %s not found and AMBERHOME is not set.\n", fname.c_str());
      return CpptrajState::ERR;
    }
    fileName = std::string(env) + "/dat/leap/cmd/" + fname;
    if (!File::Exists(fileName)) {
      mprinterr("Error: %s not found.\n", fileName.c_str());
      return CpptrajState::ERR;
    }
  }
  mprintf("\tReading %s\n", fileName.c_str());

  return CpptrajState::OK;
}
