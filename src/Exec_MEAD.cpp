#include "Exec_MEAD.h"
#include "CpptrajStdio.h"
#include "MeadInterface.h"
#include "StringRoutines.h"

// Exec_MEAD::Help()
void Exec_MEAD::Help() const
{
  mprintf("\t[ogm <ngridpoints>,<gridspacing>] ...\n");
}

// Exec_MEAD::Execute()
Exec::RetType Exec_MEAD::Execute(CpptrajState& State, ArgList& argIn)
{
  Cpptraj::MeadInterface MEAD;

  std::string ogmstr = argIn.GetStringKey("ogm");
  while (!ogmstr.empty()) {
    // Format: N,spacing TODO centering
    ArgList ogmarg( ogmstr, "," );
    if (ogmarg.Nargs() < 2) {
      mprinterr("Error: Malformed ogm key; expected <N>,<spacing>\n");
      return CpptrajState::ERR;
    }
    int ngridpts = convertToInteger( ogmarg[0] );
    double spacing = convertToDouble( ogmarg[1] );
    mprintf("\tAdding grid of %i points, spacing %g\n", ngridpts, spacing);
    if (MEAD.AddGrid(ngridpts, spacing, Vec3(0,0,0))) {
      mprinterr("Error: Adding MEAD grid.\n");
      return CpptrajState::ERR;
    }
    ogmstr = argIn.GetStringKey("ogm");
  }

  MEAD.Print();

  return CpptrajState::OK;    
}
