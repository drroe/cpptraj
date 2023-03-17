#include "Exec_MEAD.h"
#include "CpptrajStdio.h"
#include "MeadInterface.h"
#include "StringRoutines.h"

// Exec_MEAD::Help()
void Exec_MEAD::Help() const
{
  mprintf("\t[ogm <ngridpoints>,<gridspacing>] ...\n"
          "\t[crdset <COORDS set>] [radmode {gb|parse|vdw}]\n");
}

/** MEAD potential. */
int Exec_MEAD::Potential(Cpptraj::MeadInterface& MEAD, ArgList& argIn) const {
  // Sanity checks
  if (!MEAD.HasFDM()) {
    mprinterr("Error: No MEAD grid allocated.\n");
    return 1;
  }
  if (!MEAD.HasAtoms()) {
    mprinterr("Error: No MEAD atoms allocated.\n");
    return 1;
  }

  return 0;
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

  Cpptraj::MeadInterface::Radii_Mode radiiMode;
  std::string radmode = argIn.GetStringKey("radii");
  if (radmode == "gb")
    radiiMode = Cpptraj::MeadInterface::GB;
  else if (radmode == "parse")
    radiiMode = Cpptraj::MeadInterface::PARSE;
  else if (radmode == "vdw")
    radiiMode = Cpptraj::MeadInterface::VDW;
  else
    radiiMode = Cpptraj::MeadInterface::GB;

  switch (radiiMode) {
    case Cpptraj::MeadInterface::GB : mprintf("\tUsing GB radii.\n"); break;
    case Cpptraj::MeadInterface::PARSE : mprintf("\tUsing PARSE radii.\n"); break;
    case Cpptraj::MeadInterface::VDW : mprintf("\tUsing VDW radii.\n"); break;
  }

  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: No COORDS set specified.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return CpptrajState::ERR;
  }
  // TODO multiple frames
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' has no frames.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  Frame frameIn = CRD->AllocateFrame();
  CRD->GetFrame(0, frameIn);
  if (MEAD.SetupAtoms( CRD->Top(), frameIn, radiiMode )) {
    mprinterr("Error: Setting up frame/topology failed.\n");
    return CpptrajState::ERR; 
  }

  MEAD.Print();

  int err = 0;
  if (argIn.hasKey("potential")) {
    err = Potential( MEAD, argIn );
  } else {
    mprinterr("Error: No MEAD calculation keywords given.\n");
    err = 1;
  }

  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;    
}
