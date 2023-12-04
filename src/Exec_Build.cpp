#include "Exec_Build.h"
#include "CpptrajStdio.h"
#include "Structure/GenerateAngles.h"

// Exec_Build::Help()
void Exec_Build::Help() const
{
  mprintf("\tcrdset <COORDS set> [frame <#>]\n");
}

// Exec_Build::Execute()
Exec::RetType Exec_Build::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get input coords
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify input COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
  // Get frame from input coords
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);
  // Get modifiable topology
  Topology& topIn = *(coords.TopPtr());

  // Generate angles
  if (Cpptraj::Structure::GenerateAngles(topIn)) {
    mprinterr("Error: Could not generate angles/dihedrals for '%s'\n", topIn.c_str());
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
