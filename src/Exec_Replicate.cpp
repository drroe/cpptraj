#include "Exec_Replicate.h"
#include "CpptrajStdio.h"
#include "TopInfo.h" // DEBUG

// Exec_Replicate::Help()
void Exec_Replicate::Help() const
{
  mprintf("\t<mask> [%s]\n", DataSetList::TopIdxArgs);
}

// Exec_Replicate::Execute()
Exec::RetType Exec_Replicate::Execute(CpptrajState& State, ArgList& argIn)
{
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  // Check if this topology has already been used to set up an input
  // trajectory, as this will break the traj read.
  std::string fname = State.TopUsedInInputTraj(parm);
  if (!fname.empty()) {
    mprinterr("Error: Topology '%s' has already been used to set up trajectory '%s'.\n"
              "Error:   To strip this topology use the 'strip' action.\n",
              parm->c_str(), fname.c_str());
    return CpptrajState::ERR;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  if (parm->SetupIntegerMask( tempMask )) return CpptrajState::ERR;
  mprintf("\tReplicating atoms in topology '%s'\n", parm->c_str());
  tempMask.MaskInfo();

  if (parm->ReplicateAtoms( tempMask, 1 )) {
    mprinterr("Error: Replication failed.\n");
    return CpptrajState::ERR;
  }

  TopInfo tinfo( parm );
  tinfo.PrintReplicateInfo("*");
  return CpptrajState::OK;
}
