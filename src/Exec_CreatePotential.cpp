#include "Exec_CreatePotential.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Exec_CreatePotential::Exec_CreatePotential() :
  Exec(GENERAL)
{
  SetHidden(true);
}

// Exec_CreatePotential::Help()
void Exec_CreatePotential::Help() const
{

}

// Exec_CreatePotential::Execute()
Exec::RetType Exec_CreatePotential::Execute(CpptrajState& State, ArgList& argIn)
{
  int nrep = argIn.getKeyInt("nrep", 0);
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: No 'name' specified for potential function.\n");
    return CpptrajState::ERR;
  }

  DataSet* ds = State.DSL().AddSet(DataSet::POTENTIALFXN, dsname);
  if (ds == 0) {
    mprinterr("Error: Could not allocate potential function.\n");
    return CpptrajState::ERR;
  }

  mprintf("    CREATEPOTENTIAL: Created potential function '%s'\n", ds->legend());
  return CpptrajState::OK;
}
