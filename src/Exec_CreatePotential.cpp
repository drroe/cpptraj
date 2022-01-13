#include "Exec_CreatePotential.h"
#include "CpptrajStdio.h"
#include "DataSet_PotentialFxn.h"
#include "Potential/MdOpts.h"
#include "Potential/PotentialFunction.h"
#include "DataSet_Coords.h"
#include "DataSet_Topology.h"

/** CONSTRUCTOR */
Exec_CreatePotential::Exec_CreatePotential() :
  Exec(GENERAL)
{
  SetHidden(true);
}

// Exec_CreatePotential::Help()
void Exec_CreatePotential::Help() const
{
  MdOpts::PrintHelp();
}

// Exec_CreatePotential::Execute()
Exec::RetType Exec_CreatePotential::Execute(CpptrajState& State, ArgList& argIn)
{
/*  int nrep = argIn.getKeyInt("nrep", 1);
  if (nrep < 1) {
    mprinterr("Error: 'nrep' cannot be less than 1.\n");
    return CpptrajState::ERR;
  }*/
  Topology* parm = 0;
  std::string crdset = argIn.GetStringKey("crdset");
  std::string parmarg = argIn.GetStringKey("parm");
  if (!crdset.empty()) {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL().FindSetOfGroup(crdset, DataSet::COORDINATES);
    if (ds == 0) return CpptrajState::ERR;
    mprintf("\tUsing topology from data set '%s'\n", ds->legend());
    parm = ds->TopPtr();
  } else if (!parmarg.empty()) {
    DataSet_Topology* ds = (DataSet_Topology*)State.DSL().FindSetOfType(parmarg, DataSet::TOPOLOGY);
    if (ds == 0) return CpptrajState::ERR;
    mprintf("\tUsing topology '%s'\n", ds->legend());
    parm = ds->TopPtr();
  }
  int nbonds = 1;
  int nangles = 1;
  int ndihedrals = 1;
  if (parm != 0) {
    nbonds = parm->Nbonds();
    nangles = parm->Nangles();
    ndihedrals = parm->Ndihedrals();
  }
  // ----------
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: No 'name' specified for potential function.\n");
    return CpptrajState::ERR;
  }
  bool use_openmm = argIn.hasKey("openmm");
  bool useNonbond = argIn.hasKey("nonbond");
  bool use_replicate = argIn.hasKey("replicate");

  // Get potential function options.
  MdOpts opts;
  if (opts.GetOptsFromArgs(argIn)) return CpptrajState::ERR;
  opts.PrintOpts();

  // -----------------------------------
  
  DataSet* ds = State.DSL().AddSet(DataSet::POTENTIALFXN, dsname);
  if (ds == 0) {
    mprinterr("Error: Could not allocate potential function.\n");
    return CpptrajState::ERR;
  }

  DataSet_PotentialFxn& pfxn = static_cast<DataSet_PotentialFxn&>( *ds );
  // NOTE: This is a little clunky but easier right now since I don't want to
  //       implement a copy constructor for every term.
  // TODO trap if term copy constructor invoked
  PotentialFunction* potential = pfxn.AddNewFunction();
  if (use_replicate)
    potential->AddTerm( PotentialTerm::REPLICATE, opts );
  else if (use_openmm)
    potential->AddTerm( PotentialTerm::OPENMM, opts );
  else {
    if (nbonds > 0)
      potential->AddTerm( PotentialTerm::BOND, opts );
    if (nangles > 0)
      potential->AddTerm( PotentialTerm::ANGLE, opts );
    if (ndihedrals > 0)
      potential->AddTerm( PotentialTerm::DIHEDRAL, opts );
    if (useNonbond) potential->AddTerm( PotentialTerm::SIMPLE_LJ_Q, opts );
  }

  mprintf("    CREATEPOTENTIAL: Created potential function '%s'\n", ds->legend());
  potential->FnInfo();
  return CpptrajState::OK;
}
