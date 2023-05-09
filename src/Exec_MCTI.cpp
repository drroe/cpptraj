#include "Exec_MCTI.h"
#include "CpptrajStdio.h"
#include "Structure/Protonator.h"

// Exec_MCTI::Help()
void Exec_MCTI::Help() const
{
  Cpptraj::Structure::Protonator::HelpText();
}

// Exec_MCTI::Execute()
Exec::RetType Exec_MCTI::Execute(CpptrajState& State, ArgList& argIn)
{
  using namespace Cpptraj::Structure;
  Protonator protonator;
  if (protonator.SetupProtonator( State, argIn, State.Debug() )) {
    mprinterr("Error: Set up protonator failed.\n");
    return CpptrajState::ERR;
  }
  protonator.PrintOptions();
  if (protonator.CalcTitrationCurves()) {
    mprinterr("Error: Calculation of titration curves failed.\n");
    return CpptrajState::ERR;
  }
  return CpptrajState::OK;
}
