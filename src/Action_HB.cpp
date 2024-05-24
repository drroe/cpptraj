#include "Action_HB.h"
#include "CpptrajStdio.h"

// Action_HB::Help()
void Action_HB::Help() const {

}

// Action_HB::Init()
Action::RetType Action_HB::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  hbcalc_.SetDebug( debugIn );
# ifdef MPI
  hbcalc_.SetTrajComm( init.TrajComm() );
# endif
  if (hbcalc_.InitHbCalc( actionArgs, init.DslPtr(), init.DFL(), debugIn )) {
    mprinterr("Error: Could not initialize HB calc.\n");
    return Action::ERR;
  }

  hbcalc_.PrintHbCalcOpts();

  return Action::OK;
}

// Action_HB::Setup()
Action::RetType Action_HB::Setup(ActionSetup& setup)
{
  if (hbcalc_.SetupHbCalc( setup.Top(), setup.CoordInfo().TrajBox() )) {
    mprinterr("Error: Could not setup HB calc.\n");
    return Action::ERR;
  }

  return Action::OK;
}

// Action_HB::DoAction()
Action::RetType Action_HB::DoAction(int frameNum, ActionFrame& frm)
{
  if (hbcalc_.RunCalc( frm.Frm(), frameNum, frm.TrajoutNum() ))
    return Action::ERR;
  return Action::OK;
}

void Action_HB::Print() {
  hbcalc_.FinishHbCalc();
}

#ifdef MPI
int Action_HB::SyncAction() {
  return hbcalc_.SyncToMaster();
}
#endif
