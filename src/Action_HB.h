#ifndef INC_ACTION_HB_H
#define INC_ACTION_HB_H
#include "Action.h"
#include "HB/HbCalc.h"
/// <Enter description of Action_HB here>
class Action_HB : public Action {
  public:
    Action_HB() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_HB(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Cpptraj::HB::HbCalc hbcalc_;
};
#endif
