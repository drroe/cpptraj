#ifndef INC_ACTION_RAMA_H
#define INC_ACTION_RAMA_H
#include "Action.h"
#include "DihedralSearch.h"
/// <Enter description of Action_Rama here>
class Action_Rama : public Action {
  public:
    Action_Rama() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Rama(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum Type { ALPHA=0, BETA, PP2, LALPHA, NTYPES };

    static const char* TypeKeys_[];

    DihedralSearch dihSearch_;
    CharMask Mask_;
    DataFile* outfile_;
    DataFile* sumFile_;
    std::string BB_N_;
    std::string BB_H_;
    std::string BB_C_;
    std::string BB_O_;
    std::string BB_CA_;
    std::string dsetname_;
    DataSet* ds_[NTYPES];
    double phiOff_[NTYPES];
    double psiOff_[NTYPES];
    int debug_;
    int Nframe_;
};
#endif
