#ifndef INC_ACTION_RAMA_H
#define INC_ACTION_RAMA_H
#include <map>
#include "Action.h"
#include "DihedralSearch.h"
/// <Enter description of Action_Rama here>
class Action_Rama : public Action {
  public:
    Action_Rama();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Rama(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum Type { ALPHA=0, LEFT, PP2, HAIRPIN, EXTENDED, NTYPES };

    static const char* TypeKeys_[];

    class Res {
      public:
        Res() : data_(0) {}
        Res(DataSet* d, DihedralSearch::DihedralMask const& m1,
                        DihedralSearch::DihedralMask const& m2) :
          phi_(m1), psi_(m2), data_(d) {}
      private:
        DihedralSearch::DihedralMask phi_;
        DihedralSearch::DihedralMask psi_;
        DataSet* data_;
    };

    typedef std::map<int,Res> ResMapType;
    ResMapType resMap_;

    DihedralSearch dihSearch_;
    CharMask Mask_;
    DataFile* outfile_;
    DataFile* sumFile_;
    DataSetList* masterDSL_;
    std::string BB_N_;
    std::string BB_H_;
    std::string BB_C_;
    std::string BB_O_;
    std::string BB_CA_;
    std::string dsetname_;
    DataSet* ds_[NTYPES];
    double Phi_[NTYPES];
    double Psi_[NTYPES];
    double phiOff_[NTYPES];
    double psiOff_[NTYPES];
    int debug_;
    int Nframe_;
};
#endif
