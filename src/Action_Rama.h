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

    enum Type { ALPHA=0, LEFT, PP2, HAIRPIN, EXTENDED, NONE, NTYPES };

    static const char* TypeKeys_[];

    typedef DihedralSearch::DihedralMask Dmask;
    class Res {
      public:
        Res() : data_(0), isActive_(false) {}
        Res(DataSet* d, Dmask const& m1, Dmask const& m2) :
          phi_(m1), psi_(m2), data_(d) {}
        bool IsActive() const { return isActive_; }
        DataSet* Data() const { return data_; }
        void SetActive(bool a) { isActive_ = a; }
        void SetData(DataSet* d) { data_ = d; }
        void SetMasks(Dmask const& m1, Dmask const& m2) { phi_=m1; psi_=m2; }
        Dmask const& Phi() const { return phi_; }
        Dmask const& Psi() const { return psi_; }
      private:
        Dmask phi_;
        Dmask psi_;
        DataSet* data_;
        bool isActive_;
    };

//    typedef std::map<int,Res> ResMapType;
//    ResMapType resMap_;
    typedef std::vector<Res> Rarray;
    Rarray residues_;

    std::vector<int> Sum_;
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
    double phiMin_[NTYPES];
    double phiMax_[NTYPES];
    double psiMin_[NTYPES];
    double psiMax_[NTYPES];
    int debug_;
    int Nframe_;
};
#endif
