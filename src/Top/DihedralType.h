#ifndef INC_TOP_DIHEDRALTYPE_H
#define INC_TOP_DIHEDRALTYPE_H
#include <vector>
#include "NbodyTerm.h"
namespace Cpptraj {
namespace Top {
/// Hold dihedral atom indices and parameter index
/** Dihedrals can be marked normal (A1-A2-A3-A4), end (meaning 1-4 calc should
  * be skipped to avoid overcounting, e.g. for dihedrals with multiple 
  * multiplicities or certain ring dihedrals), improper, or both end and improper.
  */
class DihedralType : public NbodyTerm {
  public:
    enum Dtype { NORMAL=0, IMPROPER, END, BOTH };
    /// Set skip 1-4 (end) and improper status
    void SetFromType(Dtype t) {
      switch (t) {
        case NORMAL   : skip14_ = false; improper_ = false; break;
        case IMPROPER : skip14_ = false; improper_ = true; break;
        case END      : skip14_ = true;  improper_ = false; break;
        case BOTH     : skip14_ = true;  improper_ = true; break;
      }
    }
    /// Default constructor
    DihedralType() : skip14_(false), improper_(false) {}
    /// For use with Amber-style dihedral array; a3_ < 0 = E, a4_ < 0 = I
    DihedralType(int a1, int a2, int a3, int a4, int idx) : NbodyTerm(idx)
    {
      at_[0] = a1;
      at_[1] = a2;
      at_[2] = a3;
      at_[3] = a4;
      if (at_[2] < 0 && at_[3] < 0) { at_[2] = -at_[2]; at_[3] = -at_[3]; skip14_ = true;  improper_ = true;  }
      else if (at_[2] < 0)          { at_[2] = -a3;                       skip14_ = true;  improper_ = false; }
      else if (at_[3] < 0)          { at_[3] = -a4;                       skip14_ = false; improper_ = true;  }
      else                          {                                     skip14_ = false; improper_ = false; }
    }
    /// Takes type, no index
    DihedralType(int a1, int a2, int a3, int a4, Dtype t)
    {
      at_[0] = a1;
      at_[1] = a2;
      at_[2] = a3;
      at_[3] = a4;
      SetFromType(t);
    }
    /// Takes type and index
    DihedralType(int a1, int a2, int a3, int a4, Dtype t, int idx) : NbodyTerm(idx)
    {
      at_[0] = a1;
      at_[1] = a2;
      at_[2] = a3;
      at_[3] = a4;
      SetFromType(t);
    }

    int A1()     const { return at_[0];    }
    int A2()     const { return at_[1];    }
    int A3()     const { return at_[2];    }
    int A4()     const { return at_[3];    }

    int At(int idx) const { return at_[idx]; }
    int Nat() const { return 4; }

    void SetSkip14(bool b)    { skip14_ = b;   }
    void SetImproper(bool b)  { improper_ = b; }
    /// \return type based on skip 1-4 (end) and improper status
    Dtype Type() const {
      if (skip14_ && improper_) return BOTH;
      else if (skip14_)         return END;
      else if (improper_)       return IMPROPER;
      else                      return NORMAL;
    }
    bool Skip14()     const { return skip14_; }
    bool IsImproper() const { return improper_; }
  private:
    int at_[4];
    bool skip14_;   ///< If true the 1-4 interaction for this dihedral should be skipped.
    bool improper_; ///< If true this is an improper dihedral.
};
/// Array of dihedrals
typedef std::vector<DihedralType> DihedralArray;
}
}
#endif
