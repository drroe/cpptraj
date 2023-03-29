#ifndef INC_EXEC_MEAD_H
#define INC_EXEC_MEAD_H
#include "Exec.h"
// Fwd declares
class DataSet_Vector_Scalar;
class DataSet_3D;
namespace Cpptraj {
class MeadInterface;
}
/// Provide MEAD functionality 
class Exec_MEAD : public Exec {
  public:
    Exec_MEAD() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_MEAD(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static int CheckMead(Cpptraj::MeadInterface const&);
    int Solvate(Cpptraj::MeadInterface&, ArgList&, DataSet*, DataSet_3D*) const;
    int Potential(Cpptraj::MeadInterface&, ArgList&, DataSet_Vector_Scalar&) const;
    int MultiFlex(Cpptraj::MeadInterface&, ArgList&) const;
};
#endif
