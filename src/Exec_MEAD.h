#ifndef INC_EXEC_MEAD_H
#define INC_EXEC_MEAD_H
#include "Exec.h"
// Fwd declares
class DataSet_Vector_Scalar;
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
    int Potential(Cpptraj::MeadInterface&, ArgList&, DataSet_Vector_Scalar&) const;
};
#endif
