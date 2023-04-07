#ifndef INC_EXEC_MEAD_H
#define INC_EXEC_MEAD_H
#include "Exec.h"
#include <string>
// Fwd declares
class DataSet_Vector_Scalar;
class DataSet_3D;
namespace Cpptraj {
class MeadInterface;
class MultiFlexResults;
namespace Mead {
class MeadGrid;
}
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
    int Solvate(Cpptraj::MeadInterface&, Cpptraj::Mead::MeadGrid const&, ArgList&, DataSet*, DataSet_3D*) const;
    int Potential(Cpptraj::MeadInterface&, Cpptraj::Mead::MeadGrid const&, ArgList&, DataSet_Vector_Scalar&) const;
    int MultiFlex(Cpptraj::MeadInterface&, ArgList&, Topology const&, Frame const&, int, Cpptraj::MultiFlexResults const&) const;
    static int addGridLevel(Cpptraj::Mead::MeadGrid&, std::string const&);
};
#endif
