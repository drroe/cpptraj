#ifndef INC_EXEC_MEAD_H
#define INC_EXEC_MEAD_H
#include "Exec.h"
#include <string>
// Fwd declares
namespace Cpptraj {
namespace Mead {
class MeadGrid;
class MeadCalc;
}
}
/// Provide MEAD functionality 
class Exec_MEAD : public Exec {
  public:
    Exec_MEAD();
   ~Exec_MEAD();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_MEAD(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int CheckMead(Cpptraj::Mead::MeadGrid const&) const;
    int Solvate(CpptrajState&, ArgList&, Cpptraj::Mead::MeadGrid const&,
                std::string const&, DataFile*) const;
    int Potential(CpptrajState&, ArgList&, Cpptraj::Mead::MeadGrid const&,
                  std::string const&, DataFile*) const;
    int MultiFlex(CpptrajState&, ArgList&,
                  Cpptraj::Mead::MeadGrid const&, Cpptraj::Mead::MeadGrid const&,
                  Topology const&, Frame const&,
                  std::string const&, DataFile*) const;
    static int addGridLevel(Cpptraj::Mead::MeadGrid&, std::string const&);

    Cpptraj::Mead::MeadCalc* MEAD_;
};
#endif
