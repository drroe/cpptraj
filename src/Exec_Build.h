#ifndef INC_EXEC_BUILD_H
#define INC_EXEC_BUILD_H
#include "Exec.h"
namespace Cpptraj {
namespace Structure {
class Creator;
}
}
/// Used to build a structure 
class Exec_Build : public Exec {
  public:
    Exec_Build() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Build(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    /// Create new topology/frame using templates
    static int FillAtomsWithTemplates(Topology&, Frame&, Topology const&, Frame const&, Cpptraj::Structure::Creator const&);
    /// Map atoms in topology to template
    static std::vector<int> MapAtomsToTemplate(Topology const&, int, DataSet_Coords*);
};
#endif
