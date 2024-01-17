#ifndef INC_EXEC_BUILD_H
#define INC_EXEC_BUILD_H
#include "Exec.h"
#include "Structure/StructureEnum.h"
/// Used to build a structure 
class Exec_Build : public Exec {
  public:
    Exec_Build() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Build(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<DataSet_Coords*> Carray;

    /// Identify residue template from residue name
    static DataSet_Coords* IdTemplateFromName(Carray const&, NameType const&, Cpptraj::Structure::TerminalType);
    /// Create new topology/frame using templates
    static int FillAtomsWithTemplates(Topology&, Frame&, Carray const&, Topology const&, Frame const&, ParameterSet const&);
    /// Map atoms in topology to template
    static std::vector<int> MapAtomsToTemplate(Topology const&, int, DataSet_Coords*);
};
#endif
