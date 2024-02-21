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
    Exec_Build() : Exec(GENERAL), debug_() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Build(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;
    // For holding bonded atom pairs.
    typedef std::pair<int,int> Ipair;
    typedef std::vector<Ipair> IParray;

    /// \return true if given IParray has the given Ipair
    static inline bool hasBondingPair(IParray const&, Ipair const&);
    /// Create new topology/frame using templates
    int FillAtomsWithTemplates(Topology&, Frame&, Topology const&, Frame const&, Cpptraj::Structure::Creator const&) const;
    /// Map atoms in topology to template
    static std::vector<int> MapAtomsToTemplate(Topology const&, int, DataSet_Coords*);

    int debug_;
};
#endif
