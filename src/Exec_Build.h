#ifndef INC_EXEC_BUILD_H
#define INC_EXEC_BUILD_H
#include "Exec.h"
namespace Cpptraj {
namespace Structure {
class Creator;
class SugarBuilder;
}
}
/// Used to build a structure 
class Exec_Build : public Exec {
  public:
    Exec_Build();
    ~Exec_Build();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Build(); }
    RetType Execute(CpptrajState&, ArgList&);
    /// Stand-alone execution
    RetType BuildStructure(DataSet*, DataSetList&, int, ArgList&);
  private:
    typedef std::vector<int> Iarray;
    // Keep track of which residues are connected to each other
    typedef std::vector<Iarray> ResConnectArray;
    // For holding bonded atom pairs.
    typedef std::pair<int,int> Ipair;
    typedef std::vector<Ipair> IParray;

    /// \return true if given IParray has the given Ipair
    static inline bool hasBondingPair(IParray const&, Ipair const&);
    /// \return true if given array of residue connections has target residue index.
    static inline bool resIsConnected(Iarray const&, int);
    /// Create new topology/frame using templates
    int FillAtomsWithTemplates(Topology&, Frame&, Topology const&, Frame const&, Cpptraj::Structure::Creator const&) const;
    /// Map atoms in topology to template
    static std::vector<int> MapAtomsToTemplate(Topology const&, int, DataSet_Coords*, Cpptraj::Structure::Creator const&, std::vector<NameType>&);

    int debug_;
    int check_box_natom_;  ///< Systems larger than this will have box added so PL check can be used
    bool check_structure_; ///< If true check the resulting structure
    Cpptraj::Structure::SugarBuilder* sugarBuilder_;

    Timer t_total_;
    Timer t_hisDetect_;
    Timer t_clean_;
    Timer t_disulfide_;
    Timer t_sugar_;
    Timer t_fill_;
    Timer t_assign_;
    Timer t_check_;
};
#endif
