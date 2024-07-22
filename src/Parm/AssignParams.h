#ifndef INC_PARM_ASSIGNPARAMS_H
#define INC_PARM_ASSIGNPARAMS_H
#include <vector>
#include <string>
class AngleArray;
class AngleParmArray;
class AngleParmType;
class AngleType;
class Atom;
class AtomType;
class BondArray;
class BondParmArray;
class BondParmType;
class BondType;
class CmapArray;
class CmapGridArray;
class CmapParmHolder;
class DihedralArray;
class DihedralParmArray;
class DihedralParmHolder;
class DihedralType;
class HB_ParmType;
class ImproperParmHolder;
class NonbondType;
class ParameterSet;
class Topology;
template<typename Type> class ParmHolder;
namespace Cpptraj {
namespace Parm {
/// Used to assign parameters to a Topology
class AssignParams {
  public:
    AssignParams();
    /// Set Debug level
    void SetDebug(int);
    /// Set verbosity
    void SetVerbose(int);
    /// Replace existing parameters with those from given set
    int AssignParameters(Topology&, ParameterSet const&) const;
    /// Update existing parameters with given parameter set
    int UpdateParameters(Topology&, ParameterSet const&) const;
    // Assign nonbond parameters to topology. Used during Merge::AppendTop()
    void AssignNonbondParams(Topology&,
                             ParmHolder<AtomType> const&,
                             ParmHolder<NonbondType> const&, ParmHolder<NonbondType> const&,
                             ParmHolder<double> const&, ParmHolder<HB_ParmType> const&) const;
  private:
    typedef std::vector<Atom> AtArray;

    void AssignAtomTypeParm(AtArray&, ParmHolder<AtomType> const&) const;
    void AssignBondParm(Topology const&, ParmHolder<BondParmType> const&,
                        BondArray&, BondParmArray&, const char*) const;
    void AssignBondParams(Topology&, ParmHolder<BondParmType> const&) const;
    void AssignUBParams(Topology&, ParmHolder<BondParmType> const&) const;

    AngleArray AssignAngleParm(Topology const&, ParmHolder<AngleParmType> const&,
                               AngleArray const&, AngleParmArray&) const;
    void AssignAngleParams(Topology&, ParmHolder<AngleParmType> const&) const;

    void warn_improper_reorder(AtArray const&, DihedralType const&, DihedralType const&) const;
    void AssignImproperParm(Topology&, ImproperParmHolder const&,
                            DihedralArray&, DihedralParmArray&) const ;
    void AssignImproperParams(Topology&, ImproperParmHolder const&) const;

    DihedralArray AssignDihedralParm(Topology&,
                                     DihedralParmHolder const&, ImproperParmHolder const&,
                                     ParmHolder<AtomType> const&, DihedralArray const&,
                                     bool) const; // TODO make the bool a class var?
    DihedralArray get_unique_dihedrals(DihedralArray const&) const;
    void AssignDihedralParams(Topology&, DihedralParmHolder const&, ImproperParmHolder const&,
                              ParmHolder<AtomType> const&) const;



    static inline int cmap_anames_match(DihedralType const&, AtArray const&, std::vector<std::string> const&);
    int remap_cmap_indices(Topology const&, std::vector<int>&, CmapGridArray&, CmapArray&, CmapParmHolder const&) const;
    int AssignCmapParams(Topology const&, CmapArray&, CmapParmHolder const&, CmapGridArray&) const;
    int AssignCmapParams(Topology const&, DihedralArray const&, CmapParmHolder const&,
                         CmapGridArray&, CmapArray&) const;
    static void AddToBondArrays(Topology&, BondType const&);
    static void AddToAngleArrays(Topology&, AngleType const&);
    static void AddToDihedralArrays(Topology&, DihedralType const&);

    int debug_;
    int verbose_;
    bool deleteExtraPointAngles_; ///< If true, remove angles/torsions containing extra points.
    bool flexibleWater_;          ///< If true, allow H-O-H angle for water
};
}
}
#endif
