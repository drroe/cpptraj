#ifndef INC_PARM_ASSIGNPARAMS_H
#define INC_PARM_ASSIGNPARAMS_H
#include <vector>
class AngleArray;
class AngleParmArray;
class AngleParmType;
class Atom;
class AtomType;
class BondParmType;
class BondArray;
class BondParmArray;
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

    int AssignParameters(Topology&, ParameterSet const&) const;
    int UpdateParameters(Topology&, ParameterSet const&) const;
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

    DihedralArray AssignDihedralParm(DihedralParmHolder const&, ImproperParmHolder const&,
                                     ParmHolder<AtomType> const&, DihedralArray const&,
                                     bool) const; // TODO make the bool a class var?
    DihedralArray get_unique_dihedrals(DihedralArray const&) const;
    void AssignDihedralParams(DihedralParmHolder const&, ImproperParmHolder const&,
                              ParmHolder<AtomType> const&) const;

    void AssignNonbondParams(ParmHolder<AtomType> const&,
                             ParmHolder<NonbondType> const&, ParmHolder<NonbondType> const&,
                             ParmHolder<double> const&, ParmHolder<HB_ParmType> const&,
                             int) const; // TODO make int a class vair

    int remap_cmap_indices(std::vector<int>&, CmapGridArray&, CmapArray&, CmapParmHolder const&) const;
    int AssignCmapParams(CmapArray&, CmapParmHolder const&, CmapGridArray&) const;
    int AssignCmapParams(DihedralArray const&, CmapParmHolder const&,
                         CmapGridArray&, CmapArray&) const;



    int debug_;
};
}
}
#endif
