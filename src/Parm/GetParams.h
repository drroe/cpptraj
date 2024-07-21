#ifndef INC_PARM_GETPARAMS_H
#define INC_PARM_GETPARAMS_H
#include <vector>
class AngleArray;
class AngleParmArray;
class AngleParmType;
class Atom;
class AtomType;
class BondArray;
class BondParmArray;
class BondParmType;
class CmapArray;
class CmapGridArray;
class CmapParmHolder;
class DihedralArray;
class DihedralParmArray;
class DihedralParmHolder;
class HB_ParmType;
class ImproperParmHolder;
class NonbondParmType;
class NonbondType;
class ParameterSet;
class Residue;
class Topology;
template<typename Type> class ParmHolder;
namespace Cpptraj {
namespace Parm {
/// Used to get parameters from a Topology
class GetParams {
  public:
    /// CONSTRUCTOR
    GetParams();
    /// Set debug level
    void SetDebug(int);
    /// \return ParameterSet containing parameters from given topology.
    ParameterSet GetParameters(Topology const&) const;
    /// Get nonbonded parameters (used in AppendTop)
    void GetLJAtomTypes(ParmHolder<AtomType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<double>&,
                        ParmHolder<HB_ParmType>&,
                        std::vector<Atom> const&,
                        NonbondParmType const&) const;
  private:
    static inline void GetBondParams(ParmHolder<BondParmType>&, std::vector<Atom> const&,
                                     BondArray const&, BondParmArray const&);
    static inline void GetAngleParams(ParmHolder<AngleParmType>&, std::vector<Atom> const& atoms,
                                      AngleArray const&, AngleParmArray const&);
    static inline void GetImproperParams(ImproperParmHolder&, std::vector<Atom> const&,
                                         DihedralArray const&, DihedralParmArray const&);
    static inline void GetDihedralParams(DihedralParmHolder&, ImproperParmHolder&,
                                         std::vector<Atom> const&,
                                         DihedralArray const&, DihedralParmArray const&);
    static inline int GetCmapParams(CmapParmHolder&, CmapArray const&, CmapGridArray const&,
                                    std::vector<Atom> const&, std::vector<Residue> const&);

    int debug_;
};
}
}
#endif
