#ifndef INC_PARM_MERGE_H
#define INC_PARM_MERGE_H
#include <vector>
class Atom;
class Residue;
class BondArray;
class BondParmArray;
class AngleArray;
class AngleParmArray;
class DihedralArray;
class DihedralParmArray;
class CmapArray;
class CmapGridArray;
class Topology;
namespace Cpptraj {
namespace Parm {
/// Used to merge Topology files.
class Merge {
  public:
    /// CONSTRUCTOR
    Merge();
    /// Append topology to this one.
    int AppendTop( Topology&, Topology const& ) const;
  private:
    /// Shorthand for array of Atoms
    typedef std::vector<Atom> AtArray;
    /// Shorthand for array of residues
    typedef std::vector<Residue> ResArray;
    /// Merge two pairs of bond arrays and their parameters.
    void MergeBondArrays(BondArray&, BondArray&, BondParmArray&, AtArray const&,
                         BondArray const&, BondArray const&, BondParmArray const&, AtArray const&) const;
    /// Merge two pairs of angle arrays and their parameters.
    void MergeAngleArrays(AngleArray&, AngleArray&, AngleParmArray&, AtArray const&,
                          AngleArray const&, AngleArray const&, AngleParmArray const&, AtArray const&) const;
    /// Merge two pairs of dihedral arrays and their parameters.
    void MergeDihedralArrays(DihedralArray&, DihedralArray&, DihedralParmArray&, AtArray const&,
                             DihedralArray const&, DihedralArray const&, DihedralParmArray const&, AtArray const&) const;
    /// Merge cmap arrays and their parameters
    void MergeCmapArrays(CmapArray&, CmapGridArray&, AtArray const&, ResArray const&,
                         CmapArray const&, CmapGridArray const&, AtArray const&, ResArray const&) const;
    /// Merge two bond arrays and their parameters
    void MergeBondArray(BondArray&, BondParmArray&, AtArray const&,
                        BondArray const&, BondParmArray const&, AtArray const&) const;
    /// Merge two improper arrays and their parameters
    void MergeImproperArray(DihedralArray&, DihedralParmArray&, AtArray const&,
                            DihedralArray const&, DihedralParmArray const&, AtArray const&) const;

    int debug_;
    int verbose_;
};
}
}
#endif
