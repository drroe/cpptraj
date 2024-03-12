#ifndef INC_PARM_MERGE_H
#define INC_PARM_MERGE_H
#include <vector>
class Atom;
class BondArray;
class BondParmArray;
class AngleArray;
class AngleParmArray;
class DihedralArray;
class DihedralParmArray;
namespace Cpptraj {
namespace Parm {
/// Shorthand for array of Atoms
typedef std::vector<Atom> AtArray;
/// Merge two bond arrays and their parameters.
void MergeBondArrays(BondArray&, BondArray&, BondParmArray&, AtArray const&,
                     BondArray const&, BondArray const&, BondParmArray const&, AtArray const&);
/// Merge two angle arrays and their parameters.
void MergeAngleArrays(AngleArray&, AngleArray&, AngleParmArray&, AtArray const&,
                      AngleArray const&, AngleArray const&, AngleParmArray const&, AtArray const&);
/// Merge two dihedral arrays and their parameters.
void MergeDihedralArrays(DihedralArray&, DihedralArray&, DihedralParmArray&, AtArray const&,
                         DihedralArray const&, DihedralArray const&, DihedralParmArray const&, AtArray const&);

}
}
#endif
