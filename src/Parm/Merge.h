#ifndef INC_PARM_MERGE_H
#define INC_PARM_MERGE_H
#include <vector>
class Atom;
class BondArray;
class BondParmArray;
namespace Cpptraj {
namespace Parm {
/// Shorthand for array of Atoms
typedef std::vector<Atom> AtArray;
/// Merge two bond arrays and their parameters.
void MergeBondArrays(BondArray&, BondArray&, BondParmArray&, AtArray const&,
                     BondArray const&, BondArray const&, BondParmArray const&, AtArray const&);
}
}
#endif
