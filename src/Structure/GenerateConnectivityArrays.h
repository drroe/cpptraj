#ifndef INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#define INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#include <vector>
class BondArray;
class AngleArray;
class DihedralArray;
class Residue;
class Atom;
namespace Cpptraj {
namespace Structure {
/// Generate bond array in same order as leap
BondArray GenerateBondArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate angle array in same order as leap
AngleArray GenerateAngleArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate dihedral array in same order as leap
DihedralArray GenerateDihedralArray(std::vector<Residue> const&, std::vector<Atom> const&);
}
}
#endif
