#ifndef INC_STRUCTURE_GENERATEIMPROPERS_H
#define INC_STRUCTURE_GENERATEIMPROPERS_H
#include <vector>
class Atom;
class Residue;
class AtomType;
class DihedralArray;
template <typename Type> class ParmHolder;
namespace Cpptraj {
namespace Structure {
/// Generate improper dihedral array in same order as leap
DihedralArray GenerateImproperArray(std::vector<Residue> const&, std::vector<Atom> const&, ParmHolder<AtomType> const&);
}
}
#endif
