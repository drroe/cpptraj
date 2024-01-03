#ifndef INC_STRUCTURE_GENERATECONNECTIVITY_H
#define INC_STRUCTURE_GENERATECONNECTIVITY_H
class Topology;
class AtomType;
class BondArray;
template <typename Type> class ParmHolder;
namespace Cpptraj {
namespace Structure {
/// Generate angle/torsion arrays from bond arrays in a Topology
int GenerateAngleTorsionArrays(Topology&);
/// Generate bond/angle/torsion arrays from atom connectivity in a Topology
int GenerateBondAngleTorsionArrays(Topology&);
/// Generate impropers in a Topology
int GenerateImpropers(Topology&, ParmHolder<AtomType> const&);
/// Generate bond array in same order as leap
BondArray GenerateBondArray(Topology const&);
}
}
#endif
