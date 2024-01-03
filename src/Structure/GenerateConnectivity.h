#ifndef INC_STRUCTURE_GENERATECONNECTIVITY_H
#define INC_STRUCTURE_GENERATECONNECTIVITY_H
class Topology;
class AtomType;
template <typename Type> class ParmHolder;
namespace Cpptraj {
namespace Structure {
/// Generate bond/angle/torsion arrays from atom connectivity in a Topology
int GenerateBondAngleTorsionArrays(Topology&);
/// Generate impropers in a Topology
int GenerateImpropers(Topology&, ParmHolder<AtomType> const&);
}
}
#endif
