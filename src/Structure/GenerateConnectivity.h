#ifndef INC_STRUCTURE_GENERATECONNECTIVITY_H
#define INC_STRUCTURE_GENERATECONNECTIVITY_H
class Topology;
namespace Cpptraj {
namespace Structure {
/// Generate angle/torsion arrays from bond arrays in a Topology
int GenerateAngleTorsionArrays(Topology&);
/// Generate bond/angle/torsion arrays from atom connectivity in a Topology
int GenerateBondAngleTorsionArrays(Topology&);
}
}
#endif
