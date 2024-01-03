#ifndef INC_STRUCTURE_GENERATEANGLES_H
#define INC_STRUCTURE_GENERATEANGLES_H
class Topology;
class AtomType;
template <typename Type> class ParmHolder;
namespace Cpptraj {
namespace Structure {
/// Generate angles from bonds in a Topology
int GenerateAngles(Topology&);
/// Generate impropers
int GenerateImpropers(Topology&, ParmHolder<AtomType> const&);
}
}
#endif
