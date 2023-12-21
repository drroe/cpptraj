#ifndef INC_STRUCTURE_GENERATEANGLES_H
#define INC_STRUCTURE_GENERATEANGLES_H
class Topology;
namespace Cpptraj {
namespace Structure {
/// Generate angles from bonds in a Topology
int GenerateAngles(Topology&);
/// Generate impropers
int GenerateImpropers(Topology&);
}
}
#endif
