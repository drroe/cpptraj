#ifndef INC_STRUCTURE_CHIRALITY_H
#define INC_STRUCTURE_CHIRALITY_H
#include "StructureEnum.h"
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
/// \return Chirality at specified atom, set torsion value
BuildAtom DetermineChirality(int, Topology const&, Frame const&, int);
//ChiralType DetermineChirality(double&, int*, int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom
//ChiralType DetermineChirality(int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom, set priority
//ChiralType SetPriority(std::vector<int>&, int, Topology const&, Frame const&, int);
/// Only update priority around atom. Useful when orientation around atom is not known/correct.
int UpdatePriority(BuildAtom&, int, Topology const&, Frame const&, int);
}
}
#endif
