#ifndef INC_GUESSATOMHYBRIDIZATION_H
#define INC_GUESSATOMHYBRIDIZATION_H
#include "AtomType.h" // for HybridizationType
class Atom;
class Topology;
namespace Cpptraj {
/// Guess atom hybridization type based on element/bonding
AtomType::HybridizationType GuessAtomHybridization(Atom const&, Topology const&);
}
#endif
