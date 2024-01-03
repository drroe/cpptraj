#ifndef INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#define INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#include <vector>
class BondArray;
class Residue;
class Atom;
namespace Cpptraj {
namespace Structure {
/// Generate bond array in same order as leap
BondArray GenerateBondArray(std::vector<Residue> const&, std::vector<Atom> const&);
}
}
#endif
