#ifndef INC_STRUCTURE_RINGFINDER_H
#define INC_STRUCTURE_RINGFINDER_H
#include <vector>
class ArgList;
class AtomMask;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Used to look for rings in residues
class RingFinder {
  public:
    RingFinder();
    int InitRingFinder(ArgList&);
    int SetupRingFinder(Topology const&);
  private:
    typedef std::vector<AtomMask> Marray;

    Marray rings_; ///< Array of masks corresponding to rings
};
}
}
#endif
