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
    void SetDebug(int);
    int InitRingFinder(ArgList&);
    int SetupRingFinder(Topology const&, AtomMask const&);

    unsigned int Nrings() const { return rings_.size(); }
    void PrintRings(Topology const&) const;
    AtomMask const& operator[](int idx) const { return rings_[idx]; }
  private:
    typedef std::vector<AtomMask> Marray;

    Marray rings_; ///< Array of masks corresponding to rings
    int debug_;
};
}
}
#endif
