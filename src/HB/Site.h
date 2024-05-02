#ifndef INC_HB_SITE_H
#define INC_HB_SITE_H
#include <vector>
namespace Cpptraj {
namespace HB {
/// Potential hydrogen bond site. Can be either donor or donor/acceptor.
class Site {
  public:
    typedef std::vector<int> Iarray;
    /// CONSTRUCTOR
    Site() : idx_(-1) {}
    /// Solute site - heavy atom, hydrogen atom
    Site(int d, int h) : hlist_(1,h), idx_(d) {}
    /// Solute site - heavy atom, list of hydrogen atoms
    Site(int d, Iarray const& H) : hlist_(H), idx_(d) {}
    /// \return heavy atom index
    int Idx() const { return idx_; }
    /// \return number of hydrogen indices
    unsigned int n_hydrogens()      const { return hlist_.size(); }
    /// \return true if site is an ion (D atom == H atom)
    bool IsIon() const { return (hlist_.size()==1 && hlist_[0] == idx_); }
    /// \return iterator to beginning of hydrogen indices
    Iarray::const_iterator Hbegin() const { return hlist_.begin(); }
    /// \return iterator to end of hydrogen indices
    Iarray::const_iterator Hend()   const { return hlist_.end(); }
  private:
    Iarray hlist_; ///< List of hydrogen indices
    int idx_;      ///< Heavy atom index
};
}
}
#endif
