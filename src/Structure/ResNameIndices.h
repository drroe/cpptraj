#ifndef INC_STRUCTURE_RESNAMEINDICES_H
#define INC_STRUCTURE_RESNAMEINDICES_H
#include <vector>
#include "../NameType.h"
namespace Cpptraj {
namespace Structure {
class ResNameIndices {
    typedef std::vector<int> Iarray;
  public:
    ResNameIndices() {}
    /// CONSTRUCTOR - Residue name and index
    ResNameIndices(NameType const& n, int r) : resname_(n), resnums_(1, r) {}
    /// Add residue index
    void AddResnum(int r) { resnums_.push_back( r ); }
    /// Sort only by resname
    bool operator<(ResNameIndices const& rhs) const {
      return (resname_ < rhs.resname_);
    }
    /// \return Residue name
    NameType const& Resname() const { return resname_; }
    /// Iterator over residue numbers
    typedef Iarray::const_iterator const_iterator;
    /// \return Beginning of residue numbers
    const_iterator begin() const { return resnums_.begin(); }
    /// \return Ending of residue numbers
    const_iterator end() const { return resnums_.end(); }
  private:
    NameType resname_; ///< Residue name
    Iarray resnums_;   ///< Residue indices with this name
};
}
}
#endif
