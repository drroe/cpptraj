#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold parameters for a unique atom type
class AtomType {
  public:
    AtomType() : mass_(0.0) {}
    /// Mass only
    AtomType(double m) : mass_(m) {}
    /// Radius, well depth, mass TODO take LJparmType instead
    AtomType(double r, double d, double m) : lj_(r, d), mass_(m) {}
    /// \return default LJ parameters
    LJparmType const& LJ() const { return lj_; }
    /// \return Atom mass in amu
    double Mass()          const { return mass_;   }
    /// \return true if LJ params are less than incoming
    bool operator<(AtomType const& rhs) const { return lj_ < rhs.lj_; }
    /// \return true if LJ params are the same
    bool operator==(AtomType const& rhs) const { return lj_ == rhs.lj_; }
    /// Used to modify LJ params
    LJparmType& SetLJ() { return lj_; }
  private:
    LJparmType lj_; ///< Default Lennard-Jones parameters (always valid for self).
    double mass_;   ///< Mass in amu
};
#endif
