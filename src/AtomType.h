#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold parameters for a unique atom type
class AtomType {
  public:
    AtomType() : mass_(0.0) {}
    /// Mass only
    AtomType(double m) : mass_(m) {}
    /// Radius, well depth, mass
    AtomType(double r, double d, double m) : lj_(r, d), mass_(m) {}
    /// \return default LJ parameters
    LJparmType const& LJ()   const { return lj_;   }
    /// \return Atom mass in amu
    double Mass()            const { return mass_; }
    /// \return True if this type has same parameters as incoming type
    bool operator==(AtomType const& rhs) const {
      return ( FEQ(mass_, rhs.mass_) &&
               (lj_ == rhs.lj_) );
    }
    /// \return True if this type is at all different from incoming type
    bool operator!=(AtomType const& rhs) const {
      return ( FNE(mass_, rhs.mass_) ||
               (lj_ != rhs.lj_) );
    }
    /// \return true if LJ params are less than incoming
    bool operator<(AtomType const& rhs) const {
      if (*this != rhs) {
        if (FEQ(mass_, rhs.mass_)) {
          return lj_ < rhs.lj_;
        } else {
          return mass_ < rhs.mass_;
        }
      } else
        return false;
    }
    /// Used to modify LJ params
    LJparmType& SetLJ() { return lj_; }
    /// Combine LJ params with this and another type using Lorentz-Berthelot rules
    NonbondType Combine_LB(AtomType const&) const;
    /// \return data size
    static size_t DataSize() { return (3*sizeof(double)) + sizeof(int); }
  private:
    LJparmType lj_; ///< Default Lennard-Jones parameters (always valid for self).
    double mass_;   ///< Mass in amu
};
#endif
