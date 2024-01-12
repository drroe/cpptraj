#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold parameters for a unique atom type
class AtomType {
  public:
    /// Atom hybridization types
    enum HybridizationType { SP = 0, SP2, SP3, UNKNOWN_HYBRIDIZATION };
    /// CONSTRUCTOR
    AtomType() : mass_(0.0), polarizability_(0.0), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION) {}
    /// CONSTRUCTOR - Mass only
    AtomType(double m) : mass_(m), polarizability_(0.0), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION) {}
    /// CONSTRUCTOR - Mass, polarizability
    AtomType(double m, double p) : mass_(m), polarizability_(p), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION) {}
    /// CONSTRUCTOR - Radius, well depth, mass, polarizability
    AtomType(double r, double d, double m, double p) : lj_(r, d), mass_(m), polarizability_(p), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION) {}
    /// Set type index
    void SetTypeIdx(int i) { oidx_ = i; }
    /// \return default LJ parameters
    LJparmType const& LJ()  const { return lj_; }
    /// \return Atom mass in amu
    double Mass()           const { return mass_;   }
    /// \return Atomic polarizability in Ang^3
    double Polarizability() const { return polarizability_; }
    /// \return Original atom type index. Useful when checking for off-diagonal NB parameters.
    int OriginalIdx()       const { return oidx_; }
    /// \return Atom hybridization
    HybridizationType Hybridization() const { return hybrid_; }
    /// \return true if mass, polarizability, or LJ params are less than incoming
    bool operator<(AtomType const& rhs) const {
      if (FEQ(mass_, rhs.mass_)) {
        if (FEQ(polarizability_, rhs.polarizability_)) {
          return lj_ < rhs.lj_;
        } else {
          return polarizability_ < rhs.polarizability_;
        }
      } else {
        return mass_ < rhs.mass_;
      }
    }
    /// \return true if mass, polarizability, and LJ params are the same
    bool operator==(AtomType const& rhs) const {
      return (FEQ(mass_, rhs.mass_) &&
              FEQ(polarizability_, rhs.polarizability_) &&
              lj_ == rhs.lj_);
    }
    /// Used to modify LJ params
    LJparmType& SetLJ() { return lj_; }
    /// Set atom hybridization
    void SetHybridization(HybridizationType h) { hybrid_ = h; }
    /// \return data size  (2 double for LJparmType)
    static size_t DataSize() { return (4*sizeof(double)) + sizeof(int); }
  private:
    LJparmType lj_;         ///< Default Lennard-Jones parameters (always valid for self).
    double mass_;           ///< Mass in amu
    double polarizability_; ///< Atomic polarizability in Ang^3
    int oidx_;              ///< Original atom type index, for indexing nonbond parameters.
    HybridizationType hybrid_; ///< Atom hybridization
};
#endif
