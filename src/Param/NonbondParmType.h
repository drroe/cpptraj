#ifndef INC_PARAM_NONBONDPARMTYPE_H
#define INC_PARAM_NONBONDPARMTYPE_H
#include <vector>
#include "NonbondType.h"
#include "HB_ParmType.h"
namespace Cpptraj {
namespace Param {
/// Hold nonbonded interaction parameters
/** The nbindex array holds indices into nbarray (>=0) or hbarray (<0).
  * nbarray size should be (ntypes*(ntypes+1))/2 (half matrix).
  */
class NonbondParmType {
  public:
    /// CONSTRUCTOR
    NonbondParmType() : ntypes_(0) {}
    /// \return True if LJ params are present
    bool HasNonbond()                    const { return ntypes_ > 0; }
    /// \return Number of nonbonded types
    int Ntypes()                         const { return ntypes_;     }
    /// \return True if LJ C coefficients are present
    bool Has_C_Coeff()                   const { return !ccoef_.empty(); }
    /// \return The internal nonbonded index array
    std::vector<int> const& NBindex()    const { return nbindex_;    }
    /// \return Array of LJ 12-6 A and B parameters
    NonbondArray     const& NBarray()    const { return nbarray_;    }
    /// \return Array of LJ 10-12 (hbond) parameters
    HB_ParmArray     const& HBarray()    const { return hbarray_;    }
    /// \return Array of LJ 12-6-4 C parameters
    std::vector<double> const& LJC_Array() const { return ccoef_; }
    /// \return LJ 6-12 A and B parameter at specified index
    NonbondType const& NBarray(int i)    const { return nbarray_[i]; }
    /// \return LJ 10-12 (hbond) parameter at specified index
    HB_ParmType const& HBarray(int i)    const { return hbarray_[i]; }
    /// \return index into nonbonded array. In Amber, index < 0 means HB, otherwise LJ 6-12
    int GetLJindex(int type1, int type2) const {
      return nbindex_[ ntypes_ * type1 + type2 ];
    }

    /// Set number of types and init nonbond index array.
    void SetNtypes(unsigned int n) {
      ntypes_ = n;
      nbindex_.assign((size_t)ntypes_ * (size_t)ntypes_, -1); 
    }
    /// Set number of types, init NB index array, init LJ array.
    void SetupLJforNtypes(int n) { SetNtypes(n); nbarray_.assign((n*(n+1))/2, NonbondType()); }
    /// Set specified LJ term
    NonbondType& SetLJ(int i) { return nbarray_[i];                  }
    /// Set number of HB terms and init HB array TODO combine with SetNtypes?
    void SetNHBterms(int n)   { hbarray_.assign( n, HB_ParmType() ); }
    /// Set specified HB term
    HB_ParmType& SetHB(int i) { return hbarray_[i];                  }
    /// Add a LJ C parameter
    void AddLJC(double c) { ccoef_.push_back( c ); }
    /// Set specified nbindex location to given value.
    void SetNbIdx(int idx, int nbidx) { nbindex_[idx] = nbidx; }
    /// Add given LJ term to nonbond array and update nonbond index array.
    /** Certain routines in sander (like the 1-4 calcs) do NOT use the 
      * nonbond index array; instead they expect the nonbond arrays to be
      * indexed like '(ibig*(ibig-1)/2+isml)', where ibig is the larger atom
      * type index.
      */
    void AddLJterm(int type1, int type2, NonbondType const& LJ) {
      int ibig, isml;
      if (type1 > type2) {
        ibig = type1 + 1;
        isml = type2 + 1;
      } else {
        ibig = type2 + 1;
        isml = type1 + 1;
      }
      int ndx = (ibig*(ibig-1)/2+isml)-1;
      nbindex_[ntypes_ * type1 + type2] = ndx;
      nbindex_[ntypes_ * type2 + type1] = ndx;
      if (ndx >= (int)nbarray_.size())
        nbarray_.resize(ndx+1);
      nbarray_[ndx] = LJ;
    }
    /// Add given HB term to HB array and update the nonbond index array.
    void AddHBterm(int type1, int type2, HB_ParmType const& HB) {
      int ndx = -((int)hbarray_.size())-1;
      nbindex_[ntypes_ * type1 + type2] = ndx;
      nbindex_[ntypes_ * type2 + type1] = ndx;
      hbarray_.push_back( HB );
    }
    /// Clear all parameters
    void Clear() { ntypes_ = 0; nbindex_.clear(); nbarray_.clear(); hbarray_.clear(); }
  private:
    int ntypes_;               ///< Number of unique atom types
    std::vector<int> nbindex_; ///< Hold indices into arrays nbarray/hbarray for atom type pairs
    NonbondArray nbarray_;     ///< Hold Lennard-Jones 6-12 A and B parameters for all pairs.
    HB_ParmArray hbarray_;     ///< Hold 10-12 Amber HBond params for all pairs.
    std::vector<double> ccoef_; ///< Hold Lennard-Jones C parameters for 12-6-4 LJ potential.
};
}
}
#endif
