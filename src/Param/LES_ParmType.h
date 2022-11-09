#ifndef INC_PARAM_LES_PARMTYPE_H
#define INC_PARAM_LES_PARMTYPE_H
#include <vector>
#include "LES_AtomType.h"
namespace Cpptraj {
namespace Param {
/// Hold LES parameters
class LES_ParmType {
  public:
    /// CONSTRUCTOR
    LES_ParmType() : ntypes_(0), ncopies_(0) {}
    /// Prepare LES_ParmType to receive data based on given # atoms and # LES types.
    void Allocate(unsigned int natomsIn, unsigned int ntypesIn) {
      ntypes_ = ntypesIn;
      ncopies_ = 0;
      array_.clear();
      array_.resize( natomsIn );
      fac_.clear();
      fac_.resize( (size_t)ntypes_ * (size_t)ntypes_ );
    }
    /// \return true if there are LES parameters
    bool HasLES()                const { return ntypes_ > 0;      }
    /// \return Number of LES types
    int Ntypes()                 const { return ntypes_;          }
    /// \return max number of copies of any atom
    int Ncopies()                const { return ncopies_;         }
    /// \return LES scaling factor array
    std::vector<double> const& FAC()    const { return fac_;             }
    /// \return Array with LES info for each atom
    LES_Array           const& Array()  const { return array_;           }
    /// Set number of LES types and the LES scaling factor array
    void SetTypes(int n, std::vector<double> const& f) {
      ntypes_ = n;
      fac_ = f;
    }
    /// Add an atom to the LES array. Update max copies
    void AddLES_Atom(LES_AtomType const& lat) {
      array_.push_back( lat );
      if (array_.back().Copy() > ncopies_ )
        ncopies_ = array_.back().Copy();
    }
    /// Set LES fac at given index in scaling factor array.
    void SetFAC(int idx, double f)  { fac_[idx] = f; }
    /// Set LES type of specified atom
    void SetType(int at, int type) { array_[at].SetType( type ); }
    /// Set copy # of specified atom. Update max copies.
    void SetCopy(int at, int cnum) {
      array_[at].SetCopy( cnum );
      if (cnum > ncopies_) ncopies_ = cnum;
    }
    /// Set given LES region ID of specified atom
    void SetID(int at, int id)     { array_[at].SetID( id ); }
    /// Clear all data
    void Clear() { ntypes_ = 0; ncopies_ = 0; array_.clear(); fac_.clear(); }
  private:
    int ntypes_;              ///< Total number of LES types 
    int ncopies_;             ///< Total number of LES copies.
    LES_Array array_;         ///< LES parameters for each atom
    std::vector<double> fac_; ///< Scaling factor for typeA * typeB
};
}
}
#endif
