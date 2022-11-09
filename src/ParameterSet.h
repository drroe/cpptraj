#ifndef INC_PARAMETERSET_H
#define INC_PARAMETERSET_H
#include "Param/AtomType.h"
#include "Param/ParmHolder.h"
#include "Param/DihedralParmHolder.h"
#include "Param/HookesLawType.h"
using namespace Cpptraj::Param;
/// Hold a set of parameters for atom types, bonds, angles, etc.
class ParameterSet {
  public:
    ParameterSet() : hasLJparams_(false) {}

    ParmHolder<AtomType>& AT()         { return atomTypes_; }
    ParmHolder<NonbondType>& NB()      { return nbParm_;    }
    ParmHolder<HookesLawType>& BP()     { return bondParm_; }
    ParmHolder<HookesLawType>& AP()    { return angleParm_; }
    ParmHolder<HookesLawType>& UB()     { return ubParm_; }
    //ParmHolder<DihedralParmType>& DP() { return dihParm_; }
    ParmHolder<DihedralParmType>& IP() { return impParm_; }
    DihedralParmHolder& DP()           { return dihParm_; }

    void SetHasLJparams(bool b) { hasLJparams_ = b; }
    bool HasLJparams() const { return hasLJparams_; }

    ParmHolder<AtomType> const& AT()         const { return atomTypes_; }
    ParmHolder<NonbondType> const& NB()      const { return nbParm_;    }
    ParmHolder<HookesLawType> const& BP()     const { return bondParm_; }
    ParmHolder<HookesLawType> const& AP()    const { return angleParm_; }
    ParmHolder<HookesLawType> const& UB()     const { return ubParm_; }
    //ParmHolder<DihedralParmType> const& DP() const { return dihParm_; }
    ParmHolder<DihedralParmType> const& IP() const { return impParm_; }
    DihedralParmHolder const& DP()           const { return dihParm_; }

    void Debug(const char*) const;
    void Debug() const { return Debug(""); }

    /// Used to track what parameters were updated during UpdateParams
    class UpdateCount {
      public:
        UpdateCount() : nBondsUpdated_(0), nAnglesUpdated_(0),
                        nDihedralsUpdated_(0), nImpropersUpdated_(0),
                        nUreyBradleyUpdated_(0), nAtomTypeUpdated_(0),
                        nLJparamsUpdated_(0) {}
        unsigned int nBondsUpdated_;
        unsigned int nAnglesUpdated_;
        unsigned int nDihedralsUpdated_;
        unsigned int nImpropersUpdated_;
        unsigned int nUreyBradleyUpdated_;
        unsigned int nAtomTypeUpdated_;
        unsigned int nLJparamsUpdated_;
    };
    /// Update this set with parameters from given set
    int UpdateParamSet(ParameterSet const&, UpdateCount&, int);
    /// \return Size in memory in bytes
    size_t DataSize() const;
  private:
    //AtomTypeArray atomTypes_;
    ParmHolder<AtomType> atomTypes_;
    ParmHolder<NonbondType> nbParm_;
    ParmHolder<HookesLawType> bondParm_;
    ParmHolder<HookesLawType> angleParm_;
    ParmHolder<HookesLawType> ubParm_;
    //ParmHolder<DihedralParmType> dihParm_;
    ParmHolder<DihedralParmType> impParm_;
    DihedralParmHolder dihParm_;
    bool hasLJparams_;
};
#endif
