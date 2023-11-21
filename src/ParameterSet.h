#ifndef INC_PARAMETERSET_H
#define INC_PARAMETERSET_H
#include "ParameterTypes.h"
#include "ParameterHolders.h"
#include "AtomType.h"
class CpptrajFile;
/// Hold a set of parameters for atom types, bonds, angles, etc.
class ParameterSet {
  public:
    ParameterSet() : hasLJparams_(false) {}

    ParmHolder<AtomType>& AT()         { return atomTypes_; }
    ParmHolder<NonbondType>& NB()      { return nbParm_;    }
    ParmHolder<NonbondType>& NB14()    { return nb14Parm_;  }
    ParmHolder<BondParmType>& BP()     { return bondParm_;  }
    ParmHolder<AngleParmType>& AP()    { return angleParm_; }
    ParmHolder<BondParmType>& UB()     { return ubParm_;    }
    ParmHolder<DihedralParmType>& IP() { return impParm_;   }
    DihedralParmHolder& DP()           { return dihParm_;   }
    ParmHolder<HB_ParmType>& HB()      { return HBparm_;    }

    void SetHasLJparams(bool b) { hasLJparams_ = b; }
    bool HasLJparams() const { return hasLJparams_; }

    ParmHolder<AtomType> const& AT()         const { return atomTypes_; }
    ParmHolder<NonbondType> const& NB()      const { return nbParm_;    }
    ParmHolder<NonbondType> const& NB14()    const { return nb14Parm_;  }
    ParmHolder<BondParmType> const& BP()     const { return bondParm_;  }
    ParmHolder<AngleParmType> const& AP()    const { return angleParm_; }
    ParmHolder<BondParmType> const& UB()     const { return ubParm_;    }
    ParmHolder<DihedralParmType> const& IP() const { return impParm_;   }
    DihedralParmHolder const& DP()           const { return dihParm_;   }
    ParmHolder<HB_ParmType> const& HB()      const { return HBparm_;    }
    std::string const& NbParamName()         const { return NBname_;    }
    std::string ParamSetName()               const;
    /// Write parameters to file with given name
    void Debug(const char*) const;
    /// Write parameters to stdout
    void Debug() const { return Debug(""); }
    /// Print parameters to given file
    void Print(CpptrajFile&) const;

    /// Used to track what parameters were updated during UpdateParams
    class UpdateCount {
      public:
        UpdateCount() : nBondsUpdated_(0), nAnglesUpdated_(0),
                        nDihedralsUpdated_(0), nImpropersUpdated_(0),
                        nUreyBradleyUpdated_(0), nAtomTypeUpdated_(0),
                        nLJparamsUpdated_(0), nLJ14paramsUpdated_(0),
                        nHBparamsUpdated_(0) {}
        unsigned int nBondsUpdated_;
        unsigned int nAnglesUpdated_;
        unsigned int nDihedralsUpdated_;
        unsigned int nImpropersUpdated_;
        unsigned int nUreyBradleyUpdated_;
        unsigned int nAtomTypeUpdated_;
        unsigned int nLJparamsUpdated_;
        unsigned int nLJ14paramsUpdated_;
        unsigned int nHBparamsUpdated_;
    };
    /// Update this set with parameters from given set
    int UpdateParamSet(ParameterSet const&, UpdateCount&, int);
    /// Add hydrophilic atom type
    int AddHydrophilicAtomType(NameType const&);
    /// Set parameter set name
    void SetParamSetName(std::string const&);
    /// Set nonbond parameter set name
    void SetNbParamName(std::string const&);
    /// \return Size in memory in bytes
    size_t DataSize() const;
  private:
    typedef std::vector<NameType> NsetType;
    typedef std::vector<std::string> Sarray;

    Sarray name_;                          ///< Parameter set name(s)
    std::string NBname_;                   ///< Nonbond set name

    ParmHolder<AtomType> atomTypes_;       ///< Atom types
    ParmHolder<NonbondType> nbParm_;       ///< Lennard-Jones 6-12 A-B parameters
    ParmHolder<NonbondType> nb14Parm_;     ///< LJ 6-12 A-B parameters for 1-4 interactions
    ParmHolder<BondParmType> bondParm_;    ///< Hooke's law bond potential parameters
    ParmHolder<AngleParmType> angleParm_;  ///< Hooke's law angle potential parameters
    ParmHolder<BondParmType> ubParm_;      ///< Urey-Bradley parameters
    ParmHolder<DihedralParmType> impParm_; ///< Improper dihedral parameters
    DihedralParmHolder dihParm_;           ///< Cosine-series dihedral parameters
    ParmHolder<HB_ParmType> HBparm_;       ///< LJ 10-12 A-B parameters for hydrogen bonds
    NsetType hydrophilicAtomTypes_;        ///< Hold names of hydrophilic atom types
    bool hasLJparams_;
};
#endif
