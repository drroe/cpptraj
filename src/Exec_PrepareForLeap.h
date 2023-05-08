#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
class DataSet_Coords_CRD;
namespace Cpptraj {
class LeapInterface;
namespace Structure {
class SugarBuilder;
}
namespace Mead {
class MeadCalc_Multiflex;
}
}
/// Do common tasks to prepare a structure to be loaded into tleap 
class Exec_PrepareForLeap : public Exec {
  public:
    Exec_PrepareForLeap();
    ~Exec_PrepareForLeap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrepareForLeap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;

    /// Set PDB residue names recognized by Amber FFs
    void SetPdbResNames();
    /// Load PDB residue names recognized by Amber FFs from dat file
    int LoadPdbResNames(std::string const&);
    /// \return true if residue name is a recognized PDB name
    bool IsRecognizedPdbRes(NameType const&, Cpptraj::Structure::SugarBuilder const&) const;
    /// \return Array of residue nums with unrecognized names
    Iarray GetUnrecognizedPdbResidues(Topology const&, Cpptraj::Structure::SugarBuilder const&) const;
    /// \return Array indices of isolated unrecognized residues
    Iarray GetIsolatedUnrecognizedResidues(Topology const&, Iarray const&) const;

    /// Try to determine where TER cards should be placed based on bonds
    int FindTerByBonds(Topology&, CharMask const&) const;

    /// Remove specified atoms
    int ModifyCoords(Topology&, Frame&, bool, std::string const&, std::string const&,
                     std::string const&, Iarray const&) const;
    /// Remove hydrogen atoms
    int RemoveHydrogens(Topology&, Frame&) const;

    /// Run leap to generate topology, perform any modifications
    int RunLeap(CpptrajState&, std::string const&, std::string const&) const;
    /// Print a warning for residues that will need modification after leap
    static void LeapFxnGroupWarning(Topology const&, int);

    /// Do the protonation state calculation
    int ProtonationStateCalc(Topology&, CpptrajState&, Frame const&, Cpptraj::LeapInterface const&,
                             std::string const&, std::string const&) const;

    /// Read topology/coords written from leap
    int readLeapTop(Topology&, Frame&, std::string const&, std::string const&) const;

    // -----------------------
    typedef std::set<NameType> SetType;
    SetType pdb_res_names_; ///< PDB residue names recognized by Amber FFs

    std::string leapunitname_;
    bool errorsAreFatal_;        ///< If false, try to skip errors.
    bool doProtonationState_;    ///< If true, try to assign protonation states to titratable groups.
    int debug_;                  ///< Debug level
    double target_pH_;           ///< Target pH if assigning protonation states (doProtonationState_).
    std::string solventResName_; ///< Solvent residue name
    Cpptraj::Mead::MeadCalc_Multiflex* multiflex_; ///< For doing the protonation state calc.
    DataSet_Coords_CRD* outCoords_; ///< Hold output COORDS set
    Trajout_Single* PDB_;           ///< For writing output PDB file
};
#endif
