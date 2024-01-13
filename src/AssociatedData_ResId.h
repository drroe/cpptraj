#ifndef INC_ASSOCIATEDDATA_RESID_H
#define INC_ASSOCIATEDDATA_RESID_H
#include "AssociatedData.h"
#include "NameType.h"
#include "Structure/StructureEnum.h"
class NameType;
/// Used to identify what residues a residue template matches
class AssociatedData_ResId : public AssociatedData {
  public:
    /// CONSTRUCTOR
    AssociatedData_ResId() : AssociatedData(RESID), termType_(Cpptraj::Structure::NON_TERMINAL) {}
    /// CONSTRUCTOR - Pdb name, terminal type
    AssociatedData_ResId(NameType const& n, Cpptraj::Structure::TerminalType t) :
      AssociatedData(RESID), resName_(n), termType_(t) {}
    // ----- Inherited functions -------
    static const char* HelpText;
    int ProcessAdataArgs(ArgList&);
    AssociatedData* Copy() const { return new AssociatedData_ResId(*this); }
    void Ainfo() const;
    // ---------------------------------
    NameType const& ResName() const { return resName_; }
    Cpptraj::Structure::TerminalType TermType() const { return termType_; }
  private:
    NameType resName_;                          ///< Target residue name
    Cpptraj::Structure::TerminalType termType_; ///< Target residue terminal type
};
#endif
