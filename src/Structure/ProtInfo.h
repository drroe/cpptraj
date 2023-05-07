#ifndef INC_STRUCTURE_PROTINFO_H
#define INC_STRUCTURE_PROTINFO_H
#include <string>
#include <vector>
/** Class used to save residue name info for protonation/deprotonation. */
class ProtInfo {
  public:
    typedef std::vector<std::string> Sarray;
    enum PstateType { PROTONATED = 0, DEPROTONATED };

    ProtInfo(std::string const& pn, std::string const& dn, PstateType stateIn, Sarray const& ra) :
      defaultState_(stateIn), removeAtoms_(ra)
    {
      names_[0] = pn;
      names_[1] = dn;
    }
    /// \return Name of the residue when in the default state
    std::string const& DefaultName() const { return names_[(int)defaultState_]; }
    /// \return Protonated name
    std::string const& ProtName() const { return names_[(int)PROTONATED]; }
    /// \return Deprotonated name
    std::string const& DeprotName() const { return names_[(int)DEPROTONATED]; }
    /// \return array of atoms to remove
    Sarray const& RemoveAtoms() const { return removeAtoms_; }
  private:
    std::string names_[2];    ///< 0 is protonated name, 1 is deprotonated name
    PstateType defaultState_; ///< Which state is the default
    Sarray removeAtoms_;      ///< Optional list of atoms to remove from residue
};
#endif
