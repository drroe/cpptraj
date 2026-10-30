#ifndef INC_STRUCTURE_MODXNA_INFO_H
#define INC_STRUCTURE_MODXNA_INFO_H
#include <string>
#include <vector>
//namespace Cpptraj {
//namespace Structure {
/// Hold info used by the modXNA FF for joining fragments
class ModXNA_Info {
  public:
    /// CONSTRUCTOR
    ModXNA_Info();

    enum StatType { OK = 0, NOT_MODXNA, ERR };

    enum MaskType { HEAD = 0, HEADSTRIP, TAIL, TAILSTRIP, ANCHOR, ANCHORSTRIP };

    typedef std::vector<std::string> Sarray;

    /// Parse the given modXNA title string
    StatType ParseModxnaStr(std::string const&);
    /// \return true if there is modXNA info
    bool HasModxna() const { return (!fragmentName_.empty()); }
    /// Print to stdout
    void PrintModxna() const;
    /// Generate ModXNA metadata string for e.g. writing ModXNA mol2 files
    std::string GenMetadataString() const;

    /// \return The atom name array for specified mask type.
    Sarray const& AtomNames(MaskType) const;

    /// \return The atom mask string for specified mask type.
    std::string AtomMaskString(MaskType) const;

    std::string const& FragmentName() const { return fragmentName_; }

  private:
    /// \return an array of atom names from the given modxna header
    static Sarray atomNameList(std::string const&, int&);

    static std::string atomMaskStr(Sarray const&);

    static inline void printmxna(const char*, Sarray const&);

    std::string fragmentName_;
    Sarray head_;
    Sarray headStrip_;
    Sarray tail_;
    Sarray tailStrip_;
    Sarray anchor_;
    Sarray anchorStrip_;
};
//}
//}
#endif
