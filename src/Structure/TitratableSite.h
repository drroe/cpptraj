#ifndef INC_TITRATABLESITE_H
#define INC_TITRATABLESITE_H
#include <map>
#include <string>
class NameType;
namespace Cpptraj {
namespace Structure {
/// Hold information for titratable atoms in a residue
class TitratableSite {
  public:
    TitratableSite();
    /// Clear all data
    void Clear();
    /// Load site data from file
    int LoadSiteData(std::string const&);
  private:
    typedef std::pair<double, double> ChargePair;
    typedef std::pair<NameType, ChargePair> PairType;
    /// Map atom name to state charges
    typedef std::map<NameType, ChargePair> MapType;

    std::string resName_; ///< Site residue name
    double pKa_; ///< Site pka
    MapType nameToCharges_; ///< Map atom names to state charge pairs
};
}
}
#endif
