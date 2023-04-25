#ifndef INC_TITRATABLESITE_H
#define INC_TITRATABLESITE_H
#include <map>
#include <string>
#include "../NameType.h"
namespace Cpptraj {
namespace Structure {
/// Hold information for titratable atoms in a residue
class TitratableSite {
    typedef std::pair<double, double> ChargePair;
    /// Map atom name to state charges
    typedef std::map<NameType, ChargePair> MapType;
  public:
    TitratableSite();
    /// Clear all data
    void Clear();
    /// Load site data from file
    int LoadSiteFile(std::string const&, std::string const&);
    /// Print data to stdout
    void Print() const;

    typedef MapType::const_iterator const_iterator;
    const_iterator begin() const { return nameToCharges_.begin(); }
    const_iterator end()   const { return nameToCharges_.end(); }

    NameType const& SiteOfInterest() const { return siteOfInterest_; }
    int RefStateIdx() const { return refStateIdx_; }
    double pKa() const { return pKa_; }
    std::string const& SiteName() const { return siteName_; }
    std::string const& SiteResName() const { return resName_; }
  private:
    typedef std::pair<NameType, ChargePair> PairType;

    std::string resName_;     ///< Residue name that site belongs to
    std::string siteName_;    ///< Internal site name
    double pKa_;              ///< Site pka
    MapType nameToCharges_;   ///< Map atom names to state charge pairs
    int refStateIdx_;         ///< Reference state index
    NameType siteOfInterest_; ///< Name of the site of interest; typically first atom in site file
};
}
}
#endif
