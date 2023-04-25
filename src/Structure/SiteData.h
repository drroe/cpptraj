#ifndef INC_TITRATIONDATA_H
#define INC_TITRATIONDATA_H
#include <vector>
#include <string>
#include <map>
namespace Cpptraj {
//TODO Move to Mead namespace
namespace Structure {
// Fwd declares
class TitratableSite;
/// Hold information for all titratable sites
class SiteData {
    typedef std::pair<int, std::string> IdxNamePair;
  public:
    typedef std::vector<std::string> Sarray;

    SiteData();
    /// Load titration data from files
    int LoadSiteData(std::string const&, std::string const&);

    /// \return Array of site names for given residue number
    Sarray ResSiteNames(int) const;
    /// \return Titratable site corresponding to name
    TitratableSite const& GetSite(std::string const&) const;

    typedef std::vector<IdxNamePair> IdxNameArray;
    typedef IdxNameArray::const_iterator const_iterator;
    /// \return Array of residue indices and corresponding site names
    IdxNameArray Sites() const { return IdxNames_; }
    /// \return begin iterator
    const_iterator begin() const { return IdxNames_.begin(); }
    /// \return end iterator
    const_iterator end() const { return IdxNames_.end(); }
  private:
    typedef std::pair<int,Sarray> ResNamePair;
    typedef std::map<int,Sarray> ResNameMap;
    typedef std::pair<std::string,TitratableSite> NameSitePair;
    typedef std::map<std::string,TitratableSite> NameSiteMap;

    ResNameMap ResToSitename_;  ///< Map residue numbers to site names
    NameSiteMap NameToSite_;    ///< Map site names to titratable site data
    IdxNameArray IdxNames_;  ///< Hold res #s/site names in original order.
};
}
}
#endif
