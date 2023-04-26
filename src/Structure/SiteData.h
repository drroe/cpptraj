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
    typedef std::vector<IdxNamePair> IdxNameArray;
  public:
    /// CONSTRUCTOR
    SiteData();
    /// Load titration data and sites from a MEAD-style site file.
    int LoadMeadSiteData(std::string const&, std::string const&);
    /// Load titration data from specified directory.
    int LoadSiteDirectory(std::string const&);

    /// \return Titratable site corresponding to name
    TitratableSite const& GetSite(std::string const&) const;

    /// Iterator over defined sites
    typedef IdxNameArray::const_iterator const_iterator;
    /// \return Iterator to beginning of defined sites
    const_iterator begin() const { return IdxNames_.begin(); }
    /// \return Iterator to end of defined sites
    const_iterator end() const { return IdxNames_.end(); }

    /// \return true if No sites
    bool NoSites() const { return IdxNames_.empty(); }
  private:
    typedef std::pair<std::string,TitratableSite> NameSitePair;
    typedef std::map<std::string,TitratableSite> NameSiteMap;

    typedef std::vector<std::string> Sarray;
    typedef std::pair<std::string, Sarray> ResSitesPair;
    typedef std::map<std::string, Sarray> ResSitesMap;

    static std::string defaultSiteDir();

    NameSiteMap NameToSite_;    ///< Map site names to titratable site data
    IdxNameArray IdxNames_;  ///< Hold res #s/site names in original order.
    ResSitesMap resnameToSites_; ///< Map residue name to site names
};
}
}
#endif
