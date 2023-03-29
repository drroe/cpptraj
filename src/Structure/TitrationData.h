#ifndef INC_TITRATIONDATA_H
#define INC_TITRATIONDATA_H
#include <vector>
#include <string>
#include <map>
namespace Cpptraj {
namespace Structure {
// Fwd declares
class TitratableSite;
/// Hold information for all titratable sites
class TitrationData {
  public:
    TitrationData();
    /// Load titration data from files
    int LoadTitrationData(std::string const&, std::string const&);
  private:
    typedef std::pair<int,std::string> ResNamePair;
    typedef std::map<int,std::string> ResNameMap;
    typedef std::pair<std::string,TitratableSite> NameSitePair;
    typedef std::map<std::string,TitratableSite> NameSiteMap;

    ResNameMap ResToSitename_;  ///< Map residue numbers to site names
    NameSiteMap NameToSite_;    ///< Map site names to titratable site data
};
}
}
#endif
