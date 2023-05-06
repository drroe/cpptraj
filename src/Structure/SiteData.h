#ifndef INC_TITRATIONDATA_H
#define INC_TITRATIONDATA_H
#include <vector>
#include <string>
#include <map>
class Topology;
namespace Cpptraj {
//TODO Move to Mead namespace
namespace Structure {
// Fwd declares
class TitratableSite;
/// Hold information for all titratable sites
class SiteData {
    class Tsite;
    //typedef std::pair<int, std::string> IdxNamePair;
    typedef std::vector<Tsite> IdxNameArray;
  public:
    /// Different site types
    enum SiteType { T_NONE = 0, ///< No specific type
                    T_N,        ///< Protein N-terminal
                    T_C };      ///< Protein C-terminal

    /// CONSTRUCTOR
    SiteData();
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Load titration data and sites from a MEAD-style site file.
    int LoadMeadSiteData(std::string const&, std::string const&, Topology const&);
    /// Load titration data from specified directory.
    int LoadSiteDirectory(std::string const&);
    /// Set up sites from Topology
    int SetupSitesFromTop(Topology const&);

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
    /// \return Default directory name containing titratable site data
    static std::string DefaultSiteDir();
  private:
    typedef std::pair<std::string,TitratableSite> NameSitePair;
    typedef std::map<std::string,TitratableSite> NameSiteMap;

    typedef std::vector<std::string> Sarray;
    typedef std::pair<std::string, Sarray> ResSitesPair;
    typedef std::map<std::string, Sarray> ResSitesMap;

    class ProtInfo;
    typedef std::pair<std::string, ProtInfo> ResProtPair;
    typedef std::map<std::string, ProtInfo> ResProtMap;

    ///< Use to hold topology-related info for a site
    class Tsite {
      public:
        //Tsite() : ridx_(-1), stype_(T_NONE) {}
        Tsite(int r, std::string const& n) : ridx_(r), sname_(n) {
          // Site name start with or end with NT, N-terminal
          if ( (sname_[0] == 'N' && sname_[1] == 'T') ||
               (sname_[sname_.size()-2] == 'N' && sname_[sname_.size()-1] == 'T') )
            stype_ = T_N;
          else if ( sname_[0] == 'C' && sname_[1] == 'T' )
            stype_ = T_C;
          else
            stype_ = T_NONE;
        }
        int Ridx() const { return ridx_; }
        std::string const& Sname() const { return sname_; }
        SiteType Stype() const { return stype_; }
      private:
        int ridx_;          ///< Internal residue index that site corresponds to.
        std::string sname_; ///< Site name (in NameToSite_)
        SiteType stype_;    ///< Site type
    };

    void PrintTitrationSiteData() const;
    int loadProtInfo(std::string const&);

    NameSiteMap NameToSite_;     ///< Map site names to titratable site data
    IdxNameArray IdxNames_;      ///< Hold res #s/site names in order.
    ResSitesMap resnameToSites_; ///< Map residue name to site names
    ResProtMap resnameToProt_;   ///< Map residue names to protonation info
    int debug_;
};
}
}
#endif
