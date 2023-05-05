#include "ChangeTitratableResNames.h"
#include "SiteData.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"
#include <map>

int Cpptraj::Structure::ChangeTitratableResNames(Topology&, std::string const& sitesDirNameIn) {
  std::string sitesDirName;
  if (sitesDirNameIn.empty())
    sitesDirName = SiteData::DefaultSiteDir();
  else
    sitesDirName = sitesDirNameIn;

  // Load titratable residue name map
  std::string resMapFile = sitesDirName + "/TitratableResMap.dat";
  mprintf("\tLoading titratable residue name map from '%s'\n", resMapFile.c_str());
  std::map<std::string, std::string> resNameMap;

  

  return 0;
}
