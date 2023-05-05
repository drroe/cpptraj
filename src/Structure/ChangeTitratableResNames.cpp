#include "ChangeTitratableResNames.h"
#include "SiteData.h"
#include "StructureRoutines.h"
#include "../BufferedLine.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include <map>

int Cpptraj::Structure::ChangeTitratableResNames(Topology& topIn, std::string const& sitesDirNameIn) {
  std::string sitesDirName;
  if (sitesDirNameIn.empty())
    sitesDirName = SiteData::DefaultSiteDir();
  else
    sitesDirName = sitesDirNameIn;

  // Load titratable residue name map
  std::string resMapFile = sitesDirName + "/TitratableResMap.dat";
  mprintf("\tLoading titratable residue name map from '%s'\n", resMapFile.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( resMapFile )) {
    mprinterr("Error: Could not open titratable residue name map file '%s'\n", resMapFile.c_str());
    return 1;
  }
  typedef std::map<std::string, std::string> ResNameMapType;
  typedef std::pair<std::string, std::string> ResNamePairType;
  ResNameMapType resNameMap;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (*ptr != '#') {
      // Expected format: <original name> <new name>
      int ntokens = infile.TokenizeLine(" ");
      if (ntokens != 2) {
        mprinterr("Error: Malformed line in titratable residue name map file '%s': %s\n",
                  infile.Filename().full(), ptr);
        return 1;
      }
      std::string name0( infile.NextToken() );
      std::string name1( infile.NextToken() );
      resNameMap.insert( ResNamePairType( name0, name1 ) );
    }
    ptr = infile.Line();
  }
  infile.CloseFile();
  mprintf("\tTitratable residue name map:\n");
  for (ResNameMapType::const_iterator it = resNameMap.begin(); it != resNameMap.end(); ++it)
    mprintf("\t  %s -> %s\n", it->first.c_str(), it->second.c_str());

  // Change residue names where appropriate
  for (int ires = 0; ires < topIn.Nres(); ires++) {
    ResNameMapType::const_iterator it = resNameMap.find( topIn.Res(ires).Name().Truncated() );
    if (it != resNameMap.end()) {
      mprintf("\tChanging residue %s name to %s\n", topIn.TruncResNameNum(ires).c_str(), it->second.c_str());
      ChangeResName( topIn.SetRes(ires), it->second );
    }
  }

  return 0;
}
