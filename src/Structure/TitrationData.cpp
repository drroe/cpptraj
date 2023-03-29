#include "TitrationData.h"
#include "TitratableSite.h"
#include "../NameType.h"
#include "../CpptrajStdio.h"
#include "../BufferedLine.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
TitrationData::TitrationData()
{}

/** Load titratable sites data.
  * \param sitesFileName File name containing residue numbers and site names.
  * \param sitesDirName Directory name containins site files.
  */
int TitrationData::LoadTitrationData(std::string const& sitesFileName,
                                     std::string const& sitesDirName)
{
  if (sitesFileName.empty()) {
    mprinterr("Internal Error: sitesFileName is empty.\n");
    return 1;
  }
  BufferedLine sitesFile;
  if (sitesFile.OpenFileRead(sitesFileName)) {
    mprinterr("Error: Could not open sites file '%s'\n", sitesFileName.c_str());
    return 1;
  }
  const char* ptr = sitesFile.Line();
  while (ptr != 0) {
    // Expected format: <res#> <siteName>
    int ntokens = sitesFile.TokenizeLine(" ");
    if (ntokens != 2) {
      mprinterr("Error: Malformed line in sites file '%s': %s\n", sitesFileName.c_str(), ptr);
      mprinterr("Error: Expected 2 tokens, got %i\n", ntokens);
      return 1;
    }
    int rnum = atoi( sitesFile.NextToken() ) - 1;
    std::string sname( sitesFile.NextToken() );
    mprintf("DEBUG: Res#=%i  siteName=%s\n", rnum + 1, sname.c_str());
    ResNameMap::iterator it = ResToSitename_.lower_bound( rnum );
    if (it == ResToSitename_.end() || it->first != rnum) {
      // rnum - 1
      it = ResToSitename_.insert( it, ResNamePair(rnum, Sarray(1,sname)) );
    } else {
      // A residue number may have more than one site
      mprintf("DEBUG: Residue number %i already has an entry.\n", rnum + 1);
      // Ensure this is not a duplicate
      for (Sarray::const_iterator nm = it->second.begin(); nm != it->second.end(); ++nm) {
        if (*nm == sname) {
          mprinterr("Error: Residue number %i already has entry '%s'\n", rnum+1, sname.c_str());
          return 1;
        }
      }
      it->second.push_back( sname );
    }
    // See if there is data for this site yet.
    // Next line
    ptr = sitesFile.Line();
  }
  sitesFile.CloseFile();

  mprintf("DEBUG: Titratable sites:\n");
  for (ResNameMap::const_iterator it = ResToSitename_.begin(); it != ResToSitename_.end(); ++it)
  {
    mprintf("\t%8i :", it->first + 1);
    for (Sarray::const_iterator nm = it->second.begin(); nm != it->second.end(); ++nm)
      mprintf(" %s", nm->c_str());
    mprintf("\n");
  }
  return 0;
}
    
