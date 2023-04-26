#include "SiteData.h"
#include "TitratableSite.h"
#include "../NameType.h"
#include "../CpptrajStdio.h"
#include "../BufferedLine.h"
#include "../ArgList.h"
#include <cstdlib> // atoi, getenv

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
SiteData::SiteData()
{}

std::string SiteData::defaultSiteDir() {
  const char* env = getenv("CPPTRAJHOME");
  if (env != 0)
    return (std::string(env) + std::string("/dat/TitratableSites"));
  return std::string(".");
}
  

/** Load all sites in a given directory. */
int SiteData::LoadSiteDirectory(std::string const& sitesDirNameIn)
{
  std::string sitesDirName;
  if (sitesDirNameIn.empty())
    sitesDirName = defaultSiteDir();
  else
    sitesDirName = sitesDirNameIn;

  std::string resToSiteFile = sitesDirName + "/ResToSite.dat";
  mprintf("Loading site titration data from '%s'\n", resToSiteFile.c_str());

  BufferedLine infile;
  if (infile.OpenFileRead( resToSiteFile )) {
    mprinterr("Error: Could not open titratable site data file '%s'\n", resToSiteFile.c_str());
    return 1;
  }
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (*ptr != '#') {
      // Expected format: <resname list> <Filename list>
      int ntokens = infile.TokenizeLine(" ");
      if (ntokens != 2) {
        mprinterr("Error: Malformed line in titratable site data file '%s': %s\n", resToSiteFile.c_str(), ptr);
        mprinterr("Error: Expected 2 tokens, got %i\n", ntokens);
        return 1;
      }
      ArgList resList( infile.NextToken(), "," );
      ArgList fileList( infile.NextToken(), "," );
      mprintf("DEBUG: %i res %i files\n", resList.Nargs(), fileList.Nargs());
      // Loop over res names
      for (ArgList::const_iterator it = resList.begin(); it != resList.end(); ++it) {
        // Check if residue is already defined
        ResSitesMap::iterator rs = resnameToSites_.lower_bound( *it );
        if (rs == resnameToSites_.end() || rs->first != *it) {
          // Define residue
          resnameToSites_.insert( rs, ResSitesPair(*it, fileList.List()) );
        } else {
          mprinterr("Error: Residue '%s' already has an entry.\n", it->c_str());
        }
      }
      // Load site data
    } // END line is not comment
    ptr = infile.Line();
  } // END loop over file lines

  mprintf("DEBUG: Titration data:\n");
  for (ResSitesMap::const_iterator it = resnameToSites_.begin(); it != resnameToSites_.end(); ++it)
  {
    mprintf("\t%s :", it->first.c_str());
    for (Sarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      mprintf(" %s", jt->c_str());
    mprintf("\n");
  }
  return 0;
}


/** Load titratable sites data.
  * \param sitesFileName File name containing residue numbers and site names.
  * \param sitesDirName Directory name containins site files.
  */
int SiteData::LoadMeadSiteData(std::string const& sitesFileName,
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
    // Expect residue number to start from 1
    int rnum = atoi( sitesFile.NextToken() );
    std::string sname( sitesFile.NextToken() );
    mprintf("DEBUG: Res#=%i  siteName=%s\n", rnum, sname.c_str());
    /*ResNameMap::iterator it = ResToSitename_.lower_bound( rnum );
    if (it == ResToSitename_.end() || it->first != rnum) {
      it = ResToSitename_.insert( it, ResNamePair(rnum, Sarray(1,sname)) );
    } else {
      // A residue number may have more than one site
      mprintf("DEBUG: Residue number %i already has an entry.\n", rnum);
      // Ensure this is not a duplicate
      for (Sarray::const_iterator nm = it->second.begin(); nm != it->second.end(); ++nm) {
        if (*nm == sname) {
          mprinterr("Error: Residue number %i already has entry '%s'\n", rnum, sname.c_str());
          return 1;
        }
      }
      it->second.push_back( sname );
    }*/
    // See if there is data for this site yet.
    NameSiteMap::iterator ns = NameToSite_.lower_bound( sname );
    if (ns == NameToSite_.end() || ns->first != sname) {
      // Load site data
      std::string stFileName;
      if (sitesDirName.empty())
        stFileName = sname + ".st";
      else
        stFileName = sitesDirName + "/" + sname + ".st";
      mprintf("DEBUG: Loading data for site '%s' from '%s'\n", sname.c_str(), stFileName.c_str());
      ns = NameToSite_.insert( ns, NameSitePair(sname, TitratableSite()) );
      if (ns->second.LoadSiteFile( stFileName, sname )) {
        mprinterr("Error: Failed to load titratable site data from '%s'\n", stFileName.c_str());
        return 1;
      }
    } else {
      mprintf("DEBUG: Site '%s' already has data.\n", sname.c_str());
    }
    IdxNames_.push_back( IdxNamePair(rnum, sname) );
    // Next line
    ptr = sitesFile.Line();
  }
  sitesFile.CloseFile();

  mprintf("DEBUG: Titratable sites:\n");
  for (IdxNameArray::const_iterator it = IdxNames_.begin(); it != IdxNames_.end(); ++it)
    mprintf("\t%6i %s\n", it->first, it->second.c_str());
/*  for (ResNameMap::const_iterator it = ResToSitename_.begin(); it != ResToSitename_.end(); ++it)
  {
    mprintf("\t%8i :", it->first + 1);
    for (Sarray::const_iterator nm = it->second.begin(); nm != it->second.end(); ++nm)
      mprintf(" %s", nm->c_str());
    mprintf("\n");
  }*/
  mprintf("DEBUG: Titration site data:\n");
  for (NameSiteMap::const_iterator ns = NameToSite_.begin(); ns != NameToSite_.end(); ++ns)
  {
    mprintf("\tSite: %s\n", ns->first.c_str());
    ns->second.Print();
  }
  return 0;
}

/** \return Site corresponding to site name. */
TitratableSite const& SiteData::GetSite(std::string const& sname) const {
  NameSiteMap::const_iterator it = NameToSite_.find( sname );
  if (it == NameToSite_.end()) {
    mprinterr("Internal Error: Titratable site data '%s' not present.\n", sname.c_str());
  }
  return it->second;
}
