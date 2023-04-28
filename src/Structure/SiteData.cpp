#include "SiteData.h"
#include "TitratableSite.h"
#include "../NameType.h"
#include "../CpptrajStdio.h"
#include "../BufferedLine.h"
#include "../ArgList.h"
#include "../Topology.h"
#include <cstdlib> // atoi, getenv

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
SiteData::SiteData() :
  debug_(0)
{}

std::string SiteData::defaultSiteDir() {
  const char* env = getenv("CPPTRAJHOME");
  if (env != 0)
    return (std::string(env) + std::string("/dat/TitratableSites"));
  return std::string(".");
}

/** Print titration site data to stdout. */
void SiteData::PrintTitrationSiteData() const {
   mprintf("DEBUG: Titration site data:\n");
  for (NameSiteMap::const_iterator ns = NameToSite_.begin(); ns != NameToSite_.end(); ++ns)
  {
    mprintf("\tSite: %s\n", ns->first.c_str());
    ns->second.Print();
  }
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
  mprintf("\tLoading site titration data from '%s'\n", resToSiteFile.c_str());

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
//      mprintf("DEBUG: %i res %i files\n", resList.Nargs(), fileList.Nargs());
      // Loop over res names
      for (ArgList::const_iterator it = resList.begin(); it != resList.end(); ++it) {
        // Check if residue is already defined
        ResSitesMap::iterator rs = resnameToSites_.lower_bound( *it );
        if (rs == resnameToSites_.end() || rs->first != *it) {
          // Load site data
          Sarray siteNames;
          for (ArgList::const_iterator file = fileList.begin(); file != fileList.end(); ++file)
          {
            // Develop the site name
            FileName fname( *file );
            size_t found = fname.Base().find_last_of(".");
            std::string sname;
            if (found == std::string::npos)
              sname = fname.Base();
            else
              sname = fname.Base().substr(0, found);
            siteNames.push_back( sname );
            // Check if site is present
            NameSiteMap::iterator ns = NameToSite_.lower_bound( siteNames.back() );
            if (ns == NameToSite_.end() || ns->first != siteNames.back()) {
              // Load site data
              std::string stFileName = sitesDirName + "/" + *file;
              if (debug_ > 1) mprintf("DEBUG: Loading data for site '%s' from '%s'\n", siteNames.back().c_str(), stFileName.c_str());
              ns = NameToSite_.insert( ns, NameSitePair(siteNames.back(), TitratableSite()) );
              if (ns->second.LoadSiteFile( stFileName, siteNames.back() )) {
                mprinterr("Error: Failed to load titratable site data from '%s'\n", stFileName.c_str());
                return 1;
              }
            } else {
              mprintf("DEBUG: Site '%s' already has data.\n", siteNames.back().c_str());
            }

          } // END loop over site files for res
          // Define residue
          resnameToSites_.insert( rs, ResSitesPair(*it, siteNames) );
        } else {
          mprinterr("Error: Residue '%s' already has an entry.\n", it->c_str());
        }
      }
      // Load site data
    } // END line is not comment
    ptr = infile.Line();
  } // END loop over file lines

  mprintf("\tTitration data:\n");
  for (ResSitesMap::const_iterator it = resnameToSites_.begin(); it != resnameToSites_.end(); ++it)
  {
    mprintf("\t%s :", it->first.c_str());
    for (Sarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      mprintf(" %s", jt->c_str());
    mprintf("\n");
  }
  if (debug_ > 0) PrintTitrationSiteData();
  return 0;
}


/** Load titratable sites data.
  * \param sitesFileName File name containing residue numbers and site names.
  * \param sitesDirName Directory name containins site files.
  * \param topIn Topology for converting res #s in MEAD file to internal res #s.
  */
int SiteData::LoadMeadSiteData(std::string const& sitesFileName,
                               std::string const& sitesDirName,
                               Topology const& topIn)
{
  if (sitesFileName.empty()) {
    mprinterr("Internal Error: sitesFileName is empty.\n");
    return 1;
  }
  mprintf("\tLoading titratable sites from MEAD sites file '%s'\n", sitesFileName.c_str());
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
    // Expect residue number to start from 1 FIXME take chain ID into account 
    int rnum = atoi( sitesFile.NextToken() );
    std::string sname( sitesFile.NextToken() );
    // Convert from original residue # to internal residue index
    int ridx = -1;
    for (int ires = 0; ires < topIn.Nres(); ires++) {
      if (topIn.Res(ires).OriginalResNum() == rnum) {
        ridx = ires;
        break;
      }
    }
    if (ridx < 0 ) {
      mprinterr("Error: Residue number %i not found in topology '%s'.\n", rnum, topIn.c_str());
      return 1;
    }

    mprintf("DEBUG: Res#=%i  ResIndex=%i  siteName=%s\n", rnum, ridx, sname.c_str());
    // See if there is data for this site yet.
    NameSiteMap::iterator ns = NameToSite_.lower_bound( sname );
    if (ns == NameToSite_.end() || ns->first != sname) {
      // Load site data
      std::string stFileName;
      if (sitesDirName.empty())
        stFileName = sname + ".st";
      else
        stFileName = sitesDirName + "/" + sname + ".st";
//      mprintf("DEBUG: Loading data for site '%s' from '%s'\n", sname.c_str(), stFileName.c_str());
      ns = NameToSite_.insert( ns, NameSitePair(sname, TitratableSite()) );
      if (ns->second.LoadSiteFile( stFileName, sname )) {
        mprinterr("Error: Failed to load titratable site data from '%s'\n", stFileName.c_str());
        return 1;
      }
    }// else {
     // mprintf("DEBUG: Site '%s' already has data.\n", sname.c_str());
    //}
    IdxNames_.push_back( IdxNamePair(ridx, sname) );
    // Next line
    ptr = sitesFile.Line();
  }
  sitesFile.CloseFile();

  mprintf("\tTitratable sites:\n");
  for (IdxNameArray::const_iterator it = IdxNames_.begin(); it != IdxNames_.end(); ++it)
    mprintf("\t  %6i %s\n", it->first, it->second.c_str());

  if (debug_ > 0) PrintTitrationSiteData();

  return 0;
}

/** Set up sites array from input topology. */
int SiteData::SetupSitesFromTop(Topology const& topIn) {
  IdxNames_.clear();
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& thisRes = topIn.Res(ires);
    // See if there are sites for this res name
    ResSitesMap::const_iterator it = resnameToSites_.find( thisRes.Name().Truncated() );
    if (it != resnameToSites_.end()) {
      if (debug_ > 0) mprintf("DEBUG: %zu potential sites found for residue %s\n", it->second.size(), *(thisRes.Name()));
      // Check if this residue is terminal. TODO use residue IsTerminal?
      std::vector<int> bondedResIdxs;
      bondedResIdxs.reserve(4);
      for (int ii = thisRes.FirstAtom(); ii != thisRes.LastAtom(); ii++) {
        for (Atom::bond_iterator bat = topIn[ii].bondbegin(); bat != topIn[ii].bondend(); ++bat) {
          if (topIn[*bat].ResNum() != ires)
            bondedResIdxs.push_back( topIn[*bat].ResNum() );
        }
      }
      if (debug_ > 0) {
        mprintf("\tResidue '%s' is bonded to residues", topIn.TruncResNameNum(ires).c_str());
        for (std::vector<int>::const_iterator bri = bondedResIdxs.begin(); bri != bondedResIdxs.end(); ++bri)
          mprintf(" %s", topIn.TruncResNameNum(*bri).c_str());
        mprintf("\n");
      }
      // FIXME this is very protein-specific
      enum TerminalType { T_NONE = 0, T_N, T_C };
      TerminalType termType;
      if (bondedResIdxs.size() > 1)
        termType = T_NONE;
      else {
        if (bondedResIdxs[0] > ires)
          termType = T_N;
        else
          termType = T_C;
      }
      if (debug_ > 0) {
        if (termType == T_N)
          mprintf("\tN-TERMINAL.\n");
        else if (termType == T_C)
          mprintf("\tC-TERMINAL.\n");
      }
      // Loop over potential sites
      for (Sarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
        std::string const& sname = *jt;
        // FIXME not very general
        TerminalType siteTermType;
        // Site name must start with or end with NT
        if ( (sname[0] == 'N' && sname[1] == 'T') ||
             (sname[sname.size()-2] == 'N' && sname[sname.size()-1] == 'T') )
          siteTermType = T_N;
        else if ( sname[0] == 'C' && sname[1] == 'T' )
          siteTermType = T_C;
        else
          siteTermType = T_NONE;
        if (termType != siteTermType) continue;
            
        TitratableSite const& site = GetSite( sname );
        bool siteIsPresent = true;
        // All atoms of the site must be present
        for (TitratableSite::const_iterator at = site.begin(); at != site.end(); ++at) {
          // Get the atom index in the topology
          int aidx = topIn.FindAtomInResidue(ires, at->first);
          if (aidx < 0) {
            siteIsPresent = false;
            break;
          }
        } // END loop over site atoms
        if (siteIsPresent) {
          if (debug_ > 0) site.Print(); // DEBUG
          mprintf("\tFound Site '%s' in residue '%s'\n", jt->c_str(), topIn.TruncResNameNum(ires).c_str());
          IdxNames_.push_back( IdxNamePair(ires, sname) );
        }
      } // END loop over potential sites
    }
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
