#include "CmapParmHolder.h"
#include "CpptrajStdio.h"
#include "ParameterTypes.h"

//using namespace Cpptraj::Parm;

/** CONSTRUCTOR */
CmapParmHolder::CmapParmHolder() {}

/** Update/add to CMAP parameters.
  * CMAP terms are different in that they are applied via residue name instead
  * of atom type, hence the separate routine.
  */
ParameterHolders::RetType CmapParmHolder::AddParm(CmapGridType const& cmap1, bool allowUpdate) {
  typedef std::vector<std::string> Sarray;
  enum ResMatchType { NO_MATCH = 0, FULL_MATCH, PARTIAL_MATCH };
  // Does a CMAP for any of the residues exist?
  ResMatchType mtype = NO_MATCH;
  CmapGridArray::iterator currentCmap = CMAP_.end();
  for (CmapGridArray::iterator jt = CMAP_.begin(); jt != CMAP_.end(); ++jt)
  {
    int nMatch = 0;
    for (Sarray::const_iterator rni = cmap1.ResNames().begin(); rni != cmap1.ResNames().end(); ++rni)
    {
      for (Sarray::const_iterator rnj = jt->ResNames().begin(); rnj != jt->ResNames().end(); ++rnj)
      {
        if (*rni == *rnj) {
          nMatch++;
          break;
        }
      } // END loop over cmap j res names
    } // END loop over cmap i res names
    if (nMatch == (int)jt->ResNames().size()) {
      mtype = FULL_MATCH;
      currentCmap = jt;
      break;
    } else if (nMatch > 0) {
      mtype = PARTIAL_MATCH;
      currentCmap = jt;
      break;
    }
  } // END loop over existing CMAP terms
  if (mtype == PARTIAL_MATCH) {
    mprinterr("Internal Error: Partial match between new CMAP '%s' and existing CMAP '%s'\n"
              "Internal Error:   Not yet set up to handle this kind of parameter update.\n",
              cmap1.Title().c_str(), currentCmap->Title().c_str());
    return ParameterHolders::ERR;
  } else if (mtype == FULL_MATCH) {
    mprintf("DEBUG: Match between new CMAP '%s' and existing CMAP '%s'\n",
            cmap1.Title().c_str(), currentCmap->Title().c_str());
    if (allowUpdate) {
      *currentCmap = cmap1;
      return ParameterHolders::UPDATED;
    } else
      return ParameterHolders::ERR;
  }
  // NO_MATCH
  mprintf("DEBUG: CMAP '%s' is a new CMAP.\n", cmap1.Title().c_str());
  CMAP_.push_back( cmap1 );
  return ParameterHolders::ADDED;
}

