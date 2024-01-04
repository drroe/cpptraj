#include "GB_Params.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"

static const char* GB_RadiiTypeStr_[] = {
  "Bondi radii",
  "Amber 6 modified Bondi radii",
  "Modified Bondi radii",
  "Radii optimized for Amber charges by Huo and Kollman",
  "H(N)-modified Bondi radii",
  "PARSE radii",
  "ArgH and AspGluO modified Bondi2 radii",
  "Unknown GB radii set"
};

static const char* GB_RadiiTypeKey_[] = {
  "bondi",
  "amber6",
  "mbondi",
  "pbamber",
  "mbondi2",
  "parse",
  "mbondi3"
};

/** Assign GB radii and screening parameters based on the given radius set. */
int Cpptraj::Parm::Assign_GB_Radii(Topology& top, GB_RadiiType radiiSet)
{
  if (radiiSet == UNKNOWN_GB) {
    mprinterr("Error: Unknown GB radii set.\n");
    return 1;
  }
  mprintf("\tUsing GB radii set: %s\n", GB_RadiiTypeStr_[radiiSet]);

  return 0;
}
