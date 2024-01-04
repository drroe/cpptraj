#ifndef INC_PARM_GBPARAMS_H
#define INC_PARM_GBPARAMS_H
#include <string>
class Topology;
namespace Cpptraj {
namespace Parm {
/// Known radii types. KEEP SYNCED WITH GB_RadiiTypeStr[] in GB_Params.cpp
enum GB_RadiiType { 
  BONDI=0,          ///< 0, bondi, Bondi radii
  BONDI_AMBER6,     ///< 1, amber6, Amber6 modified Bondi radii
  MBONDI,           ///< 2, mbondi, Modified bondi radii
  PB_AMBER,         ///< 3, pbamber, PB radii from Huo and Kollman, currently unused
  MBONDI2,          ///< 6, mbondi2, H(N)-modified Bondi radii
  PARSE,            ///< 7, parse, PARSE radii
  MBONDI3,          ///< 8, ArgH and AspGluO modified Bondi2 radii
  UNKNOWN_GB
};
/// \return GB radii type corresponding to string
GB_RadiiType GbTypeFromKey(std::string const&);
/// Assign GB radii
int Assign_GB_Radii(Topology&, GB_RadiiType);
}
}
#endif 
