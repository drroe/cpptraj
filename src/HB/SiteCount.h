#ifndef INC_HB_SITECOUNT_H
#define INC_HB_SITECOUNT_H
#include "HbEnum.h"
namespace Cpptraj {
namespace HB {
/// Hold count of hydrogen bond heavy atom sites
class SiteCount {
  public:
    /// CONSTRUCTOR
    SiteCount();
    /// Add site of specified type with number of bonded hydrogens
    void AddSite(Type, unsigned int);
    /// Add ion site
    void AddIon();
  private:
    unsigned int nsites_;
    unsigned int NacceptorOnly_;
    unsigned int Nboth_;
    unsigned int NdonorOnly_;
    unsigned int NumH_;
    unsigned int NV_acceptorOnly_;
    unsigned int NV_both_;
    unsigned int NV_donorOnly_;
    unsigned int NV_H_;
    unsigned int NIons_;
};
}
}
#endif
