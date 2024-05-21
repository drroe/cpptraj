#ifndef INC_HB_SITECOUNT_H
#define INC_HB_SITECOUNT_H
#include <vector>
#include "HbEnum.h"
namespace Cpptraj {
namespace HB {
/// Hold count of hydrogen bond heavy atom sites
class SiteCount {
    typedef std::vector<int> Iarray;
  public:
    /// CONSTRUCTOR
    SiteCount();
    /// Add site of specified type
    void AddSite(Type, Iarray const&);
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
