#include "SiteCount.h"

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
SiteCount::SiteCount() :
  nsites_(0),
  NacceptorOnly_(0),
  Nboth_(0),
  NdonorOnly_(0),
  NumH_(0),
  NV_acceptorOnly_(0),
  NV_both_(0),
  NV_donorOnly_(0),
  NV_H_(0),
  NIons_(0)
{}

/** Add site of given type */
void SiteCount::AddSite(Type currentType, unsigned int nH) {
  switch (currentType) {
    case VACCEPTOR : NV_acceptorOnly_++; break;
    case ACCEPTOR  : NacceptorOnly_++; break;
    case VDONOR    : NV_donorOnly_++; NV_H_ += nH; break;
    case DONOR     : NdonorOnly_++;   NumH_ += nH; break;
    case VBOTH     : NV_both_++;      NV_H_ += nH; break;
    case BOTH      : Nboth_++;        NumH_ += nH; break;
    case UNKNOWN   : return;
  }
  nsites_++;
}

/** Add Ion site. */
void SiteCount::AddIon() {
  NIons_++;
  nsites_++;
}
