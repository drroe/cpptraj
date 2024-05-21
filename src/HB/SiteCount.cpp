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
void SiteCount::AddSite(Type currentType, Iarray const& h_atoms) {
  switch (currentType) {
    case VACCEPTOR : NV_acceptorOnly_++; break;
    case ACCEPTOR  : NacceptorOnly_++; break;
    case VDONOR    : NV_donorOnly_++; NV_H_ += h_atoms.size(); break;
    case DONOR     : NdonorOnly_++;   NumH_ += h_atoms.size(); break;
    case VBOTH     : NV_both_++;      NV_H_ += h_atoms.size(); break;
    case BOTH      : Nboth_++;        NumH_ += h_atoms.size(); break;
    case UNKNOWN   : return;
  }
  nsites_++;
}
