#include "SiteCount.h"
#include "../CpptrajStdio.h"

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

/** Print counts to stdout. */
void SiteCount::PrintCounts(bool calcSolvent) const {
  mprintf("\tTotal Number of heavy atom sites: %i\n", nsites_);
  mprintf("\t  Solute acceptor-only atoms: %u\n", NacceptorOnly_);
  mprintf("\t  Solute donor/acceptor sites: %u\n", Nboth_);
  mprintf("\t  Solute donor-only sites: %u\n", NdonorOnly_);
  mprintf("\t  %u solute hydrogens.\n", NumH_);
  if (calcSolvent) {
    mprintf("\t  Solvent acceptor-only atoms: %u\n", NV_acceptorOnly_);
    mprintf("\t  Solvent donor/acceptor sites: %u\n", NV_both_);
    mprintf("\t  Solvent donor-only sites: %u\n", NV_donorOnly_);
    mprintf("\t  %u solvent hydrogens, %u ions.\n", NV_H_, NIons_);
  }
}

/** \return Potential number of solute-solute interactions */
unsigned int SiteCount::UUsize() const {
  return (NumH_ * (Nboth_ + NacceptorOnly_));
}

/** \return Potential number of solute-solvent interactions. */
unsigned int SiteCount::UVsize() const {
  return ((NumH_ * (NV_both_ + NV_acceptorOnly_ + NIons_)) +
          (Nboth_ + NacceptorOnly_));
}
