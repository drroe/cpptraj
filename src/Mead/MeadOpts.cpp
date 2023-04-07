#include "MeadOpts.h"

Cpptraj::Mead::MeadOpts::MeadOpts() :
  epsin_(4),
  epsext_(80),
  epsvac_(1),
  solrad_(1.4),
  sterln_(2.0),
  ionicstr_(0),
  temp_(300)
{}
