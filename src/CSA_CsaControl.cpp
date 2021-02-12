#include "CSA_CsaControl.h"
#include "CpptrajStdio.h"
#include "ArgList.h"

using namespace Cpptraj::CSA;

Cpptraj::CSA::CsaControl::CsaControl() :
  num_conf_(0),
  num_iter_(0),
  num_seed_(0),
  rnd_per_iter_(0),
  nmove1_(0),
  nmove2_(0),
  nmove3_(0),
  nmove4_(0),
  move_maxnres_(0),
  move_res0_(0),
  move_res1_(0),
  idist_(NO_DIST),
  ietyp_(SCORE_NONE),
  imanip_(MANIP_NONE),
  idcut_(DCUT_NONE),
  iranseed_(0),
  nolower_cut_(0), 
  use_seeds_(false),
  free_ene_(false),
  write_banks_(false),
  min_ecut_(0),
  dcut0_fac_(0),
  dcut1_fac_(0),
  temperature_(0)
{}

void Cpptraj::CSA::CsaControl::Help() {
  mprintf("\t[nconf <#>] [niter <#>] [nseed <#>]\n");
}

/** Process options */
int Cpptraj::CSA::CsaControl::InitCsa(ArgList& argIn) {
  num_conf_ = argIn.getKeyInt("nconf", 50);
  num_iter_ = argIn.getKeyInt("niter", 99999);
  num_seed_ = argIn.getKeyInt("nseed", 30);
  rnd_per_iter_ = argIn.getKeyInt("roundsperiter", -1);

  return 0;
}
