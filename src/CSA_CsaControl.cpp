#include "CSA_CsaControl.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include <limits>

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
  idist_(DIST_NONE),
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
  nmove1_ = argIn.getKeyInt("nmove1", 10);
  nmove2_ = argIn.getKeyInt("nmove2", 10);
  nmove3_ = argIn.getKeyInt("nmove3", 10);
  nmove4_ = argIn.getKeyInt("nmove4", 10);
  move_maxnres_ = argIn.getKeyInt("movemaxnres", 5);
  move_res0_ = argIn.getKeyInt("moveres0", 4);
  move_res1_ = argIn.getKeyInt("moveres1", 8);
  std::string distarg = argIn.GetStringKey("dist");
  if (distarg == "angle")
    idist_ = DIST_ANGLE;
  else if (distarg == "rmsd")
    idist_ = DIST_RMSD;
  else if (distarg == "tm")
    idist_ = DIST_TMSCORE;
  else if (distarg == "frechet")
    idist_ = DIST_FRECHET;
  else
    idist_ = DIST_NONE;
  std::string earg = argIn.GetStringKey("etype");
  if (earg == "pe")
    ietyp_ = SCORE_PE;
  else if (earg == "om")
    ietyp_ = SCORE_OM;
  else
    ietyp_ = SCORE_NONE;
  std::string maniparg = argIn.GetStringKey("manip");
  if (maniparg == "cart")
    imanip_ = MANIP_CART;
  else if (maniparg == "ic")
    imanip_ = MANIP_IC;
  else
    imanip_ = MANIP_NONE;
  std::string dcutarg = argIn.GetStringKey("dcut");
  if (dcutarg == "linear")
    idcut_ = DCUT_LINEAR;
  else if (dcutarg == "geom")
    idcut_ = DCUT_GEOM;
  else
    idcut_ = DCUT_NONE;
  iranseed_ = argIn.getKeyInt("ranseed", -1);
  nolower_cut_ = argIn.getKeyInt("nolowercut", -1);

  use_seeds_ = argIn.hasKey("useseeds");
  free_ene_ = argIn.hasKey("freeenergy");
  write_banks_ = argIn.hasKey("writebanks");

  min_ecut_ = argIn.getKeyDouble("minecut", -std::numeric_limits<double>::max());
  dcut0_fac_ = argIn.getKeyDouble("dcut0", 2);
  dcut1_fac_ = argIn.getKeyDouble("dcut1", 5);
  temperature_ = argIn.getKeyDouble("temp", 300.0);
  

  return 0;
}
