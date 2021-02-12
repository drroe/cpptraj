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

/** Print status to STDOUT */
void Cpptraj::CSA::CsaControl::Status() const {
  //  write(outu,'(a,i6)') 'CSA> # selected atoms= ', nselected_atm 
  mprintf("\tNumber of conformations to store in bank/firstbank= %i\n", num_conf_);
  mprintf("\tNumber of CSA iterations= %i\n", num_iter_);
  mprintf("\tNumber of seeds to use each round= %i\n", num_seed_);
  if (use_seeds_)
    mprintf("\tSeeds will be used as trials.\n");
  else
    mprintf("\tSeeds will not be used as trials.\n");
  if (rnd_per_iter_ > 0) mprintf("\tRounds per iteration= %i\n", rnd_per_iter_);
  if (write_banks_) mprintf("\tWriting first bank and banks during iterations.\n");
  // TODO better description of the move types and which ranges apply to what moves
  mprintf("\tNumber of moves of type 1= %i\n", nmove1_);
  mprintf("\tNumber of moves of type 2= %i\n", nmove2_);
  mprintf("\tNumber of moves of type 3= %i\n", nmove3_);
  mprintf("\tNumber of moves of type 4= %i\n", nmove4_);
  mprintf("\tMax number of residues to manipulate= %i\n", move_maxnres_);
  mprintf("\tWhen moving blocks of residues, will move anywhere from %i to %i residues.\n", move_res0_, move_res1_);
  mprintf("\tIterations will stop if minimum conformation energy is less than %g\n", min_ecut_);
  if (nolower_cut_ > 0)
    mprintf("\tIterations will stop if lower E conf cannot be found %i times in a row.\n", nolower_cut_);
  mprintf("\tInitial value of Dcut will be <average distance> / %g\n", dcut0_fac_);
  mprintf("\tFinal value of Dcut will be <average distance> / %g", dcut1_fac_);
  mprintf("\tTemperature= %g\n", temperature_);
  //  write(outu,'(a,a)') 'CSA> Molecule type= ', trim(mtype)
  mprintf("\tRandom number generator seed= %i\n", iranseed_);
  if (free_ene_) mprintf("\tEnergy will be modified by -T*S\n");

  if (idist_ == DIST_ANGLE) mprintf("\tDistance type: cumulative angle distance.\n");
  else if (idist_ == DIST_RMSD) mprintf("\tDistance type: RMSD\n");
  else if (idist_ == DIST_TMSCORE) mprintf("\tDistance type: TM score\n");
  else if (idist_ == DIST_FRECHET) mprintf("\tDistance type: Frechet distance\n");

  if (ietyp_ == SCORE_PE) mprintf("\tEnergy type: potential energy.\n");
  else if (ietyp_ == SCORE_OM) mprintf("\tEnergy type: Onsager-Machlup action.\n");

  if (imanip_ == MANIP_IC) mprintf("\tManipulating internal coordinates.\n");
  else if (imanip_ == MANIP_CART) mprintf("\tManipulating Cartesian coordinates.\n");

  if (idcut_ == DCUT_LINEAR) mprintf("\tDcut will be decreased linearly after each update.\n");
  else if (idcut_ == DCUT_GEOM) mprintf("\tDcut will be decreased geometrically after each update.\n");

}
