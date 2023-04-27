#include "MeadInterface.h"
#include "MultiFlexResults.h"
#include "MeadOpts.h"
#include "MeadGrid.h"
#include "MeadError.h"
#include "../Vec3.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../DataSet_Vector_Scalar.h"
#include "../DataSet_3D.h"
#include "../DataSet_1D.h"
#include "../DataSet_string.h"
#include "../Structure/SiteData.h"
#include "../Structure/TitratableSite.h"
// MEAD includes
#include "../../mead/MEADexcept.h"
#include "../../mead/AtomSet.h"
#include "../../mead/ChargeDist.h"
#include "../../mead/AtomChargeSet.h"
#include "../../mead/DielectricEnvironment.h"
#include "../../mead/DielByAtoms.h"
#include "../../mead/ElectrolyteEnvironment.h"
#include "../../mead/ElectrolyteByAtoms.h"
#include "../../mead/FinDiffElstatPot.h"
#include "../../mead/Potat.h"

using namespace Cpptraj::Mead;

/** CONSTRUCTOR */
MeadInterface::MeadInterface() :
  atomset_(0),
  rmode_(GB),
  t_total_("Mead Total")
{
  t_total_.AddSubTimer(Timer("SetupAtoms"));  // 0
  t_total_.AddSubTimer(Timer("MultiFlex "));  // 1
  t_total_[1].AddSubTimer(Timer("Setup sites")); // 1,0
  t_total_[1].AddSubTimer(Timer("MAC1       ")); // 1,1
  t_total_[1].AddSubTimer(Timer("MAC2       ")); // 1,2
  t_total_[1].AddSubTimer(Timer("MOD1       ")); // 1,3
  t_total_[1].AddSubTimer(Timer("MOD2       ")); // 1,4
}

/** DESTRUCTOR */
MeadInterface::~MeadInterface() {
  if (atomset_ != 0) delete atomset_;
}

/** Set MEAD Atom from Topology Atom. */
void MeadInterface::set_at_from_top(MEAD::Atom& at, Topology const& topIn, Frame const& frameIn, int aidx)
const
{
    Atom const& thisAtom = topIn[aidx];
    at.atname.assign( thisAtom.Name().Truncated() );
    int rnum = thisAtom.ResNum();
    Residue const& thisRes = topIn.Res(rnum);
    at.resname.assign( thisRes.Name().Truncated() );
    at.resnum = thisRes.OriginalResNum();
    if (thisRes.HasChainID())
      at.chainid.assign( 1, thisRes.ChainId() );
    const double* xyz = frameIn.XYZ(aidx);
    at.coord.x = xyz[0];
    at.coord.y = xyz[1];
    at.coord.z = xyz[2];
    at.charge = thisAtom.Charge();
    switch (rmode_) {
      case MeadInterface::GB : at.rad = thisAtom.GBRadius(); break;
      case MeadInterface::PARSE : at.rad = thisAtom.ParseRadius(); break;
      case MeadInterface::VDW   : at.rad = topIn.GetVDWradius(aidx); break;
    }
}

/** Setup an AtomSet from Frame and Topology. */
int MeadInterface::SetupAtoms(Topology const& topIn, Frame const& frameIn, Radii_Mode radiiMode)
{
  t_total_[0].Start();
  rmode_ = radiiMode;
  // Sanity checking
  if (topIn.Natom() != frameIn.Natom()) {
    mprinterr("Internal Error: MeadInterface::SetupAtoms(): Top '%s' has %i atoms, frame has %i atoms.\n",
              topIn.c_str(), topIn.Natom(), frameIn.Natom());
    return 1;
  }

  if (atomset_ != 0) delete atomset_;
  atomset_ = new AtomSet();

  bool has_radii = false;

  for (int aidx = 0; aidx < topIn.Natom(); aidx++)
  {
    MEAD::Atom at;
    set_at_from_top(at, topIn, frameIn, aidx);
    if (at.rad > 0)
      has_radii = true;
    try {
      atomset_->insert( at );
    }
    catch (MEADexcept& e) {
      return ERR("SetupAtoms()", e);
    }
  }
  if (!has_radii) {
    mprinterr("Error: No radii set for topology '%s'\n", topIn.c_str());
    return 1;
  }
  t_total_[0].Stop();
  
  return 0;
}

/** Print debug info. */
void MeadInterface::Print() const {
  //std::cout << *fdm_;
  //std::cout << *mgm_;
}

/** Set verbosity of underlying MEAD library. */
void MeadInterface::MeadVerbosity(int i) const {
  blab1pt = &cnull;
  blab2pt = &cnull;
  blab3pt = &cnull;
  if (i >= 1)
    blab1pt = &std::cout;
  if (i >= 2)
    blab2pt = &std::cout;
  if (i >= 3)
    blab3pt = &std::cout;
  if (i != 0)
    mprintf("Info: MEAD verbosity set to %i\n", i);
}
// -----------------------------------------------------------------------------
/** Class for holding calculation info for a titratable site. */
class MeadInterface::TitrationCalc {
  public:
    /// CONSTRUCTOR
    TitrationCalc(Cpptraj::Structure::TitratableSite const* sd, AtomSet const& state1Atoms, AtomSet const& state2Atoms,
                  int ridxIn, const double* xyz) :
      siteData_(sd), charge_state1_(state1Atoms), charge_state2_(state2Atoms), ridx_(ridxIn)
    {
      siteOfInterest_.x = xyz[0];
      siteOfInterest_.y = xyz[1];
      siteOfInterest_.z = xyz[2];

    }

    AtomChargeSet const& ChargeState1() const { return charge_state1_; }
    AtomChargeSet const& ChargeState2() const { return charge_state2_; }
    Coord const& SiteOfInterest()       const { return siteOfInterest_; }
    int Ridx()                          const { return ridx_; }
    Cpptraj::Structure::TitratableSite const& SiteInfo() const { return *siteData_; }
    AtomChargeSet const* RefStatePtr()  const { 
      if (siteData_->RefStateIdx() == 0)
        return &charge_state1_;
      else
        return &charge_state2_;
    }
  private:
    Cpptraj::Structure::TitratableSite const* siteData_; ///< Pointer to associated TitratableSite data
    AtomChargeSet charge_state1_;   ///< Hold charges in state 1
    AtomChargeSet charge_state2_;   ///< Hold charges in state 2
    int ridx_;                      ///< Residue index in associated Topology
    Coord siteOfInterest_;          ///< Coordinates of the site of interest, used to focus grid
};

/** For debugging - print atom potential and atom charge set. */
void MeadInterface::printAtomPotentials(Topology const& topIn, Frame const& frameIn, OutPotat* outpotat, AtomChargeSet* acs)
const
{
  double sum = 0;
  for (int aidx = 0; aidx != topIn.Natom(); ++aidx) {
    MEAD::Atom at;
    set_at_from_top(at, topIn, frameIn, aidx); // what a kludge, should be easier to access potat values
    double potential_at_atom = (double) (*outpotat)[at];
    MEAD::Atom const& at_from_acs = (*acs)[at];
    //mprintf("\t  Potential at atom %6s is %f charge %f\n", *(topIn[aidx].Name()), potential_at_atom, topIn[aidx].Charge());
    mprintf("\t  Potential at atom %6s is %f charge %f\n", *(topIn[aidx].Name()), at_from_acs.charge);
    //sum += potential_at_atom * topIn[aidx].Charge();
    sum += potential_at_atom * at_from_acs.charge;
  }
  mprintf("\tProduct of potentials with charges: %g\n", sum);
}

/** \return residue index of atom bonded to given atom in different residue or -1. */
static inline int other_res_index(Topology const& topIn, Atom const& thisAtom)
{
  int idx = -1;
  for (Atom::bond_iterator bat = thisAtom.bondbegin(); bat != thisAtom.bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() != thisAtom.ResNum()) {
      idx = topIn[*bat].ResNum();
      break;
    }
  }
  return idx;
}

/** Warning for atom not found in residue. */
static inline void warn_atNotFound(const char* at, NameType const& aname, Topology const& topIn, int ridx) {
  mprintf("Warning: %s atom '%s' not found in residue %s.\n", at, *aname, topIn.TruncResNameNum(ridx).c_str());
}

/** Error for atom not found in residue. */
static inline int err_atNotFound(const char* at, NameType const& aname, Topology const& topIn, int ridx) {
  mprinterr("Error: %s atom '%s' not found in residue %s.\n", at, *aname, topIn.TruncResNameNum(ridx).c_str());
  return 1;
}

/** \return charge from atom in AtomChargeSet. */
static inline float acs_charge(AtomChargeSet const& ref_atp, Topology const& topIn, int aidx, int ridx) {
  return ref_atp[AtomID(topIn.Res(ridx).OriginalResNum(), topIn[aidx].Name().Truncated())].charge; // TODO chainid?
}

/** Create model compounds within protein.
  * The model compound contains all atoms of the residue containing the site
  * of interest, along with the peptide C=O of the previous residue and
  * the N-H and CA of the following residue (Bashford & Karplus, 1990).
  * The background term uses the neutral state charges in ref_atp, with
  * charges being for atoms in the titrating residue being zeroed after
  * this routine.
  */
int MeadInterface::createModelCompounds(AtomChargeSet& model_compound, AtomChargeSet& model_back, AtomChargeSet const& ref_atp, int ridx, Topology const& topIn, Frame const& frameIn)
const
{
  Residue const& thisRes = topIn.Res(ridx);
  // TODO make these options
  static NameType Nname("N");
  static NameType Cname("C");
  static NameType Oname("O");
  static NameType Hname("H");
  static NameType CAname("CA");

  // Get index of the previous residue
  int Nidx = topIn.FindAtomInResidue(ridx, Nname);
  if (Nidx < 0) return err_atNotFound("N", Nname, topIn, ridx);
  int prevRidx = other_res_index(topIn, topIn[Nidx]);
# ifdef DEBUG_CPPTRAJ_MEAD
  mprintf("Previous residue index = %i\n", prevRidx + 1 );
# endif
  // Get index of next residue
  int Cidx = topIn.FindAtomInResidue(ridx, Cname);
  if (Cidx < 0) return err_atNotFound("C", Cname, topIn, ridx);
  int nextRidx = other_res_index(topIn, topIn[Cidx]);
# ifdef DEBUG_CPPTRAJ_MEAD
  mprintf("Next residue index = %i\n", nextRidx + 1);
# endif
  // Insert C and O from previous residue
  if (prevRidx > -1) {
    MEAD::Atom at;
    // Get the C and O atoms of the previous residue
    int p_Cidx = topIn.FindAtomInResidue(prevRidx, Cname);
    if (p_Cidx < 0)
      warn_atNotFound("C", Cname, topIn, prevRidx);
    else {
      set_at_from_top(at, topIn, frameIn, p_Cidx);
      model_compound.insert( at );
      at.charge = acs_charge(ref_atp, topIn, p_Cidx, prevRidx);
      model_back.insert( at );
    }
    int p_Oidx = topIn.FindAtomInResidue(prevRidx, Oname);
    if (p_Oidx < 0)
      warn_atNotFound("O", Oname, topIn, prevRidx);
    else {
      set_at_from_top(at, topIn, frameIn, p_Oidx);
      model_compound.insert( at );
      at.charge = acs_charge(ref_atp, topIn, p_Oidx, prevRidx);
      model_back.insert( at );
    }
  }

  // Insert atoms from this residue
  for (int aidx = thisRes.FirstAtom(); aidx != thisRes.LastAtom(); aidx++) {
    MEAD::Atom at;
    set_at_from_top(at, topIn, frameIn, aidx);
    model_compound.insert( at );
    at.charge = acs_charge(ref_atp, topIn, aidx, ridx);
    model_back.insert( at );
  }

  // Insert N, H, and CA from next residue
  if (nextRidx > -1) {
    MEAD::Atom at;
    int n_Nidx = topIn.FindAtomInResidue(nextRidx, Nname);
    if (n_Nidx < 0)
      warn_atNotFound("N", Nname, topIn, nextRidx);
    else {
      set_at_from_top(at, topIn, frameIn, n_Nidx);
      model_compound.insert( at );
      at.charge = acs_charge(ref_atp, topIn, n_Nidx, nextRidx);
      model_back.insert( at );
    }
    int n_Hidx = topIn.FindAtomInResidue(nextRidx, Hname);
    if (n_Hidx < 0)
     warn_atNotFound("H", Hname, topIn, nextRidx);
    else {
      set_at_from_top(at, topIn, frameIn, n_Hidx);
      model_compound.insert( at );
      at.charge = acs_charge(ref_atp, topIn, n_Hidx, nextRidx);
      model_back.insert( at );
    }
    int n_CAidx = topIn.FindAtomInResidue(nextRidx, CAname);
    if (n_CAidx < 0)
      warn_atNotFound("CA", CAname, topIn, nextRidx);
    else {
      set_at_from_top(at, topIn, frameIn, n_CAidx);
      model_compound.insert( at );
      at.charge = acs_charge(ref_atp, topIn, n_CAidx, nextRidx);
      model_back.insert( at );
    }
  }
  return 0;
}

/** Set up a single site to be calculated. */
int MeadInterface::setup_titration_site_calc(std::vector<TitrationCalc>& Sites,
                                             AtomChargeSet& ref_atp,
                                             Topology const& topIn,
                                             Frame const& frameIn,
                                             Structure::TitratableSite const& site,
                                             int ridx)
const
{
  using namespace Cpptraj::Structure;
  AtomSet state1Atoms;
  AtomSet state2Atoms;
  const double* siteOfInterest_xyz = 0;
  // Set up atoms of this site for each protonation state
  for (TitratableSite::const_iterator jt = site.begin(); jt != site.end(); ++jt)
  {
    // Get the atom index in the topology
    int aidx = topIn.FindAtomInResidue(ridx, jt->first);
    if (aidx < 0) {
      mprinterr("Error: Atom '%s' not found in residue %s\n",
                *(jt->first), topIn.TruncResNameNum(ridx).c_str());
      return 1;
    }

    // Set reference state charge for this atom in ref_atp TODO chainID for AtomID?
    MEAD::Atom& mod_at = ref_atp[AtomID(topIn.Res(ridx).OriginalResNum(), topIn[aidx].Name().Truncated())];
    if (site.RefStateIdx() == 0)
      mod_at.charge = jt->second.first;
    else
      mod_at.charge = jt->second.second;

    // Is this the site of interest? Record the coordinates if so.
    if (topIn[aidx].Name() == site.SiteOfInterest()) {
      //siteOfInterest = Vec3(frameIn.XYZ(aidx));
      siteOfInterest_xyz = frameIn.XYZ(aidx);
#     ifdef DEBUG_CPPTRAJ_MEAD
      mprintf("SITE OF INTEREST: %f %f %f\n", siteOfInterest_xyz[0], siteOfInterest_xyz[1], siteOfInterest_xyz[2]);
#     endif
      //siteOfInterest.x = xyz[0];
      //siteOfInterest.y = xyz[1];
      //siteOfInterest.z = xyz[2];
    }
    MEAD::Atom at;
    set_at_from_top(at, topIn, frameIn, aidx);
    at.charge = jt->second.first;
    state1Atoms.insert( at );
    at.charge = jt->second.second;
    state2Atoms.insert( at );
#   ifdef DEBUG_CPPTRAJ_MEAD
    mprintf("DEBUG: Atom %s idx %i charge1= %f charge2= %f\n", *(jt->first), aidx+1, jt->second.first, jt->second.second);
#   endif
  } // END loop over site atoms
  if (siteOfInterest_xyz == 0) {
    mprinterr("Error: Atom of interest %s not found in residue %s\n",
              *(site.SiteOfInterest()), topIn.TruncResNameNum(ridx).c_str());
    return 1;
  }
  Sites.push_back( TitrationCalc(&site, state1Atoms, state2Atoms, ridx, siteOfInterest_xyz) );
  return 0;
}

/** Set up sites to be calculated by original site order. */
int MeadInterface::setup_titration_calcs_by_site(std::vector<TitrationCalc>& Sites,
                                         AtomChargeSet& ref_atp,
                                         Topology const& topIn, Frame const& frameIn,
                                         Cpptraj::Structure::SiteData const& titrationData)
const
{
  using namespace Cpptraj::Structure;
  Sites.clear();
  for (SiteData::const_iterator it = titrationData.begin(); it != titrationData.end(); ++it)
  {
    int oidx = it->first; // TODO this is original resnum, need ridx
    int ridx = -1;
    for (int ires = 0; ires < topIn.Nres(); ires++) {
      if (topIn.Res(ires).OriginalResNum() == oidx) {
        ridx = ires;
        break;
      }
    }
    if (ridx < 0) {
      mprinterr("Error: Residue number %i not found.\n", oidx);
      return 1;
    }
    std::string const& siteName = it->second;
    TitratableSite const& site = titrationData.GetSite( siteName );
    if (setup_titration_site_calc(Sites, ref_atp, topIn, frameIn, site, ridx))
      return 1;
  }
  return 0;
}

/** Run multiflex calc. */
int MeadInterface::MultiFlex(MultiFlexResults& results,
                             MeadOpts const& Opts,
                             MeadGrid const& ogm, MeadGrid const& mgm,
                             Topology const& topIn, Frame const& frameIn,
                             Structure::SiteData const& titrationData, int siteIdx)
{
  t_total_[1].Start();
  using namespace Cpptraj::Structure;
  typedef std::vector<double> Darray;
  typedef std::vector<Darray> Dmatrix;
  Dmatrix SiteSiteInteractionMatrix;
  std::vector<TitrationCalc> Sites;
  // Calculate the geometric center
  Vec3 vgeom_center = frameIn.VGeometricCenter(0, frameIn.Natom());
  vgeom_center.Print("Geometric center"); // DEBUG
  Coord geom_center(vgeom_center[0], vgeom_center[1], vgeom_center[2]);

  try {
    PhysCond::set_epsext(Opts.EpsExt());
    PhysCond::set_solrad(Opts.SolRad());
    PhysCond::set_sterln(Opts.SterLn());
    PhysCond::set_ionicstr(Opts.IonicStr());

    PhysCond::print();

    // Create reference atom set with charges set to reference state charges
    // for atoms in sites of interest.
    AtomChargeSet ref_atp( *atomset_ );
    // Set up sites to calc. The charge states for each site need to be set
    // up first in order to do the site-site interactions.
    t_total_[1][0].Start();
    if (setup_titration_calcs_by_site(Sites, ref_atp, topIn, frameIn, titrationData))
    {
      mprinterr("Error: Could not set up sites to titrate.\n");
      return 1;
    }
    t_total_[1][0].Stop();
    // Allocate results
    results.AllocateSets( Sites.size() );
    // Set up site-site interaction matrix.
    SiteSiteInteractionMatrix.resize( Sites.size() );
    //Dmatrix::iterator ssi_row_it = SiteSiteInteractionMatrix.begin();
    // NOTE: In this context, *atomset_ is equivalent to atlist in multiflex.cc:FD2DielEMaker
    DielectricEnvironment_lett* eps = new TwoValueDielectricByAtoms( *atomset_, Opts.EpsIn() );
    ElectrolyteEnvironment_lett* ely = new ElectrolyteByAtoms( *atomset_ );

    // Loop over titration sites
    unsigned int ibeg, iend;
    if (siteIdx == -1) {
      // All sites
      ibeg = 0;
      iend = Sites.size();
    } else {
      // Single site
      ibeg = siteIdx;
      iend = ibeg + 1;
    }
    for (unsigned int sidx = ibeg; sidx < iend; sidx++)
    {
      TitrationCalc const& tSite = Sites[sidx];
#     ifdef DEBUG_CPPTRAJ_MEAD
      mprintf("Site %s\n", tSite.SiteInfo().SiteName().c_str());
#     endif
      Darray& ssi_row = SiteSiteInteractionMatrix[sidx];
//    for (std::vector<TitrationCalc>::const_iterator tSite = Sites.begin();
//                                                    tSite != Sites.end();
//                                                  ++tSite, ++ssi_row_it)
//    {
      // Arrays use to hold site-site interaction values
      //Darray& ssi_row = *ssi_row_it;
      ssi_row.resize(Sites.size(), 0);
      // Create model compound TODO reuse
      AtomChargeSet model_compound, model_back_chrg;
      if (createModelCompounds(model_compound, model_back_chrg, ref_atp, tSite.Ridx(), topIn, frameIn)) {
        mprinterr("Error: Creating model compound failed.\n");
        return 1;
      }
      // ----- Refocus the grid --------------
      ogm.Resolve( geom_center, tSite.SiteOfInterest() );
      // Set up charges for each state and point to the reference state
      AtomChargeSet const& charge_state1 = tSite.ChargeState1();
      AtomChargeSet const& charge_state2 = tSite.ChargeState2();
      AtomChargeSet const* refstatep = tSite.RefStatePtr();
      // mackbackX is the interaction with background charges (ref_atp)
      // EXCEPT those of site X (refstatep)
      // State1
      double macself1 = 0;
      double macback1 = 0;
      t_total_[1][1].Start();
      if (charge_state1.has_charges()) {
        // TODO check for different atoms/coords
        ChargeDist rho1(new AtomChargeSet(charge_state1));
        ElstatPot phi1(ogm.FDM(), eps, rho1, ely);
        phi1.solve();
        OutPotat* state1_pot = new OutPotat(*atomset_, phi1);
        // DEBUG Print potential at atoms
        //printAtomPotentials( topIn, frameIn, state1_pot, &ref_atp );
        macself1 = (*state1_pot) * charge_state1;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MACSELF1 = %f\n", macself1);
#       endif
        macback1 = (*state1_pot) * (ref_atp) - (*state1_pot) * (*refstatep);
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MACBACK1 = %f - %f\n", (*state1_pot) * (ref_atp), (*state1_pot) * (*refstatep));
#       endif
        // Site-site interactions
        for (unsigned int is = 0; is < Sites.size(); is++)
          ssi_row[is] = (*state1_pot) * Sites[is].ChargeState1() - (*state1_pot) * Sites[is].ChargeState2();
        delete state1_pot;
      }
      t_total_[1][1].Stop();
      //mprintf("macself1= %g  macback1= %g\n", macself1, macback1);
      // State2
      double macself2 = 0;
      double macback2 = 0;
      t_total_[1][2].Start();
      if (charge_state2.has_charges()) {
        // TODO check for different atoms/coords
        ChargeDist rho2(new AtomChargeSet(charge_state2));
        ElstatPot phi2(ogm.FDM(), eps, rho2, ely);
        phi2.solve();
        OutPotat* state2_pot = new OutPotat(*atomset_, phi2);
        macself2 = (*state2_pot) * charge_state2;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MACSELF2 = %f\n", macself2);
#       endif
        macback2 = (*state2_pot) * (ref_atp) - (*state2_pot) * (*refstatep);
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MACBACK2 = %f - %f\n", (*state2_pot) * (ref_atp), (*state2_pot) * (*refstatep));
#       endif
        // Site-site interactions
        for (unsigned int is = 0; is < Sites.size(); is++)
          ssi_row[is] = ssi_row[is] - ((*state2_pot) * Sites[is].ChargeState1() - (*state2_pot) * Sites[is].ChargeState2());
        delete state2_pot;
      }
      t_total_[1][2].Stop();
      //mprintf("macself2= %g  macback2= %g\n", macself2, macback2);
      // ----- Refocus the model grid --------
      mgm.Resolve( geom_center, tSite.SiteOfInterest() );
      // Model dielectric environment
      DielectricEnvironment_lett* model_eps = new TwoValueDielectricByAtoms( model_compound, Opts.EpsIn() );
      ElectrolyteEnvironment_lett* model_ely = new ElectrolyteByAtoms( model_compound );
      // Any atoms that are "titrating" in *this* site should be zero in
      // the background set.  This requires adjustment...
      for (AtomSet::iterator b = model_back_chrg.begin(); b!=model_back_chrg.end() ; ++b) {
        if (refstatep->contains(b->first)) (b->second).charge = 0;
      }
      // Model state 1
      double modself1 = 0;
      double modback1 = 0;
      t_total_[1][3].Start();
      if (charge_state1.has_charges()) {
        ChargeDist rho1(new AtomChargeSet(charge_state1));
        ElstatPot phi1(mgm.FDM(), model_eps, rho1, model_ely);
        phi1.solve();
        OutPotat* mod1_pot = new OutPotat(model_compound, phi1);
        modself1 = (*mod1_pot) * charge_state1;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MODSELF1 = %f\n", modself1);
#       endif
        modback1 = (*mod1_pot) * model_back_chrg;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MODBACK1 = %f\n", modback1);
#       endif
        delete mod1_pot;
      }
      t_total_[1][3].Stop();
      // Model state 2
      double modself2 = 0;
      double modback2 = 0;
      t_total_[1][4].Start();
      if (charge_state2.has_charges()) {
        ChargeDist rho2(new AtomChargeSet(charge_state2));
        ElstatPot phi2(mgm.FDM(), model_eps, rho2, model_ely);
        phi2.solve();
        OutPotat* mod2_pot = new OutPotat(model_compound, phi2);
        modself2 = (*mod2_pot) * charge_state2;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MODSELF2 = %f\n", modself2);
#       endif
        modback2 = (*mod2_pot) * model_back_chrg;
#       ifdef DEBUG_CPPTRAJ_MEAD
        mprintf("DEBUG: MODBACK2 = %f\n", modback2);
#       endif
      }
      t_total_[1][4].Stop();
      // TODO use Constants instead of PhysCond
      double delta_pK_self = -(macself1-macself2-modself1+modself2)/2.0 / PhysCond::get_ln10kT();
      double delta_pK_back = -(macback1-macback2-modback1+modback2)     / PhysCond::get_ln10kT();
#     ifdef DEBUG_CPPTRAJ_MEAD
      mprintf("macself1-macself2 = %g\n", macself1 - macself2);
      mprintf("macback1-deprotback = %g\n", macback1 - macback2);
      mprintf("modself1-modself2 = %g\n", modself1 - modself2);
      mprintf("modback1-modback2 = %g\n", modback1 - modback2);
      mprintf("delta self = %g or %g pK units\n", delta_pK_self * PhysCond::get_ln10kT(), delta_pK_self);
      mprintf("delta back = %g or %g pK units\n", delta_pK_back * PhysCond::get_ln10kT(), delta_pK_back);
      mprintf("Site-site interactions:\n");
      for (Darray::const_iterator it = ssi_row.begin(); it != ssi_row.end(); ++it)
        mprintf("   %g\n", *it);
#     endif
      double pKint = tSite.SiteInfo().pKa() + (delta_pK_self + delta_pK_back);
#     ifdef DEBUG_CPPTRAJ_MEAD
      mprintf("DEBUG: pKint = %f\n", pKint);
#     endif
      // Add site result
      results.AddSiteResult(sidx,
                            tSite.SiteInfo().SiteName(),
                            topIn.Res(tSite.Ridx()).OriginalResNum(),
                            pKint,
                            tSite.SiteInfo().RefStateIdx(),
                            delta_pK_self,
                            delta_pK_back);
    } // END loop over sites

  }
  catch (MEADexcept& e) {
    return ERR("MultiFlex()", e);
  }

  // Symmetrize the interaction matrix
  for (unsigned int i = 0; i < Sites.size(); i++) {
    if (!SiteSiteInteractionMatrix[i].empty()) {
      // Zero the self interaction
      SiteSiteInteractionMatrix[i][i] = 0;
      for (unsigned int j = i + 1; j < Sites.size(); j++) {
        if (!SiteSiteInteractionMatrix[j].empty()) {
          double ave = (SiteSiteInteractionMatrix[i][j] + SiteSiteInteractionMatrix[j][i]) / 2.0;
          double dev;
          if (ave == 0) {
            mprintf("Warning: Average interaction for site %u to site %u is 0\n", i, j);
            dev = 0;
          } else {
            dev = (SiteSiteInteractionMatrix[i][j] - SiteSiteInteractionMatrix[j][i]) / ave;
            if (dev < 0) dev = -dev;
          }
#         ifdef DEBUG_CPPTRAJ_MEAD
          mprintf("%u %u   %e %e   dev = %e\n", i+1, j+1,
                  SiteSiteInteractionMatrix[i][j], SiteSiteInteractionMatrix[j][i], dev);
#         endif
          SiteSiteInteractionMatrix[i][j] = SiteSiteInteractionMatrix[j][i] = ave;
        }
      } // END loop over j
    }
  } // END loop over i
  // Add upper-triangle matrix to results
  results.AddSiteSiteMatrix( SiteSiteInteractionMatrix );

  if (!Sites.empty()) {
    // Write pKint
    static const char pkchar[] = {'A', 'C'};
    DataSet_1D const& PK = static_cast<DataSet_1D const&>( *(results.PkIntSet()) );
    DataSet_string const& SN = static_cast<DataSet_string const&>( *(results.SiteNamesSet()) );
    for (unsigned int idx = 0; idx != results.SiteIndices().size(); idx++) {
      int siteIdx = results.SiteIndices()[idx];
      results.PkIntFile()->Printf("%e %c %s\n", PK.Dval(idx), pkchar[Sites[siteIdx].SiteInfo().RefStateIdx()],
                                  SN[idx].c_str());
              //site->SiteInfo().SiteName().c_str(), topIn.Res(site->Ridx()).OriginalResNum());
    }
    // Write summ file
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *(results.Delta_pK_SelfSet()) );
    DataSet_1D const& DB = static_cast<DataSet_1D const&>( *(results.Delta_pK_BackSet()) );
    results.SummFile()->Printf("   site name           pKmod      delta self    delta back      pkint\n");
    for (unsigned int idx = 0; idx != results.SiteIndices().size(); idx++) {
      int siteIdx = results.SiteIndices()[idx];
      results.SummFile()->Printf(" %12s %13g %13g %13g %13g\n", SN[idx].c_str(), Sites[siteIdx].SiteInfo().pKa(),
                                 DS.Dval(idx), DB.Dval(idx), PK.Dval(idx));
    }
    // Write site-site interaction file
    for (unsigned int ii = 0; ii < Sites.size(); ii++)
      for (unsigned int jj = 0; jj < SiteSiteInteractionMatrix[ii].size(); jj++)
        results.Gfile()->Printf("%u %u   %e\n", ii+1, jj+1, SiteSiteInteractionMatrix[ii][jj]);
    
  }
  t_total_[1].Stop();
  
  return 0;
}
// -----------------------------------------------------------------------------

/** Run potential calc.
  * \param values Output potential values at each coordinate in fieldPoints.
  * \param epsin Internal dielectric.
  * \param epsext External dielectric.
  * \param fieldPoints Coordinates to evaluate the potential at.
  */
int MeadInterface::Potential(DataSet_Vector_Scalar& values, MeadOpts const& Opts,
                             MeadGrid const& ogm, std::vector<Vec3> const& fieldPoints)
const
{
  values.reset();
  values.Allocate( DataSet::SizeArray(1, fieldPoints.size()) );
  try {
    PhysCond::set_epsext(Opts.EpsExt());
    ChargeDist_lett* prho = new AtomChargeSet( *atomset_ );
    DielectricEnvironment_lett* peps = new TwoValueDielectricByAtoms( *atomset_, Opts.EpsIn() );
    ElectrolyteEnvironment_lett* pely = new ElectrolyteByAtoms( *atomset_ );
    //AtomChargeSet* prho = new AtomChargeSet( *atomset_ );
    //TwoValueDielectricByAtoms* peps = new TwoValueDielectricByAtoms( *atomset_, epsin );
    //ElectrolyteByAtoms* pely = new ElectrolyteByAtoms( *atomset_ );

    FinDiffElstatPot phi( ogm.FDM(), peps, prho, pely );
    /*
    if (initfield.length()) {
    string ffn = initfield;
    ffn += ".fld"; 
    const char *ffnc = ffn.c_str();
    // open it just as an existence test
    ifstream initfield_file_existence(ffnc);
    if (initfield_file_existence.good()) {
      initfield_file_existence.close();
      phi.solve_using_coarse_init(initfield);
    }
    }
    else
    */
    phi.solve();
    /*if (outfield.length() != 0)
    phi.write_coarse_field(outfield);*/
    if (!fieldPoints.empty()) {
      for (std::vector<Vec3>::const_iterator fpt = fieldPoints.begin(); 
                                             fpt != fieldPoints.end(); ++fpt)
      {
        Coord cxyz;
        cxyz.x = (*fpt)[0];
        cxyz.y = (*fpt)[1];
        cxyz.z = (*fpt)[2];
        //std::cout << phi.value(cxyz) << "\n"; // DEBUG
        values.AddElement(*fpt, phi.value(cxyz) );
        //mprintf("DEBUG: %g\n", values.LastVal());
      }
    }

    //delete prho; // FIXME these deletes cause segfaults within MEAD
    //delete peps;
    //delete pely;
  }
  catch (MEADexcept& e) {
    return ERR("Potential()", e);
  }

  return 0;
}

/** Solvate calculation. 
  * \param Output Born solvation energy in kcal/mol
  * \param epsIn Dielectric constant of molecular interior.
  * \param epsSol Dielectric constant of solvent.
  * \param epsVac Dielectric constant of vacuum.
  * \param solRad Solvent probe radius used in rolling ball procedure to determine contact surface,
  *               boundary between epsin and epsext.
  * \param sterln Ion exclusion layer thickness added to atomic radii to determine region inaccessible
  *               to salt so kappa in PB eq is zero.
  * \param ionicStr Ionic strength (mol/L)
  * \param temperature Temperature in Kelvin
  * \param rxnField If not null, calculate reaction field for given grid at all grid points.
  */
int MeadInterface::Solvate(double& Esolv, MeadOpts const& Opts,
                           MeadGrid const& ogm, DataSet_3D* rxnField)
const
{
  Esolv = 0;

  try {
    PhysCond::set_epsext(Opts.EpsExt());
    PhysCond::set_solrad(Opts.SolRad());
    PhysCond::set_sterln(Opts.SterLn());
    PhysCond::set_ionicstr(Opts.IonicStr());
    PhysCond::set_T(Opts.Temperature());

    //mprintf("DEBUG: Interior dielectric: %g\n", Opts.EpsIn());
    //mprintf("DEBUG: Physical conditions:\n");
    PhysCond::print();
    //mprintf("DEBUG: Vacuum dielectric: %g\n", Opts.EpsVac());

    ChargeDist rho(new AtomChargeSet(*atomset_));

    // Solvent
    DielectricEnvironment eps(new TwoValueDielectricByAtoms(*atomset_, Opts.EpsIn()));
    ElectrolyteEnvironment ely(new ElectrolyteByAtoms(*atomset_));
    ElstatPot phi(ogm.FDM(), eps, rho, ely);
    phi.solve();
    float prod_sol = phi * rho;
    //mprintf("DEBUG: prod_sol= %f\n", prod_sol);

    // Vacuum
    PhysCond::set_epsext(Opts.EpsVac());
    PhysCond::set_ionicstr(0.0);
    DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(*atomset_, Opts.EpsIn()));
    ElectrolyteEnvironment elyvac;  // No electrolyte is the default
    ElstatPot vac_phi(ogm.FDM(), vac_eps, rho, elyvac);
    vac_phi.solve();
    float prod_vac = vac_phi * rho;
    //mprintf("DEBUG: prod_vac= %f\n", prod_vac);

    // Total
    float solvation_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
    //mprintf("DEBUG: SOLVATION ENERGY = %f\n", solvation_energy);
    Esolv = solvation_energy;

    // Reaction field
    if (rxnField != 0) {
      for (unsigned int ix = 0; ix < rxnField->NX(); ix++) {
        for (unsigned int iy = 0; iy < rxnField->NY(); iy++) {
          for (unsigned int iz = 0; iz < rxnField->NZ(); iz++) {
            // Get xyz coord
            //Vec3 vxyz = rxnField->Bin().Center(ix, iy, iz);
            Vec3 vxyz = rxnField->Bin().Corner(ix, iy, iz);
            Coord cxyz;
            cxyz.x = vxyz[0];
            cxyz.y = vxyz[1];
            cxyz.z = vxyz[2];
            double val = (double)(phi.value(cxyz) - vac_phi.value(cxyz));
            long int idx = rxnField->CalcIndex(ix, iy, iz);
            //mprintf("DEBUG: rxnField %i %f %f %f = %f\n", idx, vxyz[0], vxyz[1], vxyz[2], val);
            //mprintf("DBG: (%f, %f, %f)\n", vxyz[0], vxyz[1], vxyz[2]);
            rxnField->UpdateVoxel(idx, val);
          }
        }
      }
    }
  }
  catch (MEADexcept& e) {
    return ERR("Solvate()", e);
  }

  return 0;
}
