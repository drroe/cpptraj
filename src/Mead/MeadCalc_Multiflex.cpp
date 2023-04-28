#include "MeadCalc_Multiflex.h"
#include "MultiFlexResults.h"
#include "MeadOpts.h"
#include "MeadGrid.h"
#include "MeadError.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_string.h"
#include "../Structure/SiteData.h"
#include "../Structure/TitratableSite.h"
// MEAD includes
#include "../../mead/PhysCond.h"
#include "../../mead/ChargeDist.h"
#include "../../mead/AtomChargeSet.h"
//#include "../../mead/DielectricEnvironment.h"
#include "../../mead/DielByAtoms.h"
//#include "../../mead/ElectrolyteEnvironment.h"
#include "../../mead/ElectrolyteByAtoms.h"
#include "../../mead/FinDiffElstatPot.h"
#include "../../mead/MEADexcept.h"
#include "../../mead/Potat.h"

using namespace Cpptraj::Mead;

/** CONSTRUCTOR */
MeadCalc_Multiflex::MeadCalc_Multiflex() :
  t_total_("Multiflex Total")
{
  t_total_.AddSubTimer(Timer("Setup sites")); // 0
  t_total_.AddSubTimer(Timer("MAC1       ")); // 1
  t_total_.AddSubTimer(Timer("MAC2       ")); // 2
  t_total_.AddSubTimer(Timer("MOD1       ")); // 3
  t_total_.AddSubTimer(Timer("MOD2       ")); // 4
}

// -----------------------------------------------------------------------------
/** Class for holding calculation info for a titratable site. */
class MeadCalc_Multiflex::TitrationCalc {
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
/*
void MeadCalc_Multiflex::printAtomPotentials(Topology const& topIn, Frame const& frameIn, OutPotat* outpotat, AtomChargeSet* acs)
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
}*/

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

/** Create model compounds within protein.
  * The model compound contains all atoms of the residue containing the site
  * of interest, along with the peptide C=O of the previous residue and
  * the N-H and CA of the following residue (Bashford & Karplus, 1990).
  * The background term uses the neutral state charges in ref_atp, with
  * charges being for atoms in the titrating residue being zeroed after
  * this routine.
  */
int MeadCalc_Multiflex::createModelCompounds(AtomChargeSet& model_compound, AtomChargeSet& model_back, AtomChargeSet const& ref_atp, int ridx, Topology const& topIn)
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
    // Get the C and O atoms of the previous residue
    int p_Cidx = topIn.FindAtomInResidue(prevRidx, Cname);
    if (p_Cidx < 0)
      warn_atNotFound("C", Cname, topIn, prevRidx);
    else {
      AtomID c_id(prevRidx, "C");
      model_compound.insert( InternalAtomset()[c_id] );
      model_back.insert( ref_atp[c_id] );
    }
    int p_Oidx = topIn.FindAtomInResidue(prevRidx, Oname);
    if (p_Oidx < 0)
      warn_atNotFound("O", Oname, topIn, prevRidx);
    else {
      AtomID o_id(prevRidx, "O"); 
      model_compound.insert( InternalAtomset()[o_id] );
      model_back.insert( ref_atp[o_id] );
    }
  }

  // Insert atoms from this residue
  for (int aidx = thisRes.FirstAtom(); aidx != thisRes.LastAtom(); aidx++) {
    AtomID atid(ridx, topIn[aidx].Name().Truncated());
    model_compound.insert( InternalAtomset()[atid] );
    model_back.insert( ref_atp[atid] );
  }

  // Insert N, H, and CA from next residue
  if (nextRidx > -1) {
    int n_Nidx = topIn.FindAtomInResidue(nextRidx, Nname);
    if (n_Nidx < 0)
      warn_atNotFound("N", Nname, topIn, nextRidx);
    else {
      AtomID n_id(nextRidx, "N");
      model_compound.insert( InternalAtomset()[n_id] );
      model_back.insert( ref_atp[n_id] );
    }
    int n_Hidx = topIn.FindAtomInResidue(nextRidx, Hname);
    if (n_Hidx < 0)
     warn_atNotFound("H", Hname, topIn, nextRidx);
    else {
      AtomID h_id(nextRidx, "H");
      model_compound.insert( InternalAtomset()[h_id] );
      model_back.insert( ref_atp[h_id] );
    }
    int n_CAidx = topIn.FindAtomInResidue(nextRidx, CAname);
    if (n_CAidx < 0)
      warn_atNotFound("CA", CAname, topIn, nextRidx);
    else {
      AtomID ca_id(nextRidx, "CA");
      model_compound.insert( InternalAtomset()[ca_id] );
      model_back.insert( ref_atp[ca_id] );
    }
  }
  return 0;
}

/** Set up a single site to be calculated. */
int MeadCalc_Multiflex::setup_titration_site_calc(std::vector<TitrationCalc>& Sites,
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
    // Set up MEAD AtomID
    AtomID atid(ridx, topIn[aidx].Name().Truncated());
    // Set reference state charge for this atom in ref_atp
    MEAD::Atom& mod_at = ref_atp[atid];
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
    MEAD::Atom at = InternalAtomset()[atid];
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
int MeadCalc_Multiflex::setup_titration_calcs_by_site(std::vector<TitrationCalc>& Sites,
                                         AtomChargeSet& ref_atp,
                                         Topology const& topIn, Frame const& frameIn,
                                         Cpptraj::Structure::SiteData const& titrationData)
const
{
  using namespace Cpptraj::Structure;
  Sites.clear();
  for (SiteData::const_iterator it = titrationData.begin(); it != titrationData.end(); ++it)
  {
    int ridx = it->first;
    if (ridx < 0 ) {
      mprinterr("Error: Residue number %i not found.\n", ridx);
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
int MeadCalc_Multiflex::MultiFlex(MultiFlexResults& results,
                             MeadOpts const& Opts,
                             MeadGrid const& ogm, MeadGrid const& mgm,
                             Topology const& topIn, Frame const& frameIn,
                             Structure::SiteData const& titrationData, int siteIdx)
{
  t_total_.Start();
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
    AtomChargeSet ref_atp( InternalAtomset() );
    // Set up sites to calc. The charge states for each site need to be set
    // up first in order to do the site-site interactions.
    t_total_[0].Start();
    if (setup_titration_calcs_by_site(Sites, ref_atp, topIn, frameIn, titrationData))
    {
      mprinterr("Error: Could not set up sites to titrate.\n");
      return 1;
    }
    t_total_[0].Stop();
    // Allocate results
    results.AllocateSets( Sites.size() );
    // Set up site-site interaction matrix.
    SiteSiteInteractionMatrix.resize( Sites.size() );
    //Dmatrix::iterator ssi_row_it = SiteSiteInteractionMatrix.begin();
    // NOTE: In this context, *atomset_ is equivalent to atlist in multiflex.cc:FD2DielEMaker
    DielectricEnvironment_lett* eps = new TwoValueDielectricByAtoms( InternalAtomset(), Opts.EpsIn() );
    ElectrolyteEnvironment_lett* ely = new ElectrolyteByAtoms( InternalAtomset() );

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
      if (createModelCompounds(model_compound, model_back_chrg, ref_atp, tSite.Ridx(), topIn)) {
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
      t_total_.Start();
      if (charge_state1.has_charges()) {
        // TODO check for different atoms/coords
        ChargeDist rho1(new AtomChargeSet(charge_state1));
        ElstatPot phi1(ogm.FDM(), eps, rho1, ely);
        phi1.solve();
        OutPotat* state1_pot = new OutPotat(InternalAtomset(), phi1);
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
      t_total_.Stop();
      //mprintf("macself1= %g  macback1= %g\n", macself1, macback1);
      // State2
      double macself2 = 0;
      double macback2 = 0;
      t_total_[2].Start();
      if (charge_state2.has_charges()) {
        // TODO check for different atoms/coords
        ChargeDist rho2(new AtomChargeSet(charge_state2));
        ElstatPot phi2(ogm.FDM(), eps, rho2, ely);
        phi2.solve();
        OutPotat* state2_pot = new OutPotat(InternalAtomset(), phi2);
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
      t_total_[2].Stop();
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
      t_total_[3].Start();
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
      t_total_[3].Stop();
      // Model state 2
      double modself2 = 0;
      double modback2 = 0;
      t_total_[4].Start();
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
      t_total_[4].Stop();
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
  t_total_.Stop();
  
  return 0;
}

