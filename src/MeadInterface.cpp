#include "MeadInterface.h"
#include "Vec3.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "DataSet_Vector_Scalar.h"
#include "DataSet_3D.h"
#include "Structure/TitrationData.h"
#include "Structure/TitratableSite.h"
// MEAD includes
#include "../mead/FinDiffMethod.h"
#include "../mead/MEADexcept.h"
#include "../mead/AtomSet.h"
#include "../mead/ChargeDist.h"
#include "../mead/AtomChargeSet.h"
#include "../mead/DielectricEnvironment.h"
#include "../mead/DielByAtoms.h"
#include "../mead/ElectrolyteEnvironment.h"
#include "../mead/ElectrolyteByAtoms.h"
#include "../mead/FinDiffElstatPot.h"
#include "../mead/Potat.h"
// FOR DEBUG
#include <iostream>

using namespace Cpptraj;

/** Corresponds to enum GridCenter_Mode */
const char* MeadInterface::GridCenter_ModeStr_[] = {
  "ON_ORIGIN", "ON_CENT_OF_INTR", "ON_GEOM_CENT"
};

/** \return Char string corresponding to GridCenter_Mode */
const char* MeadInterface::GridCenter_ModeStr(GridCenter_Mode gc) {
  return GridCenter_ModeStr_[gc];
}

/** CONSTRUCTOR */
MeadInterface::MeadInterface() :
  fdm_(0),
  atomset_(0)
{ }

/** DESTRUCTOR */
MeadInterface::~MeadInterface() {
  if (fdm_ != 0) delete fdm_;
  if (atomset_ != 0) delete atomset_;
}

/** Print MEAD error message. */
int MeadInterface::ERR(const char* fxn, MEADexcept& e) {
  mprinterr("Error: MEAD error in '%s': '%s' '%s' '%s'\n", fxn,
              e.get_error1().c_str(),
              e.get_error2().c_str(),
              e.get_error3().c_str());
  return 1;
}

/** Add a grid to the finite difference method object with explicit centering. */
int MeadInterface::AddGrid(int ngrd, float spc, Vec3 const& cntr)
{
  if (fdm_ == 0)
    fdm_ = new FinDiffMethod();

  try { 
    fdm_->add_level( ngrd, spc, Coord(cntr[0], cntr[1], cntr[2]) );
  }
  catch (MEADexcept& e) {
    return ERR("AddGrid(coord)", e);
  }
  return 0;
}

/** Add a grid to the finite difference method object with centering string. */
int MeadInterface::AddGrid(int ngrd, float spc, GridCenter_Mode ctrmode)
{
  if (fdm_ == 0)
    fdm_ = new FinDiffMethod();

  CenteringStyle censtl;
  switch (ctrmode) {
    case C_ON_ORIGIN       : censtl = ON_ORIGIN; break;
    case C_ON_CENT_OF_INTR : censtl = ON_CENT_OF_INTR; break;
    case C_ON_GEOM_CENT    : censtl = ON_GEOM_CENT; break;
  }
  try {
    fdm_->add_level( ngrd, spc, censtl );
  }
  catch (MEADexcept& e) {
    return ERR("AddGrid(style)", e);
  }
  return 0;
}

/** Set MEAD Atom from Topology Atom. */
void MeadInterface::set_at_from_top(MEAD::Atom& at, Topology const& topIn, Frame const& frameIn, int aidx, Radii_Mode radiiMode)
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
    switch (radiiMode) {
      case MeadInterface::GB : at.rad = thisAtom.GBRadius(); break;
      case MeadInterface::PARSE : at.rad = thisAtom.ParseRadius(); break;
      case MeadInterface::VDW   : at.rad = topIn.GetVDWradius(aidx); break;
    }
}

/** Setup an AtomSet from Frame and Topology. */
int MeadInterface::SetupAtoms(Topology const& topIn, Frame const& frameIn, Radii_Mode radiiMode)
{
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

    /*Atom const& thisAtom = topIn[aidx];
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
    switch (radiiMode) {
      case MeadInterface::GB : at.rad = thisAtom.GBRadius(); break;
      case MeadInterface::PARSE : at.rad = thisAtom.ParseRadius(); break;
      case MeadInterface::VDW   : at.rad = topIn.GetVDWradius(aidx); break;
    }*/
    set_at_from_top(at, topIn, frameIn, aidx, radiiMode);
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
  
  return 0;
}

/** Print debug info. */
void MeadInterface::Print() const {
  std::cout << *fdm_;
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


/** Run multiflex calc. */
int MeadInterface::MultiFlex(double epsIn, double epsSol,
                             double solRad, double sterln, double ionicStr,
                             Topology const& topIn, Frame const& frameIn,
                             Structure::TitrationData const& titrationData,
                             Radii_Mode radiiMode)
const
{
  using namespace Cpptraj::Structure;
  // Calculate the geometric center
  Vec3 vgeom_center = frameIn.VGeometricCenter(0, frameIn.Natom());
  vgeom_center.Print("Geometric center"); // DEBUG
  Coord geom_center(vgeom_center[0], vgeom_center[1], vgeom_center[2]);

  try {
    PhysCond::set_epsext(epsSol);
    PhysCond::set_solrad(solRad);
    PhysCond::set_sterln(sterln);
    PhysCond::set_ionicstr(ionicStr);

    PhysCond::print();
    // NOTE: In this context, *atomset_ is equivalent to atlist in multiflex.cc:FD2DielEMaker
    DielectricEnvironment_lett* eps = new TwoValueDielectricByAtoms( *atomset_, epsIn );
    ElectrolyteEnvironment_lett* ely = new ElectrolyteByAtoms( *atomset_ );

    // Loop over titratable sites
    for (int ridx = 0; ridx != topIn.Nres(); ridx++) {
      TitrationData::Sarray siteNames = titrationData.ResSiteNames( ridx );
      if (!siteNames.empty()) {
        mprintf("DEBUG: Residue %s site names:", topIn.TruncResNameNum(ridx).c_str());
        for (TitrationData::Sarray::const_iterator it = siteNames.begin();
                                                   it != siteNames.end(); ++it)
          mprintf(" %s", it->c_str());
        mprintf("\n");
        // Loop over the sites for this residue
        for (TitrationData::Sarray::const_iterator it = siteNames.begin();
                                                   it != siteNames.end(); ++it)
        {
          TitratableSite const& site = titrationData.GetSite( *it );
          AtomSet state1Atoms;
          AtomSet state2Atoms;
          Coord siteOfInterest;
          // Set up atoms for each state for the site of interest
          for (TitratableSite::const_iterator jt = site.begin(); jt != site.end(); ++jt)
          {
            // Get the atom index in the topology
            int aidx = topIn.FindAtomInResidue(ridx, jt->first);
            if (aidx < 0) {
              mprinterr("Error: Atom '%s' not found in residue %s\n",
                        *(jt->first), topIn.TruncResNameNum(ridx).c_str());
              return 1;
            }
            // Is this the site of interest?
            if (topIn[aidx].Name() == site.SiteOfInterest()) {
              //siteOfInterest = Vec3(frameIn.XYZ(aidx));
              const double* xyz = frameIn.XYZ(aidx);
              mprintf("SITE OF INTEREST: %f %f %f\n", xyz[0], xyz[1], xyz[2]);
              siteOfInterest.x = xyz[0];
              siteOfInterest.y = xyz[1];
              siteOfInterest.z = xyz[2];
            }
            MEAD::Atom at;
            set_at_from_top(at, topIn, frameIn, aidx, radiiMode);
            at.charge = jt->second.first;
            state1Atoms.insert( at );
            at.charge = jt->second.second;
            state2Atoms.insert( at );
            mprintf("DEBUG: Atom %s idx %i charge1= %f charge2= %f\n", *(jt->first), aidx+1, jt->second.first, jt->second.second);
          }
          // Refocus the grid
          fdm_->resolve( geom_center, siteOfInterest );
          // State1
          double macself1 = 0;
          double macback1 = 0;
          AtomChargeSet charge_state1(state1Atoms);
          if (charge_state1.has_charges()) {
            // TODO check for different atoms/coords
            ChargeDist rho1(new AtomChargeSet(charge_state1));
            ElstatPot phi1(*fdm_, eps, rho1, ely);
            phi1.solve();
            OutPotat* state1_pot = new OutPotat(*atomset_, phi1);
            macself1 = (*state1_pot) * charge_state1;
            delete state1_pot;
          }
          mprintf("macself1= %g\n", macself1);
          // State2
          double macself2 = 0;
          double macback2 = 0;
          AtomChargeSet charge_state2(state2Atoms);
          if (charge_state2.has_charges()) {
            // TODO check for different atoms/coords
            ChargeDist rho2(new AtomChargeSet(charge_state2));
            ElstatPot phi2(*fdm_, eps, rho2, ely);
            phi2.solve();
            OutPotat* state2_pot = new OutPotat(*atomset_, phi2);
            macself2 = (*state2_pot) * charge_state2;
            delete state2_pot;
          }
          mprintf("macself2= %g\n", macself2);
          mprintf("macself1-macself2 = %g\n", macself1 - macself2);
        }
      }
    }
  }
  catch (MEADexcept& e) {
    return ERR("MultiFlex()", e);
  }
  return 0;\
}

/** Run potential calc.
  * \param values Output potential values at each coordinate in fieldPoints.
  * \param epsin Internal dielectric.
  * \param epsext External dielectric.
  * \param fieldPoints Coordinates to evaluate the potential at.
  */
int MeadInterface::Potential(DataSet_Vector_Scalar& values, double epsin, double epsext, std::vector<Vec3> const& fieldPoints)
const
{
  values.reset();
  values.Allocate( DataSet::SizeArray(1, fieldPoints.size()) );
  try {
    PhysCond::set_epsext(epsext);
    ChargeDist_lett* prho = new AtomChargeSet( *atomset_ );
    DielectricEnvironment_lett* peps = new TwoValueDielectricByAtoms( *atomset_, epsin );
    ElectrolyteEnvironment_lett* pely = new ElectrolyteByAtoms( *atomset_ );
    //AtomChargeSet* prho = new AtomChargeSet( *atomset_ );
    //TwoValueDielectricByAtoms* peps = new TwoValueDielectricByAtoms( *atomset_, epsin );
    //ElectrolyteByAtoms* pely = new ElectrolyteByAtoms( *atomset_ );

    FinDiffElstatPot phi( *fdm_, peps, prho, pely );
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
        mprintf("DEBUG: %g\n", values.LastVal());
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
int MeadInterface::Solvate(double& Esolv,
                           double epsIn, double epsSol, double epsVac, double solRad, double sterln, double ionicStr,
                           double temperature, DataSet_3D* rxnField)
const
{
  Esolv = 0;

  try {
    PhysCond::set_epsext(epsSol);
    PhysCond::set_solrad(solRad);
    PhysCond::set_sterln(sterln);
    PhysCond::set_ionicstr(ionicStr);
    PhysCond::set_T(temperature);

    mprintf("DEBUG: Interior dielectric: %g\n", epsIn);
    mprintf("DEBUG: Physical conditions:\n");
    PhysCond::print();
    mprintf("DEBUG: Vacuum dielectric: %g\n", epsVac);

    ChargeDist rho(new AtomChargeSet(*atomset_));

    // Solvent
    DielectricEnvironment eps(new TwoValueDielectricByAtoms(*atomset_, epsIn));
    ElectrolyteEnvironment ely(new ElectrolyteByAtoms(*atomset_));
    ElstatPot phi(*fdm_, eps, rho, ely);
    phi.solve();
    float prod_sol = phi * rho;
    mprintf("DEBUG: prod_sol= %f\n", prod_sol);

    // Vacuum
    PhysCond::set_epsext(epsVac);
    PhysCond::set_ionicstr(0.0);
    DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(*atomset_, epsIn));
    ElectrolyteEnvironment elyvac;  // No electrolyte is the default
    ElstatPot vac_phi(*fdm_, vac_eps, rho, elyvac);
    vac_phi.solve();
    float prod_vac = vac_phi * rho;
    mprintf("DEBUG: prod_vac= %f\n", prod_vac);

    // Total
    float solvation_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
    mprintf("DEBUG: SOLVATION ENERGY = %f\n", solvation_energy);
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
            mprintf("DEBUG: rxnField %i %f %f %f = %f\n", idx, vxyz[0], vxyz[1], vxyz[2], val);
            mprintf("DBG: (%f, %f, %f)\n", vxyz[0], vxyz[1], vxyz[2]);
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
