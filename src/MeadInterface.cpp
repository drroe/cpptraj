#include "MeadInterface.h"
#include "Vec3.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "DataSet_Vector_Scalar.h"
#include "DataSet_3D.h"
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
// FOR DEBUG
#include <iostream>

using namespace Cpptraj;

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

/** Add a grid to the finite difference method object. */
int MeadInterface::AddGrid(int ngrd, float spc, Vec3 const& cntr)
{
  if (fdm_ == 0)
    fdm_ = new FinDiffMethod();

  try { 
    fdm_->add_level( ngrd, spc, Coord(cntr[0], cntr[1], cntr[2]) );
  }
  catch (MEADexcept& e) {
    return ERR("AddGrid()", e);
  }
  return 0;
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
            Vec3 vxyz = rxnField->Bin().Center(ix, iy, iz);
            Coord cxyz;
            cxyz.x = vxyz[0];
            cxyz.y = vxyz[1];
            cxyz.z = vxyz[2];
            double val = (double)(phi.value(cxyz) - vac_phi.value(cxyz));
            long int idx = rxnField->CalcIndex(ix, iy, iz);
            mprintf("DEBUG: rxnField %i %f %f %f = %f\n", idx, vxyz[0], vxyz[1], vxyz[2], val);
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
