#ifdef HAS_MEAD
#include "MeadCalc_Solvate.h"
#include "MeadError.h"
#include "MeadGrid.h"
#include "MeadOpts.h"
#include "../DataSet_3D.h"
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

using namespace Cpptraj::Mead;

MeadCalc_Solvate::MeadCalc_Solvate() {}

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
int MeadCalc_Solvate::Solvate(double& Esolv, MeadOpts const& Opts,
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

    ChargeDist rho(new AtomChargeSet(InternalAtomset()));

    // Solvent
    DielectricEnvironment eps(new TwoValueDielectricByAtoms(InternalAtomset(), Opts.EpsIn()));
    ElectrolyteEnvironment ely(new ElectrolyteByAtoms(InternalAtomset()));
    ElstatPot phi(ogm.FDM(), eps, rho, ely);
    phi.solve();
    float prod_sol = phi * rho;
    //mprintf("DEBUG: prod_sol= %f\n", prod_sol);

    // Vacuum
    PhysCond::set_epsext(Opts.EpsVac());
    PhysCond::set_ionicstr(0.0);
    DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(InternalAtomset(), Opts.EpsIn()));
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
#endif
