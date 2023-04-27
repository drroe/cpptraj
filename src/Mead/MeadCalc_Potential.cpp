#include "MeadCalc_Potential.h"
#include "MeadError.h"
#include "MeadOpts.h"
#include "MeadGrid.h"
#include "../DataSet_Vector_Scalar.h"
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

MeadCalc_Potential::MeadCalc_Potential() {}

/** Run potential calc.
  * \param values Output potential values at each coordinate in fieldPoints.
  * \param epsin Internal dielectric.
  * \param epsext External dielectric.
  * \param fieldPoints Coordinates to evaluate the potential at.
  */
int MeadCalc_Potential::Potential(DataSet_Vector_Scalar& values, MeadOpts const& Opts,
                             MeadGrid const& ogm, std::vector<Vec3> const& fieldPoints)
const
{
  values.reset();
  values.Allocate( DataSet::SizeArray(1, fieldPoints.size()) );
  try {
    PhysCond::set_epsext(Opts.EpsExt());
    ChargeDist_lett* prho = new AtomChargeSet( InternalAtomset() );
    DielectricEnvironment_lett* peps = new TwoValueDielectricByAtoms( InternalAtomset(), Opts.EpsIn() );
    ElectrolyteEnvironment_lett* pely = new ElectrolyteByAtoms( InternalAtomset() );
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
