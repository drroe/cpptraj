#include "MeadInterface.h"
#include "Vec3.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "DataSet_Vector_Scalar.h"
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
    atomset_->insert( at );
  }

  return 0;
}

/** Print debug info. */
void MeadInterface::Print() const {
  std::cout << *fdm_;
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
