#include "EnergyDecomposer.h"
#include <cmath> // sqrt for Ene_Bond etc
#include "Ene_Angle.h"
#include "Ene_Bond.h"
#include "Ene_LJ_6_12.h"
#include "Kernel_Fourier.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataFileList.h"
#include "../DataSet_Mesh.h"
#include "../DataSetList.h"
#include "../DistRoutines.h"
#include "../ParameterTypes.h"
#include "../TorsionRoutines.h"
#ifdef MPI
# include "../Stats_Reduce.h"
#endif
#include <algorithm> //std::sort

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EnergyDecomposer::EnergyDecomposer() :
  eneOut_(0),
  outfile_(0),
  debug_(0),
  currentTop_(0)
{ }

void EnergyDecomposer::HelpText() {
  mprintf("\t[<name>] [<mask>] [out <filename>]\n"
          "\t[ pme %s\n"
          "\t      %s\n"
          "\t      %s\n",
          EwaldOptions::KeywordsCommon1(),
          EwaldOptions::KeywordsCommon2(),
          EwaldOptions::KeywordsPME());
}

/** Initialize decomposer. */
int EnergyDecomposer::InitDecomposer(ArgList& argIn, DataSetList& DSLin, DataFileList& DFLin,
                                     int debugIn)
{
  debug_ = debugIn;
  // Process keywords
  outfile_ = DFLin.AddDataFile( argIn.GetStringKey("out"), argIn );
  nbcalctype_ = Ecalc_Nonbond::SIMPLE;
  if (argIn.hasKey("pme"))
    nbcalctype_ = Ecalc_Nonbond::PME;
  if (nbcalctype_ != Ecalc_Nonbond::SIMPLE) {
#   ifdef LIBPME
    ewaldOpts_.AllowLjPme(false); // TODO enable LJPME decomp
    if (ewaldOpts_.GetOptions(EwaldOptions::PME, argIn, "enedecomp"))
      return 1;
    if (ewaldOpts_.Type() == EwaldOptions::LJPME)
      nbcalctype_ = Ecalc_Nonbond::LJPME;
#   else
    mprinterr("Error: 'pme' with energy decomposition requires compilation with LIBPME.\n");
    return 1;
#   endif
  }
  // Get atom mask
  if (selectedAtoms_.SetMaskString( argIn.GetMaskNext() ))
    return 1;
  // Output DataSet
  std::string setname = argIn.GetStringNext();
  if (setname.empty())
    setname = DSLin.GenerateDefaultName("DECOMP");
  eneOut_ = DSLin.AddSet(DataSet::XYMESH, MetaData(setname));
  if (eneOut_ == 0) {
    mprinterr("Error: Could not allocate decomp. output set '%s'\n", setname.c_str());
    return 1;
  }
# ifdef MPI
  eneOut_->SetNeedsSync( false ); // Not a time series
# endif
  eneOut_->ModifyDim(Dimension::X).SetLabel("Atom");
  if (outfile_ != 0)
    outfile_->AddDataSet( eneOut_ );

  return 0;
}

/** Print options to stdout. */
void EnergyDecomposer::PrintOpts() const {
  if (eneOut_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::PrintOpts() called before initialization.\n");
    return;
  }
  mprintf("\tCalculating for atoms selected by mask: %s\n", selectedAtoms_.MaskString());
  mprintf("\tData set name: %s\n", eneOut_->legend());
  if (outfile_ != 0)
    mprintf("\tOutput file: %s\n", outfile_->DataFilename().full());
  if (nbcalctype_ != Ecalc_Nonbond::SIMPLE) {
    mprintf("\tUsing PME.\n");
    ewaldOpts_.PrintOptions();
  }
}

// -----------------------------------------------------------------------------
/** Set up bonds. */
int EnergyDecomposer::setupBonds(BndArrayType const& bondsIn) {
  for (BndArrayType::const_iterator bnd = bondsIn.begin(); bnd != bondsIn.end(); ++bnd)
  {
    if ( selectedAtoms_.AtomInCharMask( bnd->A1() ) ||
         selectedAtoms_.AtomInCharMask( bnd->A2() ) )
    {
      if (bnd->Idx() < 0) {
        mprinterr("Error: Bond %i - %i does not have parameters, cannot calculate energy.\n",
                  bnd->A1()+1, bnd->A2()+1);
        return 1;
      }
      bonds_.push_back( *bnd );
    }
  }
  return 0;
}

/** Set up angles. */
int EnergyDecomposer::setupAngles(AngArrayType const& anglesIn) {
  for (AngArrayType::const_iterator ang = anglesIn.begin(); ang != anglesIn.end(); ++ang)
  {
    if ( selectedAtoms_.AtomInCharMask( ang->A1() ) ||
         selectedAtoms_.AtomInCharMask( ang->A2() ) ||
         selectedAtoms_.AtomInCharMask( ang->A3() ) )
    {
      if (ang->Idx() < 0) {
        mprinterr("Error: Angle %i - %i - %i does not have parameters, cannot calculate energy.\n",
                  ang->A1()+1, ang->A2()+1, ang->A3()+1);
        return 1;
      }
      angles_.push_back( *ang );
    }
  }
  return 0;
}

/** Set up dihedrals. */
int EnergyDecomposer::setupDihedrals(DihArrayType const& dihedralsIn) {
  for (DihArrayType::const_iterator dih = dihedralsIn.begin(); dih != dihedralsIn.end(); ++dih)
  {
    if ( selectedAtoms_.AtomInCharMask( dih->A1() ) ||
         selectedAtoms_.AtomInCharMask( dih->A2() ) ||
         selectedAtoms_.AtomInCharMask( dih->A3() ) ||
         selectedAtoms_.AtomInCharMask( dih->A4() ) )
    {
      if (dih->Idx() < 0) {
        mprinterr("Error: Dihedral %i - %i - %i - %i does not have parameters, cannot calculate energy.\n",
                  dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
        return 1;
      }
      dihedrals_.push_back( *dih );
    }
  }
  return 0;
}

/** Topology-based setup.
  * \return 0 if setup OK, 1 if error, -1 if nothing selected.
  */
int EnergyDecomposer::SetupDecomposer(Topology const& topIn, Box const& boxIn) {
  // First set up the mask
  if (topIn.SetupCharMask( selectedAtoms_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", selectedAtoms_.MaskString());
    return 1;
  }
  selectedAtoms_.MaskInfo();
  if (selectedAtoms_.None()) {
    mprintf("Warning: Nothing selected by mask '%s'\n", selectedAtoms_.MaskString());
    return -1;
  }
  // Set up calculation arrays
  if (energies_.empty()) {
    // First time setup
//    indices_.reserve( selectedAtoms_.Nselected() );
//    for (int idx = 0; idx != topIn.Natom(); idx++)
//      if (selectedAtoms_.AtomInCharMask( idx ))
//        indices_.push_back( idx );
    energies_.resize( topIn.Natom() );
  } else {
    // Already setup. Warn if indices have changed.
    //if ((unsigned int)selectedAtoms_.Nselected() != indices_.size())
    if ((unsigned int)topIn.Natom() != energies_.size())
    {
      // FIXME implement this
      mprinterr("Error: Number of atoms has changed in topology '%s'\n", topIn.c_str());
      mprinterr("Error: Now %i atoms, expected %zu\n", topIn.Natom(), energies_.size());
      mprinterr("Error: Not yet supported by energy decomposition.\n");
      return 1;
    }
  }
  // Set up bonds
  bonds_.clear();
  if (setupBonds( topIn.Bonds() )) return 1;
  if (setupBonds( topIn.BondsH() )) return 1;
  std::sort( bonds_.begin(), bonds_.end() );
  // Set up angles
  angles_.clear();
  if (setupAngles( topIn.Angles() )) return 1;
  if (setupAngles( topIn.AnglesH() )) return 1;
  std::sort( angles_.begin(), angles_.end() );
  // Set up dihedrals
  dihedrals_.clear();
  if (setupDihedrals( topIn.Dihedrals() )) return 1;
  if (setupDihedrals( topIn.DihedralsH() )) return 1;
  std::sort( dihedrals_.begin(), dihedrals_.end() );
  // Set up nonbonds
  // Set up for all atoms FIXME
  AtomMask Imask(0, topIn.Natom());
  if (NB_.InitNonbondCalc(nbcalctype_, true, boxIn, ewaldOpts_, debug_)) {
    mprinterr("Error: Nonbond decomp init failed.\n");
    return 1;
  }
  if (NB_.SetupNonbondCalc(topIn, Imask)) {
    mprinterr("Error: Nonbond decomp setup failed.\n");
    return 1;
  }

  // DEBUG
  if (debug_ > 0) {
    mprintf("DEBUG: Saving energy for atoms:\n");
    for (int idx = 0; idx != topIn.Natom(); idx++)
      if (selectedAtoms_.AtomInCharMask( idx ))
        mprintf("\t%s\n", topIn.AtomMaskName( idx ).c_str());
    mprintf("DEBUG: Bonds:\n");
    for (BndArrayType::const_iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
      mprintf("\t%s - %s\n", topIn.AtomMaskName(bnd->A1()).c_str(), topIn.AtomMaskName(bnd->A2()).c_str());
    mprintf("DEBUG: Angles:\n");
    for (AngArrayType::const_iterator ang = angles_.begin(); ang != angles_.end(); ++ang)
      mprintf("\t%s - %s - %s\n", topIn.AtomMaskName(ang->A1()).c_str(), topIn.AtomMaskName(ang->A2()).c_str(), topIn.AtomMaskName(ang->A3()).c_str());
    mprintf("DEBUG: Dihedrals:\n");
    for (DihArrayType::const_iterator dih = dihedrals_.begin(); dih != dihedrals_.end(); ++dih)
      mprintf("\t%s - %s - %s - %s\n", topIn.AtomMaskName(dih->A1()).c_str(), topIn.AtomMaskName(dih->A2()).c_str(), topIn.AtomMaskName(dih->A3()).c_str(), topIn.AtomMaskName(dih->A4()).c_str());
  }

  currentTop_ = &topIn;

  return 0;
}

// -----------------------------------------------------------------------------
/** Save energy contribution to atom if it is selected. */
void EnergyDecomposer::saveEne(int idx, double ene_cont) {
  if (selectedAtoms_.AtomInCharMask( idx ))
    currentEne_[ idx ] += ene_cont;
}

/** Calculate bond energies. */
void EnergyDecomposer::calcBonds( Frame const& frameIn ) {
  for (BndArrayType::const_iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
  {
    BondParmType const& BP = currentTop_->BondParm()[ bnd->Idx() ];
    double ene = Ene_Bond<double>( frameIn.XYZ( bnd->A1() ),
                                   frameIn.XYZ( bnd->A2() ),
                                   BP.Req(), BP.Rk() );
#   ifdef CPPTRAJ_DEBUG_ENEDECOMP
    mprintf("DEBUG: BND %f\n", ene);
#   endif
    // Divide the energy equally between the two atoms.
    double ene_half = ene * 0.5;
    saveEne( bnd->A1(), ene_half );
    saveEne( bnd->A2(), ene_half );
  }
}

/** Calculate angle energies. */
void EnergyDecomposer::calcAngles( Frame const& frameIn ) {
  for (AngArrayType::const_iterator ang = angles_.begin(); ang != angles_.end(); ++ang)
  {
    AngleParmType const& AP = currentTop_->AngleParm()[ ang->Idx() ];
    double ene = Ene_Angle<double>( frameIn.XYZ( ang->A1() ),
                                    frameIn.XYZ( ang->A2() ),
                                    frameIn.XYZ( ang->A3() ),
                                    AP.Teq(), AP.Tk() );
#   ifdef CPPTRAJ_DEBUG_ENEDECOMP
    mprintf("DEBUG: ANG %f\n", ene);
#   endif
    // Divide the energy equally between the three atoms.
    double ene_third = ene / 3.0;
    saveEne( ang->A1(), ene_third );
    saveEne( ang->A2(), ene_third );
    saveEne( ang->A3(), ene_third );
  }
}

/** Calculate dihedral and LJ 1-4 energies. */
void EnergyDecomposer::calcDihedrals( Frame const& frameIn ) {
  for (DihArrayType::const_iterator dih = dihedrals_.begin(); dih != dihedrals_.end(); ++dih)
  {
    DihedralParmType const& DP = currentTop_->DihedralParm()[ dih->Idx() ];

    double theta = Torsion( frameIn.XYZ(dih->A1()),
                            frameIn.XYZ(dih->A2()),
                            frameIn.XYZ(dih->A3()),
                            frameIn.XYZ(dih->A4()) );
    double ene = Kernel_Fourier<double>( theta, DP.Pk(), DP.Pn(), DP.Phase() );

#   ifdef CPPTRAJ_DEBUG_ENEDECOMP
    mprintf("DEBUG: DIH %f\n", ene);
#   endif
    // Divide the energy equally between the four atoms.
    double ene_fourth = ene / 4.0;
    saveEne( dih->A1(), ene_fourth );
    saveEne( dih->A2(), ene_fourth );
    saveEne( dih->A3(), ene_fourth );
    saveEne( dih->A4(), ene_fourth );
    if (dih->Type() == DihedralType::NORMAL) {
      // 1-4 vdw energy
      double rij2 = DIST2_NoImage( frameIn.XYZ(dih->A1()), frameIn.XYZ(dih->A4()) );
      NonbondType const& LJ = currentTop_->GetLJparam(dih->A1(), dih->A4());
      double e_vdw = Ene_LJ_6_12( rij2, LJ.A(), LJ.B() );
      e_vdw /= DP.SCNB();
#     ifdef CPPTRAJ_DEBUG_ENEDECOMP
      mprintf("DEBUG: V14 %f\n", e_vdw);
#     endif
      double ene_half = e_vdw * 0.5;
      saveEne( dih->A1(), ene_half );
      saveEne( dih->A4(), ene_half );
      // 1-4 coulomb energy
      double rij = sqrt(rij2);
      double qiqj = Constants::COULOMBFACTOR * (*currentTop_)[dih->A1()].Charge() * (*currentTop_)[dih->A4()].Charge();
      double e_elec = qiqj / rij;
      e_elec /= DP.SCEE();
#     ifdef CPPTRAJ_DEBUG_ENEDECOMP
      mprintf("DEBUG: E14 %f\n", e_elec);
#     endif
      ene_half = e_elec * 0.5;
      saveEne( dih->A1(), ene_half );
      saveEne( dih->A4(), ene_half );
    }
  }
}

/** Calculate and decompose energies. */
int EnergyDecomposer::CalcEne(Frame const& frameIn) {
  t_total_.Start();
  if (currentTop_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::CalcEne() called before setup.\n");
    return 1;
  }
  currentEne_.assign( energies_.size(), 0.0 );
  // Bonds
  calcBonds(frameIn);
  // Angles
  calcAngles(frameIn);
  // Dihedrals
  calcDihedrals(frameIn);
  // Nonbonds
  double e_elec, e_vdw;
  if (NB_.DecomposedNonbondEnergy(frameIn, selectedAtoms_, e_elec, e_vdw,
                                  currentEne_, currentEne_))
  {
    mprinterr("Error: Decompose nonbond energy calc failed.\n");
    return 1;
  }

  // Accumulate the energies
  for (unsigned int idx = 0; idx != energies_.size(); idx++) {
    if (selectedAtoms_.AtomInCharMask( idx ))
      energies_[idx].accumulate( currentEne_[idx] );
  }
  t_total_.Stop();

  return 0;
}

// -----------------------------------------------------------------------------
/** Finish the calculation by putting the results into the output DataSet. */
int EnergyDecomposer::FinishCalc() {
  if (energies_.empty() || eneOut_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::FinishCalc() called before setup.\n");
    return 1;
  }
  // Only add entities that have data.
  DataSet_Mesh& set = static_cast<DataSet_Mesh&>( *eneOut_ );
  set.Clear();
  // TODO allocate?
  for (unsigned int idx = 0; idx != energies_.size(); idx++) {
    if ( energies_[idx].nData() > 0 ) {
      set.AddXY( idx+1, energies_[idx].mean() );
    }
  }
  mprintf("Timing for energy decomposition: '%s'\n", eneOut_->legend());
  t_total_.WriteTiming(0, "  Decomp total:");
  NB_.PrintTiming(t_total_.Total());
  return 0;
}

#ifdef MPI
/** Reduce the per-atom energy array to the master rank.
  * Should be called before FinishCalc().
  */
int EnergyDecomposer::ReduceToMaster(Parallel::Comm const& trajComm) {
  unsigned long maxbin;
  Stats_Reduce( trajComm, energies_, maxbin );
  return 0;
}
#endif