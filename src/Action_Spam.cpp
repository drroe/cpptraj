// Action_Spam
#include <cmath> // sqrt
#include <cstdio> // sscanf
#include <algorithm> // std::min, std::max
#include "Action_Spam.h"
#include "Box.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
#include "DistRoutines.h"
#include "StringRoutines.h" // integerToString
#include "KDE.h"
#include "OnlineVarT.h" // Stats
#include "DataSet_Mesh.h"
#include "DataSet_Vector_Scalar.h"
#include "DataIO_Peaks.h"

// ----- SolventRes class ------------------------------------------------------
/** CONSTRUCTOR */
Action_Spam::SolventRes::SolventRes() :
  at0_(-1),
  at1_(-1),
  sidx_(-1)
{}

/** Construct from res first atom, last atom, index into solvents_ */
Action_Spam::SolventRes::SolventRes(int at0, int at1, int sidx) :
  at0_(at0),
  at1_(at1),
  sidx_(sidx)
{}

/** Print solvent res info to stdout. */
void Action_Spam::SolventRes::PrintInfo() const {
  mprintf("\t%8i - %8i (%i)\n", at0_ + 1, at1_ + 1, sidx_);
}

// ----- SolventInfo class -----------------------------------------------------
/** CONSTRUCTOR */
Action_Spam::SolventInfo::SolventInfo() :
  peaksData_(0),
  site_size_(0)
{}

/** Construct from peaks data, solvent site size, solvent name */
Action_Spam::SolventInfo::SolventInfo(DataSet_Vector_Scalar const* ds,
                                      double s, std::string const& n) :
  peaksData_(ds),
  site_size_(s),
  name_(n)
{}

/** Print solvent info to stdout. */
void Action_Spam::SolventInfo::PrintInfo() const {
  mprintf("\t%s %6.2f '%s'\n", peaksData_->legend(), site_size_, name_.c_str());
}

// ----- SolventPeak class -----------------------------------------------------
/** CONSTRUCTOR */
Action_Spam::SolventPeak::SolventPeak() :
  energies_(0)
{}

/** Construct with given output energy DataSet */
Action_Spam::SolventPeak::SolventPeak(DataSet* ds) :
  energies_(ds)
{
  mprintf("DEBUG: Constructed with set '%s'\n", ds->legend());
}

// ----- PeakSite class --------------------------------------------------------
/** For adding zero energy to omitted frames. */
const double Action_Spam::PeakSite::ZERO_ = 0.0;

/** CONSTRUCTOR */
Action_Spam::PeakSite::PeakSite() {}

/** Construct from given peak location. */
Action_Spam::PeakSite::PeakSite(Vec3 const& xyzIn) :
  xyz_(xyzIn)
{}

/** Add output energy DataSets for each solvent type in given array. */
int Action_Spam::PeakSite::AddEneDataSets(std::vector<SolventInfo> const& solvents,
                                          std::string const& dsname,
                                          DataSetList& DSL, DataFile* df, unsigned int peakIdx)
{
  solvPeaks_.clear();
  MetaData meta(dsname, peakIdx);
  for (std::vector<SolventInfo>::const_iterator solv = solvents.begin();
                                                solv != solvents.end(); ++solv)
  {
    if (solvents.size() > 1)
      meta.SetAspect( solv->Name() );
    DataSet* ds = DSL.AddSet( DataSet::DOUBLE, meta );
    if (ds == 0) {
      mprinterr("Error: AddEneDataSets failed for '%s' peak %u solvent '%s'\n",
                dsname.c_str(), peakIdx, solv->Name().c_str());
      return 1;
    }
    solvPeaks_.push_back( SolventPeak(ds) );
    if (df != 0)
      df->AddDataSet( ds );
  }
  return 0;
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Action_Spam::Action_Spam() :
  debug_(0),
  DG_BULK_(-30.3), // Free energy of bulk SPCE water
  DH_BULK_(-22.2), // Enthalpy of bulk SPCE water
  temperature_(300.0),
  purewater_(false),
  reorder_(false),
  calcEnergy_(false),
  cut2_(144.0),
  onecut2_(1.0 / 144.0),
  doublecut_(24.0),
  infofile_(0),
  site_size_(2.5),
  sphere_(false),
  bulk_ene_set_(0),
  ds_dg_(0),
  ds_dh_(0),
  ds_ds_(0),
  Nframes_(0),
  overflow_(false),
  peaksData_(0)
{ }

/** Search for DataSet with peaks data. If that fails, try to load peaks
  * from a file.
  */
DataSet_Vector_Scalar* Action_Spam::GetPeaksData(std::string const& name, DataSetList const& dsl)
{
  DataSet_Vector_Scalar* pdata = 0;
  // Check for peaks DataSet.
  DataSet* ds = dsl.FindSetOfType(name, DataSet::VECTOR_SCALAR);
  if (ds == 0) {
    // No set found. See if file exists.
    FileName fname(name);
    if (!File::Exists(fname)) {
      File::ErrorMsg( fname.full() );
      mprinterr("Error: No peak data or file with name '%s'.\n", name.c_str());
      return pdata;
    }
    // Try to load peaks from file.
    DataIO_Peaks infile;
    if (infile.ReadData(fname, peaksdsl_, fname.Base())) {
      mprinterr("Error: Could not load peaks data from %s\n", fname.full());
      return pdata;
    }
    // Sanity check
    if (peaksdsl_.size() < 1 || peaksdsl_[0]->Type() != DataSet::VECTOR_SCALAR) {
      mprinterr("Error: Could not allocate peaks data set for file.\n");
      return pdata;
    }
    pdata = (DataSet_Vector_Scalar*)peaksdsl_[0];
  } else {
    pdata = (DataSet_Vector_Scalar*)ds;
  }

  // Add each peak to peakSites_ TODO check for overlaps
  for (unsigned int idx = 0; idx != pdata->Size(); idx++)
  {
    peakSites_.push_back( PeakSite( pdata->Vec(idx) ) );
  }
  return pdata;
}

void Action_Spam::Help() const {
  mprintf("\t[name <name>] [out <datafile>] [cut <cut>] [solv <solvname>]\n"
          "\t{ purewater |\n"
          "\t  <peaksname> [reorder] [info <infofile>] [summary <summary>]\n"
          "\t  [site_size <size>] [sphere] [temperature <T>]\n"
          "\t  [dgbulk <dgbulk>] [dhbulk <dhbulk>] }\n"
          "  Perform SPAM water analysis. If 'purewater' is specified calculate\n"
          "  bulk energy values for a pure water system. Otherwise determine SPAM\n"
          "  energies from peaks previously identified from the 'volmap' action.\n"
          "    <name>      : Output data set name.\n"
          "    <datafile>  : Data file with all SPAM energies for each snapshot.\n"
          "    <cut>       : Non-bonded cutoff for energy evaluation\n"
          "    <solvname>  : Name of the solvent residues\n"
          "    [purewater] : The system is pure water---used to parametrize the bulk values.\n"
          "    <peaksname> : Dataset/File (XYZ format) with the peak locations present.\n"
          "    [reorder]   : The solvent should be re-ordered so the same solvent molecule\n"
          "                  is always in the same site.\n"
          "    <infofile>  : File with stats about which sites are occupied when.\n"
          "    <summary>   : File with the summary of all SPAM results.\n"
          "    <size>      : Size of the water site around each density peak.\n"
          "    [sphere]    : Treat each site like a sphere.\n"
          "    <T>         : Temperature at which SPAM calculation was run.\n"
          "    <dgbulk>    : SPAM free energy of the bulk solvent in kcal/mol\n"
          "    <dhbulk>    : SPAM enthalpy of the bulk solvent in kcal/mol\n");
}

// Action_Spam::Init()
Action::RetType Action_Spam::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Always use imaged distances
  imageOpt_.InitImaging(true);

  // See if we're doing pure water. If so, we don't need a peak file
  purewater_ = actionArgs.hasKey("purewater");

  // Get data set name.
  std::string ds_name = actionArgs.GetStringKey("name");
  if (ds_name.empty())
    ds_name = init.DSL().GenerateDefaultName("SPAM");

  // Get output data file
  DataFile* datafile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  DataFile* summaryfile = 0;

  // Solvent residue name
  solvname_ = actionArgs.GetStringKey("solv");
  if (solvname_.empty())
    solvname_ = std::string("WAT");

  // Get energy cutoff
  double cut = actionArgs.getKeyDouble("cut", 12.0);
  if (purewater_) {
    if (pairList_.InitPairList( cut, 0.1, debug_ )) return Action::ERR;
  }
  cut2_ = cut * cut;
  doublecut_ = 2 * cut;
  onecut2_ = 1 / cut2_;

  if (purewater_) {
    // We only have one data set averaging over every water. Add it here
    bulk_ene_set_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name));
    if (bulk_ene_set_ == 0) return Action::ERR;
    if (datafile != 0) datafile->AddDataSet( bulk_ene_set_ );
    bulk_ene_set_->ModifyDim(Dimension::X).SetLabel("Index");
    DG_BULK_ = 0.0;
    DH_BULK_ = 0.0;
    // Shouldn't need any more arguments.
    if (actionArgs.NremainingArgs() > 0)
      mprintf("Warning: 'purewater' specified but more arguments remain.\n");
    // Add bulk water info; only need name if not calculating peaks
    solvents_.push_back( SolventInfo(0, 0, solvname_) );
  } else {
    // Get the file/dataset name with the peaks defined in it
    std::string peaksname = actionArgs.GetStringNext();
    if (peaksname.empty()) {
      mprinterr("Error: No Peak dataset/file specified.\n");
      return Action::ERR;
    }
    // Get the remaining optional arguments
    reorder_ = actionArgs.hasKey("reorder");
    DG_BULK_ = actionArgs.getKeyDouble("dgbulk", -30.3);
    if (!actionArgs.Contains("dgbulk"))
      mprintf("Warning: 'dgbulk' not specified; using default for SPC/E water.\n");
    DH_BULK_ = actionArgs.getKeyDouble("dhbulk", -22.2);
    if (!actionArgs.Contains("dhbulk"))
      mprintf("Warning: 'dhbulk' not specified; using default for SPC/E water.\n");
    temperature_ = actionArgs.getKeyDouble("temperature", 300.0);
    std::string infoname = actionArgs.GetStringKey("info");
    if (infoname.empty())
      infoname = std::string("spam.info");
    infofile_ = init.DFL().AddCpptrajFile(infoname, "SPAM info");
    if (infofile_ == 0) return Action::ERR;
    summaryfile = init.DFL().AddDataFile(actionArgs.GetStringKey("summary"), actionArgs);
    // Determine if energy calculation needs to happen
    calcEnergy_ = (summaryfile != 0 || datafile != 0);
    // Divide site size by 2 to make it half the edge length (or radius)
    site_size_ = actionArgs.getKeyDouble("site_size", 2.5) / 2.0;
    sphere_ = actionArgs.hasKey("sphere");
    // If it's a sphere, square the radius to compare with
    if (sphere_)
      site_size_ *= site_size_;

    // ----- END processing input args -----------


    // Get or load the peaks data for bulk
    peaksData_ = GetPeaksData(peaksname, init.DSL());
    if (peaksData_ == 0) {
      mprinterr("Error: Could not get peaks.\n");
      return Action::ERR;
    }
    // Add bulk solvent site info TODO make these non-class vars
    solvents_.push_back( SolventInfo(peaksData_, site_size_, solvname_) );

    // DEBUG print solvents
    mprintf("DEBUG: Solvents:\n");
    for (std::vector<SolventInfo>::const_iterator it = solvents_.begin(); it != solvents_.end(); ++it)
      it->PrintInfo();

    // DEBUG print peaks
    mprintf("DEBUG: Peak sites:\n");
    for (std::vector<PeakSite>::const_iterator it = peakSites_.begin(); it != peakSites_.end(); ++it)
      mprintf("DEBUG:\t%8li %8.3f %8.3f %8.3f\n", it - peakSites_.begin(),
              it->XYZ()[0], it->XYZ()[1], it->XYZ()[2]);

    // Add DataSets for each solvent to all peak sites
    if (calcEnergy_) {
      for (std::vector<PeakSite>::iterator it = peakSites_.begin(); it != peakSites_.end(); ++it)
        if (it->AddEneDataSets(solvents_, ds_name, init.DSL(), datafile, (it-peakSites_.begin())+1))
          return Action::ERR;
    }

    // Make the peak overall energy sets Mesh so we can skip unoccupied peaks
    Dimension Pdim( 1.0, 0.0, "Peak" );
    ds_dg_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"DG"));
    ds_dh_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"DH"));
    ds_ds_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"-TDS"));
    if (ds_dg_==0 || ds_dh_==0 || ds_ds_==0) return Action::ERR;
    ds_dg_->SetDim(Dimension::X, Pdim);
    ds_dh_->SetDim(Dimension::X, Pdim);
    ds_ds_->SetDim(Dimension::X, Pdim);
    if (summaryfile != 0) {
      summaryfile->AddDataSet( ds_dg_ );
      summaryfile->AddDataSet( ds_dh_ );
      summaryfile->AddDataSet( ds_ds_ );
    }
#   ifdef MPI
    ds_dg_->SetNeedsSync(false);
    ds_dh_->SetNeedsSync(false);
    ds_ds_->SetNeedsSync(false);
#   endif
  }

  if (reorder_ && solvents_.size() > 1) {
    mprinterr("Error: 'reorder' cannot be specified with multiple solvents.\n");
    return Action::ERR;
  }

  // Set function for determining if water is inside peak
  if (sphere_)
    Inside_ = &Action_Spam::inside_sphere;
  else
    Inside_ = &Action_Spam::inside_box;

  // Print info now
  mprintf("    SPAM:\n");
  if (purewater_) {
    mprintf("\tCalculating bulk value for pure solvent\n");
    if (datafile != 0)
      mprintf("\tPrinting solvent energies to %s\n", datafile->DataFilename().full());
    mprintf("\tData set '%s' index is water # * frame.\n", bulk_ene_set_->legend());
    mprintf("\tUsing a %.2f Angstrom non-bonded cutoff with shifted EEL.\n",
            sqrt(cut2_));
    if (reorder_)
      mprintf("\tWarning: Re-ordering makes no sense for pure solvent.\n");
    if (summaryfile != 0)
      mprintf("\tPrinting solvent SPAM summary to %s\n",
               summaryfile->DataFilename().full());
  } else {
    mprintf("\tSolvent [%s], %zu density peaks taken from %s.\n",
            solvname_.c_str(), peaksData_->Size(), peaksData_->legend());
    mprintf("\tOccupation information printed to %s.\n", infofile_->Filename().full());
    mprintf("\tSites are ");
    if (sphere_)
      mprintf("spheres with diameter %.3f\n", site_size_);
    else
      mprintf("boxes with edge length %.3f\n", site_size_);
    if (reorder_)
      mprintf("\tRe-ordering trajectory so each site always has the same water molecule.\n");
    if (!calcEnergy_) {
      if (!reorder_) {
        mprinterr("Error: Not re-ordering trajectory or calculating energies. Nothing to do!\n");
        return Action::ERR;
      }
      mprintf("\tNot calculating any SPAM energies\n");
    } else {
      mprintf("\tUsing a non-bonded cutoff of %.2f Ang. with a EEL shifting function.\n",
              sqrt(cut2_));
      mprintf("\tBulk solvent SPAM free energy: %.3f kcal/mol\n", DG_BULK_);
      mprintf("\tBulk solvent SPAM enthalpy: %.3f kcal/mol\n", DH_BULK_);
      mprintf("\tTemperature: %.3f K\n", temperature_);
    }
  }
  mprintf("#Citation: Cui, G.; Swails, J.M.; Manas, E.S.; \"SPAM: A Simple Approach\n"
          "#          for Profiling Bound Water Molecules\"\n"
          "#          J. Chem. Theory Comput., 2013, 9 (12), pp 5539–5549.\n");
  return Action::OK;
}

// Action_Spam::Setup()
Action::RetType Action_Spam::Setup(ActionSetup& setup) {
  // We need box info
  Box const& currentBox = setup.CoordInfo().TrajBox();
  if (!currentBox.HasBox()) {
    mprinterr("Error: SPAM: Must have explicit solvent with periodic boundaries!\n");
    return Action::ERR;
  }

  // See if our box dimensions are too small for our cutoff...
  if (currentBox.Param(Box::X) < doublecut_ ||
      currentBox.Param(Box::Y) < doublecut_ ||
      currentBox.Param(Box::Z) < doublecut_)
  {
    mprinterr("Error: SPAM: The box appears to be too small for your cutoff!\n");
    return Action::ERR;
  }
  // Set up imaging info for this parm
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  // SANITY CHECK - imaging should always be active.
  if (!imageOpt_.ImagingEnabled()) {
    mprinterr("Error: Imaging not possible for %s; required for SPAM.\n", setup.Top().c_str());
    return Action::ERR;
  }

  // Set up solvResArray_
  solvResArray_.clear();
  for (Topology::res_iterator res = setup.Top().ResStart();
                              res != setup.Top().ResEnd(); res++)
  {
    for (int sidx = 0; sidx != (int)solvents_.size(); sidx++)
    {
      if (res->Name().Truncated() == solvents_[sidx].Name()) {
        solvResArray_.push_back( SolventRes(res->FirstAtom(), res->LastAtom(), sidx) );
        break;
      }
    }
  }
  mprintf("DEBUG: Solvent residues\n");
  for (std::vector<SolventRes>::const_iterator it = solvResArray_.begin(); it != solvResArray_.end(); ++it)
    it->PrintInfo();
  if (solvResArray_.empty()) {
    mprinterr("Error: No solvent residues.\n");
    return Action::ERR;
  }

  // Set up mask_ and watidx_ 
  mask_.ResetMask();
  int idx = 0;
  watidx_.clear();
  watidx_.reserve( setup.Top().Natom() );
  for (Topology::res_iterator res = setup.Top().ResStart();
                              res != setup.Top().ResEnd(); res++)
  {
    if (res->Name().Truncated() == solvname_) {
      for (int i = res->FirstAtom(); i < res->LastAtom(); i++) {
        mask_.AddAtom( i );
        watidx_.push_back( idx ); // TODO currently purewater only - skip if not purewater?
      }
      idx++;
    }
  }
  // Reserve space to hold assigned peak for each solvent residue
  resPeakNum_.reserve( solvResArray_.size() );

  mprintf("\tFound %zu solvent residues [%s]\n", solvResArray_.size(), // FIXME fix for multiple solvents
          solvname_.c_str());

  // Set up pair list
  if (purewater_) {
    if (pairList_.SetupPairList( currentBox )) return Action::ERR;
  }

  // Set up the charge array and check that we have enough info
  if (SetupParms(setup.Top())) return Action::ERR;

  // Save topology address so we can get NB params during energy calc. 
  CurrentParm_ = setup.TopAddress();

  return Action::OK;
}

// Action_Spam::SetupParms
/** Sets the temporary charge array and makes sure that we have the necessary
  * parameters in our topology to calculate nonbonded energy terms
  */
int Action_Spam::SetupParms(Topology const& ParmIn) {
  // Store the charges, convert to Amber style so we get energies in kcal/mol
  atom_charge_.clear();
  atom_charge_.reserve( ParmIn.Natom() );
  for (Topology::atom_iterator atom = ParmIn.begin(); atom != ParmIn.end(); ++atom)
    atom_charge_.push_back( atom->Charge() * Constants::ELECTOAMBER );
  if (!ParmIn.Nonbond().HasNonbond()) {
    mprinterr("Error: SPAM: Parm does not have LJ information.\n");
    return 1;
  }
  return 0;
}

// Action_Spam::DoAction()
Action::RetType Action_Spam::DoAction(int frameNum, ActionFrame& frm) {
  Nframes_++;
  if (imageOpt_.ImagingEnabled())
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
  // Check that our box is still big enough...
  overflow_ = overflow_ || frm.Frm().BoxCrd().Param(Box::X) < doublecut_ ||
                           frm.Frm().BoxCrd().Param(Box::Y) < doublecut_ ||
                           frm.Frm().BoxCrd().Param(Box::Z) < doublecut_;
  if (purewater_)
    return DoPureWater(frameNum, frm.Frm());
  else {
    return SpamCalc(frameNum, frm.ModifyFrm());
  }
}

/** \return Energy between atoms i and j with given distance squared.
  * \param i Absolute atom index for atom i.
  * \param j Absolute atom index for atom j.
  * \param dist2 Distance squared between atoms i and j.
  */
double Action_Spam::Ecalc(int i, int j, double dist2) const {
  double qiqj = atom_charge_[i] * atom_charge_[j];
  NonbondType const& LJ = CurrentParm_->GetLJparam(i, j);
  double r2 = 1 / dist2;
  double r6 = r2 * r2 * r2;
  // Shifted electrostatics: qiqj/r * (1-r/rcut)^2 + VDW
  double shift = (1 - dist2 * onecut2_);
  double eval = (qiqj / sqrt(dist2) * shift * shift + LJ.A() * r6 * r6 - LJ.B() * r6);
  //if (i < j) {
  //  if (i > 2 && i < 6)
  //    mprintf("DEBUG: %6i %6i %8.3f %8.3f\n", i, j, sqrt(dist2), eval);
  //} else {
  //  if (j > 2 && j < 6)
  //    mprintf("DEBUG: %6i %6i %8.3f %8.3f\n", j, i, sqrt(dist2), eval);
  //}
  return eval;
}

// Action_Spam::DoPureWater
/** Carries out SPAM analysis for pure water to parametrize bulk.
  * This is relatively simple. For each frame, calculate the interaction
  * energy for every water to every other water in the system. Therefore
  * we will have NFRAMES * NWATER data points total.
  */
Action::RetType Action_Spam::DoPureWater(int frameNum, Frame const& frameIn)
{
  t_action_.Start();
  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), mask_);
  if (retVal != 0) {
    mprinterr("Error: Grid setup failed.\n");
    return Action::ERR;
  }
  int wat = 0, wat1 = 0;
  int basenum = frameNum * solvResArray_.size();
  DataSet_double& evals = static_cast<DataSet_double&>( *bulk_ene_set_ );
  // Make room for each solvent residue energy this frame.
  evals.Resize( evals.Size() + solvResArray_.size() );
  t_energy_.Start();
  // Loop over all grid cells
  for (int cidx = 0; cidx < pairList_.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = pairList_.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        wat = watidx_[it0->Idx()];
        int atomi = mask_[it0->Idx()];
        Vec3 const& xyz0 = it0->ImageCoords();
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          wat1 = watidx_[it1->Idx()];
          if ( wat != wat1 ) {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < cut2_) {
              double eval = Ecalc(atomi, mask_[it1->Idx()], D2);
              evals[basenum + wat] += eval;
              evals[basenum + wat1] += eval;
            }
          }
        } // END loop over all other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = pairList_.Cell( cellList[nidx] );
          // Translate vector for neighbor cell
          Vec3 const& tVec = pairList_.TransVec( transList[nidx] );
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            wat1 = watidx_[it1->Idx()];
            if ( wat != wat1 ) {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double D2 = dxyz.Magnitude2();
              if (D2 < cut2_) {
                double eval = Ecalc(atomi, mask_[it1->Idx()], D2);
                evals[basenum + wat] += eval;
                evals[basenum + wat1] += eval;
              }
            }
          } // END loop over atoms in neighbor cell
        } // END loop over neighbor cells
      } // END loop over atoms in thisCell
    } // END cell not empty
  } // END loop over grid cells
  t_energy_.Stop();
  t_action_.Stop();
  return Action::OK;
}

// Action_Spam::inside_box()
bool Action_Spam::inside_box(Vec3 const& gp, Vec3 const& pt, double edge) const {
  return (gp[0] + edge > pt[0] && gp[0] - edge < pt[0] &&
          gp[1] + edge > pt[1] && gp[1] - edge < pt[1] &&
          gp[2] + edge > pt[2] && gp[2] - edge < pt[2]);
}

// Action_Spam::inside_sphere()
bool Action_Spam::inside_sphere(Vec3 const& gp, Vec3 const& pt, double rad2) const {
  return ( (gp[0]-pt[0])*(gp[0]-pt[0]) + (gp[1]-pt[1])*(gp[1]-pt[1]) +
           (gp[2]-pt[2])*(gp[2]-pt[2]) < rad2 );
}

/** Do the SPAM calculation for each solvent peak site. */
Action::RetType Action_Spam::SpamCalc(int frameNum, Frame& frameIn) {
  t_action_.Start();
  // For each solvent residue, try to find a peak that is close. If two peaks
  // are close only assign the closest peak.
  t_assign_.Start();
  resPeakNum_.assign(solvResArray_.size(), -1);
  int solvResIdx;
  int solvResMax = (int)solvResArray_.size();
  for (solvResIdx = 0; solvResIdx < solvResMax; solvResIdx++)
  {
    SolventRes const& res = solvResArray_[solvResIdx];
    // Calculate solvent center of mass
    Vec3 solvCoM = frameIn.VCenterOfMass(res.At0(), res.At1());
    // Loop over each peak, determine if solvent is close enough to be "inside"
    double closest_peak_distance2 = -1.0;
    for (std::vector<PeakSite>::const_iterator peak = peakSites_.begin();
                                               peak != peakSites_.end(); ++peak)
    {
      if ( (this->*Inside_)(peak->XYZ(), solvCoM, solvents_[res.Sidx()].SiteSize()) ) {
        // Solvent is inside peak.
        if (resPeakNum_[solvResIdx] < 0) {
          // First peak that is close enough
          closest_peak_distance2 = DIST2_NoImage(solvCoM.Dptr(), peak->XYZ().Dptr());
          resPeakNum_[solvResIdx] = (int)(peak - peakSites_.begin());
        } else {
          // See if this peak is closer than previous peak
          double dist2 = DIST2_NoImage(solvCoM.Dptr(), peak->XYZ().Dptr());
          if (dist2 < closest_peak_distance2) {
            closest_peak_distance2 = dist2;
            resPeakNum_[solvResIdx] = (int)(peak - peakSites_.begin());
          }
        }
      }
    } // END loop over peaks

  } // END loop over solvent residues

  // DEBUG - print peak assignments for each solvent
  mprintf("DEBUG: Peak assignments for each solvent idx (new):\n");
  for (Iarray::const_iterator it = resPeakNum_.begin(); it != resPeakNum_.end(); ++it)
    if (*it > -1)
      mprintf("DEBUG:\t%8li %i\n", it - resPeakNum_.begin(), *it);
  t_assign_.Stop();

  t_occupy_.Start();
  // We want to make sure that each site is occupied once and only once.
  // If a site is unoccupied, add frameNum to this peak's list of
  // omitted frames. If a site is multiply-occupied, add -frameNum to
  // the list.
  // Number of times peak is associated with a solvent residue.
  std::vector<unsigned int> numTimesPeakAssigned( peakSites_.size(), 0 );
  // Hold indices into solvResArray_ for solvent occupying a peak.
  Iarray peakResIdx( peakSites_.size(), -1 );
  for (Iarray::const_iterator peak = resPeakNum_.begin();
                              peak != resPeakNum_.end(); ++peak)
  {
    if (*peak > -1) {
      numTimesPeakAssigned[*peak]++;
      peakResIdx[*peak] = (peak-resPeakNum_.begin());
    }
  }
  // Hold indices into solvResArray_ for solvent singly-occupying a peak.
  Iarray singleOccSolvResIdx;
  // Hold indices into peakSites_ for singly-occupied peaks.
  Iarray singleOccPeakIdx;
  for (unsigned int idx = 0; idx != peakSites_.size(); ++idx)
  {
    if (numTimesPeakAssigned[idx] == 1) {
      // Singly occupied peak
      singleOccSolvResIdx.push_back( peakResIdx[idx] );
      singleOccPeakIdx.push_back( idx );
    } else {
      // No occupancy or multiple occupancy - omit 
      peakSites_[idx].AddOmittedFrame( frameNum, numTimesPeakAssigned[idx] );
    }
  }
  // DEBUG - print peak assignment stats
  mprintf("DEBUG: Peak assignment stats:\n");
  for (unsigned int idx = 0; idx != peakSites_.size(); ++idx)
    if (numTimesPeakAssigned[idx] > 0)
      mprintf("DEBUG:\t%8u %u (%i)\n", idx, numTimesPeakAssigned[idx], peakResIdx[idx]);
  mprintf("DEBUG: Singly-occupied peaks:\n");
  for (Iarray::const_iterator it = singleOccSolvResIdx.begin();
                              it != singleOccSolvResIdx.end(); ++it)
    mprintf("DEBUG:\t%8i %8i - %8i\n", *it,
            solvResArray_[*it].At0()+1, solvResArray_[*it].At1()+1);
  

  t_occupy_.Stop();

  // Energy calculation
  if (calcEnergy_) {
    t_energy_.Start();
    // Energy associated with each peak
    std::vector<double> singleOccPeakEne( singleOccPeakIdx.size(), 0 );
    for (int atom0 = 0; atom0 != frameIn.Natom(); atom0++)
    {
      const double* atm1 = frameIn.XYZ(atom0);
      for (unsigned int idx = 0; idx != singleOccSolvResIdx.size(); idx++)
      {
        SolventRes const& solvRes = solvResArray_[singleOccSolvResIdx[idx]];
        for (int resat1 = solvRes.At0(); resat1 != solvRes.At1(); resat1++)
        {
          if (atom0 >= solvRes.At0() && atom0 < solvRes.At1()) continue;
          const double* atm2 = frameIn.XYZ(resat1);
          // Get imaged distance
          double dist2 = DIST2( imageOpt_.ImagingType(), atm1, atm2, frameIn.BoxCrd() );
          if (dist2 < cut2_) {
            singleOccPeakEne[idx] += Ecalc(atom0, resat1, dist2);
          }
        } // END loop over solvent residue atoms
      } // END loop over singly occupied peaks
    } // END loop over all atoms
    mprintf("DEBUG: Singly-occupied peak energies:\n");
    for (unsigned int idx = 0; idx != singleOccPeakIdx.size(); idx++)
      mprintf("%8i : %g\n", singleOccPeakIdx[idx], singleOccPeakEne[idx]);

    /// Add the energy to the singly-occupied peak sites
    for (unsigned int idx = 0; idx != singleOccPeakIdx.size(); idx++)
    {
      int peakNum = singleOccPeakIdx[idx];
      peakSites_[peakNum].AddSolventEne(frameNum, singleOccPeakEne[idx],
                                        solvResArray_[peakResIdx[peakNum]].Sidx());
    }

    t_energy_.Stop();
  } // END calcEnergy_

  Action::RetType ret = Action::OK;
  if (reorder_) {
    t_reordr_.Start();
    // For each occupied peak N, swap the location of that solvent in the frame
    // with solvent position N.
    // NOTE: This should only happen when 1 solvent is present!
    for (unsigned int peakNum = 0; peakNum != peakSites_.size(); peakNum++)
    {
      if (numTimesPeakAssigned[peakNum] == 1) {
        int solvIdxP = peakResIdx[peakNum];
        int solvIdxS = (int)peakNum;
        SolventRes const& solvP = solvResArray_[solvIdxP];
        SolventRes const& solvS = solvResArray_[solvIdxS];
        int atP = solvP.At0();
        for (int atS = solvS.At0(); atS != solvS.At1(); atS++, atP++)
          frameIn.SwapAtoms(atP, atS);
      }
    }
    ret = Action::MODIFY_COORDS;
    t_reordr_.Stop();
  }
  t_action_.Stop();

  return ret;
}

/// \return absolute value of integer
static inline int absval(int i) { if (i < 0) return -(i+1); else return i; }

/** Calculate the DELTA G of an individual water site */
int Action_Spam::Calc_G_Wat(DataSet* dsIn, int peaknum, Iarray const& SkipFrames)
{
  DataSet_1D const& dataIn = static_cast<DataSet_1D const&>( *dsIn );
  // Create energy vector containing only frames that are singly-occupied.
  // Calculate the mean (enthalpy) while doing this.
  DataSet_double enevec;
  Stats<double> Havg;
  double min = 0.0, max = 0.0;
  if (!SkipFrames.empty()) {
    Iarray::const_iterator fnum = SkipFrames.begin();
    for (int frm = 0; frm != (int)dataIn.Size(); frm++) {
      bool frameIsSkipped = (fnum != SkipFrames.end() && absval(*fnum) == frm);
      if (frameIsSkipped)
        ++fnum;
      else {
        double ene = dataIn.Dval(frm);
        if (enevec.Size() < 1) {
          min = ene; 
          max = ene;
        } else {
          min = std::min(min, ene);
          max = std::max(max, ene);
        }
        enevec.AddElement( ene );
        Havg.accumulate( ene );
      }
    }
  } else {
    min = dataIn.Dval(0);
    max = dataIn.Dval(0);
    for (unsigned int frm = 0; frm != dataIn.Size(); frm++) {
      double ene = dataIn.Dval(frm);
      min = std::min(min, ene);
      max = std::max(max, ene);
      enevec.AddElement( ene );
      Havg.accumulate( ene );
    }
  }
  if (enevec.Size() < 1)
    return 1;
  // Calculate distribution of energy values using KDE. Get the bandwidth
  // factor here since we already know the SD.
  double BWfac = KDE::BandwidthFactor( enevec.Size() );
  if (debug_ > 0)
    mprintf("DEBUG:\tNvals=%zu min=%g max=%g BWfac=%g\n", enevec.Size(), min, max, BWfac);
  // Estimate number of bins the same way spamstats.py does.
  int nbins = (int)(((max - min) / BWfac) + 0.5) + 100;
  if (nbins < 0) {
    // Probably an overflow due to extremely large energy.
    mprintf("Warning: Large magnitude energy observed for peak %i (min=%g max=%g)\n",
            peaknum+1, min, max);
    mprintf("Warning: Skipping peak.\n");
    return -1;
  }

  HistBin Xdim(nbins, min - (50*BWfac), BWfac, "P(Ewat)");
  //Xdim.CalcBinsOrStep(min - Havg.variance(), max + Havg.variance(), 0.0, nbins, "P(Ewat)");
  if (debug_ > 0) {
    mprintf("DEBUG:");
    Xdim.PrintHistBin();
  }
  DataSet_double kde1;
  KDE gkde;
  double bandwidth;
  if (enevec.Size() == 1) {
    // Special case. Juse use BWfac to avoid a zero bandwidth.
    bandwidth = BWfac;
  } else
    bandwidth = 1.06 * sqrt(Havg.variance()) * BWfac;
  if (gkde.CalcKDE( kde1, enevec, Xdim, bandwidth )) {
    mprinterr("Error: Could not calculate E KDE histogram.\n");
    return -1;
  }
  kde1.SetupFormat() = TextFormat(TextFormat::GDOUBLE, 12, 5);
  // Determine SUM[ P(Ewat) * exp(-Ewat / RT) ]
  double RT = Constants::GASK_KCAL * temperature_;
  double KB = 1.0 / RT;
  double sumQ = 0.0;
  for (unsigned int i = 0; i != kde1.Size(); i++) {
    double Ewat = kde1.Xcrd(i);
    double PEwat = kde1.Dval(i);
    sumQ += (PEwat * exp( -Ewat * KB ));
    //mprintf("DEBUG:\t\tEwat %20.10E PEwat %20.10E sumQ %20.10E\n", Ewat, PEwat, sumQ);
  }
  if (debug_ > 0)
    mprintf("DEBUG: peak %6i sumQ= %20.10E\n", peaknum+1, sumQ);
  double DG = -RT * log(BWfac * sumQ);

  double adjustedDG = DG - DG_BULK_;
  double adjustedDH = Havg.mean() - DH_BULK_;
  double ntds = adjustedDG - adjustedDH;

  if (ds_dg_ == 0) {
    mprintf("\tSPAM bulk energy values:\n"
            "\t  <G>= %g, <H>= %g +/- %g, -TdS= %g\n", adjustedDG, adjustedDH,
            sqrt(Havg.variance()), ntds);
  } else {
    ((DataSet_Mesh*)ds_dg_)->AddXY(peaknum+1, adjustedDG);
    ((DataSet_Mesh*)ds_dh_)->AddXY(peaknum+1, adjustedDH);
    ((DataSet_Mesh*)ds_ds_)->AddXY(peaknum+1, ntds);
  }

  // DEBUG
  if (debug_ > 1) {
    FileName rawname("dbgraw." + integerToString(peaknum+1) + ".dat");
    FileName kdename("dbgkde." + integerToString(peaknum+1) + ".dat");
    mprintf("DEBUG: Writing peak %u raw energy values to '%s', KDE histogram to '%s'\n",
            peaknum+1, rawname.full(), kdename.full());
    DataFile rawout;
    rawout.SetupDatafile( rawname, 0 );
    rawout.AddDataSet( &enevec );
    rawout.WriteDataOut();
    DataFile kdeout;
    kdeout.SetupDatafile( kdename, 0 );
    kdeout.AddDataSet( &kde1 );
    kdeout.WriteDataOut();
  }

  return 0;
}

#ifdef MPI
int Action_Spam::SyncAction() {
  // Get total number of frames.
  Iarray frames_on_rank( trajComm_.Size() );
  int myframes = Nframes_;
  trajComm_.GatherMaster( &myframes, 1, MPI_INT, &frames_on_rank[0] );
  if (trajComm_.Master())
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      Nframes_ += frames_on_rank[rank];

  // Sync peakFrameData_
  Iarray size_on_rank( trajComm_.Size() );
  for (unsigned int i = 0; i != peakFrameData_.size(); i++)
  {
    Iarray& Data = peakFrameData_[i];
    int mysize = (int)Data.size();
    trajComm_.GatherMaster( &mysize, 1, MPI_INT, &size_on_rank[0] );
    if (trajComm_.Master()) {
      int total = size_on_rank[0];
      for (int rank = 1; rank < trajComm_.Size(); rank++)
        total += size_on_rank[rank];
      Data.resize( total );
      int* endptr = &(Data[0]) + size_on_rank[0];
      // Receive data from each rank
      int offset = 0;
      for (int rank = 1; rank < trajComm_.Size(); rank++) {
        offset += frames_on_rank[rank-1];
        trajComm_.SendMaster( endptr, size_on_rank[rank], rank, MPI_INT );
        // Properly offset the frame numbers
        for (int j = 0; j != size_on_rank[rank]; j++, endptr++)
          if (*endptr < 0)
            *endptr -= offset;
          else
            *endptr += offset;
      }
    } else // Send data to master
      trajComm_.SendMaster( &(Data[0]), Data.size(), trajComm_.Rank(), MPI_INT );
  }
  return 0;
}
#endif

// Action_Spam::Print()
void Action_Spam::Print() {
  // Timings
  mprintf("\tSPAM timing data:\n");
  t_resCom_.WriteTiming(2, "Residue c.o.m. calc:", t_action_.Total());
  t_assign_.WriteTiming(2, "Peak assignment    :", t_action_.Total());
  t_occupy_.WriteTiming(2, "Occupancy calc.    :", t_action_.Total());
  t_energy_.WriteTiming(2, "Energy calc        :", t_action_.Total());
  if (purewater_)
    pairList_.Timing(t_energy_.Total(), 3);
  t_reordr_.WriteTiming(2, "Residue reordering :", t_action_.Total());
  t_action_.WriteTiming(1, "SPAM Action Total:");
  // Print the spam info file if we didn't do pure water
  if (!purewater_) {
    // Warn about any overflows
    if (overflow_)
      mprinterr("Warning: SPAM: Some frames had a box too small for the cutoff.\n");

    // Print information about omitted frames
    infofile_->Printf("# There are %zu density peaks and %d frames\n\n",
                      peakSites_.size(), Nframes_);
    // Loop over every peak
    for (std::vector<PeakSite>::const_iterator peak = peakSites_.begin();
                                               peak != peakSites_.end(); ++peak)
    {
      // Loop over every solvent type
      for (PeakSite::const_iterator solv = peak->begin(); solv != peak->end(); ++solv)
      {
        Iarray const& peakFrameData = solv->Omitted();
        // Skip peaks with 0 omitted frames
        if (peakFrameData.empty()) continue;
        // Find out how many multiple-occupied frames there are
        unsigned int ndouble = 0;
        for (Iarray::const_iterator it = peakFrameData.begin(); it != peakFrameData.end(); ++it)
          if (*it < 0)
            ndouble++;
        infofile_->Printf("# Peak %li has %zu omitted frames (%u double-occupied)\n", //TODO add solv idx
                          peak - peakSites_.begin() + 1, peakFrameData.size(), ndouble);
        // Print omitted frames
        for (unsigned int j = 0; j < peakFrameData.size(); j++) {
          if (j > 0 && j % 10 == 0) infofile_->Printf("\n");
          // Adjust frame number.
          int fnum = peakFrameData[j];
          if (fnum > -1)
            fnum++;
          infofile_->Printf(" %7d", fnum);
        }
        infofile_->Printf("\n\n");
      }
    }

    if (calcEnergy_) {
      int n_peaks_no_energy = 0;
      for (int p = 0; p != (int)peakSites_.size(); p++)
      {
        for (PeakSite::const_iterator solv = peakSites_[p].begin(); solv != peakSites_[p].end(); ++solv)
        {
          int err = Calc_G_Wat( solv->DS(), p, solv->Omitted() );
          if (err == 1)
            n_peaks_no_energy++;
          else if (err == -1)
            mprintf("Warning: Error calculating SPAM energies for peak %i\n", p + 1);
        }
      }
      if (n_peaks_no_energy > 0)
        mprintf("Warning: No energies for %i peaks.\n", n_peaks_no_energy);
    }
  } else
    Calc_G_Wat( bulk_ene_set_, -1, Iarray() );
}
