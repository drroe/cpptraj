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
#ifdef _OPENMP
# include <omp.h>
#endif

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
  site_size_(0),
  ds_dg_(0),
  ds_dh_(0),
  ds_ds_(0),
  sfile_(0)
{}

/** Construct with bulk solvent name. */
Action_Spam::SolventInfo::SolventInfo(std::string const& n) :
  peaksData_(0),
  site_size_(0),
  name_(n),
  ds_dg_(0),
  ds_dh_(0),
  ds_ds_(0),
  sfile_(0)
{}

/** Construct from peaks data, solvent site size, solvent name */
Action_Spam::SolventInfo::SolventInfo(DataSet_Vector_Scalar const* ds,
                                      double s, std::string const& n) :
  peaksData_(ds),
  site_size_(s),
  name_(n),
  ds_dg_(0),
  ds_dh_(0),
  ds_ds_(0),
  sfile_(0)
{}

/** Create total delta energy DataSets */
int Action_Spam::SolventInfo::CreateDeltaEneSets(std::string const& dsname, int solventIdx,
                                                 DataSetList& DSL,
                                                 DataFileList& DFL, DataFile* summaryfile)
{
  Dimension Pdim( 1.0, 0.0, "Peak" );
  // Create total DG DH TDS sets
  MetaData md(dsname);
  if (solventIdx > 0)
    md.SetIdx( solventIdx );
  md.SetAspect("DG");
  ds_dg_ = DSL.AddSet(DataSet::XYMESH, md);
  md.SetAspect("DH");
  ds_dh_ = DSL.AddSet(DataSet::XYMESH, md);
  md.SetAspect("-TDS");
  ds_ds_ = DSL.AddSet(DataSet::XYMESH, md);
  if (ds_dg_ == 0 || ds_dh_ == 0 || ds_ds_ == 0) {
    mprinterr("Error: Could not create delta energy sets for solvent '%s'\n", Name().c_str());
    return 1;
  }
  ds_dg_->SetDim(Dimension::X, Pdim);
  ds_dh_->SetDim(Dimension::X, Pdim);
  ds_ds_->SetDim(Dimension::X, Pdim);
  sfile_ = 0;
  if (summaryfile != 0) {
    if (solventIdx > 0) {
      // Beyond first solvent; prepend name with hetero solvent name
      sfile_ = DFL.AddDataFile( summaryfile->DataFilename().PrependFileName(Name() + ".") );
      if (sfile_ == 0) {
        mprinterr("Error: Could not create summary file for solvent '%s'\n", Name().c_str());
        return 1;
      }
    } else {
      // First solvent; use summaryfile
      sfile_ = summaryfile;
    }
    sfile_->AddDataSet( ds_dg_ );
    sfile_->AddDataSet( ds_dh_ );
    sfile_->AddDataSet( ds_ds_ );
  }
# ifdef MPI
  ds_dg_->SetNeedsSync(false);
  ds_dh_->SetNeedsSync(false);
  ds_ds_->SetNeedsSync(false);
# endif
  return 0;
}

/** Print solvent info to stdout, no newline. */
void Action_Spam::SolventInfo::PrintInfo() const {
  if (peaksData_ != 0)
    mprintf("'%s', Peaks= '%s', Site size= %.2f Ang.",
            name_.c_str(), peaksData_->legend(), site_size_);
  else
    mprintf("'%s'", name_.c_str());
  if (sfile_ != 0) mprintf(", output to '%s'", sfile_->DataFilename().full());
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
  //mprintf("DEBUG: Constructed with set '%s'\n", ds->legend());
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

/** Add SolventPeak with no output energy DataSet for each solvent type in given array. */
int Action_Spam::PeakSite::AddSolventPeaks(std::vector<SolventInfo> const& solvents)
{
  solvPeaks_.clear();
  for (std::vector<SolventInfo>::const_iterator solv = solvents.begin();
                                                solv != solvents.end(); ++solv)
    solvPeaks_.push_back( SolventPeak() );
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
  printFrameInfo_(true),
  cut2_(144.0),
  onecut2_(1.0 / 144.0),
  doublecut_(24.0),
  infofile_(0),
  sphere_(false),
  bulk_ene_set_(0),
  Nframes_(0),
  overflow_(false)
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
          "\t{ pure |\n"
          "\t  <peaksname> [reorder] [summary <summary>] [info <infofile>]\n"
          "\t  [noframeinfo] [site_size <size>] [sphere] [temperature <T>]\n"
          "\t  [skipE] [dgbulk <dgbulk>] [dhbulk <dhbulk>]\n"
          "\t  [hetsolvent <hname>,<hpeaks>,<hsize>] }\n"
          "  Perform SPAM solvent analysis. If 'pure' is specified calculate\n"
          "  bulk energy values for a pure solvent system. Otherwise determine SPAM\n"
          "  energies from peaks previously identified from the 'volmap' action.\n"
          "    <name>       : Output data set name.\n"
          "    <datafile>   : Data file with all SPAM energies for each snapshot.\n"
          "    <cut>        : Non-bonded cutoff for energy evaluation.\n"
          "    <solvname>   : Name of the bulk solvent residues.\n"
          "    [pure]       : The system is pure solvent. Used to parametrize the bulk values.\n"
          "    <peaksname>  : DataSet/File (XYZ format) with the peak locations for bulk solvent.\n"
          "    [reorder]    : The solvent should be re-ordered so the same solvent molecule\n"
          "                   is always in the same site.\n"
          "    <summary>    : File with the summary of all SPAM results.\n"
          "    <infofile>   : File with stats about which peak sites are occupied when.\n"
          "    [noframeinfo] : If specified do not print individual frame #s to info file.\n"
          "    <size>       : Edge length (or diameter for 'sphere') of each solvent peak site.\n"
          "    [sphere]     : Treat each site like a sphere instead of a box.\n"
          "    <T>          : Temperature at which SPAM calculation was run.\n"
          "    [skipE]      : If specified, skip SPAM energy calculation.\n"
          "    <dgbulk>     : SPAM free energy of the bulk solvent in kcal/mol.\n"
          "    <dhbulk>     : SPAM enthalpy of the bulk solvent in kcal/mol.\n"
          "    [hetsolvent] : If specified, also calculate SPAM energy for heterosolvent with.\n"
          "                   residue name <hname>, peaks DataSet/File <hpeaks>, and site size\n"
          "                   <hsize> (in Ang.).\n"
         );
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

  // See if we're doing pure solvent. If so, we don't need a peak file
  purewater_ = actionArgs.hasKey("pure") || actionArgs.hasKey("purewater");

  // Get data set name.
  std::string ds_name = actionArgs.GetStringKey("name");
  if (ds_name.empty())
    ds_name = init.DSL().GenerateDefaultName("SPAM");

  // Get energy vs frame output data file
  DataFile* datafile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);

  // Bulk solvent residue name
  std::string solvname = actionArgs.GetStringKey("solv");
  if (solvname.empty())
    solvname = std::string("WAT");

  // Get energy cutoff
  double cut = actionArgs.getKeyDouble("cut", 12.0);
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
      mprintf("Warning: 'pure' specified but more arguments remain.\n");
#   ifdef _OPENMP
    // Set up temp space for each thread
    int numthreads;
#   pragma omp parallel
    {
    if (omp_get_thread_num() == 0)
      numthreads = omp_get_num_threads();
    }
    threadResEne_.resize( numthreads );
#   endif
    // Add bulk water info; only need name if not calculating peaks
    solvents_.push_back( SolventInfo(solvname) );
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
    printFrameInfo_ = !actionArgs.hasKey("noframeinfo");
    DataFile* summaryfile = init.DFL().AddDataFile(actionArgs.GetStringKey("summary"), actionArgs);
    // Determine if energy calculation needs to happen
    calcEnergy_ = !actionArgs.hasKey("skipE");
    // Divide site size by 2 to make it half the edge length (or radius)
    double site_size = actionArgs.getKeyDouble("site_size", 2.5) / 2.0;
    sphere_ = actionArgs.hasKey("sphere");
    // If it's a sphere, square the radius to compare with
    if (sphere_)
      site_size *= site_size;
    // Get heterosolvents. Only 1 for now.
    std::string hetSolventStr = actionArgs.GetStringKey("hetsolvent");

    // ----- END processing input args -----------

    // Get or load the peaks data for bulk
    DataSet_Vector_Scalar* peaksData = GetPeaksData(peaksname, init.DSL());
    if (peaksData == 0) {
      mprinterr("Error: Could not get peaks.\n");
      return Action::ERR;
    }
    // Add bulk solvent site info
    solvents_.push_back( SolventInfo(peaksData, site_size, solvname) );
    // Process heterosolvents
    if (!hetSolventStr.empty()) {
      ArgList hetArgs(hetSolventStr, ",");
      // Expect <name>,<peaks>,<size>
      if (hetArgs.Nargs() != 3) {
        mprinterr("Error: Expected 3 comma-separated args for 'hetsolvent', got %i\n", hetArgs.Nargs());
        return Action::ERR;
      }
      std::string sname = hetArgs.GetStringNext();
      DataSet_Vector_Scalar* pdat = GetPeaksData(hetArgs.GetStringNext(), init.DSL());
      if (pdat == 0) {
        mprinterr("Error: Could not get peaks for solvent '%s'\n", hetArgs[0].c_str());
        return Action::ERR;
      }
      double ssize = hetArgs.getNextDouble(2.5) / 2.0;
      if (sphere_)
        ssize *= ssize;
      solvents_.push_back( SolventInfo(pdat, ssize, sname) );
    }

    // DEBUG print solvents
    //mprintf("DEBUG: Solvents:\n");
    //for (std::vector<SolventInfo>::const_iterator it = solvents_.begin(); it != solvents_.end(); ++it)
    // it->PrintInfo();

    // DEBUG print peaks
    if (debug_ > 0) {
      mprintf("DEBUG: Peak sites:\n");
      for (std::vector<PeakSite>::const_iterator it = peakSites_.begin(); it != peakSites_.end(); ++it)
        mprintf("DEBUG:\t%8li %8.3f %8.3f %8.3f\n", it - peakSites_.begin(),
                it->XYZ()[0], it->XYZ()[1], it->XYZ()[2]);
    }

    if (calcEnergy_) {
      // Add energy vs frame DataSets for each solvent to all peak sites
      for (std::vector<PeakSite>::iterator it = peakSites_.begin(); it != peakSites_.end(); ++it)
        if (it->AddEneDataSets(solvents_, ds_name, init.DSL(), datafile, (it-peakSites_.begin())+1))
          return Action::ERR;
      for (std::vector<SolventInfo>::iterator it = solvents_.begin(); it != solvents_.end(); ++it)
        if (it->CreateDeltaEneSets(ds_name, (it-solvents_.begin()), init.DSL(), init.DFL(), summaryfile))
          return Action::ERR;
    } else {
      // Add SolventPeak for each solvent to all peak sites.
      for (std::vector<PeakSite>::iterator it = peakSites_.begin(); it != peakSites_.end(); ++it)
        if (it->AddSolventPeaks(solvents_))
          return Action::ERR;
    }
  } // END if purewater_

  if (purewater_ || calcEnergy_) {
    if (pairList_.InitPairList( cut, 0.1, debug_ )) return Action::ERR;
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
  mprintf("\tSolvent info:\n");
  for (std::vector<SolventInfo>::const_iterator it = solvents_.begin(); it != solvents_.end(); ++it)
  {
    mprintf("\t  ");
    it->PrintInfo();
    mprintf("\n");
  }
  if (purewater_) {
    mprintf("\tCalculating bulk value for pure solvent.\n");
    if (threadResEne_.size() > 1)
      mprintf("\tParallelizing energy calculation with %zu threads.\n", threadResEne_.size());
    if (datafile != 0)
      mprintf("\tPrinting solvent energies to %s\n", datafile->DataFilename().full());
    mprintf("\tData set '%s' index is solvent # * frame.\n", bulk_ene_set_->legend());
    mprintf("\tUsing a %.2f Angstrom non-bonded cutoff with shifted EEL.\n",
            sqrt(cut2_));
    if (reorder_)
      mprintf("\tWarning: Re-ordering makes no sense for pure solvent.\n");
  } else {
    mprintf("\tOccupation information printed to %s.\n", infofile_->Filename().full());
    if (!printFrameInfo_)
      mprintf("\tNot printing individual frame #s to info file.\n");
    mprintf("\tSites are ");
    if (sphere_)
      mprintf("spheres.\n");
    else
      mprintf("boxes.\n");
    if (reorder_)
      mprintf("\tRe-ordering trajectory so each site always has the same solvent molecule.\n");
    if (!calcEnergy_) {
      mprintf("\tSkipping calculation of SPAM energies.\n");
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

  // Set up solvResArray_ and mask_ (for PairList)
  mask_.ResetMask();
  solvResArray_.clear();
  watidx_.clear();
  if (purewater_)
    watidx_.reserve( setup.Top().Natom() );
  int idx = 0;
  // Loop over all residues
  std::vector<unsigned int> solvCount(solvents_.size(), 0);
  for (Topology::res_iterator res = setup.Top().ResStart();
                              res != setup.Top().ResEnd(); res++)
  {
    for (int sidx = 0; sidx != (int)solvents_.size(); sidx++)
    {
      if (res->Name().Truncated() == solvents_[sidx].Name()) {
        solvCount[sidx]++;
        solvResArray_.push_back( SolventRes(res->FirstAtom(), res->LastAtom(), sidx) );
        break;
      }
    }
    // Set up mask/resnums only for solvent for purewater_, everything otherwise
    if (purewater_) {
      if (res->Name().Truncated() == solvents_.front().Name()) {
        for (int i = res->FirstAtom(); i < res->LastAtom(); i++) {
          mask_.AddAtom( i );
          watidx_.push_back( idx );
        }
        idx++;
      }
    } else {
      for (int i = res->FirstAtom(); i < res->LastAtom(); i++)
        mask_.AddAtom( i );
    }
  }
# ifdef _OPENMP
  for (std::vector<Darray>::iterator it = threadResEne_.begin(); it != threadResEne_.end(); ++it)
    it->resize( solvResArray_.size() );
# endif
  //mask_.MaskInfo();
  if (debug_ > 0) {
    mprintf("DEBUG: Solvent residues\n");
    for (std::vector<SolventRes>::const_iterator it = solvResArray_.begin(); it != solvResArray_.end(); ++it)
      it->PrintInfo();
  }
  if (solvResArray_.empty()) {
    mprinterr("Error: No solvent residues.\n");
    return Action::ERR;
  }

  // Reserve space to hold assigned peak for each solvent residue
  resPeakNum_.reserve( solvResArray_.size() );

  for (unsigned int sidx = 0; sidx != solvents_.size(); sidx++)
    mprintf("\tFound %u '%s' solvent residues.\n", solvCount[sidx], solvents_[sidx].Name().c_str());

  // Set up pair list
  if (purewater_ || calcEnergy_) {
    if (pairList_.SetupPairList( currentBox )) return Action::ERR;
  }

  // Set up the charge array and check that we have nonbonded params 
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

/*static inline void dbgprint(int wat, int wat1, double eval)
{
  if (wat < wat1)
    mprintf("DEBUG: %8i - %8i = %10.5f\n", wat, wat1, eval);
  else
    mprintf("DEBUG: %8i - %8i = %10.5f\n", wat, wat1, eval);
}*/
    

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
    mprinterr("Error: Grid setup for bulk energy calculation failed.\n");
    return Action::ERR;
  }
//  int wat = 0, wat1 = 0;
  int basenum = frameNum * solvResArray_.size();
  DataSet_double& evals = static_cast<DataSet_double&>( *bulk_ene_set_ );
  // Make room for each solvent residue energy this frame.
  evals.Resize( evals.Size() + solvResArray_.size() );
  t_energy_.Start();
  //std::vector<Atom> const& Atoms = CurrentParm_->Atoms();
  // Loop over all grid cells
  int cidx;
# ifdef _OPENMP
  int mythread;
  double* resEne;
# pragma omp parallel private(cidx, mythread, resEne)
  {
  mythread = omp_get_thread_num();
  threadResEne_[mythread].assign(threadResEne_[mythread].size(), 0);
  resEne = &(threadResEne_[mythread][0]);
# pragma omp for
# endif
  for (cidx = 0; cidx < pairList_.NGridMax(); cidx++)
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
        int atomi = mask_[it0->Idx()];
        int wat = watidx_[it0->Idx()]; //Atoms[atomi].ResNum();
        Vec3 const& xyz0 = it0->ImageCoords();
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          int wat1 = watidx_[it1->Idx()]; //Atoms[mask_[it1->Idx()]].ResNum();
          if ( wat != wat1 ) {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < cut2_) {
              double eval = Ecalc(atomi, mask_[it1->Idx()], D2);
#             ifdef _OPENMP
//              dbgprint(atomi, mask_[it1->Idx()], eval);
              resEne[wat] += eval;
              resEne[wat1] += eval;
#             else
              evals[basenum + wat] += eval;
              evals[basenum + wat1] += eval;
#             endif
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
            int wat1 = watidx_[it1->Idx()]; //Atoms[mask_[it1->Idx()]].ResNum();
            if ( wat != wat1 ) {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double D2 = dxyz.Magnitude2();
              if (D2 < cut2_) {
                double eval = Ecalc(atomi, mask_[it1->Idx()], D2);
#               ifdef _OPENMP
//                dbgprint(atomi, mask_[it1->Idx()], eval);
                resEne[wat] += eval;
                resEne[wat1] += eval;
#               else
                evals[basenum + wat] += eval;
                evals[basenum + wat1] += eval;
#               endif
              }
            }
          } // END loop over atoms in neighbor cell
        } // END loop over neighbor cells
      } // END loop over atoms in thisCell
    } // END cell not empty
  } // END loop over grid cells
# ifdef _OPENMP
  } // END omp parallel
  // Sum resEne into evals
  for (std::vector<Darray>::const_iterator it = threadResEne_.begin(); it != threadResEne_.end(); ++it)
  {
    int idx = basenum;
    for (Darray::const_iterator jt = it->begin(); jt != it->end(); ++jt, ++idx)
      evals[idx] += *jt;
  }
# endif
  t_energy_.Stop();
  t_action_.Stop();
  return Action::OK;
}

/// Wrap the given grid index if out of bounds; set offset to original direction
static inline int wrapidx(int nz, int MAXZ, int& offset) {
        int iz;
        if (nz < 0) {
          iz = nz + MAXZ;
          offset = -1;
        } else if (nz >= MAXZ) {
          iz = nz - MAXZ;
          offset = 1;
        } else {
          iz = nz;
          offset = 0;
        }
        return iz;
}


/** For each occupied peak, calculate the energy between the solvent
  * molecule occupying the peak and every other atom.
  * \param singleOccSolvResIdx Contain indices into solvResArray_ for solvent associated with singly-occupied peaks
  * \param singleOccPeakIdx Contain indices into peakSites_ for single-occupied peaks
  */
int Action_Spam::Peaks_Ene_Calc(Iarray const& singleOccSolvResIdx,
                                Iarray const& singleOccPeakIdx,
                                Frame const& frameIn,
                                int frameNum)
{
  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), mask_);
  if (retVal != 0) {
    mprinterr("Error: Grid setup for peaks energy calculation failed.\n");
    return Action::ERR;
  }
  std::vector<Atom> const& Atoms = CurrentParm_->Atoms();
  // Loop over singly-occupied peaks
  int idx;
  int maxidx = (int)singleOccSolvResIdx.size();
# ifdef _OPENMP
# pragma omp parallel private(idx)
  {
# pragma omp for
# endif
  for (idx = 0; idx < maxidx; idx++)
  {
    SolventRes const& solvRes = solvResArray_[singleOccSolvResIdx[idx]];
    PeakSite& currentPeak = peakSites_[singleOccPeakIdx[idx]];
    //mprintf("DEBUG: PL solvRes %i - %i peak %i\n", solvRes.At0(), solvRes.At1(), singleOccPeakIdx[idx]);
    double etot = 0;
    // Loop over solvent atoms
    for (int solvAt = solvRes.At0(); solvAt != solvRes.At1(); solvAt++)
    {
      int res0 = Atoms[solvAt].ResNum();
      // get the imaged coords corresponding to solvAt
      Vec3 xyz0 = frameIn.BoxCrd().UnitCell().TransposeMult( pairList_.FracCoords()[solvAt] );
      // Get the grid cell corresponding to solvAt
      int i1, i2, i3;
      // NOTE: This only works if the entire system is in the pair list, which
      //       should be the case when !purewater_.
      pairList_.CalcCellIdx( solvAt, i1, i2, i3 );
//      int cidx = pairList_.CalcCellIdx( solvAt, i1, i2, i3 );
//      mprintf("DEBUG: solvAt= %i i1= %i i2= %i i3= %i cidx %i %i\n", solvAt, i1, i2, i3, cidx, cidxcached);
      // Loop over cell and all surrounding cells
      int cellOffset = 3; // FIXME get from pairlist
      int minz = i3 - cellOffset;
      int maxz = i3 + cellOffset + 1;
      int miny = i2 - cellOffset;
      int maxy = i2 + cellOffset + 1;
      int minx = i1 - cellOffset;
      int maxx = i1 + cellOffset + 1;
      int nGridXY = pairList_.NX() * pairList_.NY();
      for (int nz = minz; nz != maxz; nz++) {
        int oz;
        int iz = wrapidx(nz, pairList_.NZ(), oz);
        //int Cidx3 = nz * nGridXY;
        for (int ny = miny; ny != maxy; ny++) {
          int oy;
          int iy = wrapidx(ny, pairList_.NY(), oy);
          //int Cidx2 = Cidx3 + (ny*nGridX_);
          for (int nx = minx; nx != maxx; nx++) {
            int ox;
            int ix = wrapidx(nx, pairList_.NX(), ox);
            //int thisCellIdx = Cidx2 + nx; // Absolute grid cell index
            int thisCellIdx = (iz*nGridXY)+(iy*pairList_.NX())+ix;
//            mprintf("DEBUG:\t\t cell %i abs {%i %i %i} wrapped {%i %i %i}\n", thisCellIdx, nx, ny, nz, ix, iy, iz);
            PairList::CellType const& myCell = pairList_.Cell( thisCellIdx );
            // Loop over every atom in myCell
            for (PairList::CellType::const_iterator it1 = myCell.begin();
                                                    it1 != myCell.end(); ++it1)
            {
              int res1 = Atoms[mask_[it1->Idx()]].ResNum();
//              mprintf("DEBUG: At %i res %i\n", mask_[it1->Idx()], res1);
              if ( res0 != res1 ) {
                Vec3 transVec(0.0);
                if (ox != 0 || oy != 0 || oz != 0)
                  transVec = frameIn.BoxCrd().UnitCell().TransposeMult( Vec3(ox, oy, oz) );
//                mprintf("DEBUG:\t\t\t offsets %i %i %i  transVec= %f %f %f\n",
//                        ox, oy, oz, transVec[0], transVec[1], transVec[2]);
                Vec3 const& xyz1 = it1->ImageCoords();
                Vec3 dxyz = xyz1 + transVec - xyz0;
                double D2 = dxyz.Magnitude2();
//                mprintf("DEBUG:\t\t\tDist= %g\n", sqrt(D2));
                if (D2 < cut2_) {
                  double eval = Ecalc(solvAt, mask_[it1->Idx()], D2);
//                  if (solvAt < mask_[it1->Idx()])
//                    mprintf("DEBUG1: %8i - %8i = %g\n", solvAt, mask_[it1->Idx()], eval);
//                  else
//                    mprintf("DEBUG1: %8i - %8i = %g\n", mask_[it1->Idx()], solvAt, eval);
                  etot += eval;
                }
              }
            } // END loop over every atom in myCell
          } // END nx loop
        } // END ny loop
      } // END nz loop
    } // END loop over solvent atoms
    // Add energy to singly-occupied peak site
    currentPeak.AddSolventEne(frameNum, etot, solvRes.Sidx());
    //mprintf("DEBUG: PL ene: peak %8i = %g\n", singleOccPeakIdx[idx], etot);
  } // END loop over singly-occupied peaks
# ifdef _OPENMP
  } // END omp parallel
# endif
  return 0;
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
# ifdef _OPENMP
# pragma omp parallel private(solvResIdx)
  {
# pragma omp for
# endif
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
# ifdef _OPENMP
  } // END omp parallel
# endif

  // DEBUG - print peak assignments for each solvent
/*
  mprintf("DEBUG: Peak assignments for each solvent idx (new):\n");
  for (Iarray::const_iterator it = resPeakNum_.begin(); it != resPeakNum_.end(); ++it)
    if (*it > -1)
      mprintf("DEBUG:\t%8li %i\n", it - resPeakNum_.begin(), *it);
*/
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

/*
  mprintf("DEBUG: Peak assignment stats:\n");
  for (unsigned int idx = 0; idx != peakSites_.size(); ++idx)
    if (numTimesPeakAssigned[idx] > 0)
      mprintf("DEBUG:\t%8u %u (%i)\n", idx, numTimesPeakAssigned[idx], peakResIdx[idx]);
  mprintf("DEBUG: Singly-occupied peaks:\n");
  for (Iarray::const_iterator it = singleOccSolvResIdx.begin();
                              it != singleOccSolvResIdx.end(); ++it)
    mprintf("DEBUG:\t%8i %8i - %8i\n", *it,
            solvResArray_[*it].At0(), solvResArray_[*it].At1());
*/

  t_occupy_.Stop();

  // Energy calculation
  if (calcEnergy_) {
    t_energy_.Start();
/*
    // Loop over each singly-occupied peak
    for (unsigned int idx = 0; idx != singleOccSolvResIdx.size(); idx++)
    {
      double etot = 0.0;
      SolventRes const& solvRes = solvResArray_[singleOccSolvResIdx[idx]];
      for (int resat1 = solvRes.At0(); resat1 != solvRes.At1(); resat1++)
      {
        const double* atm1 = frameIn.XYZ(resat1);
        for (int atom0 = 0; atom0 != frameIn.Natom(); atom0++)
        {
          if (atom0 >= solvRes.At0() && atom0 < solvRes.At1()) continue;
          // Get imaged distance
          double dist2 = DIST2( imageOpt_.ImagingType(), frameIn.XYZ(atom0), atm1, frameIn.BoxCrd() );
          if (dist2 < cut2_)
            etot += Ecalc(atom0, resat1, dist2);
        } // END loop over all atoms
      } // END loop over solvent residue atoms
      // Add energy to singly-occupied peak site
      peakSites_[singleOccPeakIdx[idx]].AddSolventEne(frameNum, etot, solvRes.Sidx());
    } // END loop over singly-occupied peaks
*/
/*
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
            double eval = Ecalc(atom0, resat1, dist2);
//            if (atom0 < resat1)
//              mprintf("DEBUG0: %8i - %8i = %g\n", atom0, resat1, eval);
//            else
//              mprintf("DEBUG0: %8i - %8i = %g\n", resat1, atom0, eval);
            singleOccPeakEne[idx] += eval;
            //singleOccPeakEne[idx] += Ecalc(atom0, resat1, dist2);
          }
        } // END loop over solvent residue atoms
      } // END loop over singly occupied peaks
    } // END loop over all atoms
//
//    mprintf("DEBUG: Singly-occupied peak energies:\n");
//    for (unsigned int idx = 0; idx != singleOccPeakIdx.size(); idx++)
//      mprintf("%8i : %g\n", singleOccPeakIdx[idx], singleOccPeakEne[idx]);
//
    /// Add the energy to the singly-occupied peak sites
    for (unsigned int idx = 0; idx != singleOccPeakIdx.size(); idx++)
    {
      int peakNum = singleOccPeakIdx[idx];
      peakSites_[peakNum].AddSolventEne(frameNum, singleOccPeakEne[idx],
                                        solvResArray_[peakResIdx[peakNum]].Sidx());
    }
*/
    Peaks_Ene_Calc(singleOccSolvResIdx, singleOccPeakIdx, frameIn, frameNum);

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

/** Calculate G by integrating the probability distribution constructed
  * from the input array of energy values.
  * \return 1 if energy array is empty.
  * \return -1 if number of bins cannot be calculated.
  * \return 0 if the calculation completed successfully.
  */
int Action_Spam::Calc_G(double& DG, int peaknum, double min, double max, double variance,
                        DataSet_double const& enevec)
const
{
  DG = 0;
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
  // Construct the KDE histogram
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
    bandwidth = 1.06 * sqrt(variance) * BWfac;
  if (gkde.CalcKDE( kde1, enevec, Xdim, bandwidth )) {
    mprinterr("Error: Could not calculate E KDE histogram.\n");
    return -1;
  }
  // Summation over the probability distribution 
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
  DG = -RT * log(BWfac * sumQ);

  // DEBUG
  if (debug_ > 1) {
    FileName rawname("dbgraw." + integerToString(peaknum+1) + ".dat");
    FileName kdename("dbgkde." + integerToString(peaknum+1) + ".dat");
    mprintf("DEBUG: Writing peak %u raw energy values to '%s', KDE histogram to '%s'\n",
            peaknum+1, rawname.full(), kdename.full());
    DataFile rawout;
    rawout.SetupDatafile( rawname, 0 );
    rawout.AddDataSet( (DataSet*)(&enevec) );
    rawout.WriteDataOut();
    DataFile kdeout;
    kdeout.SetupDatafile( kdename, 0 );
    kdeout.AddDataSet( &kde1 );
    kdeout.WriteDataOut();
  }

  return 0;
}

/** Calculate the G, H and TS for bulk water. */
int Action_Spam::Calc_Bulk() const {
  DataSet_double const& enevec = static_cast<DataSet_double const&>( *bulk_ene_set_ );
  Stats<double> Havg;
  double min = enevec.Dval(0);
  double max = enevec.Dval(0);
  for (unsigned int frm = 0; frm != enevec.Size(); frm++) {
    double ene = enevec.Dval(frm);
    min = std::min(min, ene);
    max = std::max(max, ene);
    Havg.accumulate( ene );
  }
  // Get the G value
  double DG = 0;
  int err = Calc_G(DG, -1, min, max, Havg.variance(), enevec);
  if (err != 0) {
    mprinterr("Error: Could not get SPAM bulk energy values.\n");
    return 1;
  }
  mprintf("\tSPAM bulk energy values:\n"
          "\t  <G>= %g, <H>= %g +/- %g, -TdS= %g\n", DG, Havg.mean(), DG - Havg.mean());
  return 0;
}

/// \return absolute value of integer
static inline int absval(int i) { if (i < 0) return -(i+1); else return i; }

/** Calculate G, H, and TS for a peak site. */
int Action_Spam::Calc_G_Peak(unsigned int peakNum, PeakSite const& peakSite) const {
  // Assume the first solvent type is to be calculated relative to bulk.
  double G_ref = DG_BULK_;
  double H_ref = DH_BULK_;

  // Loop over solvent types for this peak
  for (PeakSite::const_iterator solv = peakSite.begin(); solv != peakSite.end(); ++solv)
  {
    // Create energy vector containing only frames that are singly-occupied.
    // Calculate the mean enthalpy while doing this.
    Iarray const& SkipFrames = solv->Omitted();
    DataSet_1D const& dataIn = static_cast<DataSet_1D const&>( *(solv->DS()) );
    double min = 0.0;
    double max = 0.0;
    DataSet_double enevec;
    Stats<double> Havg;
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
    double DG = 0;
    // TODO put all peak-related messages in this function instead of Calc_G
    int err = Calc_G(DG, peakNum, min, max, Havg.variance(), enevec);
    if (err != 0) return err;

    // Calculate and record delta G, H, -TdS
    double adjustedDG = DG - G_ref;
    double adjustedDH = Havg.mean() - H_ref;
    double ntds = adjustedDG - adjustedDH;

    long int sidx = solv - peakSite.begin();

    ((DataSet_Mesh*)solvents_[sidx].DG())->AddXY(peakNum+1, adjustedDG);
    ((DataSet_Mesh*)solvents_[sidx].DH())->AddXY(peakNum+1, adjustedDH);
    ((DataSet_Mesh*)solvents_[sidx].TDS())->AddXY(peakNum+1, ntds);

    if (solv == peakSite.begin()) {
      // The reference energies for the next solvents will be this solvents
      G_ref = DG;
      H_ref = Havg.mean();
    }

  } // END loop over solvent types for peak
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
  for (std::vector<PeakSite>::iterator peak = peakSites_.begin(); peak != peakSites_.end(); ++peak)
  {
    for (unsigned int sidx = 0; sidx != peak->Nsolvs(); sidx++)
    {

      Iarray& Data = peak->ModifySolv(sidx).ModifyOmitted();
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
    } // END loop over solvents at this peak
  } // END loop over peaks
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
        if (solvents_.size() == 1)
          infofile_->Printf("# Peak %li has %zu omitted frames (%u double-occupied)\n",
                            peak - peakSites_.begin() + 1, peakFrameData.size(), ndouble);
        else
          infofile_->Printf("# Peak %li solvent %li has %zu omitted frames (%u double-occupied)\n",
                            peak - peakSites_.begin() + 1, 
                            solv - peak->begin() + 1, peakFrameData.size(), ndouble);
        // Print omitted frames
        if (printFrameInfo_) {
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
    }

    if (calcEnergy_) {
      int n_peaks_no_energy = 0;
      for (int p = 0; p != (int)peakSites_.size(); p++)
      {
        int err = Calc_G_Peak(p, peakSites_[p]);
        if (err == 1)
          n_peaks_no_energy++;
        else if (err == -1)
          mprintf("Warning: Error calculating SPAM energies for peak %i\n", p + 1);
      }
      if (n_peaks_no_energy > 0)
        mprintf("Warning: No energies for %i peaks.\n", n_peaks_no_energy);
    }
  } else
    Calc_Bulk();
}
