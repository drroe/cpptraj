#include <cmath> // sqrt
#include <algorithm> // sort
#include "Action_HydrogenBond.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_2D.h"
#include "DataSet_integer.h"
#include "DistRoutines.h"
#include "StringRoutines.h" // ByteString
#include "TorsionRoutines.h"
#ifdef _OPENMP
# include <omp.h>
#endif

// CONSTRUCTOR
Action_HydrogenBond::Action_HydrogenBond() :
  CurrentParm_(0),
  masterDSL_(0),
  NumHbonds_(0),
  NumSolvent_(0),
  NumBridge_(0),
  BridgeID_(0),
  UU_matrix_byRes_(0),
  UUseriesout_(0),
  UVseriesout_(0),
  Bseriesout_(0),
  avgout_(0),
  solvout_(0),
  bridgeout_(0),
  dcut2_(0.0),
  acut_(0.0),
  bothEnd_(0),
  Nframes_(0),
  debug_(0),
  UUmatByRes_norm_(NORM_FRAMES),
  series_(false),
  Bseries_(false),
  seriesUpdated_(false),
  useAtomNum_(false),
  noIntramol_(false),
  hasDonorMask_(false),
  hasDonorHmask_(false),
  hasAcceptorMask_(false),
  hasSolventDonor_(false),
  hasSolventAcceptor_(false),
  calcSolvent_(false),
  bridgeByAtom_(false)
{}

// void Action_HydrogenBond::Help()
void Action_HydrogenBond::Help() const {
  mprintf("\t[<dsname>] [out <filename>] [<mask>] [angle <acut>] [dist <dcut>]\n"
          "\t[donormask <dmask> [donorhmask <dhmask>]] [acceptormask <amask>]\n"
          "\t[avgout <filename>] [printatomnum] [nointramol] [image]\n"
          "\t[solventdonor <sdmask>] [solventacceptor <samask>]\n"
          "\t[solvout <filename>] [bridgeout <filename>] [bridgebyatom]\n"
          "\t[series [uuseries <filename>] [uvseries <filename>]]\n"
          "\t[bseries [bseriesfile <filename>]]\n"
          "\t[uuresmatrix [uuresmatrixnorm {none|frames}] [uuresmatrixout <file>]]\n"
          "\t[splitframe <comma-separated-list>]\n"
          "  Hydrogen bond is defined as A-HD, where A is acceptor heavy atom, H is\n"
          "  hydrogen, D is donor heavy atom. Hydrogen bond is formed when\n"
          "  A to D distance < dcut and A-H-D angle > acut; if acut < 0 it is ignored.\n"
          "  Search for hydrogen bonds using atoms in the region specified by mask.\n"
          "  If just <mask> specified donors and acceptors will be automatically searched for.\n"
          "  If donormask is specified but not acceptormask, acceptors will be\n"
          "  automatically searched for in <mask>.\n"
          "  If acceptormask is specified but not donormask, donors will be automatically\n"
          "  searched for in <mask>.\n"
          "  If both donormask and acceptor mask are specified no automatic searching will occur.\n"
          "  If donorhmask is specified atoms in that mask will be paired with atoms in\n"
          "  donormask instead of automatically searching for hydrogen atoms.\n"
          "  The 'splitframe' keyword can be used to divide the average hydrogen bond\n"
          "  analysis into parts, e.g. 'splitframe 250,500,1000' will divide analysis\n"
          "  into 1-249,250-499,500-999, etc.\n");
}

// Action_HydrogenBond::Init()
Action::RetType Action_HydrogenBond::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Get keywords
  imageOpt_.InitImaging( (actionArgs.hasKey("image")) );
  DataFile* DF = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  series_ = actionArgs.hasKey("series");
  if (series_) {
    UUseriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("uuseries"), actionArgs);
    UVseriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("uvseries"), actionArgs);
    init.DSL().SetDataSetsPending(true);
  }
  Bseries_ = actionArgs.hasKey("bseries");
  if (Bseries_) {
    Bseriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("bseriesfile"), actionArgs);
    init.DSL().SetDataSetsPending(true);
  }
  bool do_uuResMatrix = actionArgs.hasKey("uuresmatrix");
  DataFile* uuResMatrixFile = 0;
  if (do_uuResMatrix) {
    uuResMatrixFile = init.DFL().AddDataFile(actionArgs.GetStringKey("uuresmatrixout"),
                                             ArgList("nosquare2d"),
                                             actionArgs);
    std::string uuResMatrixNorm = actionArgs.GetStringKey("uuresmatrixnorm");
    if (!uuResMatrixNorm.empty()) {
      if (uuResMatrixNorm == "none")
        UUmatByRes_norm_ = NORM_NONE;
      else if (uuResMatrixNorm == "frames")
        UUmatByRes_norm_ = NORM_FRAMES;
      //else if (uuResMatrixNorm == "resmax") // DISABLE for now, needs more testing
      //  UUmatByRes_norm_ = NORM_RESMAX;
      else {
        mprinterr("Error: Invalid keyword for 'uuresmatrixnorm'.\n");
        return Action::ERR;
      }
    }
    
  }
  std::string avgname = actionArgs.GetStringKey("avgout");
  std::string solvname = actionArgs.GetStringKey("solvout");
  if (solvname.empty()) solvname = avgname;
  std::string bridgename = actionArgs.GetStringKey("bridgeout");
  if (bridgename.empty()) bridgename = solvname;
  
  useAtomNum_ = actionArgs.hasKey("printatomnum");
  acut_ = actionArgs.getKeyDouble("angle",135.0);
  noIntramol_ = actionArgs.hasKey("nointramol");
  bridgeByAtom_ = actionArgs.hasKey("bridgebyatom");
  // Hbond split analysis
  std::string splitarg = actionArgs.GetStringKey("splitframe");
  if (!splitarg.empty()) {
    ArgList splits( splitarg, "," );
    if (splits.Nargs() < 1) {
      mprinterr("Error: Invalid argument for 'splitframe': %s\n", splitarg.c_str());
      return Action::ERR;
    }
    int sf = splits.getNextInteger(-1); // User frame #s start at 1
    while (sf > 0) {
      // Since user frame #s start at 1, subtract 1 for actual frame index.
      sf--;
      // Check that frame arg is valid
      if (!splitFrames_.empty()) {
        if (sf <= splitFrames_.back()) {
          mprinterr("Error: 'splitframe' #s must be in increasing order.\n");
          return Action::ERR;
        }
      }
      splitFrames_.push_back( sf );
      sf = splits.getNextInteger(-1);
    }
    if ((int)splitFrames_.size() < splits.Nargs()) {
      mprinterr("Error: Invalid split frame arguments.\n");
      splits.CheckForMoreArgs();
      return Action::ERR;
    }
  }
  // Convert angle cutoff to radians
  acut_ *= Constants::DEGRAD;
  double dcut = actionArgs.getKeyDouble("dist",3.0);
  dcut = actionArgs.getKeyDouble("distance", dcut); // for PTRAJ compat.
  dcut2_ = dcut * dcut;
  // Get donor mask
  std::string mask = actionArgs.GetStringKey("donormask");
  if (!mask.empty()) {
    if (DonorMask_.SetMaskString(mask)) return Action::ERR;
    hasDonorMask_=true;
    // Get donorH mask (if specified)
    mask = actionArgs.GetStringKey("donorhmask");
    if (!mask.empty()) {
      if (DonorHmask_.SetMaskString(mask)) return Action::ERR;
      hasDonorHmask_=true;
    }
  }
  // Get acceptor mask
  mask = actionArgs.GetStringKey("acceptormask");
  if (!mask.empty()) {
    if (AcceptorMask_.SetMaskString(mask)) return Action::ERR;
    hasAcceptorMask_=true;
  }
  // Get solvent donor mask
  mask = actionArgs.GetStringKey("solventdonor");
  if (!mask.empty()) {
    if (SolventDonorMask_.SetMaskString(mask)) return Action::ERR;
    hasSolventDonor_ = true;
    calcSolvent_ = true;
  }
  // Get solvent acceptor mask
  mask = actionArgs.GetStringKey("solventacceptor");
  if (!mask.empty()) {
    if (SolventAcceptorMask_.SetMaskString(mask)) return Action::ERR;
    hasSolventAcceptor_ = true;
    calcSolvent_ = true;
  }
  // Get generic mask
  if (Mask_.SetMaskString(actionArgs.GetMaskNext())) return Action::ERR;

  // Setup datasets
  hbsetname_ = actionArgs.GetStringNext();
  if (hbsetname_.empty())
    hbsetname_ = init.DSL().GenerateDefaultName("HB");
  NumHbonds_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UU"));
  if (NumHbonds_==0) return Action::ERR;
  if (DF != 0) DF->AddDataSet( NumHbonds_ );
  avgout_ = init.DFL().AddCpptrajFile(avgname, "Avg. solute-solute HBonds");
  if (calcSolvent_) {
    NumSolvent_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UV"));
    if (NumSolvent_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( NumSolvent_ );
    NumBridge_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "Bridge"));
    if (NumBridge_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( NumBridge_ );
    BridgeID_ = init.DSL().AddSet(DataSet::STRING, MetaData(hbsetname_, "ID"));
    if (BridgeID_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( BridgeID_ );
    solvout_ = init.DFL().AddCpptrajFile(solvname,"Avg. solute-solvent HBonds");
    bridgeout_ = init.DFL().AddCpptrajFile(bridgename,"Solvent bridging info");
  }
  if (do_uuResMatrix) {
    UU_matrix_byRes_ = (DataSet_2D*)
                       init.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(hbsetname_, "UUresmat"));
    if (UU_matrix_byRes_ == 0) return Action::ERR;
    UU_matrix_byRes_->ModifyDim(Dimension::X).SetLabel("Res");
    UU_matrix_byRes_->ModifyDim(Dimension::Y).SetLabel("Res");
    if (uuResMatrixFile != 0)
      uuResMatrixFile->AddDataSet( UU_matrix_byRes_ );
  }
    
# ifdef _OPENMP
  // Each thread needs temp. space to store found hbonds every frame
  // to avoid memory clashes when adding/updating in map.
# pragma omp parallel
  {
# pragma omp master
  {
  thread_HBs_.resize( omp_get_num_threads() );
  }
  }
# endif

  mprintf( "  HBOND: ");
  if (!hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Searching for Hbond donors/acceptors in region specified by %s\n",
            Mask_.MaskString());
  else if (hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Donor mask is %s, acceptors will be searched for in region specified by %s\n",
            DonorMask_.MaskString(), Mask_.MaskString());
  else if (hasAcceptorMask_ && !hasDonorMask_)
    mprintf("Acceptor mask is %s, donors will be searched for in a region specified by %s\n",
            AcceptorMask_.MaskString(), Mask_.MaskString());
  else
    mprintf("Donor mask is %s, Acceptor mask is %s\n",
            DonorMask_.MaskString(), AcceptorMask_.MaskString());
  if (hasDonorHmask_)
    mprintf("\tSeparate donor H mask is %s\n", DonorHmask_.MaskString() );
# ifdef _OPENMP
  if (thread_HBs_.size() > 1)
    mprintf("\tParallelizing calculation with %zu threads.\n", thread_HBs_.size());
# endif
  if (noIntramol_)
    mprintf("\tOnly looking for intermolecular hydrogen bonds.\n");
  if (hasSolventDonor_)
    mprintf("\tWill search for hbonds between solute and solvent donors in [%s]\n",
            SolventDonorMask_.MaskString());
  if (hasSolventAcceptor_)
    mprintf("\tWill search for hbonds between solute and solvent acceptors in [%s]\n",
            SolventAcceptorMask_.MaskString());
  mprintf("\tDistance cutoff = %.3f, Angle Cutoff = %.3f\n",dcut,acut_*Constants::RADDEG);
  if (DF != 0) 
    mprintf("\tWriting # Hbond v time results to %s\n", DF->DataFilename().full());
  if (avgout_ != 0)
    mprintf("\tWriting Hbond avgs to %s\n",avgout_->Filename().full());
  if (!splitFrames_.empty()) {
    mprintf("\tWill split analysis at frames:");
    for (Iarray::const_iterator it = splitFrames_.begin(); it != splitFrames_.end(); ++it)
      mprintf(" %i", *it + 1);
    mprintf("\n");
  }
  if (calcSolvent_) {
    if (solvout_ != 0)
      mprintf("\tWriting solute-solvent hbond avgs to %s\n", solvout_->Filename().full());
    if (bridgeout_ != 0)
      mprintf("\tWriting solvent bridging info to %s\n", bridgeout_->Filename().full());
    if (bridgeByAtom_)
      mprintf("\tSolvent bridges will be determined between solute atoms.\n");
    else
      mprintf("\tSolvent bridges will be determined between solute residues.\n");
  }
  if (useAtomNum_)
    mprintf("\tAtom numbers will be written to output.\n");
  if (series_) {
    mprintf("\tTime series data for each hbond will be saved for analysis.\n");
    if (UUseriesout_ != 0) mprintf("\tWriting solute-solute time series to %s\n",
                                   UUseriesout_->DataFilename().full());
    if (UVseriesout_ != 0) mprintf("\tWriting solute-solvent time series to %s\n",
                                   UVseriesout_->DataFilename().full());
  }
  if (Bseries_) {
    mprintf("\tTime series data for each bridge will be saved for analysis.\n");
    if (Bseriesout_ != 0)
      mprintf("\tWriting bridge time series to '%s'\n", Bseriesout_->DataFilename().full());
  }
  if (UU_matrix_byRes_ != 0) {
    mprintf("\tCalculating solute-solute residue matrix: %s\n", UU_matrix_byRes_->legend());
    if (uuResMatrixFile != 0)
      mprintf("\tWriting solute-solute residue matrix to '%s'\n", uuResMatrixFile->DataFilename().full());
    if (UUmatByRes_norm_ == NORM_NONE)
      mprintf("\tNot normalizing solute-solute residue matrix.\n");
    else if (UUmatByRes_norm_ == NORM_FRAMES)
      mprintf("\tNormalizing solute-solute residue matrix by frames.\n");
    else if (UUmatByRes_norm_ == NORM_RESMAX)
      mprintf("\tNormalizing solute-solute residue matrix by max possible h. bonds between residues.\n");
  }
  if (imageOpt_.UseImage())
    mprintf("\tImaging enabled.\n");
  masterDSL_ = init.DslPtr();

  return Action::OK;
}

// IsFON()
/** Default criterion for being a hydrogen bond donor/acceptor. */
inline bool IsFON(Atom const& atm) {
  return (atm.Element() == Atom::FLUORINE ||
          atm.Element() == Atom::OXYGEN ||
          atm.Element() == Atom::NITROGEN);
}

/** Create legend for hydrogen bond based on given atoms. */
static inline std::string CreateHBlegend(Topology const& topIn, int a_atom, int h_atom, int d_atom)
{
  if (a_atom == -1)
    return (topIn.TruncResAtomName(h_atom) + "-V");
  else if (d_atom == -1)
    return (topIn.TruncResAtomName(a_atom) + "-V");
  else
    return (topIn.TruncResAtomName(a_atom) + "-" +
            topIn.TruncResAtomName(d_atom) + "-" +
            topIn[h_atom].Name().Truncated());
}

/** Create legend for bridge based on given indices. */
static inline std::string CreateBridgeLegend(std::string const& prefix, std::set<int> indices)
{
  std::string blegend(prefix);
  for (std::set<int>::const_iterator brs = indices.begin(); brs != indices.end(); ++brs)
    blegend.append("_" + integerToString(*brs + 1));
  return blegend;
}

// Action_HydrogenBond::Setup()
Action::RetType Action_HydrogenBond::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );

  // Set up generic mask
  if (!hasDonorMask_ || !hasAcceptorMask_) {
    if ( setup.Top().SetupIntegerMask( Mask_ ) ) return Action::ERR;
    if ( Mask_.None() ) {
      mprintf("Warning: Mask has no atoms.\n");
      return Action::SKIP;
    }
  }

  // ACCEPTOR MASK SETUP
  if (hasAcceptorMask_) {
    // Acceptor mask specified
    if ( setup.Top().SetupIntegerMask( AcceptorMask_ ) ) return Action::ERR;
    if (AcceptorMask_.None()) {
      mprintf("Warning: AcceptorMask has no atoms.\n");
      return Action::SKIP;
    }
  } else {
    // No specified acceptor mask; search generic mask.
    AcceptorMask_.ResetMask();
    for (AtomMask::const_iterator at = Mask_.begin(); at != Mask_.end(); ++at) {
      // Since an acceptor mask was not specified ignore solvent.
      int molnum = setup.Top()[*at].MolNum();
      if (!setup.Top().Mol(molnum).IsSolvent() && IsFON( setup.Top()[*at] ))
        AcceptorMask_.AddSelectedAtom( *at );
    }
    AcceptorMask_.SetNatoms( Mask_.NmaskAtoms() );
  }
  int lastAcceptor = -1;
  if (!AcceptorMask_.None())
    lastAcceptor = AcceptorMask_.back();

  // Clear existing sites
  Both_.clear();
  Acceptor_.clear();
  SolventSites_.clear(); 

  // SOLUTE DONOR/ACCEPTOR SITE SETUP
  Sarray donorOnly;
  if (hasDonorMask_) {
    // Donor heavy atom mask specified
    if ( setup.Top().SetupIntegerMask( DonorMask_ ) ) return Action::ERR;
    if (DonorMask_.None()) {
      mprintf("Warning: DonorMask has no atoms.\n");
      return Action::SKIP;
    }
    if ( hasDonorHmask_ ) {
      // Donor hydrogen mask also specified
      if ( setup.Top().SetupIntegerMask( DonorHmask_ ) ) return Action::ERR;
      if ( DonorHmask_.None() ) {
        mprintf("Warning: Donor H mask has no atoms.\n");
        return Action::SKIP;
      }
      if ( DonorHmask_.Nselected() != DonorMask_.Nselected() ) {
        mprinterr("Error: There is not a 1 to 1 correspondance between donor and donorH masks.\n");
        mprinterr("Error: donor (%i atoms), donorH (%i atoms).\n",DonorMask_.Nselected(),
                  DonorHmask_.Nselected());
        return Action::ERR;
      }
      int maxAtom = std::max( lastAcceptor, DonorMask_.back() ) + 1;
      AtomMask::const_iterator a_atom = AcceptorMask_.begin();
      AtomMask::const_iterator d_atom = DonorMask_.begin();
      AtomMask::const_iterator h_atom = DonorHmask_.begin();
      bool isDonor, isAcceptor;
      int H_at = -1;
      for (int at = 0; at != maxAtom; at++)
      {
        isDonor = false;
        isAcceptor = false;
        if ( d_atom != DonorMask_.end() && *d_atom == at ) {
          isDonor = true;
          ++d_atom;
          H_at = *(h_atom++);
        } 
        if ( a_atom != AcceptorMask_.end() && *a_atom == at ) {
          isAcceptor = true;
          ++a_atom;
        }
        if (isDonor && isAcceptor)
          Both_.push_back( Site(at, H_at) );
        else if (isDonor)
          donorOnly.push_back( Site(at, H_at) );
        else if (isAcceptor)
          Acceptor_.push_back( at );
      }
    } else {
      // No donor hydrogen mask; use any hydrogens bonded to donor heavy atoms.
      int maxAtom = std::max( lastAcceptor, DonorMask_.back() ) + 1;
      AtomMask::const_iterator a_atom = AcceptorMask_.begin();
      AtomMask::const_iterator d_atom = DonorMask_.begin();
      bool isDonor, isAcceptor;
      for (int at = 0; at != maxAtom; at++)
      {
        Iarray Hatoms;
        isDonor = false;
        isAcceptor = false;
        if ( d_atom != DonorMask_.end() && *d_atom == at ) {
          ++d_atom;
          for (Atom::bond_iterator H_at = setup.Top()[at].bondbegin();
                                   H_at != setup.Top()[at].bondend(); ++H_at)
            if (setup.Top()[*H_at].Element() == Atom::HYDROGEN)
              Hatoms.push_back( *H_at );
          isDonor = !Hatoms.empty();
        }
        if ( a_atom != AcceptorMask_.end() && *a_atom == at ) {
          isAcceptor = true;
          ++a_atom;
        }
        if (isDonor && isAcceptor)
          Both_.push_back( Site(at, Hatoms) );
        else if (isDonor)
          donorOnly.push_back( Site(at, Hatoms) );
        else if (isAcceptor)
          Acceptor_.push_back( at );
      }
    }
  } else {
    // No specified donor mask; search generic mask.
    int maxAtom = std::max( lastAcceptor, Mask_.back() ) + 1;
    AtomMask::const_iterator a_atom = AcceptorMask_.begin();
    AtomMask::const_iterator d_atom = Mask_.begin();
    bool isDonor, isAcceptor;
    for (int at = 0; at != maxAtom; at++)
    {
      // Since an acceptor mask was not specified ignore solvent.
      Atom const& atom = setup.Top()[at];
      int molnum = atom.MolNum();
      if (!setup.Top().Mol(molnum).IsSolvent())
      {
        Iarray Hatoms;
        isDonor = false;
        isAcceptor = false;
        if ( d_atom != Mask_.end() && *d_atom == at) {
          ++d_atom;
          if ( IsFON( atom ) ) {
            for (Atom::bond_iterator H_at = atom.bondbegin();
                                     H_at != atom.bondend(); ++H_at)
              if (setup.Top()[*H_at].Element() == Atom::HYDROGEN)
                Hatoms.push_back( *H_at );
            isDonor = !Hatoms.empty();
          }
        }
        if ( a_atom != AcceptorMask_.end() && *a_atom == at ) {
          isAcceptor = true;
          ++a_atom;
        }
        if (isDonor && isAcceptor)
          Both_.push_back( Site(at, Hatoms) );
        else if (isDonor)
          donorOnly.push_back( Site(at, Hatoms) );
        else if (isAcceptor)
          Acceptor_.push_back( at );
      }
    }
  }
  // Place donor-only sites at the end of Both_
  bothEnd_ = Both_.size();
  for (Sarray::const_iterator site = donorOnly.begin(); site != donorOnly.end(); ++site)
    Both_.push_back( *site );

  // Print solute site stats
  mprintf("\tAcceptor-only atoms: %zu\n", Acceptor_.size());
  if (debug_ > 0)
    for (Iarray::const_iterator at = Acceptor_.begin(); at != Acceptor_.end(); ++at)
      mprintf("\t  %20s %8i\n", setup.Top().TruncResAtomName(*at).c_str(), *at+1);
  unsigned int hcount = 0;
  mprintf("\tDonor/acceptor sites: %u\n", bothEnd_);
  Sarray::const_iterator END = Both_.begin() + bothEnd_;
  for (Sarray::const_iterator si = Both_.begin(); si != END; ++si) {
    hcount += si->n_hydrogens();
    if (debug_ > 0) {
      mprintf("\t  %20s %8i", setup.Top().TruncResAtomName(si->Idx()).c_str(), si->Idx()+1);
      for (Iarray::const_iterator at = si->Hbegin(); at != si->Hend(); ++at)
        mprintf(" %s", setup.Top()[*at].c_str());
      mprintf("\n");
    }
  }
  mprintf("\tDonor-only sites: %zu\n", Both_.size() - bothEnd_);
  for (Sarray::const_iterator si = END; si != Both_.end(); ++si) {
    hcount += si->n_hydrogens();
    if (debug_ > 0) {
      mprintf("\t  %20s %8i", setup.Top().TruncResAtomName(si->Idx()).c_str(), si->Idx()+1);
      for (Iarray::const_iterator at = si->Hbegin(); at != si->Hend(); ++at)
        mprintf(" %s", setup.Top()[*at].c_str());
      mprintf("\n");
    }
  }
  mprintf("\t%u solute hydrogens.\n", hcount);

  // For backwards compat. store donor/acceptor indices for determining data set index.
  if (series_) {
    // Donor hydrogen indices
    Iarray at_array;
    at_array.reserve( hcount );
    for (Sarray::const_iterator site = Both_.begin(); site != Both_.end(); ++site)
      for (Iarray::const_iterator at = site->Hbegin(); at != site->Hend(); ++at)
        at_array.push_back( *at );
    std::sort(at_array.begin(), at_array.end());
    int didx = 0;
    for (Iarray::const_iterator at = at_array.begin(); at != at_array.end(); ++at)
      DidxMap_[*at] = didx++;
    // Acceptor indices
    at_array.clear();
    at_array.reserve( Both_.size() + Acceptor_.size() );
    for (Sarray::const_iterator site = Both_.begin(); site != Both_.end(); ++site)
      at_array.push_back( site->Idx() );
    for (Iarray::const_iterator at = Acceptor_.begin(); at != Acceptor_.end(); ++at)
      at_array.push_back( *at );
    std::sort(at_array.begin(), at_array.end());
    int aidx = 0;
    for (Iarray::const_iterator at = at_array.begin(); at != at_array.end(); ++at)
      AidxMap_[*at] = aidx++;
  }

  // SOLVENT SITE SETUP
  if (calcSolvent_) {
    int at_beg = 0;
    int at_end = 0;
    // Set up solvent donor/acceptor masks
    if (hasSolventDonor_) {
      if (setup.Top().SetupIntegerMask( SolventDonorMask_ )) return Action::ERR;
      if (SolventDonorMask_.None()) {
        mprintf("Warning: SolventDonorMask has no atoms.\n");
        return Action::SKIP;
      }
      at_beg = SolventDonorMask_[0];
      at_end = SolventDonorMask_.back() + 1;
    }
    if (hasSolventAcceptor_) {
      if (setup.Top().SetupIntegerMask( SolventAcceptorMask_ )) return Action::ERR;
      if (SolventAcceptorMask_.None()) {
        mprintf("Warning: SolventAcceptorMask has no atoms.\n");
        return Action::SKIP;
      }
      if (!hasSolventDonor_) {
        at_beg = SolventAcceptorMask_[0];
        at_end = SolventAcceptorMask_.back() + 1;
      } else {
        at_beg = std::min( SolventDonorMask_[0], SolventAcceptorMask_[0] );
        at_end = std::max( SolventDonorMask_.back(), SolventAcceptorMask_.back() ) + 1;
      }
    }
    AtomMask::const_iterator a_atom = SolventAcceptorMask_.begin();
    AtomMask::const_iterator d_atom = SolventDonorMask_.begin();
    bool isDonor, isAcceptor;
    for (int at = at_beg; at != at_end; at++)
    {
      Iarray Hatoms;
      isDonor = false;
      isAcceptor = false;
      if ( d_atom != SolventDonorMask_.end() && *d_atom == at ) {
        ++d_atom;
        Atom const& atom = setup.Top()[at];
        if ( IsFON( atom ) ) {
          for (Atom::bond_iterator H_at = atom.bondbegin();
                                   H_at != atom.bondend(); ++H_at)
            if (setup.Top()[*H_at].Element() == Atom::HYDROGEN)
              Hatoms.push_back( *H_at );
        } else if ( atom.Nbonds() == 0 ) {
          // If no bonds to this atom assume it is an ion. Set the H atom
          // to be the same as D atom; this will skip the angle calc.
          Hatoms.push_back( at );
        }
        isDonor = !Hatoms.empty();
      }
      if ( a_atom != SolventAcceptorMask_.end() && *a_atom == at ) {
        isAcceptor = true;
        ++a_atom;
      }
      if (isDonor || isAcceptor)
        SolventSites_.push_back( Site(at, Hatoms) );
    }

    hcount = 0;
    unsigned int icount = 0;
    mprintf("\tSolvent sites (%zu)\n", SolventSites_.size());
    for (Sarray::const_iterator si = SolventSites_.begin(); si != SolventSites_.end(); ++si) {
      if (si->IsIon())
        icount++;
      else
        hcount += si->n_hydrogens();
      if (debug_ > 0) {
        mprintf("\t  %20s %8i", setup.Top().TruncResAtomName(si->Idx()).c_str(), si->Idx()+1);
        for (Iarray::const_iterator at = si->Hbegin(); at != si->Hend(); ++at)
          mprintf(" %s", setup.Top()[*at].c_str());
        mprintf("\n");
      }
    }
    mprintf("\t%u solvent hydrogens, %u ions.\n", hcount, icount);
  }

  mprintf("\tEstimated max potential memory usage: %s\n",
          MemoryUsage(Both_.size()         * bothEnd_+Acceptor_.size(),
                      SolventSites_.size() * bothEnd_+Acceptor_.size() +
                      Both_.size()         * SolventSites_.size(),
                      masterDSL_->MaxFrames()).c_str());

  if (Both_.empty() && Acceptor_.empty() && SolventSites_.empty()) {
    mprintf("Warning: Nothing selected for hydrogen bond analysis.\n");
    return Action::SKIP;
  }

  // Hbond matrix setup
  if (UU_matrix_byRes_ != 0) {
    // Find highest solute index
    int highest_U_idx = 0;
    for (int rnum = 0; rnum < setup.Top().Nres(); rnum++)
    {
      // Find molecule number
      int at0 = setup.Top().Res(rnum).FirstAtom();
      int mnum = setup.Top()[at0].MolNum();
      if (!setup.Top().Mol(mnum).IsSolvent())
        highest_U_idx = rnum;
    }
    mprintf("\tHighest solute residue # = %i\n", highest_U_idx+1);
    if (UU_matrix_byRes_->AllocateHalf(highest_U_idx+1)) {
      mprinterr("Error: Allocating solute-solute hbond matrix failed.\n");
      return Action::ERR;
    }
    // Set up normalization matrix if necesary
    if (UUmatByRes_norm_ == NORM_RESMAX) {
      UU_norm_byRes_.AllocateHalf( UU_matrix_byRes_->Nrows() );
      for (unsigned int sidx0 = 0; sidx0 < Both_.size(); sidx0++)
      {
        Site const& Site0 = Both_[sidx0];
        int mol0 = (*CurrentParm_)[Site0.Idx()].MolNum();
        int rnum0 = (*CurrentParm_)[Site0.Idx()].ResNum();
        // Loop over solute sites that can be both donor and acceptor
        for (unsigned int sidx1 = sidx0 + 1; sidx1 < bothEnd_; sidx1++)
        {
          Site const& Site1 = Both_[sidx1];
          if (!noIntramol_ || mol0 != (*CurrentParm_)[Site1.Idx()].MolNum()) {
            int rnum1 = (*CurrentParm_)[Site1.Idx()].ResNum();
            // The max # of hbonds between sites depends on # hydrogens.
            // However, given H1-X-H2 Y, you tend to have either Y..H1-X or
            // Y..H2-X, not both, so do based on sites only.
            //UU_norm_byRes_.UpdateElement(rnum0, rnum1, Site0.n_hydrogens() + Site1.n_hydrogens());
            UU_norm_byRes_.UpdateElement(rnum0, rnum1, 2.0);
          }
        } // END loop over solute sites that are D/A
        // Loop over solute acceptor-only
        for (Iarray::const_iterator a_atom = Acceptor_.begin(); a_atom != Acceptor_.end(); ++a_atom)
        {
          if (!noIntramol_ || mol0 != (*CurrentParm_)[*a_atom].MolNum()) {
            int rnum1 = (*CurrentParm_)[*a_atom].ResNum();
            //UU_norm_byRes_.UpdateElement(rnum0, rnum1, Site0.n_hydrogens());
            UU_norm_byRes_.UpdateElement(rnum0, rnum1, 1.0);
          }
        } // END loop over acceptor atoms
      } // END loop over all D/A sites
      mprintf("DEBUG: Residue normalization matrix:\n");
      for (unsigned int rnum0 = 0; rnum0 < UU_norm_byRes_.Nrows(); rnum0++) {
        for (unsigned int rnum1 = rnum0; rnum1 < UU_norm_byRes_.Ncols(); rnum1++) {
          mprintf("\t%6i %6i %f\n", rnum0+1, rnum1+1, UU_norm_byRes_.GetElement(rnum1, rnum0));
        }
      }
    }
  }

  return Action::OK;
}

// Action_HydrogenBond::Angle()
/** Calculate angle between 3 atoms, optionally with imaging. */
double Action_HydrogenBond::Angle(const double* XA, const double* XH, const double* XD, Box const& boxIn) const
{
  if (imageOpt_.ImagingType() == ImageOption::NO_IMAGE)
    return (CalcAngle(XA, XH, XD));
  else {
    double angle;
    Vec3 VH = Vec3(XH);
    Vec3 H_A = MinImagedVec(VH, Vec3(XA), boxIn.UnitCell(), boxIn.FracCell());
    Vec3 H_D = Vec3(XD) - VH;
    double rha = H_A.Magnitude2();
    double rhd = H_D.Magnitude2();
    if (rha > Constants::SMALL && rhd > Constants::SMALL) {
      angle = (H_A * H_D) / sqrt(rha * rhd);
      if      (angle >  1.0) angle =  1.0;
      else if (angle < -1.0) angle = -1.0;
      angle = acos(angle);
    } else
      angle = 0.0;
    return angle;
  }
}

// Action_HydrogenBond::AddUV()
void Action_HydrogenBond::AddUV(double dist, double angle,int fnum, int a_atom,int h_atom,int d_atom,bool udonor, int onum)
{
  int hbidx, solventres, soluteres;
  // TODO: Option to use solvent mol num?
  if (udonor) {
    // Index U-H .. V hydrogen bonds by solute H atom.
    hbidx = h_atom;
    if (bridgeByAtom_)
      soluteres = h_atom;
    else
      soluteres = (*CurrentParm_)[d_atom].ResNum();
    solventres = (*CurrentParm_)[a_atom].ResNum();
  } else {
    // Index U .. H-V hydrogen bonds by solute A atom.
    hbidx = a_atom;
    if (bridgeByAtom_)
      soluteres = a_atom;
    else
      soluteres = (*CurrentParm_)[a_atom].ResNum();
    solventres = (*CurrentParm_)[d_atom].ResNum();
  }
  solvent2solute_[solventres].insert( soluteres );
  UVmapType::iterator it = UV_Map_.lower_bound( hbidx );
  if (it == UV_Map_.end() || it->first != hbidx)
  {
//      mprintf("DBG1: NEW hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
    DataSet_integer* ds = 0;
    if (series_) {
      ds = (DataSet_integer*)
           masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,"solventhb",hbidx));
      if (UVseriesout_ != 0) UVseriesout_->AddDataSet( ds );
    }
    Hbond hb;
    if (udonor) { // Do not care about which solvent acceptor
      if (ds != 0) ds->SetLegend( CreateHBlegend(*CurrentParm_, -1, h_atom, d_atom) );
      hb = Hbond(ds, -1, h_atom, d_atom, splitFrames_);
    } else {           // Do not care about which solvent donor
      if (ds != 0) ds->SetLegend( CreateHBlegend(*CurrentParm_, a_atom, -1, -1) );
      hb = Hbond(ds, a_atom, -1, -1, splitFrames_);
    }
    it = UV_Map_.insert(it, std::pair<int,Hbond>(hbidx,hb));
  } else {
//      mprintf("DBG1: OLD hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
  }
  it->second.Update(dist, angle, fnum, splitFrames_, onum);
}

//  Action_HydrogenBond::CalcSolvHbonds()
/** Calculate hydrogen bonds between solute site and solvent acceptor,
  * or solvent site and solute acceptor.
  */
void Action_HydrogenBond::CalcSolvHbonds(int frameNum, double dist2,
                                         Site const& SiteD, const double* XYZD,
                                         int a_atom,        const double* XYZA,
                                         Frame const& frmIn, int& numHB, bool soluteDonor,
                                         int trajoutNum)
{
  // The two sites are close enough to hydrogen bond.
  int d_atom = SiteD.Idx();
  // Determine if angle cutoff is satisfied
  for (Iarray::const_iterator h_atom = SiteD.Hbegin(); h_atom != SiteD.Hend(); ++h_atom)
  {
    double angle = 0.0;
    bool angleSatisfied = true;
    // For ions, donor atom will be same as h atom so no angle needed.
    if (d_atom != *h_atom) {
      angle = Angle(XYZA, frmIn.XYZ(*h_atom), XYZD, frmIn.BoxCrd());
      angleSatisfied = !(angle < acut_);
    }
    if (angleSatisfied)
    {
#     ifdef _OPENMP
      // numHB holds thread number, will be counted later on.
      thread_HBs_[numHB].push_back( Hbond(sqrt(dist2), angle, a_atom, *h_atom, d_atom, (int)soluteDonor) );
#     else
      ++numHB;
      AddUV(sqrt(dist2), angle, frameNum, a_atom, *h_atom, d_atom, soluteDonor, trajoutNum);
#     endif
    }
  }
}

// Action_HydrogenBond::UU_Set_Idx()
/** Determine solute-solute hbond index for backwards compatibility:
  *   hbidx = (donorIndex * #acceptors) + acceptorIndex
  */
int Action_HydrogenBond::UU_Set_Idx(int a_atom, int h_atom) const {
  IdxMapType::const_iterator it = DidxMap_.find( h_atom );
  int didx = it->second;
  it = AidxMap_.find( a_atom );
  int aidx = it->second;
  return (didx * AidxMap_.size()) + aidx;
}

//  Action_HydrogenBond::UUset()
/** \return solute-solute hydrogen bond time series set with legend set. */
DataSet_integer* Action_HydrogenBond::UUset(int a_atom, int h_atom, int d_atom) {
  DataSet_integer* ds = (DataSet_integer*)
    masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,"solutehb",UU_Set_Idx(a_atom,h_atom)));
  if (UUseriesout_ != 0) UUseriesout_->AddDataSet( ds );
  ds->SetLegend( CreateHBlegend(*CurrentParm_, a_atom, h_atom, d_atom) );
  return ds;
}

// Action_HydrogenBond::AddUU()
/** Add or update a solute-solute hydrogen bond with given angle/distance. */
void Action_HydrogenBond::AddUU(double dist, double angle, int fnum, int a_atom, int h_atom, int d_atom, int onum)
{
  // Index UU hydrogen bonds by DonorH-Acceptor
  Hpair hbidx(h_atom, a_atom);
  UUmapType::iterator it = UU_Map_.lower_bound( hbidx );
  if (it == UU_Map_.end() || it->first != hbidx)
  {
//      mprintf("DBG1: NEW hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
    DataSet_integer* ds = 0;
    if (series_) {
      ds = UUset(a_atom, h_atom, d_atom);
    }
    it = UU_Map_.insert(it, std::pair<Hpair,Hbond>(hbidx,Hbond(ds, a_atom, h_atom, d_atom, splitFrames_)));
  } else {
//      mprintf("DBG1: OLD hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
  }
  it->second.Update(dist, angle, fnum, splitFrames_, onum);
  if (UU_matrix_byRes_ != 0) {
    int a_res = (*CurrentParm_)[a_atom].ResNum();
    int d_res = (*CurrentParm_)[d_atom].ResNum();
    UU_matrix_byRes_->UpdateElement(a_res, d_res, 1.0);
  }
}

// Action_HydrogenBond::CalcSiteHbonds()
/** Calculate hydrogen bonds between given solute donor site and 
  * solute acceptor atom.
  * The distance cutoff should already be satisfied between donor and
  * acceptor heavy atoms.
  */
void Action_HydrogenBond::CalcSiteHbonds(int frameNum, double dist2,
                                         Site const& SiteD, const double* XYZD,
                                         int a_atom,        const double* XYZA,
                                         Frame const& frmIn, int& numHB,
                                         int trajoutNum)
{
  // The two sites are close enough to hydrogen bond.
  int d_atom = SiteD.Idx();
  // Determine if angle cutoff is satisfied
  for (Iarray::const_iterator h_atom = SiteD.Hbegin(); h_atom != SiteD.Hend(); ++h_atom)
  {
    double angle = Angle(XYZA, frmIn.XYZ(*h_atom), XYZD, frmIn.BoxCrd());
    if ( !(angle < acut_) )
    {
#     ifdef _OPENMP
      // numHB holds thread number, will be counted later on.
      thread_HBs_[numHB].push_back( Hbond(sqrt(dist2), angle, a_atom, *h_atom, d_atom) );
#     else
      ++numHB;
      AddUU(sqrt(dist2), angle, frameNum, a_atom, *h_atom, d_atom, trajoutNum);
#     endif
    }
  }
}

// Action_HydrogenBond::DoAction()
Action::RetType Action_HydrogenBond::DoAction(int frameNum, ActionFrame& frm) {
# ifdef TIMER
  t_action_.Start();
# endif
  if (imageOpt_.ImagingEnabled()) {
    //frm.Frm().BoxCrd().PrintDebug("hbond");
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
  }
  // Loop over all solute donor sites
# ifdef TIMER
  t_uu_.Start();
# endif
  int numHB = 0;
  int sidx0;
  int sidx0end = (int)Both_.size();
  if (noIntramol_) {
    // Ignore intramolecular hydrogen bonds
    int mol0 = -1;
#   ifdef _OPENMP
    // Use numHB to track thread. Will be actually counted after the parallel section.
#   pragma omp parallel private(sidx0, numHB) firstprivate(mol0)
    {
    numHB = omp_get_thread_num();
#   pragma omp for
#   endif
    for (sidx0 = 0; sidx0 < sidx0end; sidx0++)
    {
      Site const& Site0 = Both_[sidx0];
      const double* XYZ0 = frm.Frm().XYZ( Site0.Idx() );
      mol0 = (*CurrentParm_)[Site0.Idx()].MolNum(); 
      // Loop over solute sites that can be both donor and acceptor
      for (unsigned int sidx1 = sidx0 + 1; sidx1 < bothEnd_; sidx1++)
      {
        Site const& Site1 = Both_[sidx1];
        if (mol0 != (*CurrentParm_)[Site1.Idx()].MolNum()) {
          const double* XYZ1 = frm.Frm().XYZ( Site1.Idx() );
          double dist2 = DIST2( imageOpt_.ImagingType(), XYZ0, XYZ1, frm.Frm().BoxCrd() );
          if ( !(dist2 > dcut2_) )
          {
            // Site 0 donor, Site 1 acceptor
            CalcSiteHbonds(frameNum, dist2, Site0, XYZ0, Site1.Idx(), XYZ1, frm.Frm(), numHB, frm.TrajoutNum());
            // Site 1 donor, Site 0 acceptor
            CalcSiteHbonds(frameNum, dist2, Site1, XYZ1, Site0.Idx(), XYZ0, frm.Frm(), numHB, frm.TrajoutNum());
          }
        }
      }
      // Loop over solute acceptor-only
      for (Iarray::const_iterator a_atom = Acceptor_.begin(); a_atom != Acceptor_.end(); ++a_atom)
      {
        if (mol0 != (*CurrentParm_)[*a_atom].MolNum()) {
          const double* XYZ1 = frm.Frm().XYZ( *a_atom );
          double dist2 = DIST2( imageOpt_.ImagingType(), XYZ0, XYZ1, frm.Frm().BoxCrd() );
          if ( !(dist2 > dcut2_) )
            CalcSiteHbonds(frameNum, dist2, Site0, XYZ0, *a_atom, XYZ1, frm.Frm(), numHB, frm.TrajoutNum());
        }
      }
    }
#   ifdef _OPENMP
    } // END pragma omp parallel, nointramol
#   endif
  } else {
    // All hydrogen bonds
#   ifdef _OPENMP
    // Use numHB to track thread. Will be actually counted after the parallel section.
#   pragma omp parallel private(sidx0, numHB)
    {
    numHB = omp_get_thread_num();
#   pragma omp for
#   endif
    for (sidx0 = 0; sidx0 < sidx0end; sidx0++)
    {
      Site const& Site0 = Both_[sidx0];
      const double* XYZ0 = frm.Frm().XYZ( Site0.Idx() );
      // Loop over solute sites that can be both donor and acceptor
      for (unsigned int sidx1 = sidx0 + 1; sidx1 < bothEnd_; sidx1++)
      {
        Site const& Site1 = Both_[sidx1];
        const double* XYZ1 = frm.Frm().XYZ( Site1.Idx() );
        double dist2 = DIST2( imageOpt_.ImagingType(), XYZ0, XYZ1, frm.Frm().BoxCrd() );
        if ( !(dist2 > dcut2_) )
        {
          // Site 0 donor, Site 1 acceptor
          CalcSiteHbonds(frameNum, dist2, Site0, XYZ0, Site1.Idx(), XYZ1, frm.Frm(), numHB, frm.TrajoutNum());
          // Site 1 donor, Site 0 acceptor
          CalcSiteHbonds(frameNum, dist2, Site1, XYZ1, Site0.Idx(), XYZ0, frm.Frm(), numHB, frm.TrajoutNum());
        }
      }
      // Loop over solute acceptor-only
      for (Iarray::const_iterator a_atom = Acceptor_.begin(); a_atom != Acceptor_.end(); ++a_atom)
      {
        const double* XYZ1 = frm.Frm().XYZ( *a_atom );
        double dist2 = DIST2( imageOpt_.ImagingType(), XYZ0, XYZ1, frm.Frm().BoxCrd() );
        if ( !(dist2 > dcut2_) )
          CalcSiteHbonds(frameNum, dist2, Site0, XYZ0, *a_atom, XYZ1, frm.Frm(), numHB, frm.TrajoutNum());
      }
    }
#   ifdef _OPENMP
    } // END pragma omp parallel
#    endif
  } // END if nointramol

# ifdef _OPENMP
  // Add all found hydrogen bonds
  numHB = 0; 
  for (std::vector<Harray>::iterator it = thread_HBs_.begin(); it != thread_HBs_.end(); ++it) {
    numHB += (int)it->size();
    for (Harray::const_iterator hb = it->begin(); hb != it->end(); ++hb)
      AddUU(hb->Dist(), hb->Angle(), frameNum, hb->A(), hb->H(), hb->D(), frm.TrajoutNum());
    it->clear();
  }
# endif
  NumHbonds_->Add(frameNum, &numHB);
# ifdef TIMER
  t_uu_.Stop();
# endif
  // Loop over all solvent sites
  if (calcSolvent_) {
#   ifdef TIMER
    t_uv_.Start();
#   endif
    solvent2solute_.clear();
    numHB = 0;
    int vidx;
    int vidxend = (int)SolventSites_.size();
#   ifdef _OPENMP
    // Use numHB to track thread. Will be actually counted after the parallel section.
#   pragma omp parallel private(vidx, numHB)
    {
    numHB = omp_get_thread_num();
#   pragma omp for
#   endif
    for (vidx = 0; vidx < vidxend; vidx++)
    {
      Site const& Vsite = SolventSites_[vidx];
      const double* VXYZ = frm.Frm().XYZ( Vsite.Idx() );
      // Loop over solute sites that can be both donor and acceptor
      for (unsigned int sidx = 0; sidx < bothEnd_; sidx++)
      {
        const double* UXYZ = frm.Frm().XYZ( Both_[sidx].Idx() );
        double dist2 = DIST2( imageOpt_.ImagingType(), VXYZ, UXYZ, frm.Frm().BoxCrd() );
        if ( !(dist2 > dcut2_) )
        {
          // Solvent site donor, solute site acceptor
          CalcSolvHbonds(frameNum, dist2, Vsite, VXYZ, Both_[sidx].Idx(), UXYZ, frm.Frm(), numHB, false, frm.TrajoutNum());
          // Solvent site acceptor, solute site donor
          CalcSolvHbonds(frameNum, dist2, Both_[sidx], UXYZ, Vsite.Idx(), VXYZ, frm.Frm(), numHB, true, frm.TrajoutNum());
        }
      }
      // Loop over solute sites that are donor only
      for (unsigned int sidx = bothEnd_; sidx < Both_.size(); sidx++)
      {
        const double* UXYZ = frm.Frm().XYZ( Both_[sidx].Idx() );
        double dist2 = DIST2( imageOpt_.ImagingType(), VXYZ, UXYZ, frm.Frm().BoxCrd() );
        if ( !(dist2 > dcut2_) )
          // Solvent site acceptor, solute site donor
          CalcSolvHbonds(frameNum, dist2, Both_[sidx], UXYZ, Vsite.Idx(), VXYZ, frm.Frm(), numHB, true, frm.TrajoutNum());
      }
      // Loop over solute sites that are acceptor only
      for (Iarray::const_iterator a_atom = Acceptor_.begin(); a_atom != Acceptor_.end(); ++a_atom)
      {
        const double* UXYZ = frm.Frm().XYZ( *a_atom );
        double dist2 = DIST2( imageOpt_.ImagingType(), VXYZ, UXYZ, frm.Frm().BoxCrd() );
        if ( !(dist2 > dcut2_) )
          // Solvent site donor, solute site acceptor
          CalcSolvHbonds(frameNum, dist2, Vsite, VXYZ, *a_atom, UXYZ, frm.Frm(), numHB, false, frm.TrajoutNum());
      }
    } // END loop over solvent sites
#   ifdef _OPENMP
    } // END pragma omp parallel
    // Add all found hydrogen bonds
    numHB = 0; 
    for (std::vector<Harray>::iterator it = thread_HBs_.begin(); it != thread_HBs_.end(); ++it) {
      numHB += (int)it->size();
      for (Harray::const_iterator hb = it->begin(); hb != it->end(); ++hb)
        AddUV(hb->Dist(), hb->Angle(), frameNum, hb->A(), hb->H(), hb->D(), (bool)hb->Frames(), frm.TrajoutNum());
      it->clear();
    }
#   endif
    NumSolvent_->Add(frameNum, &numHB);
#   ifdef TIMER
    t_uv_.Stop();
#   endif
    // Determine number of bridging waters
#   ifdef TIMER
    t_bridge_.Start();
#   endif
    numHB = 0;
    std::string bridgeID;
    for (RmapType::const_iterator bridge = solvent2solute_.begin();
                                  bridge != solvent2solute_.end(); ++bridge)
    {
      // bridge->first is solvent residue number.
      // bridge->second is a set of solute residue numbers the solvent
      // residue is bound to.
      // If solvent molecule is bound to 2 or more different solute residues,
      // it is bridging. 
      if ( bridge->second.size() > 1) {
        bool isBridge = true;
        if (noIntramol_) {
          // If all residues belong to the same molecule and 'nointramol',
          // do not consider this bridge.
          int firstmol = -1;
          unsigned int nequal = 1;
          for (std::set<int>::const_iterator res = bridge->second.begin();
                                             res != bridge->second.end(); ++res)
          {
            int currentMol;
            if (bridgeByAtom_)
              currentMol = (*CurrentParm_)[*res].MolNum();
            else
              currentMol = (*CurrentParm_)[CurrentParm_->Res(*res).FirstAtom()].MolNum();
            if ( firstmol == -1 )
              firstmol = currentMol;
            else if (currentMol == firstmol)
              ++nequal;
          }
          isBridge = (nequal < bridge->second.size());
        }
        if (isBridge) {
          // numHB is used to track the number of bridges
          ++numHB;
          // Bridging Solvent residue number
          bridgeID.append(integerToString( bridge->first+1 ) + "(");
          // Loop over solute residues this solvent is bound to.
          for (std::set<int>::const_iterator res = bridge->second.begin();
                                             res != bridge->second.end(); ++res)
            // Solute residue number being bridged
            bridgeID.append( integerToString( *res+1 ) + "+" );
          bridgeID.append("),");
          // Find bridge in map based on this combo of residues (bridge->second)
          BmapType::iterator b_it = BridgeMap_.lower_bound( bridge->second );
          if (b_it == BridgeMap_.end() || b_it->first != bridge->second) {
            // New Bridge
            DataSet_integer* bds = 0; 
            if (Bseries_) {
              bds = (DataSet_integer*)
                masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,CreateBridgeLegend("bridge",bridge->second),BridgeMap_.size()));
              // Create a legend from the indices.
              bds->SetLegend( CreateBridgeLegend( "B", bridge->second ) );
              if (Bseriesout_ != 0) Bseriesout_->AddDataSet( bds );
            }
            b_it = BridgeMap_.insert( b_it, std::pair<std::set<int>,Bridge>(bridge->second, Bridge(bds, splitFrames_)) );
          }
          // Increment bridge #frames
          b_it->second.Update(frameNum, splitFrames_, frm.TrajoutNum());
        }
      }
    } // END LOOP OVER solvent2solute_
    if (bridgeID.empty())
      bridgeID.assign("None");
    NumBridge_->Add(frameNum, &numHB);
    BridgeID_->Add(frameNum, bridgeID.c_str());
#   ifdef TIMER
    t_bridge_.Stop();
#   endif
  } // END if calcSolvent_

  Nframes_++;
# ifdef TIMER
  t_action_.Stop();
# endif
  return Action::OK;
}

/** Ensure DataSet time series has at least N frames, fill out if necessary. */
void Action_HydrogenBond::FinishSeries(DataSet_integer* data, unsigned int N) {
  static const int ZERO = 0;
  if (data != 0 && N > 0) {
    if ( data->Size() < N ) {
      data->Add( N-1, &ZERO );
#     ifdef MPI
      data->SetNeedsSync( false );
#     endif
    }
  }
}

#ifdef MPI
// Action_HydrogenBond::GetRankNhbonds()
/** Determine how many hydrogen bonds are on each rank. */
std::vector<int> Action_HydrogenBond::GetRankNhbonds( int num_hb, Parallel::Comm const& commIn )
{
  std::vector<int> nhb_on_rank;
  if (commIn.Master())
    nhb_on_rank.resize( commIn.Size() );
  commIn.GatherMaster( &num_hb, 1, MPI_INT, &nhb_on_rank[0] );
  return nhb_on_rank;
}

/** Add Hbond class to flat arrays. */
void Action_HydrogenBond::HbondToArray(std::vector<double>& Dvals, std::vector<int>& Ivals, Hbond const& hb)
{
  Dvals.push_back( hb.Dist() );
  Dvals.push_back( hb.Angle() );
  for (unsigned int idx = 0; idx != hb.Nparts(); idx++) {
    Dvals.push_back( hb.PartDist(idx).nData() );
    Dvals.push_back( hb.PartDist(idx).mean() );
    Dvals.push_back( hb.PartDist(idx).M2() );
    Dvals.push_back( hb.PartAngle(idx).nData() );
    Dvals.push_back( hb.PartAngle(idx).mean() );
    Dvals.push_back( hb.PartAngle(idx).M2() );
  }
  Ivals.push_back( hb.A() );
  Ivals.push_back( hb.H() );
  Ivals.push_back( hb.D() );
  Ivals.push_back( hb.Frames() );
}

/** PARALLEL NOTES:
  * The following tags are used for MPI send/receive:
  *   1300  : Array containing hbond double info on rank.
  *   1301  : Array containing hbond integer info on rank.
  *   1302  : Number of bridges to expect from rank.
  *   1303  : Array containing bridge integer info on rank.
  *   1304+X: Array of hbond X series info from rank.
  */
int Action_HydrogenBond::SyncAction() {
  // Make sure all time series are updated at this point.
  UpdateSeries();
  // TODO consolidate # frames / offset calc code with Action_NAstruct
  // Get total number of frames
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &Nframes_, 1, MPI_INT, &rank_frames[0] );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      Nframes_ += rank_frames[rank];
  }
  // Convert rank frames to offsets.
  std::vector<int> rank_offsets( trajComm_.Size(), 0 );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      rank_offsets[rank] = rank_offsets[rank-1] + rank_frames[rank-1];
  }

  // Need to send hbond data from all ranks to master.
  std::vector<double> Dvals;           // Hold dist_ and angle_ for each hbond, (as well as n_/mean_/M2_ for dist/angle each part)
  std::vector<int> Ivals;              // Hold A_, H_, D_, and frames_ for each hbond
  unsigned int dvalsPerHbond;
  unsigned int nParts;
  if (!splitFrames_.empty()) {
    nParts = splitFrames_.size() + 1;
    dvalsPerHbond = 2 + (nParts * 6);
  } else {
    nParts = 0;
    dvalsPerHbond = 2;
  }
  // Need to know how many hbonds on each process.
  std::vector<int> nuu_on_rank = GetRankNhbonds( UU_Map_.size(), trajComm_ );
  std::vector<int> nuv_on_rank = GetRankNhbonds( UV_Map_.size(), trajComm_ );
  if (trajComm_.Master()) {
    // MASTER RANK
    for (int rank = 1; rank < trajComm_.Size(); rank++)
    {
      int n_total_on_rank = nuu_on_rank[rank] + nuv_on_rank[rank];
      if (n_total_on_rank > 0)
      {
        std::vector<DataSet_integer*> Svals;
        if (series_) Svals.reserve( n_total_on_rank );
        Dvals.resize( dvalsPerHbond * n_total_on_rank );
        Ivals.resize( 4             * n_total_on_rank );
        trajComm_.Recv( &(Dvals[0]), Dvals.size(), MPI_DOUBLE, rank, 1300 );
        trajComm_.Recv( &(Ivals[0]), Ivals.size(), MPI_INT,    rank, 1301 );
        // Loop over all received hydrogen bonds
        // Dvals = dist, angle
        // Avals = A, H, D, frames
        const int* IV = &Ivals[0];
        const double* DV = &Dvals[0];
        // UU Hbonds
        for (int in = 0; in != nuu_on_rank[rank]; in++, IV += 4, DV += dvalsPerHbond)
        {
          Hpair hbidx(IV[1], IV[0]);
          UUmapType::iterator it = UU_Map_.lower_bound( hbidx );
          DataSet_integer* ds = 0;
          if (it == UU_Map_.end() || it->first != hbidx)
          {
            // Hbond on rank that has not been found on master
            if (series_)
              ds = UUset(IV[0], IV[1], IV[2]);
            it = UU_Map_.insert(it, std::pair<Hpair,Hbond>(hbidx,Hbond(DV[0],DV[1],ds,IV[0],IV[1],IV[2],IV[3])));
            it->second.SetupParts(nParts, DV+2);
          } else {
            // Hbond on rank and master. Update on master.
            it->second.Combine(DV[0], DV[1], IV[3]);
            it->second.CombineParts(nParts, DV+2);
            ds = it->second.Data();
          }
          Svals.push_back( ds );
        }
        // UV Hbonds
        for (int in = 0; in != nuv_on_rank[rank]; in++, IV += 4, DV += dvalsPerHbond)
        {
          int hbidx;
          if (IV[1] != -1)
            // Index U-H .. V hydrogen bonds by solute H atom.
            hbidx = IV[1];
          else
            // Index U .. H-V hydrogen bonds by solute A atom.
            hbidx = IV[0];
          DataSet_integer* ds = 0;
          UVmapType::iterator it = UV_Map_.lower_bound( hbidx );
          if (it == UV_Map_.end() || it->first != hbidx)
          {
            // Hbond on rank that has not been found on master
            if (series_) {
              ds = (DataSet_integer*)
                   masterDSL_->AddSet(DataSet::INTEGER, MetaData(hbsetname_,"solventhb",hbidx));
              ds->SetLegend( CreateHBlegend(*CurrentParm_, IV[0], IV[1], IV[2]) );
              if (UVseriesout_ != 0) UVseriesout_->AddDataSet( ds );
            }
            it = UV_Map_.insert(it, std::pair<int,Hbond>(hbidx,Hbond(DV[0],DV[1],ds,IV[0],IV[1],IV[2],IV[3])));
            it->second.SetupParts(nParts, DV+2);
          } else {
            // Hbond on rank and master. Update on master.
            it->second.Combine(DV[0], DV[1], IV[3]);
            it->second.CombineParts(nParts, DV+2);
            ds = it->second.Data();
          }
          Svals.push_back( ds );
        }
        // Update all time series
        if (series_) {
          for (int in = 0; in != n_total_on_rank; in++) {
            DataSet_integer* ds = Svals[in]; 
            //ds->Resize( Nframes_ );
            //int* d_beg = ds->Ptr() + rank_offsets[ rank ];
            //rprintf("Resizing hbond series data to %i, starting frame %i, # frames %i from rank %i (%i)\n",
            //        Nframes_, rank_offsets[rank], rank_frames[rank], rank, 1304 + in);
            ds->Recv(Nframes_, rank_offsets[ rank ], rank_frames[ rank ],
                     rank, 1304 + in, trajComm_);
            //trajComm_.Recv( d_beg, rank_frames[ rank ], MPI_INT, rank, 1304 + in );
            ds->SetNeedsSync( false );
          }
        }
      }
    } // END master loop over ranks
    // At this point we have all hbond sets from all ranks. Mark all HB sets
    // smaller than Nframes_ as synced and ensure the time series has been
    // updated to reflect overall # frames.
    if (series_) {
      for (UUmapType::iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb)
        FinishSeries( hb->second.Data(), Nframes_ );
      for (UVmapType::iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb)
        FinishSeries( hb->second.Data(), Nframes_ );
    }
  } else {
    // NON-MASTER RANK
    if (!UU_Map_.empty()) {
      unsigned int ntotal = UU_Map_.size() + UV_Map_.size();
      Dvals.reserve( ntotal * dvalsPerHbond );
      Ivals.reserve( ntotal * 4 );
      // Store UU bonds in flat arrays.
      for (UUmapType::const_iterator it = UU_Map_.begin(); it != UU_Map_.end(); ++it) {
        HbondToArray(Dvals, Ivals, it->second);
      }
      // Store UV bonds in flat arrays
      for (UVmapType::const_iterator it = UV_Map_.begin(); it != UV_Map_.end(); ++it) {
        HbondToArray(Dvals, Ivals, it->second);
      }
      trajComm_.Send( &(Dvals[0]), Dvals.size(), MPI_DOUBLE, 0, 1300 );
      trajComm_.Send( &(Ivals[0]), Ivals.size(), MPI_INT,    0, 1301 );
      // Send series data to master
      if (series_) {
        int in = 0; // For tag
        for (UUmapType::const_iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb, in++) {
          //rprintf("Sending %zu frames to master (%i).\n", hb->second.Data()->Size(), 1304+in);
          hb->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          hb->second.Data()->SetNeedsSync( false );
        }
        for (UVmapType::const_iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb, in++) {
          //rprintf("Sending %zu frames to master (%i).\n", hb->second.Data()->Size(), 1304+in);
          hb->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          hb->second.Data()->SetNeedsSync( false );
        }
      }
    }
  } // END COMMUNICATING HBOND DATA TO MASTER

  if (calcSolvent_) {
    // Sync bridging data
    // iArray will contain for each bridge: Nres, res1, ..., resN, Frames[, Npart1, ..., NpartN]
    std::vector<int> iArray;
    int iSize;
    if (trajComm_.Master()) {
      // MASTER RANK
      for (int rank = 1; rank < trajComm_.Size(); rank++)
      {
        std::vector<DataSet_integer*> Svals;
        // Receive size of iArray
        trajComm_.Recv( &iSize,           1, MPI_INT, rank, 1302 );
        //mprintf("DEBUG: Receiving %i bridges from rank %i\n", iSize, rank);
        iArray.resize( iSize );
        trajComm_.Recv( &(iArray[0]), iSize, MPI_INT, rank, 1303 );
        unsigned int idx = 0;
        while (idx < iArray.size()) {
          std::set<int> residues;
          unsigned int i2 = idx + 1;
          for (int ir = 0; ir != iArray[idx]; ir++, i2++)
            residues.insert( iArray[i2] );
          BmapType::iterator b_it = BridgeMap_.lower_bound( residues );
          DataSet_integer* bds = 0;
          if (b_it == BridgeMap_.end() || b_it->first != residues ) {
            // Bridge not found on master. Create new Bridge.
            if (Bseries_) {
              bds = (DataSet_integer*)
                masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,CreateBridgeLegend("bridge",residues),BridgeMap_.size()));
                // Create a legend from the indices.
                bds->SetLegend( CreateBridgeLegend( "B", residues ) );
            }
            b_it = BridgeMap_.insert( b_it, std::pair<std::set<int>,Bridge>(residues, Bridge(bds, iArray[i2])) );
            b_it->second.SetupParts(nParts, &iArray[0] + i2 + 1);
            //if (bds != 0) mprintf("DEBUG: '%s' was not on master.\n", bds->legend());
            if (Bseriesout_ != 0) Bseriesout_->AddDataSet( bds );
          } else {
            // Bridge on master and rank. Increment bridge #frames.
            bds = b_it->second.Data();
            b_it->second.Combine( iArray[i2] );
            b_it->second.CombineParts(nParts, &iArray[0] + i2 + 1);
            //if (bds != 0) mprintf("DEBUG: '%s' was already on master.\n", bds->legend());
          }
          Svals.push_back( bds );
          idx = i2 + 1 + nParts;
        }
        // Update all time series
        if (Bseries_) {
          for (unsigned int in = 0; in != Svals.size(); in++) {
            DataSet_integer* ds = Svals[in]; 
            //ds->Resize( Nframes_ );
            //int* d_beg = ds->Ptr() + rank_offsets[ rank ];
            //rprintf("Receiving %i frames of bridge series data for %s, starting frame %i, # frames %i from rank %i (%i)\n",
            //        Nframes_, ds->legend(), rank_offsets[rank], rank_frames[rank], rank, 1304 + in);
            ds->Recv(Nframes_, rank_offsets[ rank ], rank_frames[ rank ],
                     rank, 1304 + in, trajComm_);
            //trajComm_.Recv( d_beg, rank_frames[ rank ], MPI_INT, rank, 1304 + in );
            ds->SetNeedsSync( false );
          }
        }
      } // END LOOP OVER MASTER RANKS
      // At this point we have all bridges from all ranks. Mark all bridge sets
      // smaller than Nframes_ as synced and ensure the time series has been
      // updated to reflect overall # frames.
      if (Bseries_) {
        for (BmapType::iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
          FinishSeries( b->second.Data(), Nframes_ );
      }
    } else {
       // NON-MASTER
       // Construct bridge info array.
       for (BmapType::const_iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
       {
         iArray.push_back( b->first.size() ); // # of bridging res
         for ( std::set<int>::const_iterator r = b->first.begin(); r != b->first.end(); ++r)
           iArray.push_back( *r ); // Bridging res
         iArray.push_back( b->second.Frames() ); // # frames
         if (nParts > 0) {
           for (unsigned int part = 0; part != nParts; part++)
             iArray.push_back( b->second.PartFrames(part) );
         }
      }
      // Since the size of each bridge can be different (i.e. differing #s of
      // residues may be bridged), first send size of the transport array.
      iSize = (int)iArray.size();
      trajComm_.Send( &iSize,           1, MPI_INT, 0, 1302 );
      trajComm_.Send( &(iArray[0]), iSize, MPI_INT, 0, 1303 );
      // Send series data to master
      if (Bseries_) {
        int in = 0; // For tag
        for (BmapType::const_iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b, in++) {
          //rprintf("Sending %zu frames of %s to master (%i).\n", b->second.Data()->Size(), b->second.Data()->legend(), 1304+in);
          b->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          b->second.Data()->SetNeedsSync( false );
        }
      }
    }
  } // END COMMUNICATING BRIDGE DATA TO MASTER
  return 0;
}
#endif /* MPI */

// Action_HydrogenBond::UpdateSeries()
/** Ensure all time series data is up-to-date with Nframes.
  * Should only be called once.
  */
void Action_HydrogenBond::UpdateSeries() {
  if (seriesUpdated_) return;
  if (series_ && Nframes_ > 0) {
    for (UUmapType::iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb)
      FinishSeries(hb->second.Data(), Nframes_);
    for (UVmapType::iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb)
      FinishSeries(hb->second.Data(), Nframes_);
  }
  if (Bseries_ && Nframes_ > 0) {
    for (BmapType::iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
      FinishSeries( b->second.Data(), Nframes_ );
  }
  seriesUpdated_ = true;
}

// Action_Hbond::MemoryUsage()
/** Estimate the memory usage of the hbond command. */
std::string Action_HydrogenBond::MemoryUsage(size_t n_uu_pairs, size_t n_uv_pairs, size_t nFrames) const
{
  static const size_t sizeHbond = sizeof(Hbond);
  // NOTE: Assuming an overhead of 32 bytes per map element.
  static const size_t sizeElt = 32;
  static const size_t sizeUUmapElt = sizeElt + sizeof(Hpair) + sizeHbond;
  static const size_t sizeUVmapElt = sizeElt + sizeof(int) + sizeHbond;
  static const size_t sizeBRmapElt = sizeElt + sizeof(std::set<int>) + sizeof(Bridge);
  // Solute-solute pairs
  size_t memTotal = sizeof(UUmapType) + (n_uu_pairs * sizeUUmapElt);
  // Solute-solvent pairs
  memTotal += sizeof(UVmapType) + (n_uv_pairs * sizeUVmapElt);
  // Time series TODO bridge series
  if (series_ && nFrames > 0) {
    size_t seriesSet = (nFrames * sizeof(int)) + sizeof(DataSet_integer);
    memTotal += (seriesSet * (n_uu_pairs + n_uv_pairs));
  }
  // Solute-solvent bridges
  // Cannot really estimate bridging, so always base it on BridgeMap_
  memTotal += sizeof(BmapType);
  for (BmapType::const_iterator it = BridgeMap_.begin(); it != BridgeMap_.end(); ++it)
    memTotal += (sizeBRmapElt + it->first.size()*sizeof(int));
  // Matrices
  if (UU_matrix_byRes_ != 0)
    memTotal += UU_matrix_byRes_->MemUsageInBytes();
 
  return ByteString( memTotal, BYTE_DECIMAL );
}

/** Print header for summary by parts. */
void Action_HydrogenBond::summary_Parts_header(CpptrajFile* avgout, unsigned int nParts)
{
  if (nParts < 1) return;
  for (unsigned int idx = 0; idx != nParts; idx++) {
    std::string spart(integerToString(idx+1));
    std::string frames( "Frames"  + spart);
    std::string frac(   "Frac"    + spart);
    std::string avgdist("AvgDist" + spart);
    std::string avgang( "AvgAng"  + spart);
    avgout->Printf(" %8s %12s %12s %12s", frames.c_str(), frac.c_str(), avgdist.c_str(), avgang.c_str());
  }
}

/** Print summary by parts for given hbond. */
void Action_HydrogenBond::summary_Parts(CpptrajFile* avgout, Hbond const& hb) const {
  for (unsigned int idx = 0; idx != hb.Nparts(); idx++)
    avgout->Printf(" %8i %12.4f %12.4f %12.4f",
                   hb.PartFrames(idx), hb.PartFrac(idx, Nframes_),
                   hb.PartDist(idx).mean(), hb.PartAngle(idx).mean()*Constants::RADDEG);
}

/** Used to associate Bridge with solute atoms/residues and sort Bridges by frames. */
class Action_HydrogenBond::bridgeSorter {
  public:
    /// CONSTRUCTOR - list of solute atoms/residues, Bridge info
    bridgeSorter(std::set<int> const& uIdx, Bridge const& bridge) :
      uIdx_(uIdx), bridge_(bridge) {}
    /// Used to sort by bridge # frames
    bool operator<(bridgeSorter const& rhs) const {
      if (bridge_.Frames() == rhs.bridge_.Frames())
        return (uIdx_ < rhs.uIdx_);
      else
        return (bridge_.Frames() > rhs.bridge_.Frames());
    }
    /// \return List of solute atom/residue #s
    std::set<int> const& Uidx() const { return uIdx_; }
    /// \return Bridging info
    Bridge const& Binfo()       const { return bridge_; }
  private:
    std::set<int> uIdx_; ///< Hold solute atom/residue #s
    Bridge bridge_;      ///< Hold bridging info.
};

// Action_HydrogenBond::Print()
/** Print average occupancies over all frames for all detected Hbonds. */
void Action_HydrogenBond::Print() {
  Harray HbondList; // For sorting
  std::string Aname, Hname, Dname;

  // Final memory usage
  mprintf("    HBOND: Actual memory usage is %s\n",
          MemoryUsage(UU_Map_.size(), UV_Map_.size(), Nframes_).c_str());
  mprintf("\t%zu solute-solute hydrogen bonds.\n", UU_Map_.size());
  if (calcSolvent_) {
   mprintf("\t%zu solute-solvent hydrogen bonds.\n", UV_Map_.size());
   mprintf("\t%zu unique solute-solvent bridging interactions.\n", BridgeMap_.size());
  }
# ifdef TIMER
  t_uu_.WriteTiming(      2,"Solute-Solute   :",t_action_.Total());
  if (calcSolvent_) {
    t_uv_.WriteTiming(    2,"Solute-Solvent  :",t_uv_.Total());
    t_bridge_.WriteTiming(2,"Bridging waters :",t_action_.Total());
  }
  t_action_.WriteTiming(1,"Total:");
# endif
  // Ensure all series have been updated for all frames.
  UpdateSeries();
  // Matrix normalization
  if (UU_matrix_byRes_ != 0) {
    if (UUmatByRes_norm_ == NORM_FRAMES) {
      double norm = 1.0 / ((double)Nframes_);
      for (unsigned int r = 0; r != UU_matrix_byRes_->Nrows(); r++)
        for (unsigned int c = 0; c != UU_matrix_byRes_->Ncols(); c++)
          UU_matrix_byRes_->SetElement(c, r, UU_matrix_byRes_->GetElement(c, r) * norm);
    }
  }

  if (CurrentParm_ == 0) return;
  // Calculate necessary column width for strings based on how many residues.
  // ResName+'_'+ResNum+'@'+AtomName | NUM = 4+1+R+1+4 = R+10
  int NUM = DigitWidth( CurrentParm_->Nres() ) + 10;
  // If useAtomNum_ +'_'+AtomNum += 1+A
  if (useAtomNum_) NUM += ( DigitWidth( CurrentParm_->Natom() ) + 1 );

  // Solute Hbonds 
  if (avgout_ != 0) { 
    // Place all detected Hbonds in a list and sort.
    for (UUmapType::const_iterator it = UU_Map_.begin(); it != UU_Map_.end(); ++it) {
      HbondList.push_back( it->second );
      // Calculate average distance and angle for this hbond.
      HbondList.back().CalcAvg();
    }
    UU_Map_.clear();
    // Sort and Print 
    sort( HbondList.begin(), HbondList.end() );
    avgout_->Printf("%-*s %*s %*s %8s %12s %12s %12s", NUM, "#Acceptor", 
                    NUM, "DonorH", NUM, "Donor", "Frames", "Frac", "AvgDist", "AvgAng");
    if (!splitFrames_.empty())
      summary_Parts_header(avgout_, splitFrames_.size()+1);
    avgout_->Printf("\n");
    for (Harray::const_iterator hbond = HbondList.begin(); hbond != HbondList.end(); ++hbond ) 
    {
      double avg = ((double)hbond->Frames()) / ((double) Nframes_);
      Aname = CurrentParm_->TruncResAtomName(hbond->A());
      Hname = CurrentParm_->TruncResAtomName(hbond->H());
      Dname = CurrentParm_->TruncResAtomName(hbond->D());
      if (useAtomNum_) {
        Aname.append("_" + integerToString(hbond->A()+1));
        Hname.append("_" + integerToString(hbond->H()+1));
        Dname.append("_" + integerToString(hbond->D()+1));
      }
      avgout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     hbond->Frames(), avg, hbond->Dist(), hbond->Angle());
      if (!splitFrames_.empty())
        summary_Parts(avgout_, *hbond);
      avgout_->Printf("\n");
    }
  }

  // Solute-solvent Hbonds 
  if (solvout_ != 0 && calcSolvent_) {
    HbondList.clear();
    for (UVmapType::const_iterator it = UV_Map_.begin(); it != UV_Map_.end(); ++it) {
      HbondList.push_back( it->second );
      // Calculate average distance and angle for this hbond.
      HbondList.back().CalcAvg();
    }
    UV_Map_.clear();
    sort( HbondList.begin(), HbondList.end() );
    // Calc averages and print
    solvout_->Printf("#Solute-Solvent Hbonds:\n");
    solvout_->Printf("%-*s %*s %*s %8s %12s %12s %12s", NUM, "#Acceptor", 
                     NUM, "DonorH", NUM, "Donor", "Count", "Frac", "AvgDist", "AvgAng");
    if (!splitFrames_.empty())
      summary_Parts_header(solvout_, splitFrames_.size()+1);
    solvout_->Printf("\n");
    for (Harray::const_iterator hbond = HbondList.begin(); hbond != HbondList.end(); ++hbond )
    {
      // Average has slightly diff meaning since for any given frame multiple
      // solvent can bond to the same solute.
      double avg = ((double)hbond->Frames()) / ((double) Nframes_);
      if (hbond->A()==-1) // Solvent acceptor
        Aname = "SolventAcc";
      else {
        Aname = CurrentParm_->TruncResAtomName(hbond->A());
        if (useAtomNum_) Aname.append("_" + integerToString(hbond->A()+1));
      }
      if (hbond->D()==-1) { // Solvent donor
        Dname = "SolventDnr";
        Hname = "SolventH";
      } else {
        Dname = CurrentParm_->TruncResAtomName(hbond->D());
        Hname = CurrentParm_->TruncResAtomName(hbond->H());
        if (useAtomNum_) {
          Dname.append("_" + integerToString(hbond->D()+1));
          Hname.append("_" + integerToString(hbond->H()+1));
        }
      }
      solvout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     hbond->Frames(), avg, hbond->Dist(), hbond->Angle());
      if (!splitFrames_.empty())
        summary_Parts(solvout_, *hbond);
      solvout_->Printf("\n");
    }
    HbondList.clear();
  }

  // BRIDGING INFO
  if (bridgeout_ != 0 && calcSolvent_) {
    if (bridgeByAtom_)
      bridgeout_->Printf("#Bridging Solute Atoms:\n");
    else
      bridgeout_->Printf("#Bridging Solute Residues:\n");
    // Place bridging values in a vector for sorting
    typedef std::vector<bridgeSorter> Bvec;
    Bvec bridgevector;
    bridgevector.reserve( BridgeMap_.size() );
    for (BmapType::const_iterator it = BridgeMap_.begin();
                                  it != BridgeMap_.end(); ++it)
      bridgevector.push_back( bridgeSorter(it->first, it->second) );
    std::sort( bridgevector.begin(), bridgevector.end() );
    for (Bvec::const_iterator bv = bridgevector.begin(); bv != bridgevector.end(); ++bv)
    {
      if (bridgeByAtom_) {
        bridgeout_->Printf("Bridge Atm");
        for (std::set<int>::const_iterator atm = bv->Uidx().begin();
                                           atm != bv->Uidx().end(); ++atm)
          bridgeout_->Printf(" %s", CurrentParm_->TruncAtomNameNum(*atm).c_str());
      } else {
        bridgeout_->Printf("Bridge Res");
        for (std::set<int>::const_iterator res = bv->Uidx().begin();
                                           res != bv->Uidx().end(); ++res)
          bridgeout_->Printf(" %i:%s", *res+1, CurrentParm_->Res( *res ).Name().Formatted(4).c_str());
      }
      bridgeout_->Printf(", %i frames.", bv->Binfo().Frames());
      if (!splitFrames_.empty()) {
        bridgeout_->Printf(" Parts:");
        for (unsigned int idx = 0; idx != bv->Binfo().Nparts(); idx++)
          bridgeout_->Printf(" %i", bv->Binfo().PartFrames(idx));
      }
      bridgeout_->Printf("\n");
    } 
  }
}

// ===== Action_HydrogenBond::Hbond ============================================
/** Calculate average distance and angle for hbond. */
void Action_HydrogenBond::Hbond::CalcAvg() {
  double dFrames = (double)frames_;
  dist_ /= dFrames;
  angle_ /= dFrames;
  angle_ *= Constants::RADDEG;
}

