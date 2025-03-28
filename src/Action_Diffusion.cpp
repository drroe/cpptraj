#include <cmath> // sqrt, round
#include "Action_Diffusion.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // validDouble
#include "DataSet_1D.h" // LinearRegression
#include "DataSet_Mat3x3.h"
#ifdef TIMER
# include "Timer.h"
#endif
#ifdef _OPENMP
# include <omp.h>
#endif
#ifdef MPI
#include "DataSet_double.h"
#endif

// CONSTRUCTOR
Action_Diffusion::Action_Diffusion() :
  avg_x_(0), avg_y_(0), avg_z_(0), avg_r_(0), avg_a_(0),
  printIndividual_(false),
  calcDiffConst_(false),
  time_(1.0),
  debug_(0),
  outputx_(0), outputy_(0), outputz_(0), outputr_(0), outputa_(0),
  masterDSL_(0),
  avgucell_(0),
  allowMultipleTimeOrigins_(false)
# ifdef MPI
  ,multipleTimeOrigins_(false)
# endif
{}

static inline void ShortHelp() {
  mprintf("\t[{out <filename> | separateout <suffix>}] [time <time per frame>] [noimage]\n"
          "\t[<mask>] [<set name>] [individual] [diffout <filename>] [nocalc]\n"
          "\t[avgucell <avg ucell set>]\n");
#ifdef MPI
  mprintf("\t[allowmultipleorigins]\n");
#endif
}

void Action_Diffusion::Help() const {
  ShortHelp();
  mprintf("  Compute a mean square displacement plot for the atoms in the <mask>.\n"
          "  By default the average displacements are calculated unless 'individual'\n"
          "  is specified. Diffusion constants will be calculated from best-fit linear\n"
          "  regression lines of MSDs vs time unless 'nocalc' is specified.\n");
}

static inline int CheckTimeArg(double dt) {
  if (dt <= 0.0) {
    mprinterr("Error: Diffusion time per frame incorrectly specified, must be > 0.0.\n");
    return 1;
  }
  return 0;
}

// Action_Diffusion::Init()
Action::RetType Action_Diffusion::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
  multipleTimeOrigins_ = false;
# endif
  debug_ = debugIn;
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")) );
  // Determine if this is old syntax or new.
  if (actionArgs.Nargs() > 2 && actionArgs.ArgIsMask(1) && validDouble(actionArgs[2]))
  {
    // Old syntax: <mask> <time per frame> [average] [<prefix>]
    printIndividual_ = !(actionArgs.hasKey("average"));
    calcDiffConst_ = false;
    mprintf("Warning: Deprecated syntax for 'diffusion'. Consider using new syntax:\n");
    ShortHelp();
    mask_.SetMaskString( actionArgs.GetMaskNext() );
    time_ = actionArgs.getNextDouble(1.0);
    if (CheckTimeArg(time_)) return Action::ERR;
    std::string outputNameRoot = actionArgs.GetStringNext();
    // Default filename: 'diffusion'
    if (outputNameRoot.empty())
      outputNameRoot.assign("diffusion");
    // Open output files
    ArgList oldArgs("prec 8.3 noheader");
    DataFile::DataFormatType dft = DataFile::DATAFILE;
    outputx_ = init.DFL().AddDataFile(outputNameRoot+"_x.xmgr", dft, oldArgs);
    outputy_ = init.DFL().AddDataFile(outputNameRoot+"_y.xmgr", dft, oldArgs);
    outputz_ = init.DFL().AddDataFile(outputNameRoot+"_z.xmgr", dft, oldArgs);
    outputr_ = init.DFL().AddDataFile(outputNameRoot+"_r.xmgr", dft, oldArgs);
    outputa_ = init.DFL().AddDataFile(outputNameRoot+"_a.xmgr", dft, oldArgs);
  } else {
    // New syntax: [{separateout <suffix> | out <filename>}] [time <time per frame>]
    //             [<mask>] [<set name>] [individual] [diffout <filename>] [nocalc]
    printIndividual_ = actionArgs.hasKey("individual");
    calcDiffConst_ = !(actionArgs.hasKey("nocalc"));
    std::string suffix = actionArgs.GetStringKey("separateout");
    std::string outname = actionArgs.GetStringKey("out");
    if (!outname.empty() && !suffix.empty()) {
      mprinterr("Error: Specify either 'out' or 'separateout', not both.\n");
      return Action::ERR;
    }
    results_.AddDiffOut(init.DFL(), actionArgs.GetStringKey("diffout"));
    time_ = actionArgs.getKeyDouble("time", 1.0);
    if (CheckTimeArg(time_)) return Action::ERR;
    mask_.SetMaskString( actionArgs.GetMaskNext() );
    // Open output files.
    if (!suffix.empty()) {
      FileName FName( suffix );
      outputx_ = init.DFL().AddDataFile(FName.PrependFileName("x_"), actionArgs);
      outputy_ = init.DFL().AddDataFile(FName.PrependFileName("y_"), actionArgs);
      outputz_ = init.DFL().AddDataFile(FName.PrependFileName("z_"), actionArgs);
      outputr_ = init.DFL().AddDataFile(FName.PrependFileName("r_"), actionArgs);
      outputa_ = init.DFL().AddDataFile(FName.PrependFileName("a_"), actionArgs);
      if (!results_.HasDiffOut() && calcDiffConst_)
        results_.AddDiffOut(init.DFL(), FName.PrependFileName("diff_").Full());
    } else if (!outname.empty()) {
      outputr_ = init.DFL().AddDataFile( outname, actionArgs );
      outputx_ = outputy_ = outputz_ = outputa_ = outputr_;
    }
  }
  if (results_.HasDiffOut()) calcDiffConst_ = true;
  allowMultipleTimeOrigins_ = actionArgs.hasKey("allowmultipleorigins");
  // Get average box
  std::string avgucellstr = actionArgs.GetStringKey("avgucell");
  if (!avgucellstr.empty()) {
    avgucell_ = init.DSL().GetDataSet( avgucellstr );
    if (avgucell_ == 0) {
      mprinterr("Error: No data set selected by '%s'\n", avgucellstr.c_str());
      return Action::ERR;
    }
    if (avgucell_->Type() != DataSet::MAT3X3) {
      mprinterr("Error: Average unit cell set '%s' is not a 3x3 matrix set.\n", avgucell_->legend());
      return Action::ERR;
    }
    if (avgucell_->Size() < 1) {
      mprinterr("Error: Average unit cell set '%s' is empty.\n", avgucell_->legend());
      return Action::ERR;
    }
    DataSet_Mat3x3 const& matset = static_cast<DataSet_Mat3x3 const&>( *avgucell_ );
#   ifdef MPI
    Matrix_3x3 ucell = matset[0];
    // Broadcast the master unit cell
    ucell.BroadcastMatrix( trajComm_ );
#   else
    Matrix_3x3 const& ucell = matset[0];
#   endif
    if (avgbox_.SetupFromUcell( ucell )) {
      mprinterr("Error: Could not set up box from unit cell parameters in '%s'\n", avgucell_->legend());
      return Action::ERR;
    }
  }
  // Add DataSets
  dsname_ = actionArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = init.DSL().GenerateDefaultName("Diff");
  avg_x_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "X"));
  avg_y_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Y"));
  avg_z_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Z"));
  avg_r_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "R"));
  avg_a_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "A"));
  if (avg_x_ == 0 || avg_y_ == 0 || avg_z_ == 0 || avg_r_ == 0 || avg_a_ == 0) {
    mprinterr("Error: Could not allocate one or more average sets.\n");
    return Action::ERR;
  }
  if (outputr_ != 0) outputr_->AddDataSet( avg_r_ );
  if (outputx_ != 0) outputx_->AddDataSet( avg_x_ );
  if (outputy_ != 0) outputy_->AddDataSet( avg_y_ );
  if (outputz_ != 0) outputz_->AddDataSet( avg_z_ );
  if (outputa_ != 0) outputa_->AddDataSet( avg_a_ );
  // Set X dim
  Xdim_ = Dimension(0.0, time_, "Time");
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);
  // Add DataSets for diffusion constant calc
  if (calcDiffConst_) {
    if (results_.CreateDiffusionSets(init.DSL(), dsname_)) return Action::ERR;
  }
  // Save master data set list, needed when printIndividual_
  masterDSL_ = init.DslPtr();

  mprintf("    DIFFUSION:\n");
  mprintf("\tAtom Mask is [%s]\n", mask_.MaskString());
  if (printIndividual_)
    mprintf("\tBoth average and individual diffusion will be calculated.\n");
  else
    mprintf("\tOnly average diffusion will be calculated.\n");
  mprintf("\tData set base name: %s\n", avg_x_->Meta().Name().c_str());
  if (imageOpt_.UseImage())
    mprintf("\tCorrections for imaging enabled.\n");
  else
    mprintf("\tCorrections for imaging disabled.\n");
  // If one file defined, assume all are.
  if (outputx_ != 0) {
    mprintf("\tOutput files:\n"
            "\t  %s: (x) Mean square displacement(s) in the X direction (in Ang.^2).\n"
            "\t  %s: (y) Mean square displacement(s) in the Y direction (in Ang.^2).\n"
            "\t  %s: (z) Mean square displacement(s) in the Z direction (in Ang.^2).\n"
            "\t  %s: (r) Overall mean square displacement(s) (in Ang.^2).\n"
            "\t  %s: (a) Total distance travelled (in Ang.).\n",
            outputx_->DataFilename().full(), outputy_->DataFilename().full(),
            outputz_->DataFilename().full(), outputr_->DataFilename().full(),
            outputa_->DataFilename().full());
  }
  mprintf("\tThe time between frames is %g ps.\n", time_);
  if (calcDiffConst_) {
    results_.Info();
  } else
    mprintf("\tTo calculate diffusion constant from mean squared displacement plots,\n"
            "\t  calculate the slope of MSD vs time and multiply by 10.0/2*N (where N\n"
            "\t  is # of dimensions); this will give units of 1x10^-5 cm^2/s.\n");
  if (avgucell_ != 0) {
    mprintf("\tUsing average unit cell vectors from set '%s' to remove box fluctuations.\n", avgucell_->legend());
    avgbox_.PrintInfo();
  }
  if (allowMultipleTimeOrigins_) {
# ifdef MPI
    if (imageOpt_.UseImage()) {
      if (trajComm_.Size() > 1)
        mprintf("\tTrajectories that require unit cell imaging will be averaged over %i time origins.\n", trajComm_.Size());
    } else
      mprintf("\tImaging is disabled, ignoring 'allowmultipleorigins'\n");
# else
    mprintf("\tThe 'allowmultipleorigins' keyword is only relevant in parallel, ignoring.\n");
# endif /* MPI */
  }
# ifdef _OPENMP
# pragma omp parallel
  {
# pragma omp master
  {
  mprintf("\tParallelizing calculation with %i threads.\n", omp_get_num_threads());
  }
  }
# endif
  return Action::OK;
}

// Action_Diffusion::Setup()
Action::RetType Action_Diffusion::Setup(ActionSetup& setup) {
# ifdef TIMER
  Timer time_setup, time_addsets;
  time_setup.Start();
# endif
  // Setup atom mask
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }

  // Set up imaging info for this parm
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  if (imageOpt_.ImagingEnabled()) {
    mprintf("\tImaging enabled.\n");
#   ifdef MPI
    if (trajComm_.Size() > 1) {
      if (!allowMultipleTimeOrigins_) {
        mprinterr("Error: Imaging for 'diffusion' is not supported in parallel as there is\n"
                  "Error:   no way to correct for imaging that has taken place on preceding\n"
                  "Error:   MPI ranks. To use 'diffusion' in parallel, the trajectory should\n"
                  "Error:   be unwrapped. If this trajectory has already been unwrapped please\n"
                  "Error:   specify the 'noimage' keyword.\n");
        mprinterr("Error: To calculate diffusion in parallel by dividing the trajectory among\n"
                  "Error:   multiple time origins, specify the 'allowmultipleorigins' keyword.\n");
        return Action::ERR;
      } else {
        // Indicate data sets will have multiple time origins
        multipleTimeOrigins_ = true;
      }
    }
#   endif
  } else {
    mprintf("\tImaging disabled.\n");
    if (avgucell_ != 0)
      mprintf("Warning: 'avgucell' specified but trajectory has no unit cell.\n");
  }

  // If imaging, reserve space for the previous fractional coordinates array
  if (imageOpt_.ImagingEnabled())
    previousFrac_.reserve( mask_.Nselected() );

# ifdef MPI
  if (multipleTimeOrigins_ && initial_.empty()) {
    mprintf("Warning: Calculating diffusion in parallel with imaging turned on.\n"
            "Warning:   Mean-squared distance calculation will be averaged starting\n"
            "Warning:   from %i independent time origins.\n", trajComm_.Size());
  }
# endif

  // If initial frame already set and current # atoms > # atoms in initial
  // frame this will probably cause an error.
  if (!initial_.empty() && setup.Top().Natom() > initial_.Natom()) {
    mprintf("Warning: # atoms in current parm (%s, %i) > # atoms in initial frame (%i)\n",
             setup.Top().c_str(), setup.Top().Natom(), initial_.Natom());
    mprintf("Warning: This may lead to segmentation faults.\n");
  }

  // Set up sets for individual atoms if necessary
  if (printIndividual_) {
#   ifdef MPI
    if (multipleTimeOrigins_ && trajComm_.Size() > 1) {
      mprinterr("Error: Cannot track individual atom diffusion in parallel with imaging.\n");
      return Action::ERR;
    }
#   endif
    // Create as many spots for sets as needed. All do not have to be used.
    if (mask_.back() >= (int)atom_x_.size()) {
      int newSize = mask_.back() + 1;
      atom_x_.resize( newSize, 0 );
      atom_y_.resize( newSize, 0 );
      atom_z_.resize( newSize, 0 );
      atom_r_.resize( newSize, 0 );
      atom_a_.resize( newSize, 0 );
    }
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); at++)
    {
      if (atom_x_[*at] == 0) {
#       ifdef TIMER
        time_addsets.Start();
#       endif
        atom_x_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aX", *at+1));
        atom_y_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aY", *at+1));
        atom_z_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aZ", *at+1));
        atom_r_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aR", *at+1));
        atom_a_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aA", *at+1));
#       ifdef TIMER
        time_addsets.Stop();
#       endif
        if (outputx_ != 0) outputx_->AddDataSet(atom_x_[*at]);
        if (outputy_ != 0) outputy_->AddDataSet(atom_y_[*at]);
        if (outputz_ != 0) outputz_->AddDataSet(atom_z_[*at]);
        if (outputr_ != 0) outputr_->AddDataSet(atom_r_[*at]);
        if (outputa_ != 0) outputa_->AddDataSet(atom_a_[*at]);
        atom_x_[*at]->SetDim(Dimension::X, Xdim_);
        atom_y_[*at]->SetDim(Dimension::X, Xdim_);
        atom_z_[*at]->SetDim(Dimension::X, Xdim_);
        atom_r_[*at]->SetDim(Dimension::X, Xdim_);
        atom_a_[*at]->SetDim(Dimension::X, Xdim_);
      }
    }
  }
# ifdef TIMER
  time_setup.Stop();
  time_addsets.WriteTiming(3, "Diffusion Add Sets", time_setup.Total());
  time_setup.WriteTiming(2, "Diffusion Setup");
# endif
  return Action::OK;
}

// Action_Diffusion::DoAction()
Action::RetType Action_Diffusion::DoAction(int frameNum, ActionFrame& frm) {
  //rprintf("DEBUG: DIFFUSION FRAME %i\n", frameNum);
  Matrix_3x3 const* ucell = 0;
  if (imageOpt_.ImagingEnabled()) {
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
    if (avgucell_ == 0)
      ucell = &(frm.Frm().BoxCrd().UnitCell());
    else
      ucell = &(avgbox_.UnitCell());
  }

  // ----- Load initial frame if necessary -------
  if (initial_.empty()) {
    initial_ = frm.Frm();
#   ifdef MPI
    if (imageOpt_.ImagingEnabled()) {

    } else if (trajComm_.Size() > 1) {
      if (trajComm_.Master())
        for (int rank = 1; rank < trajComm_.Size(); ++rank)
          initial_.SendFrame( rank, trajComm_ );
      else
        initial_.RecvFrame( 0, trajComm_ );
    }
#   endif
    if (imageOpt_.ImagingEnabled()) {
      // Imaging. Store initial fractional coordinates.
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        previousFrac_.push_back( frm.Frm().BoxCrd().FracCell() * Vec3( initial_.XYZ(*atom) ) );
      }
      if (avgucell_ != 0) {
        // If using the average unit cell, we do not want it on the 
        // initial frame since that can shift coordinates and the
        // first frame should have a distance of 0.
        ucell = &(frm.Frm().BoxCrd().UnitCell());
      }
    }
  } // END initial frame load
  // ----- Diffusion calculation -----------------
  // For averaging over selected atoms
  double average2 = 0.0;
  double avgx = 0.0;
  double avgy = 0.0;
  double avgz = 0.0;
  int imask;
# ifdef _OPENMP
# pragma omp parallel private(imask) reduction(+ : average2, avgx, avgy, avgz)
  {
# pragma omp for
# endif
  for (imask = 0; imask < mask_.Nselected(); imask++)
  {
    int at = mask_[imask];
    // Get current coords for this atom.
    const double* XYZ = frm.Frm().XYZ(at);
    // Get initial coords for this atom
    const double* initialXYZ = initial_.XYZ(at);
    // Calculate distance from initial position. 
    double delx, dely, delz;
    if ( imageOpt_.ImagingEnabled() ) {
      // Convert current Cartesian coords to fractional coords
      Vec3 xyz_frac = frm.Frm().BoxCrd().FracCell() * Vec3(XYZ);
      // Correct current frac coords
      Vec3 ixyz = xyz_frac - previousFrac_[imask];
      ixyz[0] = ixyz[0] - round(ixyz[0]);
      ixyz[1] = ixyz[1] - round(ixyz[1]);
      ixyz[2] = ixyz[2] - round(ixyz[2]);
      xyz_frac = previousFrac_[imask] + ixyz;
      // Back to Cartesian
      Vec3 xyz_cart1 = ucell->TransposeMult( xyz_frac );
      // Update reference frac coords
      previousFrac_[imask] = xyz_frac;
      // Calculate the distance between the fixed coordinates
      // and reference (initial) frame coordinates.
      delx = xyz_cart1[0] - initialXYZ[0];
      dely = xyz_cart1[1] - initialXYZ[1];
      delz = xyz_cart1[2] - initialXYZ[2];
    } else {
      // No imaging. Calculate distance from current position to initial position.
      delx = XYZ[0] - initialXYZ[0];
      dely = XYZ[1] - initialXYZ[1];
      delz = XYZ[2] - initialXYZ[2];
    }
    // Calc distances for this atom
    double distx = delx * delx;
    double disty = dely * dely;
    double distz = delz * delz;
    double dist2 = distx + disty + distz;
    // Accumulate averages
    avgx += distx;
    avgy += disty;
    avgz += distz;
    average2 += dist2;
    // Store distances for this atom
    if (printIndividual_) {
      float fval = (float)distx;
      atom_x_[at]->Add(frameNum, &fval);
      fval = (float)disty;
      atom_y_[at]->Add(frameNum, &fval);
      fval = (float)distz;
      atom_z_[at]->Add(frameNum, &fval);
      fval = (float)dist2;
      atom_r_[at]->Add(frameNum, &fval);
      dist2 = sqrt(dist2);
      fval = (float)dist2;
      atom_a_[at]->Add(frameNum, &fval);
    }
  } // END loop over selected atoms
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  // Calc averages
  double dNselected = 1.0 / (double)mask_.Nselected();
  avgx *= dNselected;
  avgy *= dNselected;
  avgz *= dNselected;
  average2 *= dNselected;
  // Save averages
  avg_x_->Add(frameNum, &avgx);
  avg_y_->Add(frameNum, &avgy);
  avg_z_->Add(frameNum, &avgz);
  avg_r_->Add(frameNum, &average2);
  average2 = sqrt(average2);
  avg_a_->Add(frameNum, &average2);
  return Action::OK;
}

#ifdef MPI
void Action_Diffusion::average_multiple_time_origins(DataSet* set,
                                                     std::vector<int> const& nvals_on_rank,
                                                     int max)
const
{
  // Sanity check - underlying set should be double
  if (set->Type() != DataSet::DOUBLE) {
    mprinterr("Internal Error: Action_Diffusion::average_multiple_time_origins(): Set %s is not double.\n", set->legend());
    return;
  }
  DataSet_double const& dset = static_cast<DataSet_double const&>( *set );
  // Average over all the values
  if (trajComm_.Master()) {
    std::vector<double> Sum_all( max, 0.0 );
    std::vector<double> Ntotal( max, 0.0 );
    std::vector<double> fromChild( max, 0.0 );
    // Sum master rank
    for (unsigned int idx = 0; idx < set->Size(); idx++) {
      Sum_all[idx] = dset[idx];
      Ntotal[idx] = 1.0;
    }
    for (int rank = 1; rank < trajComm_.Size(); rank++) {
      // Receive from child
      trajComm_.Recv( &fromChild[0], nvals_on_rank[rank], MPI_DOUBLE, rank, 3000 );
      // Sum child rank
      for (int idx = 0; idx < nvals_on_rank[rank]; idx++) {
        Sum_all[idx] += fromChild[idx];
        Ntotal[idx] += 1.0;
      }
    }
    //for (int idx = 0; idx < max; idx++) 
    //  mprintf("DEBUG sum=%g total=%g avg=%g\n", Sum_all[idx], Ntotal[idx], Sum_all[idx] / Ntotal[idx]);
    // Ensure master set is big enough to hold result
    DataSet_double& updateSet = static_cast<DataSet_double&>( *set );
    if (dset.Size() < (unsigned int)max)
      updateSet.Resize( max );
    // Update master set
    for (int idx = 0; idx < max; idx++)
      updateSet[idx] = Sum_all[idx] / Ntotal[idx];
  } else {
    // Send the values to master
    trajComm_.Send(dset.DvalPtr(), dset.Size(), MPI_DOUBLE, 0, 3000);
  }
  // Set no longer needs to be synced
  set->SetNeedsSync( false );
}

/** See if we need to average over multiple time origins. */
int Action_Diffusion::SyncAction() {
  if (multipleTimeOrigins_) {
    mprintf("    DIFFUSION: Calculating diffusion by averaging over %i time origins.\n", trajComm_.Size());
    int max = 0;
    std::vector<int> nvals_on_rank;
    if (trajComm_.Master()) {
      // Get how many values each rank has
      nvals_on_rank.assign( trajComm_.Size(), 0 );
      int nvals = (int)avg_a_->Size();
      trajComm_.GatherMaster( &nvals, 1, MPI_INT, &nvals_on_rank[0] );
      //mprintf("DEBUG: %s array sizes:", avg_a_->legend());
      //for (std::vector<int>::const_iterator it = nvals_on_rank.begin(); it != nvals_on_rank.end(); ++it)
      //  mprintf(" %i", *it);
      //mprintf("\n");
      max = nvals;
      for (unsigned int idx = 1; idx < nvals_on_rank.size(); idx++)
        if (nvals_on_rank[idx] > max) max = nvals_on_rank[idx];
      mprintf("\tMax time is %g ps.\n", (double)(max-1) * time_);

    } else {
      // Tell master how many values we have
      int nvals = (int)avg_a_->Size();
      trajComm_.GatherMaster( &nvals, 1, MPI_INT, 0);
    }
    average_multiple_time_origins( avg_x_, nvals_on_rank, max );
    average_multiple_time_origins( avg_y_, nvals_on_rank, max );
    average_multiple_time_origins( avg_z_, nvals_on_rank, max );
    average_multiple_time_origins( avg_r_, nvals_on_rank, max );
    average_multiple_time_origins( avg_a_, nvals_on_rank, max );
  }
  return 0;
}
#endif /* MPI */

// Action_Diffusion::Print()
void Action_Diffusion::Print() {
  if (!calcDiffConst_) return;
  mprintf("    DIFFUSION: Calculating diffusion constants from slopes.\n");
  std::string const& name = avg_r_->Meta().Name();
  unsigned int set = 0;
  results_.CalcDiffusionConst( set, avg_r_, 3, name + "_AvgDr" );
  results_.CalcDiffusionConst( set, avg_x_, 1, name + "_AvgDx" );
  results_.CalcDiffusionConst( set, avg_y_, 1, name + "_AvgDy" );
  results_.CalcDiffusionConst( set, avg_z_, 1, name + "_AvgDz" );
  if (printIndividual_) {
    CalcDiffForSet( set, atom_r_, 3, name + "_dr" );
    CalcDiffForSet( set, atom_x_, 3, name + "_dx" );
    CalcDiffForSet( set, atom_y_, 3, name + "_dy" );
    CalcDiffForSet( set, atom_z_, 3, name + "_dz" );
  }
}

// Action_Diffusion::CalcDiffForSet()
void Action_Diffusion::CalcDiffForSet(unsigned int& set, Dlist const& Sets, int Ndim,
                                      std::string const& label) const
{
  for (Dlist::const_iterator ds = Sets.begin(); ds != Sets.end(); ds++)
    if (*ds != 0)
      results_.CalcDiffusionConst(set, *ds, Ndim, label + "_" + integerToString( (*ds)->Meta().Idx() ));
}
