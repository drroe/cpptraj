#include <cmath> // log
#include <limits> // Minimum double val for checking zero
#include "Analysis_RemLog.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "Analysis_Lifetime.h"
#include "StringRoutines.h"    // integerToString
#include "DataSet_Mesh.h"      // slope Regression
#include "DataSet_MatrixFlt.h" // replica time matrix
#include "Trajout_Single.h"    // trajout
#include "TrajectoryFile.h"    // trajout
#include "ParmFile.h"          // trajout

Analysis_RemLog::Analysis_RemLog() :
  debug_(0),
  calculateStats_(false),
  calculateLifetimes_(false),
  printIndividualTrips_(false),
  calc_JSD_(true),
  masterDSL_(0),
  remlog_(0),
  repTimeMatrix_(0),
  mode_(NONE),
  lifetimes_(0),
  statsout_(0),
  reptime_(0),
  JSD_file_(0),
  calcRepFracSlope_(0),
  repFracSlope_(0)
{}

void Analysis_RemLog::Help() const {
  mprintf("\t{<remlog dataset> | <remlog filename>} [out <filename>] [crdidx | repidx]\n"
          "\t[stats [statsout <file>] [printtrips] [reptime <file>]] [lifetime <file>]\n"
          "\t[reptimeslope <n> reptimeslopeout <file>] [acceptout <file>] [name <setname>]\n"
          "\t[trajout <file> [parmout <file>] [jsd <file>]\n"
          "    crdidx: Print coordinate index vs exchange; output sets contain replica indices.\n"
          "    repidx: Print replica index vs exchange; output sets contain coordinate indices.\n"
          "  Analyze previously read in replica log data.\n");
}

// Analysis_RemLog::Setup()
Analysis::RetType Analysis_RemLog::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  // Get remlog dataset
  std::string remlogName = analyzeArgs.GetStringNext();
  if (remlogName.empty()) {
    mprinterr("Error: no remlog data set or file name specified.\n");
    return Analysis::ERR;
  }
  // Check if data set exists
  remlog_ = (DataSet_RemLog*)setup.DSL().FindSetOfType( remlogName, DataSet::REMLOG );
  if (remlog_ == 0) {
    mprinterr("Error: remlog data with name %s not found.\n", remlogName.c_str());
    return Analysis::ERR;
  }
  if (remlog_->Size() < 1 || remlog_->NumExchange() < 1) {
    mprinterr("Error: remlog data set appears to be empty.\n");
    return Analysis::ERR;
  }
  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("RL");
  acceptout_ = setup.DFL().AddCpptrajFile( analyzeArgs.GetStringKey("acceptout"), "replica acceptance",
                                           DataFileList::TEXT, true );
  if (acceptout_ == 0) return Analysis::ERR;
  calculateStats_ = analyzeArgs.hasKey("stats");
  JSD_file_ = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("jsd"), analyzeArgs );
  calc_JSD_ = (JSD_file_ != 0);
  if (calc_JSD_) calculateStats_ = true;
  trajoutName_.SetFileName( analyzeArgs.GetStringKey("trajout") );
  std::string pout = analyzeArgs.GetStringKey("parmout");
  if (!trajoutName_.empty() && pout.empty())
    parmoutName_.SetFileName_NoExpansion( trajoutName_.DirPrefix() + "/" +
                                          trajoutName_.Base() + ".parm7" );
  lifetimes_ = setup.DFL().AddCpptrajFile( analyzeArgs.GetStringKey("lifetime"), "remlog lifetimes" );
  calculateLifetimes_ = (lifetimes_ != 0);
  if (calculateStats_) {
    statsout_ = setup.DFL().AddCpptrajFile( analyzeArgs.GetStringKey("statsout"), "remlog stats",
                                       DataFileList::TEXT, true );
    reptime_ = setup.DFL().AddCpptrajFile( analyzeArgs.GetStringKey("reptime"), "replica times",
                                      DataFileList::TEXT, true );
    if (statsout_ == 0 || reptime_ == 0) return Analysis::ERR;
    repTimeMatrix_ = setup.DSL().AddSet( DataSet::MATRIX_FLT, MetaData(dsname_, "reptime") );
    if (repTimeMatrix_ == 0) return Analysis::ERR;
  }
  calcRepFracSlope_ = analyzeArgs.getKeyInt("reptimeslope", 0);
  std::string rfs_name = analyzeArgs.GetStringKey("reptimeslopeout");
  if (!calculateStats_) {
    calcRepFracSlope_ = 0;
    rfs_name.clear();
  }
  if ( (calcRepFracSlope_ > 0) != (!rfs_name.empty()) ) {
    mprinterr("Error: Both reptimeslope and reptimeslopeout must be specified.\n");
    return Analysis::ERR;
  }
  repFracSlope_ = setup.DFL().AddCpptrajFile( rfs_name, "replica fraction slope" );
  printIndividualTrips_ = analyzeArgs.hasKey("printtrips");
  // Get mode
  if (analyzeArgs.hasKey("crdidx"))
    mode_ = CRDIDX;
  else if (analyzeArgs.hasKey("repidx"))
    mode_ = REPIDX;
  else
    mode_ = NONE;
  std::string aspect;
  const char* yaxis = 0;
  if (mode_ == CRDIDX) {
    aspect = "repidx";
    yaxis = "ylabel CrdIdx";
  } else if (mode_ == REPIDX) {
    aspect = "crdidx";
    yaxis = "ylabel RepIdx";
  }
  // Set up an output set for each replica // FIXME should happen in Analyze()?
  DataFile* dfout = 0;
  if (mode_ != NONE) {
    // Get output filename
    std::string outname = analyzeArgs.GetStringKey("out");
    if (!outname.empty()) {
      dfout = setup.DFL().AddDataFile( outname, analyzeArgs );
      if (dfout == 0 ) return Analysis::ERR;
      if (yaxis != 0 ) dfout->ProcessArgs(yaxis);
    }
    MetaData md(dsname_, aspect);
    for (int i = 0; i < (int)remlog_->Size(); i++) {
      md.SetIdx(i+1);
      DataSet_integer* ds = (DataSet_integer*)setup.DSL().AddSet(DataSet::INTEGER, md);
      if (ds == 0) return Analysis::ERR;
      outputDsets_.push_back( (DataSet*)ds );
      if (dfout != 0) dfout->AddDataSet( (DataSet*)ds );
      ds->Resize( remlog_->NumExchange() ); 
    }
  }
  mprintf("   REMLOG: %s, %i replicas, %i exchanges\n", remlog_->legend(),
          remlog_->Size(), remlog_->NumExchange());
  if (mode_ == CRDIDX)
    mprintf("\tGetting coordinate index vs exchange.\n");
  else if (mode_ == REPIDX)
    mprintf("\tGetting replica index vs exchange.\n");
  if (mode_ != NONE && dfout != 0)
    mprintf("\tOutput is to %s\n", dfout->DataFilename().base());
  if (calculateStats_) {
    mprintf("\tGetting replica exchange stats, output to %s\n", statsout_->Filename().full());
    if (printIndividualTrips_)
      mprintf("\tIndividual round trips will be printed.\n");
    mprintf("\tWriting time spent at each replica to %s\n", reptime_->Filename().full());
  }
  if (calculateLifetimes_)
    mprintf("\tThe lifetime of each crd at each replica will be calculated.\n");
  if (acceptout_ != 0)
    mprintf("\tOverall exchange acceptance % will be written to %s\n",
            acceptout_->Filename().full());
  if (!trajoutName_.empty()) {
    mprintf("\tPseudo trajectory for crdidx will be written to '%s'\n", trajoutName_.full());
    mprintf("\tPseudo topology for crdidx will be written to '%s'\n", parmoutName_.full());
  }
  if (calc_JSD_) {
    mprintf("\tJensen-Shannon divergence of crd replica residence time from flat"
            " will be written to '%s'\n", JSD_file_->DataFilename().full());
  }
  masterDSL_ = setup.DslPtr();
  return Analysis::OK;
}

// Analysis_RemLog::Analyze()
Analysis::RetType Analysis_RemLog::Analyze() {
  if (remlog_->Size() < 1) {
    mprinterr("Error: No replicas in remlog data '%s'\n", remlog_->legend());
    return Analysis::ERR;
  }
  int Ndims = remlog_->DimTypes().Ndims();
  mprintf("\t'%s' %i replicas, %i exchanges, %i dims.\n", remlog_->legend(),
         remlog_->Size(), remlog_->NumExchange(), Ndims);
  // Set up output sets for JSD
  std::vector<DataSet*> JSD_sets;
  if (calc_JSD_) {
    MetaData md(dsname_, "JSD");
    for (int crdidx = 0; crdidx != (int)remlog_->Size(); crdidx++) {
      md.SetIdx( crdidx + 1 );
      DataSet* ds = masterDSL_->AddSet( DataSet::DOUBLE, md );
      if (ds == 0) return Analysis::ERR;
      JSD_file_->AddDataSet( ds );
      JSD_sets.push_back( ds );
    }
  }
  // Set up arrays for tracking replica stats.
  std::vector<RepStats> DimStats;
  std::vector<TripStats> DimTrips;
  for (int i = 0; i != Ndims; i++) {
    DimStats.push_back( RepStats(remlog_->Size()) );
    if (calculateStats_)
      DimTrips.push_back( TripStats(remlog_->Size()) );
  }
  std::vector<Iarray> replicaFrac;
  if (calculateStats_) {
    replicaFrac.resize( remlog_->Size() ); // [replica][crdidx]
    for (std::vector<Iarray>::iterator it = replicaFrac.begin();
                                       it != replicaFrac.end(); ++it)
      it->resize( remlog_->Size(), 0 );
  }
  // Variables for calculating replica lifetimes
  Analysis_Lifetime Lifetime;
  Array1D dsLifetime;
  std::vector< std::vector<DataSet_integer> > series; // 2D - repidx, crdidx
  if (calculateLifetimes_) {
    mprintf("\tData size used for lifetime analysis= %zu bytes.\n",
            remlog_->Size() * remlog_->Size() * remlog_->NumExchange() * sizeof(int));
    series.resize( remlog_->Size() );
    for (unsigned int i = 0; i < remlog_->Size(); i++) {
      series[i].resize( remlog_->Size() );
      for (unsigned int j = 0; j < remlog_->Size(); j++) {
        series[i][j].Resize( remlog_->NumExchange() );
        series[i][j].SetLegend("Rep"+integerToString(i+1)+",Crd"+integerToString(j+1));
        dsLifetime.push_back( (DataSet_1D*)&(series[i][j]) );
      }
    }
    if (Lifetime.ExternalSetup( dsLifetime, lifetimes_ ) == Analysis::ERR) {
      mprinterr("Error: Could not set up remlog lifetime analysis.\n");
      return Analysis::ERR;
    }
  }

  DataSet_Mesh mesh;
  if ( calcRepFracSlope_ > 0 ) {
    mesh.CalculateMeshX( remlog_->Size(), 1, remlog_->Size() );
    repFracSlope_->Printf("%-8s", "#Exchg");
    for (int crdidx = 0; crdidx < (int)remlog_->Size(); crdidx++)
      repFracSlope_->Printf("  C%07i_slope C%07i_corel", crdidx + 1, crdidx + 1);
    repFracSlope_->Printf("\n");
  }
  // Pseudo traj for tracking replica motion
  Topology RG_top;
  Trajout_Single repGroup;
  Frame RG_frame;
  if (!trajoutName_.empty()) {
    // Only good up to 3D
    if (Ndims > 3) {
      mprinterr("Error: crdidx pseudo trajectory generation only valid up to 3 dims.\n");
      return Analysis::ERR;
    }
    // Create an atom for each crdidx
    for (int atom = 0; atom != (int)remlog_->Size(); atom++)
      RG_top.AddTopAtom( Atom( "C" + integerToString(atom+1), "H "),
                         Residue("CRD", atom+1, ' ', ' ') );
    ParmFile RG_top_out;
    RG_top_out.WriteTopology( RG_top, parmoutName_, ArgList(), ParmFile::UNKNOWN_PARM, debug_ );
    if (repGroup.PrepareTrajWrite(trajoutName_, ArgList(), &RG_top, CoordinateInfo(),
                                 remlog_->NumExchange(), TrajectoryFile::UNKNOWN_TRAJ))
      return Analysis::ERR;
    RG_frame.SetupFrame( remlog_->Size() );
    RG_frame.ZeroCoords();
  }
  ProgressBar progress( remlog_->NumExchange() );
  for (int frame = 0; frame < remlog_->NumExchange(); frame++) {
    progress.Update( frame );
    for (int replica = 0; replica < (int)remlog_->Size(); replica++) {
      DataSet_RemLog::ReplicaFrame const& frm = remlog_->RepFrame( frame, replica );
      int crdidx = frm.CoordsIdx() - 1;
      int repidx = frm.ReplicaIdx() - 1;
      int dim = frm.Dim();
      // Exchange acceptance.
      // NOTE: Because currently the direction of the attempt is not always
      //       known unless the attempt succeeds for certain remlog types,
      //       the results will be skewed if dimension size is 2 since in that
      //       case the left partner is the right partner.
      if (replica == 0) DimStats[dim].attempts_++; // Assume same # attempts for every rep in dim
      if (frm.Success()) {
        if (frm.PartnerIdx() - 1 == remlog_->ReplicaInfo()[replica][dim].RightID())
          DimStats[dim].acceptUp_[replica]++;
        else // Assume down
          DimStats[dim].acceptDown_[replica]++;
      }
      if (mode_ == CRDIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[repidx]) );
        ds[frame] = frm.CoordsIdx();
      } else if (mode_ == REPIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[crdidx]) );
        ds[frame] = frm.ReplicaIdx();
      }
      // Pseudo traj for tracking replica motion
      if (!trajoutName_.empty()) {
        double* XYZ = RG_frame.xAddress() + (crdidx * 3);
        for (int didx = 0; didx != Ndims; didx++)
          XYZ[didx] = remlog_->ReplicaInfo()[replica][didx].GroupID();
      }
      if (calculateLifetimes_)
        series[repidx][crdidx][frame] = 1;
      if (calculateStats_) {
        TripStats& trip = static_cast<TripStats&>( DimTrips[dim] );
        // Fraction spent at each replica
        replicaFrac[repidx][crdidx]++;
        // Replica round-trip calculation
        if (trip.status_[crdidx] == UNKNOWN) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::BOTTOM) {
            trip.status_[crdidx] = HIT_BOTTOM;
            trip.bottom_[crdidx] = frame;
          }
        } else if (trip.status_[crdidx] == HIT_BOTTOM) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::TOP)
            trip.status_[crdidx] = HIT_TOP;
        } else if (trip.status_[crdidx] == HIT_TOP) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::BOTTOM) {
            int rtrip = frame - trip.bottom_[crdidx];
            if (printIndividualTrips_) {
              if (Ndims < 2)
                statsout_->Printf("[%i] CRDIDX %i took %i exchanges to travel"
                                 " up and down (exch %i to %i)\n",
                                 replica, crdidx+1, rtrip, trip.bottom_[crdidx]+1, frame+1);
              else
                statsout_->Printf("[%i] CRDIDX %i took %i exchanges to travel"
                                 " up and down in dim %i (exch %i to %i)\n",
                                 replica, crdidx+1, rtrip, dim+1, trip.bottom_[crdidx]+1, frame+1);
            }
            trip.roundTrip_[crdidx].AddElement( rtrip );
            trip.status_[crdidx] = HIT_BOTTOM;
            trip.bottom_[crdidx] = frame;
          }
        }
      }
    } // END loop over replicas
    // Pseudo traj for tracking replica motion
    if (!trajoutName_.empty())
      repGroup.WriteSingle( frame, RG_frame );
    // Replica "slope" calculation.
    if (calcRepFracSlope_ > 0 && frame > 0 && (frame % calcRepFracSlope_) == 0) {
      repFracSlope_->Printf("%8i", frame+1);
      for (int crdidx = 0; crdidx < (int)remlog_->Size(); crdidx++) {
        for (int replica = 0; replica < (int)remlog_->Size(); replica++)
          mesh.SetY(replica, (double)replicaFrac[replica][crdidx] / (double)frame);
        double slope, intercept, correl;
        mesh.LinearRegression(slope, intercept, correl, 0);
        repFracSlope_->Printf("  %14.7g %14.7g", slope * 100.0, correl);
                //frame+1, crdidx, slope * 100.0, intercept * 100.0, correl
      }
      repFracSlope_->Printf("\n");
    }
    // DEBUG
    // Jensen-Shannon divergence of each crdidx probability from "ideal"
    if (calc_JSD_) {
      double dFrame = (double)frame + 1.0; // Already processed frame.
      double normQ = 1.0 / (double)remlog_->Size(); // ideal dist, already normalized to 1.0
      //std::vector<double> normP( remlog_->Size(), 0.0 );
      //std::vector<double> normM( remlog_->Size(), 0.0 );
      for (int crdidx = 0; crdidx != (int)remlog_->Size(); crdidx++) {
        //double sumQ = 0.0; // DEBUG
        //double sumP = 0.0; // DEBUG
        //double sumM = 0.0; // DEBUG
        double div_PM = 0.0;
        double div_QM = 0.0;
        bool PM_valid = true;
        bool QM_valid = true;
        //mprintf("\nFrame %i CRDIDX %i\n", frame+1, crdidx+1);
        for (int replica = 0; replica != (int)remlog_->Size(); replica++) {
          double normP = (double)replicaFrac[replica][crdidx] / dFrame;
          double normM = 0.5 * (normP + normQ);
          //mprintf("P= %g  Q= %g  M= %g\n", normP, normQ, normM);
 
          bool Pzero = (normP <= std::numeric_limits<double>::min());
          bool Qzero = (normQ <= std::numeric_limits<double>::min());
          bool Mzero = (normM <= std::numeric_limits<double>::min());

          if (PM_valid) {
            if (!Pzero && !Mzero)
              div_PM += (log( normP / normM) * normP);
            else if (Pzero != Mzero) {
              PM_valid = false;
              div_PM = 0.0;
            }
          }

          if (QM_valid) {
            if (!Qzero && !Mzero)
              div_QM += (log( normQ / normM) * normQ);
            else if (Qzero != Mzero) {
              QM_valid = false;
              div_QM = 0.0;
            }
          }
          // DEBUG
          //sumQ += normQ;
          //sumP += normP[replica];
          //sumM += normM[replica];
        }
        double divergence = (0.5 * div_PM) + (0.5 * div_QM);
        //mprintf("PM= %g  QM= %g  JSD= %g\n", div_PM, div_QM, divergence);
        JSD_sets[crdidx]->Add(frame, &divergence);
        //mprintf("DEBUG: crdidx %i: JS_div = %g\n", crdidx, divergence);
        //mprintf("DEBUG: crdidx %i: sumQ= %g  sumP= %g  sumM= %g\n", crdidx, sumQ, sumP, sumM);
      }
    }
  } // END loop over exchanges
  if (!trajoutName_.empty())
    repGroup.EndTraj();
  // Exchange acceptance calc.
  for (int dim = 0; dim != Ndims; dim++) {
    // Assume number of exchange attempts is actually /2 since in Amber
    // attempts alternate up/down.
    if (Ndims > 1)
      acceptout_->Printf("#DIMENSION %i\n", dim+1);
    if (debug_ > 0) {
      for (int replica = 0; replica != (int)remlog_->Size(); replica++)
        mprintf("Rep %i attempts %i up %i down %i\n", replica, DimStats[dim].attempts_,
                DimStats[dim].acceptUp_[replica], DimStats[dim].acceptDown_[replica]);
    }
    acceptout_->Printf("%-8s %8s %8s\n", "#Replica", "%UP", "%DOWN");
    double exchangeAttempts = (double)DimStats[dim].attempts_ / 2.0;
    // Write according to group order
    for (DataSet_RemLog::GroupDimType::const_iterator Group = remlog_->GroupDims()[dim].begin();
                                                      Group != remlog_->GroupDims()[dim].end();
                                                    ++Group)
    {
      for (DataSet_RemLog::GroupArray::const_iterator Rep = Group->begin();
                                                      Rep != Group->end(); ++Rep)
        acceptout_->Printf("%8i %8.3f %8.3f\n", Rep->Me(),
          ((double)DimStats[dim].acceptUp_[Rep->Me()-1] / exchangeAttempts) * 100.0,
          ((double)DimStats[dim].acceptDown_[Rep->Me()-1] / exchangeAttempts) * 100.0);
    }
  }
  if (calculateStats_) {
    statsout_->Printf("# %i replicas, %i exchanges.\n", remlog_->Size(), remlog_->NumExchange());
    for (int dim = 0; dim != Ndims; dim++) {
      if (Ndims > 1)
        statsout_->Printf("#Dim%i Round-trip stats:\n", dim+1);
      else
        statsout_->Printf("#Round-trip stats:\n");
      statsout_->Printf("#%-8s %8s %12s %12s %12s %12s\n", "CRDIDX", "RndTrips", 
                       "AvgExch.", "SD_Exch.", "Min", "Max");
      unsigned int idx = 1;
      for (DSI_array::const_iterator rt = DimTrips[dim].roundTrip_.begin();
                                     rt != DimTrips[dim].roundTrip_.end(); ++rt)
      {
        double stdev = 0.0;
        double avg = rt->Avg( stdev );
        statsout_->Printf("%-8u %8i %12.4f %12.4f %12.0f %12.0f\n", 
                          idx++, rt->Size(), avg, stdev, rt->Min(), rt->Max());
      }
    }
    // Time spent at each replica
    DataSet_MatrixFlt& RTM = static_cast<DataSet_MatrixFlt&>( *repTimeMatrix_ );
    RTM.Allocate2D( remlog_->Size(), remlog_->Size() );
    reptime_->Printf("#Percent time spent at each replica:\n%-8s", "#Replica");
    for (int crd = 0; crd < (int)remlog_->Size(); crd++)
      reptime_->Printf(" CRD_%04i", crd + 1);
    reptime_->Printf("\n");
    double dframes = (double)remlog_->NumExchange();
    for (int replica = 0; replica < (int)remlog_->Size(); replica++) {
      reptime_->Printf("%8i", replica+1);
      for (int crd = 0; crd < (int)remlog_->Size(); crd++) {
        double frac = ((double)replicaFrac[replica][crd] / dframes) * 100.0;
        reptime_->Printf(" %8.3f", frac);
        RTM.AddElement((float)frac);
      }
      reptime_->Printf("\n");
    }
  }
  if (calculateLifetimes_) {
    mprintf("\tCalculating remlog lifetimes:\n");
    Lifetime.Analyze();
  }
  return Analysis::OK;
}
