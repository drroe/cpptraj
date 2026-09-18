#include "Analysis_MDANCE.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "ProgressBar.h"
#include "StringRoutines.h" // integerToString
#include "Trajout_Single.h"
#ifdef HAS_EIGEN
# include "Mdance/KMeans.h"
# include "Mdance/helm.h"
using namespace Cpptraj::Mdance;
#endif

/** CONSTRUCTOR */
Analysis_MDANCE::Analysis_MDANCE() :
  debug_(0),
  coords_(0),
  kClusters_(0),
  percentage_(0),
  vthresh_(0),
  metric_(ExtendedSimilarity::NO_METRIC),
  kinit_(MD::KinitType::StratAll),
  cnumvtime_(0),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ)
{}

/** DESTRUCTOR */
Analysis_MDANCE::~Analysis_MDANCE() {
}

const char* Analysis_MDANCE::kinitKeys_[] = {
  "all",
  "reduced",
  "compsim",
  "divselect",
  "kmeanspp",
  "random",
  "vanillakmpp",
  0
};

const char* Analysis_MDANCE::kinitStr_[] = {
  "All",
  "Reduced",
  "CompSim",
  "DivSelect",
  "Kmeans++",
  "Random",
  "Vanilla Kmeans++"
};

const Cpptraj::Mdance::MD::KinitType Analysis_MDANCE::kinitTypes_[] = {
  MD::KinitType::StratAll,
  MD::KinitType::StratReduced,
  MD::KinitType::CompSim,
  MD::KinitType::DivSelect,
  MD::KinitType::KmeansPP,
  MD::KinitType::Random,
  MD::KinitType::VanillaKmeansPP
};

// Analysis_MDANCE::Help()
void Analysis_MDANCE::Help() const {
# ifdef HAS_EIGEN
  mprintf("\tcrdset <COORDS set> clusters <#> [mask <mask>]\n"
          "\t[metric <metric>] [vthresh <vectthreshold>]\n"
          "\t[kinit <init>] [pct <percentage>]\n"
          "\t[name <set name>] [out <cnumvtime file>]\n"
          "\t[clusterout <trajfileprefix> [clusterfmt <trajformat>]]\n");
  mprintf("  <metric> = %s\n", ExtendedSimilarity::MetricKeys().c_str());
  mprintf("  <init>   =");
  for (int i = 0; kinitKeys_[i] != 0; i++)
    mprintf(" %s", kinitKeys_[i]);
  mprintf("\n");
# else
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
# endif
}


/** Get arguments related to writing cluster data to trajectories.
  * Copied from Cluster/Results_Coords.
  */
void Analysis_MDANCE::getClusterTrajArgs(ArgList& argIn,
                                         const char* trajKey, const char* fmtKey,
                                         std::string& trajName,
                                         TrajectoryFile::TrajFormatType& fmt) const
{
  trajName = argIn.GetStringKey( trajKey );
  fmt = TrajectoryFile::WriteFormatFromString( argIn.GetStringKey(fmtKey), fmt );
  // If file name specified but not format, try to guess from name
  if (!trajName.empty() && fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromFname( trajName, TrajectoryFile::AMBERTRAJ );
}


// Analysis_MDANCE::Setup()
Analysis::RetType Analysis_MDANCE::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
# ifdef HAS_EIGEN
  debug_ = debugIn;
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Target # of clusters, percentage, threshold
  kClusters_ = analyzeArgs.getKeyInt("clusters", 0);
  percentage_ = analyzeArgs.getKeyInt("pct", 10);
  vthresh_ = analyzeArgs.getKeyInt("vthresh", 16);
  // Metric
  std::string mstr = analyzeArgs.GetStringKey("metric");
  metric_ = ExtendedSimilarity::NO_METRIC;
  if (!mstr.empty()) {
    metric_ = ExtendedSimilarity::TypeFromKeyword( mstr );
    if (metric_ == ExtendedSimilarity::NO_METRIC) {
      mprinterr("Error: Metric '%s' not recognized.\n", mstr.c_str());
      return Analysis::ERR;
    }
  } else {
    metric_ = ExtendedSimilarity::MSD;
  }
  // Init strategy
  std::string kstr = analyzeArgs.GetStringKey("kinit");
  int iKinit = -1;
  if (!kstr.empty()) {
    for (int i = 0; kinitKeys_[i] != 0; i++) {
      const char* key = kinitKeys_[i];
      if (key != 0 && kstr == std::string(key)) {
        iKinit = i;
        break;
      }
    }
    if (iKinit < 0) {
      mprinterr("Error: Unrecognized keyword for 'kinit': %s\n", kstr.c_str());
      return Analysis::ERR;
    }
    kinit_ = kinitTypes_[iKinit];
  } else {
    kinit_ = MD::KinitType::StratAll;
    iKinit = 0;
  }
  // Atom mask
  std::string maskstr = analyzeArgs.GetStringKey("mask");
  if (maskstr.empty())
    maskstr = "*";
  if (mask_.SetMaskString(maskstr)) {
    mprinterr("Error: Could not set mask string '%s'\n", maskstr.c_str());
    return Analysis::ERR;
  }
  // Set up results that depend on COORDS DataSet
  getClusterTrajArgs(analyzeArgs, "clusterout",   "clusterfmt",   clusterfile_,  clusterfmt_);
  // Output files/data
  DataFile* cnumvtimefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Overall set name extracted here. All other arguments should already be processed. 
  std::string dsname = analyzeArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = setup.DSL().GenerateDefaultName("MDANCE");
  // ---------------------------------------------
    
  // Cluster number vs time data set
  cnumvtime_ = setup.DSL().AddSet(DataSet::INTEGER, dsname, "Cnum");
  if (cnumvtime_ == 0) return Analysis::ERR;
  if (cnumvtimefile != 0) cnumvtimefile->AddDataSet( cnumvtime_ );

  // Check input
  if (kClusters_ < 1) {
    mprinterr("Error: 'clusters' must be > 0.\n");
    return Analysis::ERR;
  }

  // Info
  mprintf("    MDANCE:\n");
  mprintf("\tCOORDS set       : %s\n", coords_->legend());
  mprintf("\t# clusters       : %i\n", kClusters_);
  mprintf("\tMetric           : %s\n", ExtendedSimilarity::metricStr(metric_));
  mprintf("\tInit. Strat.     : %s\n", kinitStr_[iKinit]);
  mprintf("\tPercentage       : %i%%\n", percentage_);
  mprintf("\tVect. threshhold : %i\n", vthresh_);
  mprintf("\tAtom selection   : %s\n", mask_.MaskString());
  mprintf("\tData set name          : %s\n", dsname.c_str());
  mprintf("\tCluster # vs time set  : %s\n", cnumvtime_->legend());
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time file : %s\n", cnumvtimefile->DataFilename().full());
  if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s.cX, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));

  return Analysis::OK;
# else /* HAS_EIGEN */
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
  return Analysis::ERR;
# endif /* HAS_EIGEN */
}

/** Write frames in each cluster to a trajectory file.  */
void Analysis_MDANCE::writeClusterTraj(ClusterArray const& Clusters) const {
  Topology* clusterparm = coords_->TopPtr();
  // Loop over all clusters
  for (unsigned int cidx = 0; cidx != Clusters.size(); cidx++)
  {
    Iarray const& cluster = Clusters[cidx];
    // Create filename based on cluster number.
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cidx );
    // Set up trajectory file 
    Trajout_Single clusterout;
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), clusterparm,
                                    coords_->CoordsInfo(), cluster.size(),
                                    clusterfmt_))
    {
      mprinterr("Error: Could not set up cluster trajectory %s for write.\n",
                cfilename.c_str());  
      return;
    } 
    // Loop over all frames in cluster
    unsigned int set = 0;
    Frame clusterframe = coords_->AllocateFrame();
    for (Iarray::const_iterator fnum = cluster.begin(); fnum != cluster.end(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterout.WriteSingle(set++, clusterframe);
    }
    // Close traj
    clusterout.EndTraj();
  }
}

// Analysis_MDANCE::Analyze()
Analysis::RetType Analysis_MDANCE::Analyze() {
# ifdef HAS_EIGEN
  mprintf("    MDANCE: Starting MDANCE.\n");
  if (coords_ == 0) {
    mprinterr("Error: COORDS are null.\n");
    return Analysis::ERR;
  }
  DataSet_Coords& CRD = static_cast<DataSet_Coords&>( *coords_ );
  // Get the atom selection
  if (CRD.Top().SetupIntegerMask( mask_ )) {
    mprinterr("Error: Could not set up atom mask.\n");
    return Analysis::ERR;
  }
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: No atoms selected.\n");
    return Analysis::ERR;
  }
  // This is an Eigen matrix. Each row is a frame, each column is a coordinate.
  unsigned int nSelectedAtoms = mask_.Nselected();
  unsigned int ncoords = nSelectedAtoms * 3;
  ArrayXXd data( CRD.Size(), ncoords );
  mprintf("\tSaving Eigen matrix (%zd rows/frames, %zd cols/coords).\n", data.rows(), data.cols());
  ProgressBar progress(CRD.Size());
  Frame frmIn = CRD.AllocateFrame();
  for (unsigned int idx = 0; idx != CRD.Size(); idx++)
  {
    progress.Update(idx);
    CRD.GetFrame(idx, frmIn);
    unsigned int icrd = 0;
    for (unsigned int iat = 0; iat < nSelectedAtoms; iat++)
    {
      const double* XYZ = frmIn.XYZ(mask_[iat]);
      data( idx, icrd   ) = XYZ[0];
      data( idx, icrd+1 ) = XYZ[1];
      data( idx, icrd+2 ) = XYZ[2];
      icrd += 3;
    }
  }

  // Convert metric to internal MDANCE
  MD::Metric mt = MD::Metric::MSD;
  switch(metric_) {
    case ExtendedSimilarity::MSD : mt = MD::Metric::MSD; break;
    case ExtendedSimilarity::BUB : mt = MD::Metric::BUB; break;
    case ExtendedSimilarity::FAI : mt = MD::Metric::Fai; break;
    case ExtendedSimilarity::GLE : mt = MD::Metric::Gle; break;
    case ExtendedSimilarity::JA  : mt = MD::Metric::Ja; break;
    case ExtendedSimilarity::JT  : mt = MD::Metric::JT; break;
    case ExtendedSimilarity::RT  : mt = MD::Metric::RT; break;
    case ExtendedSimilarity::RR  : mt = MD::Metric::RR; break;
    case ExtendedSimilarity::SM  : mt = MD::Metric::SM; break;
    case ExtendedSimilarity::SS1 : mt = MD::Metric::SS1; break;
    case ExtendedSimilarity::SS2 : mt = MD::Metric::SS2; break;
    case ExtendedSimilarity::NO_METRIC :
      mprinterr("Internal Error: Analysis_MDANCE::Analyze(): No metric.\n");
      return Analysis::ERR;
  }
  // Initialize Kmeans
  KmeansNANI kmeans(data, kClusters_, mt, kinit_, nSelectedAtoms, percentage_, vthresh_);
  // Results
  // First check the clustering assignments.
  // MDANCE labels each frame with the cluster number
  ClusterArray Clusters;
  Clusters.resize( kClusters_ );
  Veci cluster_of_frame = kmeans.getLabels();
  cnumvtime_->Allocate(DataSet::SizeArray(1, cluster_of_frame.size()));
  for (int ifrm = 0; ifrm < cluster_of_frame.size(); ifrm++) {
    int cnum = cluster_of_frame[ifrm];
    if (cnum < 0 || cnum >= kClusters_) {
      mprinterr("Error: Cluster of frame %i is out of bounds: %i\n", ifrm+1, cnum);
    }
    Clusters[cnum].push_back( ifrm );
    cnumvtime_->Add(ifrm, &cnum);
    mprintf("DEBUG: Frame %8i Cluster %8i\n", ifrm+1, cnum);
  }
  // Get the pseudo-F (Calinski-Harabasz) and DBI scores
  std::pair<double,double> scores = kmeans.computeScores();
  mprintf("\tDBI      : %f\n", scores.second);
  mprintf("\tpseudo-F : %f\n", scores.first);
  // Get centers
  Mat centers = kmeans.getCenters();
  mprintf("DEBUG: centers rows %zd, cols %zd\n", centers.rows(), centers.cols());
  // Write cluster trajectories
  if (!clusterfile_.empty())
    writeClusterTraj( Clusters );

  return Analysis::OK; // DEBUG
# else /* HAS_EIGEN */
  return Anlysis::ERR;
# endif /* HAS_EIGEN */
}
