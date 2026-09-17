#include "Analysis_MDANCE.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "ProgressBar.h"
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
  cnumvtime_(0)
{}

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
  mprintf("\tcrdset <COORDS set> clusters <#>\n"
          "\t[metric <metric>] [out <file>] [vthresh <vectthreshold>]\n"
          "\t[kinit <init>] [pct <percentage>]\n");
  mprintf("  <metric> = %s\n", ExtendedSimilarity::MetricKeys().c_str());
  mprintf("  <init>   =");
  for (int i = 0; kinitKeys_[i] != 0; i++)
    mprintf(" %s", kinitKeys_[i]);
  mprintf("\n");
# else
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
# endif
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
  // Output files/data
  DataFile* cnumvtimefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Overall set name extracted here. All other arguments should already be processed. 
  std::string dsname = analyzeArgs.GetStringNext();
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
  mprintf("\tData set name          : %s\n", dsname.c_str());
  mprintf("\tCluster # vs time set  : %s\n", cnumvtime_->legend());
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time file : %s\n", cnumvtimefile->DataFilename().full());

  return Analysis::OK;
# else /* HAS_EIGEN */
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
  return Analysis::ERR;
# endif /* HAS_EIGEN */
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
  // This is an Eigen matrix. Each row is a frame, each column is a coordinate.
  ArrayXXd data( CRD.Size(), CRD.Top().Natom()*3 );
  mprintf("\tSaving Eigen matrix (%zd rows/frames, %zd cols/coords).\n", data.rows(), data.cols());
  ProgressBar progress(CRD.Size());
  Frame frmIn = CRD.AllocateFrame();
  for (unsigned int idx = 0; idx != CRD.Size(); idx++)
  {
    progress.Update(idx);
    CRD.GetFrame(idx, frmIn);
    unsigned int icrd = 0;
    for (int iat = 0; iat < CRD.Top().Natom(); iat++)
    {
      const double* XYZ = frmIn.XYZ(iat);
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
  KmeansNANI kmeans(data, kClusters_, mt, kinit_, CRD.Top().Natom(), percentage_, vthresh_);
  // Results
  // First check the clustering assignments.
  // MDANCE labels each frame with the cluster number
  Veci cluster_of_frame = kmeans.getLabels();
  cnumvtime_->Allocate(DataSet::SizeArray(1, cluster_of_frame.size()));
  for (int i = 0; i < cluster_of_frame.size(); i++) {
    int cnum = cluster_of_frame[i];
    if (cnum < 0 || cnum >= kClusters_) {
      mprinterr("Error: Cluster of frame %i is out of bounds: %i\n", i+1, cnum);
    }
    cnumvtime_->Add(i, &cnum);
    mprintf("DEBUG: Frame %8i Cluster %8i\n", i+1, cnum);
  }
  // Get the pseudo-F (Calinski-Harabasz) and DBI scores
  std::pair<double,double> scores = kmeans.computeScores();
  mprintf("\tDBI      : %f\n", scores.second);
  mprintf("\tpseudo-F : %f\n", scores.first);
  // Get centers
  Mat centers = kmeans.getCenters();
  mprintf("DEBUG: centers rows %zd, cols %zd\n", centers.rows(), centers.cols());

  return Analysis::OK; // DEBUG
# else /* HAS_EIGEN */
  return Anlysis::ERR;
# endif /* HAS_EIGEN */
}
