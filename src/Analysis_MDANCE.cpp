#include "Analysis_MDANCE.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
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
  metric_(ExtendedSimilarity::NO_METRIC)
{}

// Analysis_MDANCE::Help()
void Analysis_MDANCE::Help() const {
# ifdef HAS_EIGEN
  mprintf("\tcrdset <COORDS set> clusters <#>\n"
          "\t[metric <metric>] [out <file>]\n");
  mprintf("  <metric> = %s\n", ExtendedSimilarity::MetricKeys().c_str());
# else
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
# endif
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

const Cpptraj::Mdance::MD::KinitType Analysis_MDANCE::kinitTypes_[] = {
  MD::KinitType::StratAll,
  MD::KinitType::StratReduced,
  MD::KinitType::CompSim,
  MD::KinitType::DivSelect,
  MD::KinitType::KmeansPP,
  MD::KinitType::Random,
  MD::KinitType::VanillaKmeansPP
};

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
  // Target # of clusters
  kClusters_ = analyzeArgs.getKeyInt("clusters", 0);
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
  if (!kstr.empty()) {
    int iKinit = -1;
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
  } else
    kinit_ = MD::KinitType::StratAll;

  // Check input
  if (kClusters_ < 1) {
    mprinterr("Error: 'clusters' must be > 0.\n");
    return Analysis::ERR;
  }

  // Info
  mprintf("    MDANCE:\n");
  mprintf("\tCOORDS set: %s\n", coords_->legend());

  return Analysis::OK;
# else /* HAS_EIGEN */
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
  return Analysis::ERR;
# endif /* HAS_EIGEN */
}

// Analysis_MDANCE::Analyze()
Analysis::RetType Analysis_MDANCE::Analyze() {
# ifdef HAS_EIGEN
  if (coords_ == 0) {
    mprinterr("Error: COORDS are null.\n");
    return Analysis::ERR;
  }
  DataSet_Coords& CRD = static_cast<DataSet_Coords&>( *coords_ );
  // This is an Eigen matrix. Each row is a frame, each column is a coordinate.
  ArrayXXd data( CRD.Size(), CRD.Top().Natom()*3 );


  return Analysis::OK; // DEBUG
# else /* HAS_EIGEN */
  return Anlysis::ERR;
# endif /* HAS_EIGEN */
}
