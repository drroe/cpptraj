#include "Analysis_MDANCE.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"

/** CONSTRUCTOR */
Analysis_MDANCE::Analysis_MDANCE() :
  debug_(0),
  coords_(0)
{}

// Analysis_MDANCE::Help()
void Analysis_MDANCE::Help() const {
# ifdef HAS_EIGEN
  mprintf("\tcrdset <COORDS set>\n");
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
  using namespace Cpptraj::Mdance;
  if (coords_ == 0) {
    mprinterr("Error: COORDS are null.\n");
    return Analysis::ERR;
  }
  DataSet_Coords& CRD = static_cast<DataSet_Coords&>( *coords_ );
  // This is an Eigen matrix. Each row is a frame, each column is a coordinate.
  ArrayXXd data( CRD.Size(), CRD.Top().Natom()*3 );


  return Analysis::OK; // DEBUG
}
