#include <cmath> // ceil
#include "Analysis_TwoParticleDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"

// Analysis_TwoParticleDiffusion::Help()
void Analysis_TwoParticleDiffusion::Help() const {
  mprintf("\t[crdset <crd set>] [<name>] [<mask>] [out <filename>]\n"
          "\t[stop <maxlag>] rmax <max> rstep <step>\n");
}

// Analysis_TwoParticleDiffusion::Setup()
Analysis::RetType Analysis_TwoParticleDiffusion::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: rmsavgcorr: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  // Get Keywords
  maxlag_ = analyzeArgs.getKeyInt("stop", -1);
  rmax_ = analyzeArgs.getKeyDouble("rmax", 10.0);
  rstep_ = analyzeArgs.getKeyDouble("rstep", 1.0);
  // Some sanity checking
  if (rmax_ <= 0.0) {
    mprinterr("Error: 'rmax' must be > 0.0\n");
    return Analysis::ERR;
  }
  if (rstep_ <= 0.0) {
    mprinterr("Error: 'rstep' must be > 0.0\n");
    return Analysis::ERR;
  }
  DataFile* df = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Get mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );
  // Set up output data set
  out_ = setup.DSL().AddSet(DataSet::MATRIX_DBL, analyzeArgs.GetStringNext(), "D_RR");
  if (out_ == 0) return Analysis::ERR;
  if (df != 0) df->AddDataSet( out_ );

  mprintf("    TWOPARTICLEDIFFUSION: COORDS set '%s', mask [%s]\n",
          coords_->legend(), mask_.MaskString());
  mprintf("\tOutput set: %s\n", out_->legend());
  if (df != 0) mprintf("\tOutput set written to: %s\n", df->DataFilename().full());
  if (maxlag_ < 1)
    mprintf("\tMax lag will be half the number of frames in input COORDS set.\n");
  else
    mprintf("\tMax lag is %i\n", maxlag_);
  mprintf("\tRmax is %g, Rstep is %g\n", rmax_, rstep_);
  return Analysis::OK;
}

// Analysis_TwoParticleDiffusion::Analyze()
Analysis::RetType Analysis_TwoParticleDiffusion::Analyze() {
  // Check that there is data
  if (coords_->Size() < 1) {
    mprinterr("Error: No data in COORDS set '%s'\n", coords_->legend());
    return Analysis::ERR;
  }
  // Check number of frames
  if (maxlag_ < 1)
    maxlag_ = (int)coords_->Size() / 2;
  if (maxlag_ >= (int)coords_->Size()) {
    maxlag_ = (int)coords_->Size() - 1;
    mprintf("Warning: Max lag is > # frames, setting to %i\n", maxlag_);
  }
  mprintf("\tMax lag: %i\n", maxlag_);
  // Set up atom mask
  if (coords_->Top().SetupIntegerMask( mask_ )) return Analysis::ERR;
  mask_.MaskInfo();
  if (mask_.None()) return Analysis::ERR;
  // Set up input frames
  Frame frame0;
  frame0.SetupFrameFromMask( mask_, coords_->Top().Atoms() );
  Frame frame1 = frame0;
  // Determine size of the 'R' dimension
  int numRbins = (int)ceil(rmax_ / rstep_);
  mprintf("\tNumber of bins in the 'R' dimension: %i\n", numRbins);
  // Allocate matrix
  DataSet_MatrixDbl& mat = static_cast<DataSet_MatrixDbl&>( *out_ );
  
  return Analysis::OK;
}
