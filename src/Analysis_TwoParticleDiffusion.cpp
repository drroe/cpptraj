#include "Analysis_TwoParticleDiffusion.h"
#include "CpptrajStdio.h"

// Analysis_TwoParticleDiffusion::Help()
void Analysis_TwoParticleDiffusion::Help() const {
  mprintf("\t[crdset <crd set>] [<name>] [<mask>] [out <filename>]\n"
          "\t[stop <maxlag>]\n");
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
  DataFile* df = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Get mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );
  // Set up output data set
  out_ = setup.DSL().AddSet(DataSet::GRID_DBL, analyzeArgs.GetStringNext(), "D_RR");
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
  return Analysis::OK;
}

// Analysis_TwoParticleDiffusion::Analyze()
Analysis::RetType Analysis_TwoParticleDiffusion::Analyze() {
  return Analysis::ERR;
}
