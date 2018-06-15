#include <cmath> // ceil
#include "Analysis_TwoParticleDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"
#include "DistRoutines.h"

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
  if (maxlag_ != -1 && maxlag_ < 2) {
    mprinterr("Error: maxlag must be > 1\n");
    return Analysis::ERR;
  }
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
  // FIXME: Currently assume unwrapped coords!
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
  // Allocate matrix: cols=lag, rows=R
  DataSet_MatrixDbl& mat = static_cast<DataSet_MatrixDbl&>( *out_ );
  mat.Allocate2D( maxlag_-1, numRbins );
  // Store atom pair Rbin indices via calculating initial pair distances.
  Matrix<int> PairBins;
  // col == 0 means Triangle matrix
  PairBins.resize( 0, mask_.Nselected() );
  coords_->GetFrame( 0, frame0, mask_ );
  double cut2 = rmax_ * rmax_;
  double one_over_spacing = 1.0 / rstep_;
  for (int at0 = 0; at0 < frame0.Natom(); at0++)
    for (int at1 = at0+1; at1 < frame0.Natom(); at1++)
    {
      double dist2 = DIST2_NoImage( frame0.XYZ(at0), frame0.XYZ(at1) );
      if (dist2 < cut2) {
        double dist = sqrt(dist2);
        int idx = (int)(dist * one_over_spacing);
        // No bounds checking since dist already less than cutoff
        // TODO use addElement
        PairBins.setElement( at0, at1, idx );
      } else
        PairBins.setElement( at0, at1, -1 );
    }
 
  // Loop over frames and lag times 
  int startFrame = 0;
  int endFrame = startFrame + maxlag_;
  int offset = 1;

  for (int frm = startFrame; frm < endFrame; frm += offset)
  {
    coords_->GetFrame( frm, frame0, mask_ );
    for (int lag = 1; lag < maxlag_; lag++)
    {
      int frm1 = frm + lag;
      coords_->GetFrame( frm1, frame1, mask_ );
      // Loop over atom pairs
      for (int at0 = 0; at0 < frame0.Natom(); at0++)
      {
        const double* xyz00 = frame0.XYZ(at0);
        const double* xyz10 = frame1.XYZ(at0);
        Vec3 vec0( xyz10[0] - xyz00[0],
                   xyz10[1] - xyz00[1],
                   xyz10[2] - xyz00[2] );
        for (int at1 = at0+1; at1 < frame0.Natom(); at1++)
        {
          int idx = PairBins.element(at0, at1);
          if (idx != -1) {
            const double* xyz01 = frame0.XYZ(at1);
            const double* xyz11 = frame1.XYZ(at1);
            Vec3 vec1( xyz11[0] - xyz01[0],
                       xyz11[1] - xyz01[1],
                       xyz11[2] - xyz01[2] );
            double dot = vec0 * vec1;
            mat.Element(lag-1, idx) += dot;
          }
        }
      }
    }
  }

  return Analysis::OK;
}
