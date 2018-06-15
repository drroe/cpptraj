#include <cmath> // ceil, sqrt, fabs
#include <algorithm> // std::max
#include "Analysis_TwoParticleDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"
#include "DistRoutines.h"
#include "Constants.h"
#include "OnlineVarT.h"

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
  if (df != 0) {
    df->ProcessArgs("xlabel lag ylabel r");
    df->AddDataSet( out_ );
  }

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
  double one_over_spacing = 1.0 / rstep_;
  int numRbins = (int)ceil(rmax_ / rstep_);
  mprintf("\tNumber of bins in the 'R' dimension: %i\n", numRbins);
  // Allocate matrix: cols=lag, rows=R
  DataSet_MatrixDbl& outputMatrix = static_cast<DataSet_MatrixDbl&>( *out_ );
  outputMatrix.Allocate2D( maxlag_-1, numRbins );
  // Use Stats for actual accumulation
  Matrix< Stats<double> > mat;
  mat.resize( maxlag_-1, numRbins );
/*
  // Store atom pair Rbin indices via calculating initial pair distances.
  Matrix<int> PairBins;
  // col == 0 means Triangle matrix
  PairBins.resize( 0, mask_.Nselected() );
  coords_->GetFrame( 0, frame0, mask_ );
  double cut2 = rmax_ * rmax_;
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
*/

  // Loop over frames and lag times 
  int startFrame = 0;
  int endFrame = startFrame + maxlag_;
  int offset = 1;

  unsigned int skipOutOfRange = 0;
  double maxD = 0.0;
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
        // vec0 is displacement of atom0 from frame0 to frame1
        Vec3 vec0( xyz10[0] - xyz00[0],
                   xyz10[1] - xyz00[1],
                   xyz10[2] - xyz00[2] );
        for (int at1 = at0+1; at1 < frame0.Natom(); at1++)
        {
          const double* xyz01 = frame0.XYZ(at1);
          /// Vector connecting atom pair at time frm TODO precalculate
          Vec3 pairVec( xyz01[0] - xyz00[0],
                        xyz01[1] - xyz00[1],
                        xyz01[2] - xyz00[2] );
          // Atom pair distance at time frm TODO precalculate
          double d0 = pairVec.Normalize(); 
          maxD = std::max( maxD, d0 );
          // Calculate bin index based on distance at time frm
          // TODO idx should never be negative, make certain
          int idx = (int)(d0 * one_over_spacing);
          if (idx < numRbins) {
            const double* xyz11 = frame1.XYZ(at1);
            // vec1 is displacement of atom1 from frame 0 to frame1
            Vec3 vec1( xyz11[0] - xyz01[0],
                       xyz11[1] - xyz01[1],
                       xyz11[2] - xyz01[2] );
            // Longitudinal part (Drr), XY
            double ddl = (vec0[0]*pairVec[0] + vec0[1]*pairVec[1]) *
                         (vec1[0]*pairVec[0] + vec1[1]*pairVec[1]);
            // Orthogonal unit vector, switching Y to X and X to -Y, perpendicular.
            double px = pairVec[1];
            double py = -pairVec[0];
            // Transverse part (theta, Dtt) XY
            double ddt = (vec0[0]*px + vec0[1]*py) *
                         (vec1[0]*px + vec1[1]*py);
            mat.element(lag-1, idx).accumulate( ddl );
          } else
            skipOutOfRange++;
        } // END inner loop over atoms
      } // END outer loop over atoms
    } // END loop over lag values
  } // END loop over frames

  mprintf("\t%u pair calculations skipped because R out of range.\n", skipOutOfRange);
  mprintf("\tMax observed distance: %g Ang.\n", maxD);
  // TODO const_iterator?
  for (Matrix< Stats<double> >::iterator it = mat.begin(); 
                                               it != mat.end(); ++it)
    outputMatrix.AddElement( it->mean() );
/*
  // Normalize: number of values into each lag time is endFrame - lag
  for (int lag = 1; lag < maxlag_; lag++) {
    double norm = 1.0 / (double)(endFrame - lag);
    for (int idx = 0; idx != numRbins; idx++)
      mat.Element(lag-1, idx) *= norm;
  }*/
  return Analysis::OK;
}
