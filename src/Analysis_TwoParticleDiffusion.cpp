#include <cmath> // ceil, sqrt, fabs
#include <algorithm> // std::max
#include "Analysis_TwoParticleDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"
#include "DistRoutines.h"
#include "Constants.h"
#include "Timer.h"
#include "StringRoutines.h"
#include "ProgressBar.h"

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
  // Output DataFile setup.
  DataFile* dfr = 0;
  DataFile* dft = 0;
  std::string dfname = analyzeArgs.GetStringKey("out");
  if (!dfname.empty()) {
    FileName fname(dfname);
    dfr = setup.DFL().AddDataFile( fname.PrependExt(".drr").Full(), analyzeArgs );
    dft = setup.DFL().AddDataFile( fname.PrependExt(".dtt").Full(), analyzeArgs );
  }
  // Get mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );
  // Set up output data set
  std::string dsname = analyzeArgs.GetStringNext();
  if (dsname.empty())
    dsname = setup.DSL().GenerateDefaultName( "TPD" );
  outDrr_ = setup.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(dsname, "Drr"));
  if (outDrr_ == 0) return Analysis::ERR;
  outDtt_ = setup.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(dsname, "Dtt"));
  if (outDtt_ == 0) return Analysis::ERR;
  Dimension rdim(0.0, rstep_, "R");
  outDrr_->SetDim(1, rdim);
  outDtt_->SetDim(1, rdim);
  Dimension tdim(1.0, 1.0, "lag");
  outDrr_->SetDim(0, tdim);
  outDtt_->SetDim(0, tdim);
  const std::string dfargs = ("xlabel Lag ylabel R");
  if (dfr != 0) {
    dfr->ProcessArgs(dfargs);
    dfr->AddDataSet( outDrr_ );
  }
  if (dft != 0) {
    dft->ProcessArgs(dfargs);
    dft->AddDataSet( outDtt_ );
  }

  mprintf("    TWOPARTICLEDIFFUSION: COORDS set '%s', mask [%s]\n",
          coords_->legend(), mask_.MaskString());
  mprintf("\tOutput Drr set: %s\n", outDrr_->legend());
  mprintf("\tOutput Dtt set: %s\n", outDtt_->legend());
  if (dfr != 0) mprintf("\tWrite Drr set to: %s\n", dfr->DataFilename().full());
  if (dft != 0) mprintf("\tWrite Dtt set to: %s\n", dft->DataFilename().full());
  if (maxlag_ < 1)
    mprintf("\tMax lag will be half the number of frames in input COORDS set.\n");
  else
    mprintf("\tMax lag is %i\n", maxlag_);
  mprintf("\tRmax is %g, Rstep is %g\n", rmax_, rstep_);
  return Analysis::OK;
}

/// Class for accumulating average Drr and Dtt values.
template <class Float> class Accumulator {
  public:
    Accumulator() : n_(0.0), drrMean_(0.0), drrM2_(0.0), dttMean_(0.0), dttM2_(0.0) {} 
    void accumulate(const Float x, const Float y)
    {
      Float delta;

      n_++;
      delta = x - drrMean_;
      drrMean_ += delta / n_;
      drrM2_ += delta * (x - drrMean_);
      delta = y - dttMean_;
      dttMean_ += delta / n_;
      dttM2_ += delta * (y - dttMean_);
    }

    Float DrrMean() const { return drrMean_; };
    Float DrrVariance() const { 
      if (n_ < 2) return 0.0;
      return drrM2_ / (n_ - 1.0); 
    };
    Float DttMean() const { return dttMean_; };
    Float DttVariance() const { 
      if (n_ < 2) return 0.0;
      return dttM2_ / (n_ - 1.0); 
    };

    Float nData() const { return n_; };
  private:
    Float n_;
    Float drrMean_;
    Float drrM2_;
    Float dttMean_;
    Float dttM2_;
};

// Analysis_TwoParticleDiffusion::Analyze()
Analysis::RetType Analysis_TwoParticleDiffusion::Analyze() {
  Timer t_total;
  Timer t_pairloop;
# ifdef TIMER
  Timer t_frame;
  Timer t_precalc;
  Timer t_calc;
# endif
  t_total.Start();
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
  // Allocate matrices: cols=lag, rows=R
  DataSet_MatrixDbl& outputDrr = static_cast<DataSet_MatrixDbl&>( *outDrr_ );
  outputDrr.Allocate2D( maxlag_, numRbins );
  DataSet_MatrixDbl& outputDtt = static_cast<DataSet_MatrixDbl&>( *outDtt_ );
  outputDtt.Allocate2D( maxlag_, numRbins );
  // Use Stats-like class for actual accumulation
  typedef Accumulator<double> Dstats;
  Matrix< Dstats > mat;
  mat.resize( maxlag_, numRbins );

  // For storing atom pair vector and Rbin indices
  Matrix<int> Frame0Idxs;
  Matrix<Vec3> Frame0Vecs;
  // col == 0 means Triangle matrix
  Frame0Idxs.resize( 0, mask_.Nselected() );
  mprintf("\tSize of index matrix: %s\n", ByteString(Frame0Idxs.DataSize(), BYTE_DECIMAL).c_str());
  Frame0Vecs.resize( 0, mask_.Nselected() );
  mprintf("\tSize of vector matrix: %s\n", ByteString(Frame0Vecs.DataSize(), BYTE_DECIMAL).c_str());

  // Loop over frames and lag times 
  int startFrame = 0;
  int endFrame = startFrame + maxlag_;
  int endLag = maxlag_ + 1; // because lag starts at 1
  int offset = 1;

  t_pairloop.Start();
  unsigned int skipOutOfRange = 0;
  double maxD = 0.0;
  ParallelProgress progress( endFrame );
  for (int frm = startFrame; frm < endFrame; frm += offset)
  {
    progress.Update( frm );
#   ifdef TIMER
    t_frame.Start();
#   endif
    coords_->GetFrame( frm, frame0, mask_ );
#   ifdef TIMER
    t_frame.Stop();
#   endif
    // Precalculate atom pair vectors
#   ifdef TIMER
    t_precalc.Start();
#   endif
    int pidx = 0;
    for (int at0 = 0; at0 < frame0.Natom(); at0++) {
      const double* xyz00 = frame0.XYZ(at0);
      for (int at1 = at0+1; at1 < frame0.Natom(); at1++, pidx++)
      {
        const double* xyz01 = frame0.XYZ(at1);
        /// Vector connecting atom pair at time frm
        Vec3 pairVec( xyz01[0] - xyz00[0],
                      xyz01[1] - xyz00[1],
                      xyz01[2] - xyz00[2] );
        // Atom pair distance at time frm
        double d0 = pairVec.Normalize(); 
        maxD = std::max( maxD, d0 );
        // Calculate bin index based on distance at time frm
        // TODO ridx should never be negative, make certain
        // TODO use cutoff^2
        int ridx = (int)(d0 * one_over_spacing);
        // Store
        Frame0Vecs[pidx] = pairVec;
        if (ridx < numRbins)
          Frame0Idxs[pidx] = ridx;
        else
          Frame0Idxs[pidx] = -1;
      }
    }
#   ifdef TIMER
    t_precalc.Stop();
#   endif
    // Inner loop over lag values
    for (int lag = 1; lag < endLag; lag++)
    {
      int frm1 = frm + lag;
#     ifdef TIMER
      t_frame.Start();
#     endif
      coords_->GetFrame( frm1, frame1, mask_ );
#     ifdef TIMER
      t_frame.Stop();
#     endif
      // Loop over atom pairs
      pidx = 0;
      for (int at0 = 0; at0 < frame0.Natom(); at0++)
      {
        const double* xyz00 = frame0.XYZ(at0);
        const double* xyz10 = frame1.XYZ(at0);
        // vec0 is displacement of atom0 from frame0 to frame1
        Vec3 vec0( xyz10[0] - xyz00[0],
                   xyz10[1] - xyz00[1],
                   xyz10[2] - xyz00[2] );
        for (int at1 = at0+1; at1 < frame0.Natom(); at1++, pidx++)
        {
#         ifdef TIMER
          t_calc.Start();
#         endif
          const double* xyz01 = frame0.XYZ(at1);
/*
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
*/
          int ridx = Frame0Idxs[pidx];
          if (ridx != -1) {
            Vec3 const& pairVec = Frame0Vecs[pidx];
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
            mat.element(lag-1, ridx).accumulate( ddl, ddt );
          } else
            skipOutOfRange++;
#         ifdef TIMER
          t_calc.Stop();
#         endif
        } // END inner loop over atoms
      } // END outer loop over atoms
    } // END loop over lag values
  } // END loop over frames
  progress.Finish();
  t_pairloop.Stop();

  mprintf("\t%u pair calculations skipped because R out of range.\n", skipOutOfRange);
  mprintf("\tMax observed distance: %g Ang.\n", maxD);
  // TODO const_iterator?
  for (Matrix< Dstats >::iterator it = mat.begin(); it != mat.end(); ++it)
  {
    outputDrr.AddElement( it->DrrMean() );
    outputDtt.AddElement( it->DttMean() );
  }
/*
  // Normalize: number of values into each lag time is endFrame - lag
  for (int lag = 1; lag < maxlag_; lag++) {
    double norm = 1.0 / (double)(endFrame - lag);
    for (int idx = 0; idx != numRbins; idx++)
      mat.Element(lag-1, idx) *= norm;
  }*/
  t_total.Stop();

  t_total.WriteTiming(1, "Total TwoParticleDiff time:");
  t_pairloop.WriteTiming(2, "Pair loop time:", t_total.Total());
# ifdef TIMER
  t_frame.WriteTiming(  3, "Loop frame time   :", t_pairloop.Total());
  t_precalc.WriteTiming(3, "Loop precalc time :", t_pairloop.Total());
  t_calc.WriteTiming(   3, "Loop Calc time    :", t_pairloop.Total());
# endif
  return Analysis::OK;
}
