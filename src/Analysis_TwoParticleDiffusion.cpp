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
#ifdef _OPENMP
# include <omp.h>
#endif

/** CONSTRUCTOR */
Analysis_TwoParticleDiffusion::Analysis_TwoParticleDiffusion() :
  coords_(0),
  outDrr_(0),
  outDtt_(0),
  rmax_(0),
  rstep_(0),
  maxlag_(0)
{}

// Analysis_TwoParticleDiffusion::Help()
void Analysis_TwoParticleDiffusion::Help() const {
  mprintf("\t[crdset <crd set>] [<name>] [<mask>] [out <filename>]\n"
          "\t[stop <maxlag>] rmax <max> rstep <step> [noimage]\n");
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
  imageOpt_.InitImaging( !analyzeArgs.hasKey("noimage") );
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
  if (imageOpt_.UseImage())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tDistamces will not be imaged.\n");
  return Analysis::OK;
}

/// Class for accumulating average Drr and Dtt values.
template <class Float> class Accumulator {
  public:
    //Accumulator() : n_(0), drrMean_(0.0), drrM2_(0.0), dttMean_(0.0), dttM2_(0.0) {}
    Accumulator() : n_(0), drrMean_(0.0), dttMean_(0.0) {}
    // NOTE: The below summation is slower but much more accurate and stable
/* 
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
*/
    void accumulate(const Float x, const Float y) {
      n_++;
      drrMean_ += x;
      dttMean_ += y;
    }
    Float DrrMean() const { if (n_ < 1) return 0.0; else return drrMean_ / n_; }
    Float DttMean() const { if (n_ < 1) return 0.0; else return dttMean_ / n_; }

    Float nData() const { return n_; };
#   ifdef _OPENMP
    void operator+=(Accumulator const& rhs) {
      n_ += rhs.n_;
      drrMean_ += rhs.drrMean_;
      dttMean_ += rhs.dttMean_;
    }
#   endif
  private:
    Float n_;
    Float drrMean_;
//    Float drrM2_;
    Float dttMean_;
//    Float dttM2_;
};

// Analysis_TwoParticleDiffusion::Analyze()
Analysis::RetType Analysis_TwoParticleDiffusion::Analyze() {
  Timer t_total;
  Timer t_pairloop;
# ifdef TIMER
  Timer t_frame;
  Timer t_precalc;
  Timer t_calc;
  Timer t_mat;
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
  // Determine imaging type
  imageOpt_.SetupImaging( coords_->CoordsInfo().TrajBox().HasBox() );
  if (imageOpt_.ImagingEnabled())
    mprintf("\tImaging is enabled.\n");
  else
    mprintf("\tImaging is disabled.\n");

  // Determine size of the 'R' dimension
  double one_over_spacing = 1.0 / rstep_;
  int numRbins = (int)ceil(rmax_ / rstep_);
  mprintf("\tNumber of bins in the 'R' dimension: %i\n", numRbins);
  // Allocate matrices: cols=lag, rows=R
  DataSet_MatrixDbl& outputDrr = static_cast<DataSet_MatrixDbl&>( *outDrr_ );
  outputDrr.Allocate2D( maxlag_, numRbins );
  mprintf("\tSize of output Drr matrix: %s\n",
          ByteString(outDrr_->Size()*sizeof(double), BYTE_DECIMAL).c_str());
  DataSet_MatrixDbl& outputDtt = static_cast<DataSet_MatrixDbl&>( *outDtt_ );
  outputDtt.Allocate2D( maxlag_, numRbins );
  mprintf("\tSize of output Dtt matrix: %s\n",
          ByteString(outDtt_->Size()*sizeof(double), BYTE_DECIMAL).c_str());
  // Use Stats-like class for actual accumulation
  typedef Accumulator<double> Dstats;
# ifdef _OPENMP
  // Each thread needs their own accumulator to avoid memory clashes
  typedef std::vector< Matrix< Dstats > > DMarrayType;
  DMarrayType DMarray;
# pragma omp parallel
  {
#   pragma omp master
    {
      DMarray.resize( omp_get_num_threads() );
    }
  }
  mprintf("\tParallelizing calculation with %zu threads.\n", DMarray.size());
  for (DMarrayType::iterator dm = DMarray.begin(); dm != DMarray.end(); ++dm)
    dm->resize( maxlag_, numRbins );
# else /* not _OPENMP */
  Matrix< Dstats > mat;
  // TODO - test reversing order of indices to see if memory pattern improved.
  mat.resize( maxlag_, numRbins );
# endif /* _OPENMP */
# ifdef USE_TPD_MEMCACHE
  mprintf("\tUsing memory to cache frame(t) related values.\n");
  // For storing atom pair vector and Rbin indices
  Matrix<int> Frame0Idxs;
  Matrix<Vec3> Frame0Vecs;
  // col == 0 means Triangle matrix
  mprintf("\tSize of index matrix: %s\n",
          ByteString(Matrix<int>::sizeInBytes(0, mask_.Nselected()), BYTE_DECIMAL).c_str());
  mprintf("\tSize of vector matrix: %s\n",
          ByteString(Matrix<Vec3>::sizeInBytes(0, mask_.Nselected()), BYTE_DECIMAL).c_str());
  Frame0Idxs.resize( 0, mask_.Nselected() );
  Frame0Vecs.resize( 0, mask_.Nselected() );
# else
  mprintf("\tNot using memory to cache frame(t) related values.\n");
# endif
  // Loop over frames and lag times 
  int startFrame = 0;
  int endFrame = startFrame + maxlag_;
  int offset = 1;
  int endLag = maxlag_ + 1; // because lag starts at 1
  if ( ((endFrame-offset)+(endLag-1)) > (int)coords_->Size()) {
    mprintf("Warning: End lag (%i) > number of frames %zu;",
            endLag, coords_->Size());
    endLag = (int)coords_->Size() - endFrame + offset + 1;
    mprintf(" setting to %i\n", endLag);
  }
  double cut2 = rmax_ * rmax_;

  t_pairloop.Start();
  unsigned int skipOutOfRange = 0;
  double maxD = 0.0;
  double maxD2 = 0.0;
  ParallelProgress progress( endFrame );
  int frm;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(frm, mythread) firstprivate(progress, frame0, frame1) reduction(max: maxD, maxD2) reduction(+: skipOutOfRange)
  {
  mythread = omp_get_thread_num();
  progress.SetThread(mythread);
  Matrix< Dstats >& mat = DMarray[mythread];
# pragma omp for
# endif
  for (frm = startFrame; frm < endFrame; frm += offset)
  {
    progress.Update( frm );
#   ifdef TIMER
    t_frame.Start();
#   endif
    coords_->GetFrame( frm, frame0, mask_ );
    // Determine imaging
    if (imageOpt_.ImagingEnabled())
      imageOpt_.SetImageType( frame0.BoxCrd().Is_X_Aligned_Ortho() );
    //mprintf("DEBUG: Frame %i box:", frm+1);
    //frame0.BoxCrd().PrintInfo();
#   ifdef TIMER
    t_frame.Stop();
#   endif
#   ifdef USE_TPD_MEMCACHE
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
        Vec3 pairVec;
        switch (imageOpt_.ImagingType()) {
          case ImageOption::NO_IMAGE  : pairVec = Vec3(xyz01) - Vec3(xyz00); break;
          case ImageOption::ORTHO     : pairVec = MinImagedVec(xyz01, xyz00, frame0.BoxCrd()); break;
          case ImageOption::NONORTHO  : pairVec = MinImagedVec(xyz01, xyz00, frame0.BoxCrd().UnitCell(), frame0.BoxCrd().FracCell()); break;
        }
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
#   endif /* USE_TPD_MEMCACHE */
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
#     ifdef USE_TPD_MEMCACHE
      pidx = 0;
#     endif
      for (int at0 = 0; at0 < frame0.Natom(); at0++)
      {
        const double* xyz00 = frame0.XYZ(at0);
        const double* xyz10 = frame1.XYZ(at0);
        // vec0 is displacement of atom0 from frame0 to frame1
        Vec3 vec0( xyz10[0] - xyz00[0],
                   xyz10[1] - xyz00[1],
                   xyz10[2] - xyz00[2] );
#       ifdef USE_TPD_MEMCACHE
        for (int at1 = at0+1; at1 < frame0.Natom(); at1++, pidx++)
#       else
        for (int at1 = at0+1; at1 < frame0.Natom(); at1++)
#       endif
        {
#         ifdef TIMER
          t_calc.Start();
#         endif
#         ifdef USE_TPD_MEMCACHE
          int ridx = Frame0Idxs[pidx];
#         else
          /// Vector connecting atom pair at time frm
          const double* xyz01 = frame0.XYZ(at1);
          Vec3 pairVec(0.0);
          switch (imageOpt_.ImagingType()) {
            case ImageOption::NO_IMAGE  : pairVec = Vec3(xyz01) - Vec3(xyz00); break;
            case ImageOption::ORTHO     : pairVec = MinImagedVec(xyz01, xyz00, frame0.BoxCrd()); break;
            case ImageOption::NONORTHO  : pairVec = MinImagedVec(xyz01, xyz00, frame0.BoxCrd().UnitCell(), frame0.BoxCrd().FracCell()); break;
          }
/*
          Vec3 pairVec( xyz01[0] - xyz00[0],
                        xyz01[1] - xyz00[1],
                        xyz01[2] - xyz00[2] );
*/
          // Atom pair distance at time frm
          double dist2 = pairVec.Magnitude2();
          int ridx;
          if (dist2 < cut2) {
            double d0 = sqrt(dist2);
            double one_over_d0 = 1.0 / d0;
            pairVec *= one_over_d0;
            maxD = std::max( maxD, d0 );
            // Calculate bin index based on distance at time frm
            // TODO ridx should never be negative, make certain
            ridx = (int)(d0 * one_over_spacing);
          } else {
            ridx = -1;
            maxD2 = std::max( maxD2, dist2 );
          }
#         endif /* USE_TPD_MEMCACHE */
          if (ridx != -1) {
#           ifdef USE_TPD_MEMCACHE
            Vec3 const& pairVec = Frame0Vecs[pidx];
#           endif
            const double* xyz01 = frame0.XYZ(at1);
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
#           ifdef TIMER
            t_mat.Start();
#           endif
            mat.element(lag-1, ridx).accumulate( ddl, ddt );
#           ifdef TIMER
            t_mat.Stop();
#           endif
          } else
            skipOutOfRange++;
#         ifdef TIMER
          t_calc.Stop();
#         endif
        } // END inner loop over atoms
      } // END outer loop over atoms
    } // END loop over lag values
  } // END loop over frames
# ifdef _OPENMP
  } // END pragma omp parallel
  // Consolidate threads
  for (unsigned int thread = 1; thread < DMarray.size(); thread++)
    for (unsigned int idx = 0; idx < DMarray[thread].size(); idx++)
      DMarray[0][idx] += DMarray[thread][idx];
# endif
  progress.Finish();
  t_pairloop.Stop();

  mprintf("\t%u pair calculations skipped because R out of range.\n", skipOutOfRange);
# ifndef USE_TPD_MEMCACHE
  maxD = std::max( maxD, sqrt(maxD2) );
# endif
  mprintf("\tMax observed distance: %g Ang.\n", maxD);
  // TODO const_iterator?
# ifdef _OPENMP
  Matrix< Dstats > const& mat = DMarray[0];
# endif
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
  t_mat.WriteTiming(4, "Matrix Store time:", t_calc.Total());
# endif
  return Analysis::OK;
}
