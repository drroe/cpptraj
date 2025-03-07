#include "Exec_CrdTransform.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include <algorithm> // std::max,min,sort
#include <cmath> // floor

/// Get the minimum and maximum coordinates in a given frame, store in min and max
static inline void getMaxMin(Frame const& frmIn, Vec3& max, Vec3& min) {
  for (int at = 0; at != frmIn.Natom(); at++) {
    const double* xyz = frmIn.XYZ(at);
    for (int i = 0; i != 3; i++) {
      max[i] = std::max(max[i], xyz[i]);
      min[i] = std::min(min[i], xyz[i]);
    }
  }
}

/** Normalize coordinates between 0 and 1. */
int Exec_CrdTransform::normalizeCoords(DataSet_Coords* crdIn,
                                       DataSet_Coords* crdOut)
const
{
  mprintf("\tNormalize coordinates between 0 and 1.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdOut->legend());
  // Get max and min X Y and Z
  Frame frmIn = crdIn->AllocateFrame();
  crdIn->GetFrame(0, frmIn);
  Vec3 xyzmax( frmIn.XYZ(0) );
  Vec3 xyzmin( frmIn.XYZ(0) );
  getMaxMin(frmIn, xyzmax, xyzmin);
  for (unsigned int idx = 1; idx < crdIn->Size(); idx++) {
    crdIn->GetFrame(idx, frmIn);
    getMaxMin(frmIn, xyzmax, xyzmin);
  }
  mprintf("\tMax: %g %g %g\n", xyzmax[0], xyzmax[1], xyzmax[2]);
  mprintf("\tMin: %g %g %g\n", xyzmin[0], xyzmin[0], xyzmin[0]);
  // Choose the overall max and min
  double max = xyzmax[0];
  max = std::max(max, xyzmax[1]);
  max = std::max(max, xyzmax[2]);
  double min = xyzmin[0];
  min = std::min(min, xyzmin[1]);
  min = std::min(min, xyzmin[2]);

  Vec3 norm( max - min );
  // Protect against bad values
  bool hasBadValues = false;
  static const char dirStr[] = { 'X', 'Y', 'Z' };
  for (int ii = 0; ii != 3; ii++) {
    if (norm[ii] < 0) {
      mprinterr("Error: Min value > max value for %c coordinate.\n", dirStr[ii]);
      hasBadValues = true;
    }
    if (norm[ii] < Constants::SMALL) {
      mprinterr("Error: Min value == max value for %c coordinate.\n", dirStr[ii]);
      hasBadValues = true;
    }
  }
  if (hasBadValues) return 1;
  // Transform coords between 0 and 1
  for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
    crdIn->GetFrame(idx, frmIn);
    for (int crdidx = 0; crdidx < frmIn.size(); crdidx+=3) {
      frmIn[crdidx  ] = (frmIn[crdidx  ] - min) / norm[0];
      frmIn[crdidx+1] = (frmIn[crdidx+1] - min) / norm[1];
      frmIn[crdidx+2] = (frmIn[crdidx+2] - min) / norm[2];
    }
    crdOut->SetCRD(idx, frmIn);
  }

  return 0;
}

/** Transform coordinates by RMS-fitting to an average structure, calculating
  * a new average, then RMS-fitting to that average and so on until a
  * tolerance is reached. Essentially the procedure described by 
  * Klem et al. J. Chem. Theory Comput. 2022, 18, 3218−3230.
  */
int Exec_CrdTransform::iterativeRmsRefinement(AtomMask const& maskIn,
                                              bool useMass,
                                              double tolIn,
                                              DataSet_Coords* crdIn,
                                              DataSet_Coords* crdOut)
const
{
  mprintf("\tRMS iterative refinement.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdOut->legend());
  mprintf("\tAtom mask: %s\n", maskIn.MaskString());
  mprintf("\tRMS Tolerance: %g Ang.\n", tolIn);
  if (useMass)
    mprintf("\tMass-weighting on.\n");
  else
    mprintf("\tMass-weighting off.\n");
  // Do the initial fit to the first frame.
  Frame frmIn = crdIn->AllocateFrame();
  crdIn->GetFrame(0, frmIn);
  Frame selectedRef;
  selectedRef.SetupFrameFromMask( maskIn, crdIn->Top().Atoms() );
  selectedRef.SetCoordinates( frmIn, maskIn );
  // Ensure reference is centered on the origin
  Vec3 refTrans = selectedRef.CenterOnOrigin( useMass );
  // Set up frame for selected incoming atoms
  Frame selectedTgt = selectedRef;
  // Set up frame to hold average 
  Frame avgFrm = selectedTgt;

  double currentTol = tolIn + 9999.0;

  unsigned int iteration = 0;
  while (currentTol > tolIn) {
    avgFrm.ZeroCoords();
    Vec3 tgtTrans(0.0);
    Matrix_3x3 rot(0.0);
    for (unsigned int idx = 0; idx != crdIn->Size(); idx++) {
      crdIn->GetFrame(idx, frmIn);
      selectedTgt.SetCoordinates( frmIn, maskIn );
      selectedTgt.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
      frmIn.Trans_Rot_Trans(tgtTrans, rot, refTrans);
      crdOut->SetCRD(idx, frmIn);
      avgFrm.AddByMask( frmIn, maskIn );
    }
    avgFrm.Divide( (double)crdIn->Size() );
    // Calc RMS of current average to current reference
    currentTol = avgFrm.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
    mprintf("\t%8u %12.4f\n", iteration+1, currentTol);
    // Fit the current average TODO is this necessary?
    avgFrm.Trans_Rot_Trans(tgtTrans, rot, refTrans);
    // Set current average to be new reference
    selectedRef = avgFrm;
    iteration++;
  }

  return 0;
}

/** Strings corresponding to CriterionType */
const char* Exec_CrdTransform::CriterionStr_[] = {
  "comp_sim", "sim_to_medoid", "No criterion"
};


/** Trim a desired percentage of outliers (most dissimilar) from the COORDS
  * data set by calculating the largest complement similarity.
  */
int Exec_CrdTransform::trimOutliers(int n_trimmed, double cutoffIn,
                                    ExtendedSimilarity::MetricType metric,
                                    CriterionType criterion,
                                    DataSet_Coords* crdIn,
                                    DataSet_Coords* crdOut)
const
{
  mprintf("\tTrimming outliers.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdOut->legend());
  mprintf("\tUsing metric: %s\n", ExtendedSimilarity::metricStr(metric));
  mprintf("\tCriterion: %s\n", CriterionStr_[criterion]);
  unsigned int Ncoords = crdIn->Top().Natom() * 3;
  unsigned int Nframes = crdIn->Size();
  Frame frmIn = crdIn->AllocateFrame();
  mprintf("\t'%s' has %u coordinates, %u frames.\n", crdIn->legend(), Ncoords, Nframes);
  if (Nframes < 2) {
    mprintf("Warning: Less than 2 frames, nothing to trim.\n"); // TODO something with crdOut?
    return 0;
  }
  // Specify n_trimmed or cutoff, but not both.
  if (n_trimmed < 0 && cutoffIn < 0) {
    mprinterr("Internal Error: Must specify either number to trim or cutoff.\n");
    return 1;
  }
  if (n_trimmed >= 0 && cutoffIn > 0) {
    mprinterr("Error: Must specify either number to trim or cutoff, but not both.\n");
    return 1;
  }
  // Set cutoff value
  unsigned int cutoff;
  if (n_trimmed >= 0) {
    cutoff = (unsigned int)n_trimmed;
    mprintf("\t# to trim: %i\n", n_trimmed);
  } else {
    cutoff = (unsigned int)(floor(Nframes * cutoffIn));
    mprintf("\tFraction of outliers to remove: %f\n", cutoffIn);
  }
  mprintf("\tUsing cutoff value: %u\n", cutoff);
  // Do extended similarity calculation for each frame
  ExtendedSimilarity ExtSim;
  if (ExtSim.SetOpts( metric, crdIn->Size(), Ncoords )) {
    mprinterr("Error: Extended similarity setup for trim failed.\n");
    return 1;
  }
  ExtendedSimilarity::Darray csimvals = ExtSim.CalculateCompSim( *crdIn );
  if (csimvals.empty()) {
    mprinterr("Error: No comparitive similarity values for trim calculated.\n");
    return 1;
  }
  // DEBUG
/*  CpptrajFile dbg;
  dbg.OpenWrite("test.cpptraj.out");
  for (ExtendedSimilarity::Darray::const_iterator it = csimvals.begin(); it != csimvals.end(); ++it)
    dbg.Printf("%8li %16.8f\n", it - csimvals.begin(), *it);
  dbg.CloseFile();*/
  // Array type for sorting by sim. value while preserving indices
  typedef std::pair<unsigned int, double> IdxValPairType;
  std::vector<IdxValPairType> comp_sims;
  comp_sims.reserve( crdIn->Size() );
  // For sorting IdxValPairType by values
  struct IdxValPairCmp {
    inline bool operator()(IdxValPairType const& first, IdxValPairType const& second) {
      return (first.second < second.second);
    }
  };
  struct ReverseIdxValPairCmp {
    inline bool operator()(IdxValPairType const& first, IdxValPairType const& second) {
      return (second.second < first.second);
    }
  };

  if (criterion == COMP_SIM) {
    // ----- Comp sim ------------------
    // Place values in sortable array and sort lowest to highest.
    for (unsigned int idx = 0; idx < crdIn->Size(); idx++)
      comp_sims.push_back( IdxValPairType(idx, csimvals[idx]) );
    std::sort(comp_sims.begin(), comp_sims.end(), IdxValPairCmp());
  } else if (criterion == SIM_TO_MEDOID) {
    // ----- SIM TO MEDOID -------------
    mprintf("\tMedoid index= %li\n", ExtSim.MedoidIndex() + 1);
    if (ExtSim.MedoidIndex() < 0) {
      mprinterr("Error: Medoid index is < 0\n");
      return 1;
    }
    // Get comp sim values to medoid (excluding medoid)
    unsigned int medoid_index = (unsigned int)ExtSim.MedoidIndex();
    Frame medoid = crdIn->AllocateFrame();
    crdIn->GetFrame(medoid_index, medoid);
    //mprintf("[");
    for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
      if (idx != medoid_index) {
        //mprintf("DEBUG\n");
        crdIn->GetFrame(idx, frmIn);
        comp_sims.push_back( IdxValPairType(idx, ExtSim.CalculateCompSim(frmIn, medoid)) );
        //mprintf(" %10.8g", comp_sims.back().second);
      }
    }
    //mprintf("]\n");
    // Sort highest to lowest
    std::sort(comp_sims.begin(), comp_sims.end(), ReverseIdxValPairCmp());
  }
  // Remove frames up to the cutoff
  std::vector<bool> keepFrame(comp_sims.size(), true);
  if (debug_ > 0) mprintf("[");
  unsigned int nToRemove = 0;
  for (unsigned int idx = 0; idx < cutoff; idx++) {
    if (debug_ > 0) mprintf(" %u", comp_sims[idx].first+1);
    keepFrame[comp_sims[idx].first] = false;
    nToRemove++;
  }
  if (debug_ > 0) mprintf("]\n");
  mprintf("\tRemoving %u frames.\n", nToRemove);

  // Populate the output trajectory
  for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
    if (keepFrame[idx]) {   
      crdIn->GetFrame(idx, frmIn);
      crdOut->AddFrame( frmIn );
    }
  }  
  return 0;
}


// Exec_CrdTransform::Help()
void Exec_CrdTransform::Help() const
{
  mprintf("\t<input crd set> [name <output crd set>]\n"
          "\t{ rmsrefine [mask <mask>] [mass] [rmstol <tolerance>] |\n"
          "\t  normcoords |\n"
          "\t  trim [metric <metric>] [{ntrimmed <#>|cutoff <val>}]\n"
          "\t       [criterion {comp|medoid}]]\n"
          "\t}\n");
  mprintf("  <metric> = %s\n", ExtendedSimilarity::MetricKeys().c_str());
  mprintf("  Transform a trajectory in one of several ways:\n"
          "  - rmsrefine  : Do an iterative RMS refinement of all frames.\n"
          "  - normcoords : Normalize coordinates between 0.0 and 1.0.\n"
          "  - trim       : Remove trajectory frames using extended similarity metrics.\n");
}

// Exec_CrdTransform::Execute()
Exec::RetType Exec_CrdTransform::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  AtomMask mask;
  bool useMass = false;
  double rmsTol = -1.0;

  int n_trimmed = -1;
  double cutoff = -1.0;
  ExtendedSimilarity::MetricType metric = ExtendedSimilarity::NO_METRIC;
  CriterionType criterion = NO_CRITERION;
  // Determine mode
  bool lengthWillBeModified = false;
  enum ModeType { RMSREFINE = 0, NORMCOORDS, TRIM, UNSPECIFIED };
  ModeType mode = UNSPECIFIED;
  if (argIn.hasKey("rmsrefine")) {
    mode = RMSREFINE;
    mask.SetMaskString( argIn.GetStringKey("mask") );
    useMass = argIn.hasKey("mass");
    rmsTol = argIn.getKeyDouble("rmstol", 0.0001);
  } else if (argIn.hasKey("normcoords")) {
    mode = NORMCOORDS;
  } else if (argIn.hasKey("trim")) {
    lengthWillBeModified = true;
    mode = TRIM;
    n_trimmed = argIn.getKeyInt("ntrimmed", -1);
    cutoff = argIn.getKeyDouble("cutoff", -1.0);
    std::string mstr = argIn.GetStringKey("metric");
    if (!mstr.empty()) {
      metric = ExtendedSimilarity::TypeFromKeyword( mstr );
      if (metric == ExtendedSimilarity::NO_METRIC) {
        mprinterr("Error: Metric '%s' not recognized.\n", mstr.c_str());
        return CpptrajState::ERR;
      }
    } else {
      metric = ExtendedSimilarity::MSD;
    }
    std::string cstr = argIn.GetStringKey("criterion");
    if (!cstr.empty()) {
      if (cstr == "comp")
        criterion = COMP_SIM;
      else if (cstr == "medoid")
        criterion = SIM_TO_MEDOID;
      else {
        mprinterr("Error: Unrecognized criterion: %s\n", cstr.c_str());
        return CpptrajState::ERR;
      }
    } else {
      criterion = COMP_SIM;
    }
  } else {
    mprinterr("Error: Expected 'trim', 'rmsrefine', or 'normcoords'\n");
    return CpptrajState::ERR;
  }
 
  // Get COORDS set
  std::string outname = argIn.GetStringKey("name");
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify input COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' has no frames.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  if (CRD->Type() == DataSet::TRAJ) {
    mprinterr("Error: TRAJ sets not yet supported.\n"); // FIXME
    return CpptrajState::ERR;
  }
  // ----- No more input args processed below here ---------

  // Set up output set if needed TODO FRAMES set
  DataSet_Coords* OUT = 0;
  if (!outname.empty()) {
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, MetaData(outname) );
    if (OUT == 0) {
      mprinterr("Error: Could not create output coords set %s\n", outname.c_str());
      return CpptrajState::ERR;
    }
    if (OUT->CoordsSetup( CRD->Top(), CRD->CoordsInfo() )) {
      mprinterr("Error: Could not set up output coords set %s\n", OUT->legend());
      return CpptrajState::ERR;
    }
    if (!lengthWillBeModified) OUT->Allocate( DataSet::SizeArray(1, CRD->Size()) );
  }
  bool needToDeleteCRD = false;
  if (OUT == 0) {
    if (!lengthWillBeModified)
      OUT = CRD;
    else {
      // Need to replace CRD with a new set. Remove CRD from master DSL.
      needToDeleteCRD = true;
      State.DSL().PopSet( CRD );
      OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, CRD->Meta() );
      if (OUT == 0) {
        mprinterr("Error: Could not replace coords set %s\n", CRD->legend());
        return CpptrajState::ERR;
      }
      if (OUT->CoordsSetup( CRD->Top(), CRD->CoordsInfo() )) {
        mprinterr("Error: Could not set up replacement coords set %s\n", OUT->legend());
        return CpptrajState::ERR;
      }
    }
  }

  // Set up mask
  if (mask.MaskStringSet()) {
    if (CRD->Top().SetupIntegerMask( mask )) {
      mprinterr("Error: Could not set up mask.\n");
      return CpptrajState::ERR;
    }
    mask.MaskInfo();
  }

  int err = 0;
  switch (mode) {
    case RMSREFINE  : err = iterativeRmsRefinement(mask, useMass, rmsTol, CRD, OUT); break;
    case NORMCOORDS : err = normalizeCoords(CRD, OUT); break;
    // TODO pass in criterion
    case TRIM       : err = trimOutliers(n_trimmed, cutoff, metric, criterion, CRD, OUT); break;
    default         : err = 1; break;
  }
  if (needToDeleteCRD) delete CRD;
  if (err != 0) return CpptrajState::ERR;

  return CpptrajState::OK;
}
