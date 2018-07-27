#include <cmath> // fabs, sqrt
#include <algorithm> // sort, max
#include "MolSurf.h"
#include "DistRoutines.h"

int MolSurf::GetNeighbors() {
  // Find max radius out of all atoms
  double maxrad = 0.0;
  for (AtomArray::const_iterator at = atoms_.begin(); at != atoms_.end(); ++at)
    maxrad = std::max( at->Radius(), maxrad );

  double probe_diam = 2 * probe_radius_;
  double cutoff = (2 * maxrad + probe_diam);

  for (AtomArray::iterator at0 = atoms_.begin(); at0 != atoms_.end(); ++at0)
  {
    at0->SetBuried(false);
    double d_ext = probe_diam + at0->Radius();
    // Find all atoms near this atom
    for (AtomArray::iterator at1 = at0 + 1; at1 != atoms_.end(); ++at1)
    {
      at1->SetBuried(false);
      bool isFar = (fabs(at0->XYZ()[0] - at1->XYZ()[0]) > cutoff ||
                    fabs(at0->XYZ()[1] - at1->XYZ()[1]) > cutoff ||
                    fabs(at0->XYZ()[2] - at1->XYZ()[2]) > cutoff);
      if (!isFar) {
        double dist = sqrt(DIST2_NoImage(at0->XYZ(), at1->XYZ()));
        if (dist < (d_ext + at1->Radius())) {
          // at0 and at1 are neighbors
          if (at0->Radius() + dist < at1->Radius())
            at0->SetBuried(true);
          if (at1->Radius() + dist < at0->Radius())
            at1->SetBuried(true);
          at0->AddNeighbor( NbrIdx(at1-atoms_.begin(), dist) );
          at1->AddNeighbor( NbrIdx(at0-atoms_.begin(), dist) );
        }
      }
    } // END loop over all other atoms
    at0->SortNeighbors();
  }

  return 0;
}

void MolSurf::UpdatePositions(Frame const& frmIn) {
  AtomArray::iterator at = atoms_.begin();
  for (AtomMask::const_iterator it = mask_.begin(); it != mask_.end(); ++it, ++at)
    at->UpdatePosition( frmIn.XYZ(*it) );
}

double MolSurf::CalcSurface(Frame const& frmIn) {

  UpdatePositions(frmIn);

  GetNeighbors();

  return 0.0;
}
// -----------------------------------------------------------------------------
void MolSurf::atom::SortNeighbors() {
  std::sort(neighbors_.begin(), neighbors_.end());
}
