#include <cmath>
#include "Accumulator.h"
#include "Constants.h"

double Accumulator_Scalar::StdDev() const {
  return sqrt(popVariance());
}

// =============================================================================
/** Accumulate torsion in degrees. */
void Accumulator_Torsion::accumulate(double thetaDeg) {
  n_++;
  double thetaRad = thetaDeg * Constants::DEGRAD;
  onlineAccumulate( sin(thetaRad), n_, meanSp_, MS2_ );
  onlineAccumulate( cos(thetaRad), n_, meanCp_, MC2_ );
}

double Accumulator_Torsion::mean() const {
  return atan2(meanSp_, meanCp_) * Constants::RADDEG;
}

double Accumulator_Torsion::popVariance() const {
  double meanRp = sqrt((meanCp_*meanCp_) + (meanSp_*meanSp_));
  return (1 - meanRp) * Constants::RADDEG;
}

double Accumulator_Torsion::StdDev() const {
  double meanRp = sqrt((meanCp_*meanCp_) + (meanSp_*meanSp_));
  return sqrt(-2 * log(meanRp)) * Constants::RADDEG;
}
