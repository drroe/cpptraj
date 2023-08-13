#include "Hbond.h"
#include "../Constants.h"
#include "../DataSet_integer.h"

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
Hbond::Hbond() : dist_(0.0), angle_(0.0), data_(0), A_(-1), H_(-1), D_(-1), frames_(0) {}

/** CONSTRUCTOR - New hydrogen bond */
Hbond::Hbond(DataSet_integer* s, int ia, int ih, int id, Iarray const& splits) :
  dist_(0), angle_(0), data_(s), A_(ia), H_(ih), D_(id), frames_(0)
{
  if (!splits.empty()) {
    partsDist_.resize(splits.size()+1);
    partsAng_.resize(splits.size()+1);
  }
}

/** Update time series */
void Hbond::Update(double distIn, double angIn, int fnum, Iarray const& splitFrames, int onum) {
  dist_ += distIn;
  angle_ += angIn;
  ++frames_;
  if (data_ != 0) data_->AddVal(fnum, 1);
  if (!splitFrames.empty()) {
    //mprintf("DEBUG: SPLIT: onum= %i partsDistSize= %zu splitFramesSize= %zu\n", onum, partsDist_.size(), splitFrames.size());
    // Find the correct part NOTE assumes onum never out of range
    unsigned int part = 0;
    while (part < splitFrames.size() && onum >= splitFrames[part]) part++;
    partsDist_[part].accumulate( distIn );
    partsAng_[part].accumulate( angIn );
  }
}

/** Calculate average distance and angle for hbond. */
void Hbond::CalcAvg() {
  double dFrames = (double)frames_;
  dist_ /= dFrames;
  angle_ /= dFrames;
  angle_ *= Constants::RADDEG;
}

