#include "Bridge.h"
#include "../DataSet_integer.h"

using namespace Cpptraj::HB;

/// CONSTRUCTOR - new bridge
Bridge::Bridge(DataSet_integer* bds, Iarray const& splits) : data_(bds), frames_(0) {
  if (!splits.empty())
    partsFrames_.assign(splits.size()+1, 0);
}

    /// Update frames/time series
    void Bridge::Update(int fnum, Iarray const& splitFrames, int onum) {
      ++frames_;
      if (data_ != 0) data_->AddVal(fnum, 1);
      if (!splitFrames.empty()) {
        // Find the correct part NOTE assumes onum never out of range
        unsigned int part = 0;
        while (part < splitFrames.size() && onum >= splitFrames[part]) part++;
        partsFrames_[part]++;
      }
    }

