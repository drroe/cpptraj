#ifndef INC_HB_BRIDGE_H
#define INC_HB_BRIDGE_H
#include <vector>
class DataSet_integer;
namespace Cpptraj {
namespace HB {
/// Track solvent bridge between 2 or more solute residues.
class Bridge {
  public:
    typedef std::vector<int> Iarray;
    /// CONSTRUCTOR - new bridge
    Bridge(DataSet_integer*, Iarray const&);
#   ifdef MPI
    /// Constructor - new bridge with given # frames
    Bridge(DataSet_integer* bds, int f) : data_(bds), frames_(f) {}
    /// Increment number of frames
    void Combine(int n) { frames_ += n; }
    /// Set up parts with given int array containing # frames for each part
    void SetupParts(unsigned int nparts, const int* ivals) {
      if (nparts == 0) return;
      partsFrames_.clear();
      partsFrames_.reserve( nparts );
      for (unsigned int idx = 0; idx != nparts; idx++)
        partsFrames_.push_back( ivals[idx] );
    }
    /// Update parts with given int array containing # frames for each part
    void CombineParts(unsigned int nparts, const int* ivals) {
      for (unsigned int idx = 0; idx != nparts; idx++)
        partsFrames_[idx] += ivals[idx];
    }
#   endif
    /// \return internal data set
    DataSet_integer* Data() const { return data_; }
    int Frames()            const { return frames_; }
    /// Update frames/time series
    void Update(int, Iarray const&, int);
    /// \return true if bridge has more frames than rhs.
    bool operator()(Bridge const& rhs) const {
      if (frames_ > rhs.frames_)
        return true;
      else if (frames_ < rhs.frames_)
        return false;
      else
        return false;
    }
    // Summary by parts
    unsigned int Nparts() const { return partsFrames_.size(); }
    unsigned int PartFrames(unsigned int idx) const { return (unsigned int)partsFrames_[idx]; }
  private:
    DataSet_integer* data_; ///< Hold time series data TODO
    int frames_;            ///< # frames this bridge has been present.
    Iarray partsFrames_;    ///< Hold # frames bridge present for each part.

};
}
}
#endif
