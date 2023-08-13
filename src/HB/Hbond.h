#ifndef INC_HB_HBOND_H
#define INC_HB_HBOND_H
#include <cstddef> // size_t
#include <vector>
#include "../OnlineVarT.h"
class DataSet_integer;
namespace Cpptraj {
namespace HB {
/// Track specific hydrogen bond.
class Hbond {
  public:
    typedef std::vector<int> Iarray;
    /// CONSTRUCTOR
    Hbond();
    /// New hydrogen bond
    Hbond(DataSet_integer*, int, int, int, Iarray const&);
#   ifdef _OPENMP
    /// Just record that hbond exists
    Hbond(double d, double a, int ia, int ih, int id) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(0) {}
    /// This version is for UV hbonds; a 1 in frames_ indicates soluteDonor
    Hbond(double d, double a, int ia, int ih, int id, int sd) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(sd) {}
#   endif
    double Dist()           const { return dist_;   }
    double Angle()          const { return angle_;  }
    int Frames()            const { return frames_; }
    int A()                 const { return A_;      }
    int H()                 const { return H_;      }
    int D()                 const { return D_;      }
    DataSet_integer* Data() const { return data_; }
    // Summary by parts
    unsigned int Nparts() const { return partsDist_.size(); }
    unsigned int PartFrames(unsigned int idx) const { return (unsigned int)partsDist_[idx].nData(); }
    double PartFrac(unsigned int idx, unsigned int Nframes) const { return partsDist_[idx].nData() / (double)Nframes; }
    Stats<double> const& PartDist(unsigned int idx)  const { return partsDist_[idx]; }
    Stats<double> const& PartAngle(unsigned int idx) const { return partsAng_[idx]; }
#   ifdef MPI
    /// CONSTRUCTOR - New hydrogen bond with given # frames
    Hbond(double d, double a, DataSet_integer* s, int ia, int ih, int id, int n) :
      dist_(d), angle_(a), data_(s), A_(ia), H_(ih), D_(id), frames_(n) {}
    /// Update distance/angle/number frames
    void Combine(double d, double a, int n) {
      dist_ += d;
      angle_ += a;
      frames_ += n;
    }
    /// Set up parts with given double array containing N, mean, and M2 for each part
    void SetupParts(unsigned int nparts, const double* dvals) {
      if (nparts == 0) return;
      partsDist_.clear();
      partsAng_.clear();
      partsDist_.reserve( nparts );
      partsAng_.reserve( nparts );
      const double* dptr = dvals;
      for (unsigned int idx = 0; idx != nparts; idx++, dptr += 6) {
        partsDist_.push_back( Stats<double>(dptr[0], dptr[1], dptr[2]) );
        partsAng_.push_back( Stats<double>(dptr[3], dptr[4], dptr[5]) );
      }
    }
    /// Update parts with given double array containing N, mean, and M2 for each part
    void CombineParts(unsigned int nparts, const double* dvals) {
      const double* dptr = dvals;
      for (unsigned int idx = 0; idx != nparts; idx++, dptr += 6) {
        partsDist_[idx].Combine( Stats<double>(dptr[0], dptr[1], dptr[2]) );
        partsAng_[idx].Combine( Stats<double>(dptr[3], dptr[4], dptr[5]) );
      }
    }
#   endif /* MPI */
    /// First sort by frames (descending), then distance (ascending).
    bool operator<(const Hbond& rhs) const {
      if (frames_ == rhs.frames_)
        return (dist_ < rhs.dist_);
      else
        return (frames_ > rhs.frames_);
    }
    /// Update distance/angle/time series
    void Update(double, double, int, Iarray const&, int);
    /// Calculate averages
    void CalcAvg();
  private:
    double dist_;  ///< Used to calculate average distance of hydrogen bond
    double angle_; ///< Used to calculate average angle of hydrogen bond
    DataSet_integer* data_; ///< Hold time series data
    int A_; ///< Acceptor atom index
    int H_; ///< Hydrogen atom index
    int D_; ///< Donor atom index
    int frames_; ///< # frames this hydrogen bond has been present
    std::vector<Stats<double> > partsDist_; ///< Hold avg. distance, split by parts
    std::vector<Stats<double> > partsAng_;  ///< Hold avg. angle, split by parts
};
}
}
#endif
