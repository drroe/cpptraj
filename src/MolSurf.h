#ifndef INC_MOLSURF_H
#define INC_MOLSURF_H
#include <vector>
#include "Frame.h"
class MolSurf {
  public:
    MolSurf();

    double CalcSurface(Frame const&);
  private:
    void UpdatePositions(Frame const&);
    int GetNeighbors();

    /// Hold neighbor index/distance pair
    class NbrIdx {
      public:
        NbrIdx() : idx_(-1), dist_(0.0) {}
        NbrIdx(int i, double d) : idx_(i), dist_(d) {}
        bool operator<(NbrIdx const& rhs) const { return dist_ < rhs.dist_; }
        int Idx() const { return idx_; }
        double Dist() const { return dist_; }
      private:
        int idx_;     ///< Neighbor atom index
        double dist_; ///< Neighbor atom distance
    };

    class atom {
      public:
        atom() {}
        double Radius() const { return radius_; }
        bool IsBuried() const { return buried_; }
        Vec3 XYZ()      const { return xyz_;    }

        void SetBuried(bool b) { buried_ = b; }
        void AddNeighbor(NbrIdx const& n) { neighbors_.push_back( n ); }
        void SortNeighbors();
        void UpdatePosition(const double* xyzIn) {
          xyz_[0] = xyzIn[0];
          xyz_[1] = xyzIn[1];
          xyz_[2] = xyzIn[2];
        }
      private:
        Vec3 xyz_;      ///< Atom position
        double radius_; ///< Atom radius
        bool buried_;   ///< True if atom is buried

        std::vector<NbrIdx> neighbors_; ///< Indices of neighbor atoms
    };

    typedef std::vector<atom> AtomArray;
    AtomArray atoms_; ///< Hold info for all selected atoms
    double probe_radius_;
    AtomMask mask_; ///< Selected atoms from input frame
};
#endif
