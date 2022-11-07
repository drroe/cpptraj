#ifndef INC_PARAM_NONBONDTYPE_H
#define INC_PARAM_NONBONDTYPE_H
#include <cmath> // pow, fabs
#include <vector>
namespace Cpptraj {
namespace Param {
/// Hold Lennard-Jones 6-12 interaction A and B parameters
class NonbondType {
    /** Tolerance for comparison. A little larger than SMALL because A
      * and B tend to be large.
      */
    // NOTE: Probably should check __cpluscplus here instead of using a
    //       define, but this is guaranteed to be portable.
#     define tol_ 0.00000001
      //static const double tol_ = 0.00000001;
  public:
    NonbondType() : A_(0), B_(0) {}
    NonbondType(double a, double b) : A_(a), B_(b) {}
    double A() const { return A_; }
    double B() const { return B_; }
    void SetA(double a) { A_ = a; }
    void SetB(double b) { B_ = b; }
    double Radius() const {
      if (B_ > 0.0)
        return (0.5 * pow(2.0 * A_ / B_, (1.0/6.0)));
      else
        return 0.0;
    }
    double Depth() const {
      if (A_ > 0.0)
        return ( (B_ * B_) / (4.0 * A_) );
      else
        return 0.0;
    }
    /// \return True if A and B match
    bool operator==(NonbondType const& rhs) const {
      return ( (fabs(A_ - rhs.A_) < tol_) &&
               (fabs(B_ - rhs.B_) < tol_) );
    }
    /// \return True if A and B do not match
    bool operator!=(NonbondType const& rhs) const {
      return ( (fabs(A_ - rhs.A_) > tol_) ||
               (fabs(B_ - rhs.B_) > tol_) );
    }
    /// \return True if A less than zero, or B if A is equal.
    bool operator<(NonbondType const& rhs) const {
      if (*this != rhs) {
        if ( (fabs(A_ - rhs.A_) < tol_) )
          return (B_ < rhs.B_);
        else
          return (A_ < rhs.A_);
      } else
        return false;
    }
  private:
    double A_;
    double B_;
#   undef tol_
};
typedef std::vector<NonbondType> NonbondArray;

}
}
#endif
