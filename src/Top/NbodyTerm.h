#ifndef INC_TOP_NBODYTERM_H
#define INC_TOP_NBODYTERM_H
namespace Cpptraj {
namespace Top {
/// Abstract base class for an N-body (connected via bonds) term.
class NbodyTerm {
  public:
    /// CONSTRUCTOR
    NbodyTerm() : idx_(-1) {}
    /// CONSTRUCTOR - parameter index 
    NbodyTerm(int idx) : idx_(idx) {}
    /// DESTRUCTOR - virtual since inherited
    virtual ~NbodyTerm() {}
    /// COPY CONSTRUCTOR
    NbodyTerm(NbodyTerm const& rhs) : idx_(rhs.idx_) {}
    /// ASSIGNMENT
    NbodyTerm& operator=(NbodyTerm const& rhs) {
      if (this == &rhs) return *this;
      idx_ = rhs.idx_;
      return *this;
    }

    /// \return Index into parameter set for this term.
    int Idx() const { return idx_; }
    /// Set index into parameter set for this term.
    void SetIdx(int idx) { idx_ = idx; }
    /// \return Atom # at specified position
    virtual int At(int idx) const = 0;
    /// \return # atoms in this term
    virtual int Nat() const = 0;
    /// \return true if rhs is less than this
    bool operator<(NbodyTerm const& rhs) {
      if (Nat() == rhs.Nat()) {
        for (int i = 0; i < Nat(); i++) {
          if (At(i) != rhs.At(i)) {
            return (At(i) < rhs.At(i));
          }
        }
      }
      return (Nat() < rhs.Nat());
    }
  private:
    int idx_; ///< Index into corresponding parameter set.
};
}
}
#endif
