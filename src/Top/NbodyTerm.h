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
  private:
    int idx_; ///< Index into corresponding parameter set.
};
}
}
#endif
