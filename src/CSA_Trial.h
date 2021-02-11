#ifndef INC_CSA_TRIAL_H
#define INC_CSA_TRIAL_H
namespace Cpptraj {
namespace CSA {

/** A CSA trial is a potential structure or path that is being optimized. */
class Trial {
  public:
    Trial() {}
    // Virtual since inherited
    virtual ~Trial() {}

    double Score() const { return score_; }
  private:
    double score_; ///< Current score of trial
};

}
}
#endif
