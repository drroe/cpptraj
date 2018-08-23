#ifndef INC_PROGRESSTIMERBAR_H
#define INC_PROGRESSTIMERBAR_H
#include "Timer.h"
/// Print progress/elapsed time to screen.
/** Keeps track of progress via a progress bar (printing percentage). After a
  * certain time interval, print estimate of time remaining as well.
  * Assumes that iterations start at 0.
  */
class ProgressTimerBar {
  public:
    ProgressTimerBar() {}
    /// Prepare progress bar with # iterations, time interval, and pct interval
    void SetupProgress(int, double, double);
    /// Update progress
    void Update(int);
    /// Finish
    void Finish() const;
  private:
    Timer time_;         ///< Keep track of elapsed time.
    double tgtTime_;     ///< Next target time to print estimate time remaining.
    double tInterval_;   ///< Interval between updates in seconds.
    double pctInterval_; ///< Progress bar percentage interval.
    double tgtPct_;      ///< Next target percentage.
    int total_it_;       ///< Total number of iterations process should take.
    int tgt_it_;         ///< Next target iteration to print progress bar.
    int iInterval_;      ///< Progress interval when total not known
    /// total_it_ will be set to this when total # frames not known.
    static const int UNKNOWN_FRAMES_;
    /// Character to print when total # iterations unknown.
    static const char PROGCHAR_;
};
#endif
