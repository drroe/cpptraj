#ifndef INC_TRAJFRAMECOUNTER_H
#define INC_TRAJFRAMECOUNTER_H
#include "Counter.h"
class ArgList;
namespace Cpptraj {
/// Used to keep track of input frames in trajectory
class TrajFrameCounter {
  public:
    /// CONSTRUCTOR
    TrajFrameCounter();
    /// DESTRUCTOR
    ~TrajFrameCounter();
    /// Set up from arguments
    int CheckFrameArgs(int, ArgList&);
    /// Print counter info to stdout
    void PrintInfoLine(const char*) const;

    /// Check if the counter is finished
    bool CheckFinished() const { return counter_->IsFinished(); }
    /// /return Current frame number
    int Current() const { return counter_->CurrentNumber(); }
    /// Update internal counter
    void UpdateCounters() { counter_->UpdateCounter(); }
  private:
    /// Regular start/stop/offset
    int startStopOffset(ArgList&);

    Counter* counter_;      ///< Used to keep track of frames.
    int total_frames_;      ///< Total number of frames in the trajectory.
    int total_read_frames_; ///< Total number of frames that will be read.
};
} // END namespace cpptraj
#endif
