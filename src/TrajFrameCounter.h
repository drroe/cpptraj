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

    // NOTE: Even though Begin() and UpdateCounters() _can_ be const since
    //       they are only accessing a pointer, make them non-const to
    //       reflect the fact that they are modifying counter_.
    /// Prepare counter for use
    void Begin()          { counter_->StartCounter(); }
    /// Update internal counter
    void UpdateCounters() { counter_->UpdateCounter(); }

    /// Check if the counter is finished
    bool CheckFinished()      const { return counter_->IsFinished(); }
    /// /return Current frame number
    int Current()             const { return counter_->CurrentNumber(); }
    /// \return Previous frame number (before UpdateCounters() was called)
    int PreviousFrameNumber() const { return counter_->PreviousNumber(); }
    /// \return How many times UpdateCounters has been called
    int NumFramesProcessed()  const { return (int)counter_->CurrentIdx(); }
    /// \return Total number of frames that will be read by the counter.
    int TotalReadFrames()     const { return counter_->CounterTotal(); }
    /// \return Frame number corresponding to given index
    int IdxToFrame(int idx)   const { return counter_->NumberAtIdx(idx); }
  private:
    /// Regular start/stop/offset
    int startStopOffset(ArgList&);

    Counter* counter_;      ///< Used to keep track of frames.
    int total_frames_;      ///< Total number of frames in the trajectory.
    int total_read_frames_; ///< Total number of frames that will be read.
};
} // END namespace cpptraj
#endif
