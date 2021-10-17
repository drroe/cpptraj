#ifndef INC_TRAJFRAMECOUNTER_H
#define INC_TRAJFRAMECOUNTER_H
class ArgList;
namespace Cpptraj {
class Counter;
/// Used to keep track of input frames in trajectory
class TrajFrameCounter {
  public:
    /// CONSTRUCTOR
    TrajFrameCounter();
    /// DESTRUCTOR
    ~TrajFrameCounter();
    /// Set up from arguments
    int CheckFrameArgs(int, ArgList&);
  private:
    /// Regular start/stop/offset
    int startStopOffset(ArgList&);

    Counter* counter_;      ///< Used to keep track of frames.
    int total_frames_;      ///< Total number of frames in the trajectory.
    int total_read_frames_; ///< Total number of frames that will be read.
};
} // END namespace cpptraj
#endif
