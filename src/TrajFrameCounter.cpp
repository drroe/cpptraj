#include "TrajFrameCounter.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

#include "Counter_Regular.h"

using namespace Cpptraj;

/** CONSTRUCTOR */
TrajFrameCounter::TrajFrameCounter() :
  counter_(0),
  total_frames_(-1),
  total_read_frames_(-1)
{}

/** DESTRUCTOR */
TrajFrameCounter::~TrajFrameCounter() {
  if (counter_ != 0)
    delete counter_;
}

/** Parse argument list for trajectory-related frame args. */
int TrajFrameCounter::CheckFrameArgs(int nframes, ArgList &argIn) {
  total_frames_ = nframes;
  if (total_frames_==0) {
    mprinterr("Error: trajectory contains no frames.\n");
    return 1;
  }

  if (startStopOffset(argIn)) {
    mprinterr("Error: Could not process trajectory start/stop/offset arguments.\n");
    return 1;
  }
  // Sanity check
  if (counter_ == 0) {
    mprinterr("Internal error: Counter was not allocated.\n");
    return 1;
  }
  return 0;
}

/** Parse regular start/stop/offset frame args. Frame args start at
  * 1, internal frame #s start at 0. So for a traj with 10 frames:
  * - Internal #: 0 1 2 3 4 5 6 7 8 9
  * - Frame Arg#: 1 2 3 4 5 6 7 8 9 10
  * - Defaults: start=1, stop=-1, offset=1
  */

int TrajFrameCounter::startStopOffset(ArgList& argIn) {
  // ----- Regular start/stop/offset -------------
  int start = 0;
  int stop = -1;
  int offset = 1;

  if (argIn.hasKey("lastframe")) {
    // lastframe is a special case where only the last frame will be selected
    if (total_frames_ > 0) {
      start = total_frames_;
      stop = total_frames_;
      offset = 1;
    } else {
      mprinterr("Error: lastframe specified but # frames could not be determined.\n");
      return 1;
    }
  } else {
    start = argIn.getNextInteger(1);
    // Last explicitly selects final frame as stop arg.
    if (argIn.hasKey("last"))
      stop = -1;
    else
      stop = argIn.getNextInteger(-1);
    offset = argIn.getNextInteger(1);
  }
  // Check that start argument is valid.
  if (start != 1) {
    if (start == 0) {
      mprintf("Warning: start argument is 0, setting to 1.\n", start);
      start = 1; //start = 0;
    } else if (start < 0) {
      // Negative start means we want that many frames before stop
      if (stop == -1) {
        if (total_frames_ >=0)
          stop = total_frames_;
        else {
          mprinterr("Error: For start < 0, stop argument must be specified when # frames unknown.\n");
          return 1;
        }
      }
      mprintf("\tStarting %i frames before frame %i\n", -start, stop);
      start = stop + start;
      if (start < 1) {
        mprintf("Warning: would start before frame 1, setting start to 1.\n");
        start = 1;
      }
    } else if (total_frames_ >= 0 && start > total_frames_) {
      // start==stop and greater than # frames, archaic 'lastframe'.
      if (start == stop) {
        mprintf("Warning: start %i > #Frames (%i), setting to last frame.\n",
                start, total_frames_);
        start = total_frames_; //start = total_frames_ - 1;
      } else {
        mprinterr("Error: start %i > #Frames (%i), no frames will be processed.\n",
                  start, total_frames_);
        //start=start - 1;
        return 1;
      }
    }
  }
  start--; // Internal frame nums start from 0.
  // Check that stop argument is valid
  if (stop != -1) {
    if ( (stop - 1) < start) { // Internal frame nums start from 0.
      mprinterr("Error: stop %i < start, no frames will be processed.\n", stop);
      //stop = start;
      return 1;
    } else if (total_frames_ >= 0 && stop > total_frames_) {
      mprintf("Warning: stop %i > #Frames (%i), setting to max.\n", stop, total_frames_);
      stop = total_frames_;
    } 
  } else if (total_frames_ >= 0) // -1 means use last frame
    stop = total_frames_;
  // Check that offset argument is valid.
  if (offset != 1) {
    if (offset < 1) {
      mprintf("Warning: offset %i < 1, setting to 1.\n", offset);
      offset = 1;
    } else if (stop != -1 && offset >= (stop - start)) {
      mprintf("Warning: offset %i is so large that only 1 set will be processed.\n",
              offset);
    }
  }
  //mprintf("DEBUG SetArgs: Start %i Stop %i  Offset %i\n", start, stop, offset);
  // Calculate actual number of frames that will be read based on start,
  // stop, and offset.
  total_read_frames_ = -1;
  if (stop != -1) {
    int Nframes = stop - start;
    total_read_frames_ = Nframes / offset;
    // Round up
    if ( (Nframes % offset) > 0 )
      ++total_read_frames_;
    if (total_read_frames_ == 0) {
      mprinterr("Error: No frames will be read based on start, stop, "
                "and offset values (%i, %i, %i)\n", start+1, stop, offset);
      return 1;
    }
  }

  counter_ = new Counter_Regular(start, stop, offset);
  if (counter_ == 0) {
    mprinterr("Internal Error: Could not allocate regular counter.\n");
    return 1;
  }
  return 0;
}

