#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "ArgList.h"
#include "TrajFrameCounter.h"
#include "TrajFrameIndex.h"

using namespace Cpptraj;

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

//const double SMALL        = 0.00000000000001;

//bool f_not_equals(double d1, double d2) {
//  return (fabs(d1 - d2) > SMALL);
//}

// Fake InputTrajCommon class for testing
class Traj0 {
  public:
    Traj0() {}

    int Set(int nframes, ArgList& argIn) {
      return counter_.CheckFrameArgs(nframes, argIn);
    }

    TrajFrameCounter const& Counter() const { return counter_; }
    TrajFrameCounter& SetCounter() { return counter_; }

  private:
    TrajFrameCounter counter_;
};

// Fake Trajin class for testing
class Fake {
  public:
    Fake() {}

    Traj0 const& Traj() const { return traj_; }
    int FakeSetup(std::string const& nameIn, int nframes, ArgList const& argIn) {
      name_ = nameIn;
      ArgList tmp = argIn;
      return traj_.Set(nframes, tmp);
    }

    void Print() const {
      traj_.Counter().PrintInfoLine(name_.c_str());
    }

    int BeginTraj() {
      traj_.SetCounter().Begin();
      return 0;
    }

    int GetNextFrame(int& frameIn) {
      if (traj_.Counter().CheckFinished()) return 0;
      frameIn = traj_.Counter().Current();
      //if (ReadTrajFrame( traj_.Counter().Current(), frameIn )) return 0;
      traj_.SetCounter().UpdateCounters();
      return 1;
    }

  private:
    Traj0 traj_;
    std::string name_;
};

typedef std::vector<Fake*> ListType;

int main() {
  ListType trajectories;

  trajectories.push_back( new Fake() );
  if (trajectories.back()->FakeSetup("trajectory1", 100, ArgList("5 22 3")))
    return Err("Could not set up fake trajectory 1.");
  trajectories.push_back( new Fake() );
  if (trajectories.back()->FakeSetup("trajectory2", 20, ArgList("7 14 2")))
    return Err("Could not set up fake trajectory 2.");
  trajectories.push_back( new Fake() );
  if (trajectories.back()->FakeSetup("trajectory3", 46, ArgList("42 46")))
    return Err("Could not set up fake trajectory 3.");

  int idx = 0;
  for (ListType::const_iterator it = trajectories.begin(); it != trajectories.end(); ++it)
  {
    (*it)->Print();
    (*it)->BeginTraj();
    int fnum;
    while ( (*it)->GetNextFrame(fnum) ) {
      printf("[%8i] %8i\n", idx++, fnum+1);
    }
  }

  TrajFrameIndex<ListType> TFI;
  if (TFI.MaxFrames( trajectories ) != 15)
    return Err("Number of frames covered by TrajFrameIndex is not 15.");

  if (TFI.FindIndex(0, trajectories) != 4) return Err("FindIndex(0) failed.");
  if (TFI.CurrentTrajNum() != 0) return Err("Current traj is not 0.");
  if (TFI.FindIndex(4, trajectories) != 16) return Err("FindIndex(4) failed.");
  if (TFI.CurrentTrajNum() != 0) return Err("Current traj is not 0.");

  if (TFI.FindIndex(7, trajectories) != 8) return Err("FindIndex(7) failed.");
  if (TFI.CurrentTrajNum() != 1) return Err("Current traj is not 1.");
  if (!TFI.TrajHasChanged()) return Err("Traj change detection failed.");
  if (TFI.FindIndex(9, trajectories) != 12) return Err("FindIndex(9) failed.");
  if (TFI.CurrentTrajNum() != 1) return Err("Current traj is not 1.");
  if (TFI.TrajHasChanged()) return Err("No Traj change detection failed.");

  if (TFI.FindIndex(11, trajectories) != 42) return Err("FindIndex(11) failed.");
  if (TFI.CurrentTrajNum() != 2) return Err("Current traj is not 2.");
  if (TFI.FindIndex(14, trajectories) != 45) return Err("FindIndex(14) failed.");
  if (TFI.CurrentTrajNum() != 2) return Err("Current traj is not 2.");

  for (ListType::iterator it = trajectories.begin(); it != trajectories.end(); ++it)
    delete *it;

  return 0;
}
