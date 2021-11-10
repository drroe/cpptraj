#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "ArgList.h"
//#include "Counter_Regular.h"
//#include "Counter_Array.h"
#include "TrajFrameCounter.h"

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

  for (ListType::iterator it = trajectories.begin(); it != trajectories.end(); ++it)
    delete *it;

  return 0;
}
