#include <cstdio>
#include <cmath>
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

  private:
    TrajFrameCounter counter_;
};

// Fake Trajin class for testing
class Fake {
  public:
    Fake() {}

    Traj0 const& Traj() const { return traj_; }
    int FakeSetup(int nframes, ArgList& argIn) {
      return traj_.Set(nframes, argIn);
    }
      
  private:
    Traj0 traj_;
};


int main() {
  ArgList args("5 23 3");

  Fake fakeTraj;
  if (fakeTraj.FakeSetup(100, args))
    return Err("Could not set up fake trajectory.");


  return 0;
}
