#include "MeadOpts.h"
#include "../CpptrajStdio.h"
#include "../../mead/globals.h"

Cpptraj::Mead::MeadOpts::MeadOpts() :
  epsin_(4),
  epsext_(80),
  epsvac_(1),
  solrad_(1.4),
  sterln_(2.0),
  ionicstr_(0),
  temp_(300)
{}

/** Set verbosity of underlying MEAD library. */
void Cpptraj::Mead::MeadOpts::MeadVerbosity(int i) {
  blab1pt = &cnull;
  blab2pt = &cnull;
  blab3pt = &cnull;
  if (i >= 1)
    blab1pt = &std::cout;
  if (i >= 2)
    blab2pt = &std::cout;
  if (i >= 3)
    blab3pt = &std::cout;
  if (i != 0)
    mprintf("Info: MEAD verbosity set to %i\n", i);
}

