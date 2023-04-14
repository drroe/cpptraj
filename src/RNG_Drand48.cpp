#include "RNG_Drand48.h"
#include <cstdlib>

Cpptraj::RNG_Drand48::RNG_Drand48() {}

int Cpptraj::RNG_Drand48::SetupRng() {
  srand48( Seed() );
  return 0;
}

/** Generate a random number between 0 and RAND_MAX. */
unsigned int Cpptraj::RNG_Drand48::Number() {
  return (int)(drand48() * (double)RAND_MAX);
}

/** Generate a random number between 0 and 1. 
  */
double Cpptraj::RNG_Drand48::Generate() {
  return drand48();
}
