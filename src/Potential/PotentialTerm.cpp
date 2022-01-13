#include "PotentialTerm.h"

const char* PotentialTerm::typeStr_[] = {
  "Bond",
  "Angle",
  "Dihedral",
  "SimpleNonbond",
  "OpenMM",
  "Replicate",
  0
};

const char* PotentialTerm::TypeStr(Type typeIn) {
  return typeStr_[typeIn];
}
