// Unit test for NameType class
#include <cstdio>
#include "Frame.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

static void printAtoms(Frame const& frm1, const char* msg) {
  frm1.Info(msg);
  for (int iat = 0; iat < frm1.Natom(); iat++) {
    frm1.printAtomCoord(iat);
  }
}

int main() {
  Frame frm1;

  frm1.SetupFrame( 3 );
  frm1.ClearAtoms();
  for (int iat = 1; iat <= 3; iat++)
    frm1.AddVec3( Vec3((double)iat, (double)iat, (double)iat) );
  printAtoms(frm1, "frm1");

  Frame frm2 = frm1;

  frm1.AppendFrame( frm2 );
  printAtoms(frm1, "frm1 after append");

  return 0;
}
