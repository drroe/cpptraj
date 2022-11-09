#include <cstdio>
#include <algorithm>
#include "Param/HookesLawType.h"
#include "Param/DihedralParmType.h"
#include "Param/HB_ParmType.h"
#include "Param/NonbondType.h"
#include "Param/LJparmType.h"
#include "Param/CmapGridType.h"
#include "Param/NonbondParmType.h"
#include "Param/LES_AtomType.h"
#include "Param/LES_ParmType.h"
#include "Param/CapParmType.h"
#include "Top/BondType.h"
#include "Top/AngleType.h"
#include "Top/DihedralType.h"
#include "Top/CmapType.h"
using namespace Cpptraj::Param;
using namespace Cpptraj::Top;

void PrintBonds(BondArray const& bonds) {
  printf("%zu bonds.\n", bonds.size());
  for (BondArray::const_iterator it = bonds.begin(); it != bonds.end(); ++it)
    printf("\t%8i -- %8i (%i)\n", it->A1(), it->A2(), it->Idx());
}

void TestBonds() {
  BondArray bonds;
  bonds.push_back( BondType(1, 2) );
  BondType b1;
  b1 = BondType(0, 2, 1);
  bonds.push_back( b1 );
  BondType b0( 0, 1, 0 );
  bonds.push_back( b0 );
  PrintBonds(bonds);
  std::sort(bonds.begin(), bonds.end());
  PrintBonds(bonds);
}

int main (int argc, char** argv) {

  HookesLawType BP;
  DihedralParmType DP;
  HB_ParmType HP;
  NonbondType NB;
  LJparmType LJ;
  CmapGridType cmgrid;
  NonbondParmType NBparam;
  LES_AtomType lat;
  LES_ParmType lesparm;
  CapParmType capparm;

  BondType bond;
  AngleType angle;
  DihedralType dih;
  CmapType cmap;

  TestBonds();

  return 0;
}
