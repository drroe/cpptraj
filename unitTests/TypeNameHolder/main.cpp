// Unit test for TypeNameHolder class
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include "TypeNameHolder.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

/*struct reverseSort {
  bool operator() (const NameType& lhs, const NameType& rhs) const {
    return ( lhs > rhs );
  }
} reverseSortObj;*/

int main() {
  TypeNameHolder type1;
  type1.AddName("CT");
  TypeNameHolder type2("O");
  bool ismatch = type1.Match_NoWC( type2 );
  if (ismatch) return Err("TypeNameHolder single type no match failed.");
  TypeNameHolder type3("CT");
  ismatch = type1.Match_NoWC( type3 );
  if (!ismatch) return Err("TypeNameHolder single type match failed.");

  TypeNameHolder bp1(2);
  bp1.AddName("CT");
  bp1.AddName("O");
  TypeNameHolder bp2(2);
  bp2.AddName("O");
  bp2.AddName("CT");
  ismatch = bp1.Match_NoWC( bp2 );
  if (!ismatch) return Err("TypeNameHolder two type reverse match failed.");

  TypeNameHolder ap1(3);
  ap1.AddName("CT");
  ap1.AddName("O");
  ap1.AddName("HO");
  TypeNameHolder ap2(3);
  ap2.AddName("HO");
  ap2.AddName("O");
  ap2.AddName("CT");
  ismatch = ap2.Match_NoWC( ap1 );
  if (!ismatch) return Err("TypeNameHolder three type reverse match failed.");

  TypeNameHolder dp1(4);
  dp1.AddName("X");
  dp1.AddName("CX");
  dp1.AddName("O");
  dp1.AddName("X");
  TypeNameHolder dp2(4);
  dp2.AddName("CT");
  dp2.AddName("CX");
  dp2.AddName("O");
  dp2.AddName("HO");
  ismatch = dp1.Match_NoWC( dp2 );
  if (ismatch) return Err("TypeNameHolder four type no WC match failed.");
  ismatch = dp1.Match_WC( dp2, "X" );
  if (!ismatch) return Err("TypeNameHolder four type WC match failed.");

  TypeNameHolder ip1(4);
  ip1.AddName("O");
  ip1.AddName("HO");
  ip1.AddName("CT");
  ip1.AddName("N");
  ip1.SortImproperByAlpha("X");
  if (ip1[0] != "HO") return Err("TypeNameHolder improper alpha sort failed (pos 0)");
  if (ip1[1] != "N") return Err("TypeNameHolder improper alpha sort failed (pos 1)");
  if (ip1[2] != "CT") return Err("TypeNameHolder improper alpha sort failed (pos 2)");
  if (ip1[3] != "O") return Err("TypeNameHolder improper alpha sort failed (pos 3)");

  TypeNameHolder ip2(4);
  ip2.AddName("O");
  ip2.AddName("X");
  ip2.AddName("CT");
  ip2.AddName("X");
  ip2.SortImproperByAlpha("X");
  if (ip2[0] != "X") return Err("TypeNameHolder improper2 alpha sort failed (pos 0)");
  if (ip2[1] != "X") return Err("TypeNameHolder improper2 alpha sort failed (pos 1)");
  if (ip2[2] != "CT") return Err("TypeNameHolder improper2 alpha sort failed (pos 2)");
  if (ip2[3] != "O") return Err("TypeNameHolder improper2 alpha sort failed (pos 3)");

  //printf("ip1 %s %s %s %s\n", *ip1[0], *ip1[1], *ip1[2], *ip1[3]);
  //printf("ip2 %s %s %s %s\n", *ip2[0], *ip2[1], *ip2[2], *ip2[3]);

  ismatch = ip2.Match_NoWC( ip1 );
  if (ismatch) return Err("TypeNameHolder improper no WC match failed.");
  ismatch = ip2.Match_WC( ip1, "X" );
  if (!ismatch) return Err("TypeNameHolder improper WC match failed.");

  return 0;
}
