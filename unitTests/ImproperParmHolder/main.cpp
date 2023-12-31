// Unit test for ImproperParmHolder class
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include "ParameterHolders.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
  ImproperParmHolder IP;
  IP.SetWildcard('X');

  TypeNameHolder ip0(4);
  ip0.AddName("N*");
  ip0.AddName("CX");
  ip0.AddName("CT");
  ip0.AddName("HN");

  TypeNameHolder ip1(4);
  ip1.AddName("O");
  ip1.AddName("HO");
  ip1.AddName("CT");
  ip1.AddName("N");
  //ip1.SortImproperByAlpha("X");

  TypeNameHolder ip2(4);
  ip2.AddName("X");
  ip2.AddName("X");
  ip2.AddName("CT");
  ip2.AddName("O");
  ParameterHolders::RetType ret = IP.AddParm( ip2, DihedralParmType( 2.0, 1.0, 3.14159/2.0 ), false );
  if (ret == ParameterHolders::ERR) return Err("Could not add improper parameter");
  //ip2.SortImproperByAlpha("X");
  
  //printf("ip1 %s %s %s %s\n", *ip1[0], *ip1[1], *ip1[2], *ip1[3]);
  //printf("ip2 %s %s %s %s\n", *ip2[0], *ip2[1], *ip2[2], *ip2[3]);

  bool found;
  DihedralParmArray impropers = IP.FindParam( ip1, found );
  if (!found) {
    return Err("Improper parameter search with wildcard match failed.");
  }
  impropers = IP.FindParam( ip0, found );
  if (found) {
    return Err("Improper parameter search found something when it should not have.");
  }

  return 0;
}
