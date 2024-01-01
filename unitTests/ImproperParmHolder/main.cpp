// Unit test for ImproperParmHolder class
#include <cstdio>
#include "ParameterHolders.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
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

  TypeNameHolder ip1a(4);
  ip1a.AddName("N");
  ip1a.AddName("O");
  ip1a.AddName("CT");
  ip1a.AddName("HO");

  ImproperParmHolder IP0;
  ParameterHolders::RetType ret = IP0.AddParm( ip1, DihedralParmType( 3.0, 1.0, 0.0 ), false );
  if (ret == ParameterHolders::ERR) return Err("Could not add improper parameter");
  bool found;
  DihedralParmArray impropers = IP0.FindParam( ip1a, found );
  if (!found) return Err("Could not find improper parameter (no wildcards).");

  TypeNameHolder ip2(4);
  ip2.AddName("X");
  ip2.AddName("X");
  ip2.AddName("CT");
  ip2.AddName("O");

  ImproperParmHolder IP;
  ret = IP.AddParm( ip2, DihedralParmType( 2.0, 1.0, 3.14159/2.0 ), false );
  if (ret == ParameterHolders::ERR) return Err("Could not add improper parameter");
  //ip2.SortImproperByAlpha("X");
  
  impropers = IP.FindParam( ip1, found);
  if (found) return Err("Improper parameter search before wildcard added failed.");

  IP.SetWildcard('X');
  impropers = IP.FindParam( ip1, found );
  if (!found) return Err("Improper parameter search with wildcard match failed.");

  impropers = IP.FindParam( ip0, found );
  if (found) return Err("Improper parameter search found something when it should not have.");

  return 0;
}
