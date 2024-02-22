// Unit test for TextFormat class
#include <cstdio>
#include "TextFormat.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
  TextFormat f0;
  printf("Default format: %s\n", f0.fmt());
  if (f0.Fmt() != "%8.3f") return Err("Default format string is not %%8.3f\n");

  return 0;
}
