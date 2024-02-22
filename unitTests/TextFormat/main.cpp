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

  TextFormat f1(TextFormat::SCIENTIFIC, 16, 8, 5);
  f1.SetFormatAlign(TextFormat::NO_SPACES);
  f1.SetLongFormat( true );
  printf("5FE16.8 format: '%s' width=%i\n", f1.fmt(), f1.ColumnWidth());

  return 0;
}
