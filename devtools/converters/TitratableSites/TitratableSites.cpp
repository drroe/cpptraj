#include <cstdio>
#include <string>
#include "../../../src/FileName.h"

/** Used to create the ResToSite.dat file in dat/TitratableSites */

int main() {
  static const std::string EXT = ".st";

  File::NameArray stfiles = File::ExpandToFilenames("*" + EXT);

  for (File::NameArray::const_iterator file = stfiles.begin(); file != stfiles.end(); ++file)
  {
    printf("\t%s\n", file->full());
  }
}
