#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include "../../../src/FileName.h"

typedef std::vector<std::string> Sarray;
typedef std::pair<std::string, Sarray> ResFilesPair;
typedef std::map<std::string, Sarray> ResFilesMap;

static inline int FileErr() {
  fprintf(stderr,"Error reading file.\n");
  return 1;
}

/** Used to create the ResToSite.dat file in dat/TitratableSites */
int main() {
  char buffer[1024];
  int bufsize = 1023;
  static const std::string EXT = ".st";
  FILE* outfile = stdout;

  ResFilesMap resnameToFiles;

  File::NameArray stfiles = File::ExpandToFilenames("*" + EXT);

  for (File::NameArray::const_iterator file = stfiles.begin(); file != stfiles.end(); ++file)
  {
    //printf("\t%s\n", file->full());
    // Open the file to get the residue name
    FILE* infile = fopen( file->full(), "rb" );
    if (infile == 0) {
      fprintf(stderr,"Error opening %s\n", file->full());
      return 1;
    }
    if (fgets( buffer, bufsize, infile) == 0) return FileErr(); // pka
    if (fgets( buffer, bufsize, infile) == 0) return FileErr(); // first line
    char resnamebuf[32];
    sscanf(buffer, "%s", resnamebuf);
    //printf("\t%s\n", resnamebuf);
    fclose(infile);
    std::string resname(resnamebuf);
    ResFilesMap::iterator it = resnameToFiles.lower_bound( resname );
    if (it == resnameToFiles.end() || it->first != resname) {
      // New res
      resnameToFiles.insert( it, ResFilesPair(resname, Sarray(1, file->Full())) );
    } else {
      // Existing res
      it->second.push_back( file->Full() );
    }
  }

  fprintf(outfile, "#<resname list> <Filename list>\n");
  for (ResFilesMap::const_iterator it = resnameToFiles.begin(); it != resnameToFiles.end(); ++it)
  {
    fprintf(outfile, "%s ", it->first.c_str());
    for (Sarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
      if (jt == it->second.begin())
        fprintf(outfile, "%s", jt->c_str());
      else
        fprintf(outfile, ",%s", jt->c_str());
    }
    fprintf(outfile,"\n");
  }
  return 0;
}
