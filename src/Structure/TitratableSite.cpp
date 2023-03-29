#include "TitratableSite.h"
#include "../NameType.h"
#include "../CpptrajStdio.h"
#include "../BufferedLine.h"
#include <cstdlib> // atof

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
TitratableSite::TitratableSite() :
  pKa_(0)
{}

/** Clear all data. */
void TitratableSite::Clear() {
  pKa_ = 0;
  resName_.clear();
  nameToCharges_.clear();
}

/** Load site titration data from a file. */
int TitratableSite::LoadSiteData(std::string const& fname)
{
  if (fname.empty()) {
    mprinterr("Internal Error: TitratableSite::LoadSiteData(): No file name given.\n");
    return 1;
  }
  // Clear existing data
  Clear();
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open titratable site file '%s'\n", fname.c_str());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read pKa from site file '%s'\n", fname.c_str());
    return 1;
  }
  pKa_ = atof( ptr ); // TODO check only 1 token?
  ptr = infile.Line();
  while (ptr != 0) {
    // Expected format: <resname> <atomname> <Qstate1> <Qstate2>
    int ntokens = infile.TokenizeLine(" ");
    if (ntokens != 4) {
      mprinterr("Error: Malformed line in '%s' : %s\n", fname.c_str(), ptr);
      mprinterr("Error: Expected 4 tokens, got %i\n", ntokens);
      return 1;
    }
    const char* rname = infile.NextToken();
    if (resName_.empty())
      resName_.assign( rname );
    NameType aname( infile.NextToken() );
    double q1 = atof( infile.NextToken() );
    double q2 = atof( infile.NextToken() );

    MapType::iterator it = nameToCharges_.lower_bound( aname );
    if (it == nameToCharges_.end() || it->first != aname ) {
      nameToCharges_.insert( it, PairType( aname, ChargePair(q1, q2) ) );
    } else {
      mprinterr("Error: Atom '%s' is duplicated in site file '%s'\n", *aname, fname.c_str());
      return 1;
    }
    ptr = infile.Line();
  }
  infile.CloseFile();
  return 0;
}
