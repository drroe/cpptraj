#include "DataSet_LeapOpts.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_LeapOpts::DataSet_LeapOpts() :
  // 0 dim indicates DataSet-specific write
  DataSet(LEAPOPTS, GENERIC, TextFormat(TextFormat::STRING, 12, 0), 0),
  pbradii_(Cpptraj::Parm::MBONDI)
{

}

/** Set default GB radii from keyword. */
int DataSet_LeapOpts::SetGbRadii(std::string const& keyword) {
  Cpptraj::Parm::GB_RadiiType radType = Cpptraj::Parm::GbTypeFromKey( keyword );
  if (radType == Cpptraj::Parm::UNKNOWN_GB) {
    mprinterr("Error: Unrecognized GB radii type: %s\n", keyword.c_str());
    return 1;
  }
  pbradii_ = radType;
  mprintf("\tSet default GB radii: %s\n", Cpptraj::Parm::GbTypeStr(pbradii_).c_str());
  return 0;
}
