#include "Remote.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;

/** CONSTRUCTOR */
Remote::Remote() {}

/** CONSTRUCTOR - base url */
Remote::Remote(std::string const& baseUrl) :
  url_(baseUrl)
{}

/** Download file from base URL */
int Remote::DownloadFile(std::string const& fname, std::string const& outputFname)
const
{
  return 0;
}
  
