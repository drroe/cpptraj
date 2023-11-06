#include "Remote.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include <cstdlib> // system

using namespace Cpptraj;

/** The remote download command. */
std::string Remote::cmd_ = "";

std::string Remote::oflag_ = "";

/** CONSTRUCTOR */
Remote::Remote() {
  setRemoteDownloadCommand();
}

/** CONSTRUCTOR - base url */
Remote::Remote(std::string const& baseUrl) :
  url_(baseUrl)
{
  setRemoteDownloadCommand();
}

/** Set the remote download command. */
int Remote::setRemoteDownloadCommand() {
  if (!cmd_.empty()) return 0;
  // First try curl
  int err = system("curl --version");
  if (err == 0) {
    mprintf("\tcurl found.\n");
    cmd_.assign("curl -L ");
    oflag_.assign("-o ");
  } else {
    err = system("wget --version");
    if (err == 0) {
      mprintf("\twget found.\n");
      cmd_.assign("wget ");
      oflag_.assign("-O ");
    } else {
      mprinterr("Error: No working remote command found.\n");
      return 1;
    }
  }
  return 0;
}

/** Download file from base URL */
int Remote::DownloadFile(std::string const& fname, std::string const& outputFnameIn)
const
{
  if (cmd_.empty()) {
    mprinterr("Error: No remote download command set.\n");
    return 1;
  }
  if (fname.empty()) {
    mprinterr("Error: No file name to download specified.\n");
    return 1;
  }
  FileName fileName(fname);
  // Set output file name if necessary
  FileName outputFname;
  if (outputFnameIn.empty())
    outputFname = fileName;
  else
    outputFname.SetFileName( outputFnameIn );

  std::string remoteUrl = url_ + "/" + fileName.Full();
  mprintf("\t %s => %s\n", remoteUrl.c_str(), outputFname.full());

  std::string remoteCmd = cmd_ + remoteUrl + " " + oflag_ + outputFname.Full();
  mprintf("DEBUG: %s\n", remoteCmd.c_str());

  return 0;
}
  
