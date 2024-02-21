#include "PdbCleaner.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
PdbCleaner::PdbCleaner() :
  remove_water_(false),
  remove_h_(false)
{}

/** Initialize */
int PdbCleaner::InitPdbCleaner(ArgList& argIn, std::string const& solventResNameIn,
                               Iarray const& pdbResToRemoveIn)
{
  resnumsToRemove_ = pdbResToRemoveIn;
  remove_water_ = argIn.hasKey("nowat");
  if (solventResNameIn.empty())
    waterMask_ = ":HOH";
  else
    waterMask_ = argIn.GetStringKey("watermask", ":" + solventResNameIn);
  remove_h_ = argIn.hasKey("noh");
  altLocArg_ = argIn.GetStringKey("keepaltloc");
  if (!altLocArg_.empty()) {
    if (altLocArg_ != "highestocc" &&
        altLocArg_.size() > 1)
    {
      mprinterr("Error: Invalid keyword for 'keepaltloc' '%s'; must be 'highestocc' or 1 character.\n",
                altLocArg_.c_str());
      return 1;
    }
  }
  stripMask_ = argIn.GetStringKey("stripmask");

  return 0;
}
