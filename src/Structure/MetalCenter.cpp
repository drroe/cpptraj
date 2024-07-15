#include "MetalCenter.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
MetalCenter::MetalCenter() {}

/** Init with args */
int MetalCenter::InitMetalCenters(ArgList& argIn)
{
  std::string metalMaskStr = argIn.GetStringKey("metalmask");
  if (metalMaskStr.empty()) {
    metalMaskStr.assign("@ZN,MG");
    mprintf("\tUsing default metal mask string.\n");
  }
  if (metalMask_.SetMaskString( metalMaskStr )) {
    mprinterr("Error: Invalid mask '%s' given for 'metalmask'\n", metalMaskStr.c_str());
    return 1;
  }

  std::string coordAtomMaskStr = argIn.GetStringKey("coordatommask");
  if (coordAtomMaskStr.empty()) {
    coordAtomMaskStr.assign("@/O,S");
    mprintf("\tUsing default coordinating atom mask string.\n");
  }
  if (coordAtomMask_.SetMaskString( coordAtomMaskStr )) {
    mprinterr("Error: Invalid mask '%s' given for 'coordatommask'\n", coordAtomMaskStr.c_str());
    return 1;
  }

  return 0;
}

/** Print info to stdout */
void MetalCenter::PrintMetalCenterInfo() const {
  mprintf("\tMetal center mask: %s\n", metalMask_.MaskString());
  mprintf("\tCoordinating atom mask: %s\n", coordAtomMask_.MaskString());
}
