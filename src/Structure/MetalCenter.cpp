#include "MetalCenter.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

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

/** Find metal centers. */
int MetalCenter::FindMetalCenters(Topology const& topIn, Frame const& frameIn)
{
  if (topIn.SetupIntegerMask( metalMask_, frameIn )) {
    mprinterr("Error: Could not set up metal center mask '%s'\n", metalMask_.MaskString());
    return 1;
  }
  if (metalMask_.None()) {
    mprintf("Warning: Nothing selected by metal center mask '%s'\n", metalMask_.MaskString());
    return 0;
  }
  metalMask_.MaskInfo();

  if (topIn.SetupIntegerMask( coordAtomMask_, frameIn )) {
    mprinterr("Error: Could not set up metal center mask '%s'\n", coordAtomMask_.MaskString());
    return 1;
  }
  if (coordAtomMask_.None()) {
    mprintf("Warning: Nothing selected by coordinating atom mask '%s'\n", coordAtomMask_.MaskString());
    return 0;
  }
  coordAtomMask_.MaskInfo();

  return 0;
}
