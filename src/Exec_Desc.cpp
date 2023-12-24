#include "Exec_Desc.h"
#include "CpptrajStdio.h"

int Exec_Desc::desc_atom(Topology const& topIn, int iat) {
  Atom const& AT = topIn[iat];

  const int bondorder = 1; // FIXME
  mprintf("ATOM\n");
  mprintf("Name:         %-5s\n", *(AT.Name()));
  mprintf("Type:         %-5s\n", *(AT.Type()));
  mprintf("Charge:       %6.4f\n", AT.Charge());
  mprintf("Polarization: %6.4f\n", 0.0); // FIXME
  mprintf("Element:      %-5s\n", AT.ElementName());
  if (AT.Nbonds() == 0)
    mprintf("  NO BONDS\n");
  else {
    for (Atom::bond_iterator bat = AT.bondbegin(); bat != AT.bondend(); ++bat) {
      mprintf("  Bonded to %s by a", topIn.AtomMaskName(*bat).c_str());
      switch (bondorder) {
        case 1 : mprintf(" single"); break;
        case 2 : mprintf(" double"); break;
        case 3 : mprintf(" triple"); break;
        case 4 : mprintf(" n aromatic"); break;
        default : mprintf(" ????"); break;
      }
      mprintf(" bond.\n");
    }
  }
  return 0;
}
 

// Exec_Desc::Help()
void Exec_Desc::Help() const
{

}

// Exec_Desc::Execute()
Exec::RetType Exec_Desc::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get Topology
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) {
    mprinterr("Error: No topologies loaded.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tUsing topology: %s\n", parm->c_str());
  // Get output file
  std::string outname = argIn.GetStringKey("out");
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(outname, "Atom description",
                                                    DataFileList::TEXT, true);
  if (outfile == 0) return CpptrajState::ERR;
  mprintf("\tOutput to '%s'\n", outfile->Filename().full());
  // Get mask
  std::string maskstr = argIn.GetMaskNext();
  AtomMask mask;
  if (mask.SetMaskString(maskstr)) {
    mprinterr("Error: Invalid atom mask: %s\n", maskstr.c_str());
    return CpptrajState::ERR;
  }

  // Set up mask
  if (parm->SetupIntegerMask( mask )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }
  mask.MaskInfo();
  if (mask.None()) {
    mprintf("Warning: '%s' selects no atoms.\n", mask.MaskString());
    return CpptrajState::OK;
  }

  // Loop over atoms
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
    desc_atom( *parm, *at );

  return CpptrajState::OK;
}
