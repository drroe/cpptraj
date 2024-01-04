#include "GB_Params.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"

static const char* GB_RadiiTypeStr_[] = {
  "Bondi radii",
  "Amber 6 modified Bondi radii",
  "Modified Bondi radii",
  "Radii optimized for Amber charges by Huo and Kollman",
  "H(N)-modified Bondi radii",
  "PARSE radii",
  "ArgH and AspGluO modified Bondi2 radii",
  "Unknown GB radii set"
};

static const char* GB_RadiiTypeKey_[] = {
  "bondi",   // 0
  "amber6",  // 1
  "mbondi",  // 2
  "pbamber", // 3
  "mbondi2", // 6
  "parse",   // 7
  "mbondi3", // 8
  0
};

/// These are the numeric iGBparm values used in LEaP
static const int GB_RadiiTypeIGB_[] = {
  0, // bondi
  1, // amber6
  2, // mbondi
  3, // pbamber
  6, // mbondi2
  7, // parse
  8, // mbondi3
  999999
};

/** Assign GB radii. */
static void assign_gb_radii(Topology& top, Cpptraj::Parm::GB_RadiiType radiiSet)
{
  int iGBparm = GB_RadiiTypeIGB_[radiiSet];
  for (int iat = 0; iat != top.Natom(); iat++)
  {
    Atom const& currentAtom = top[iat];
    Residue const& currentRes = top.Res(currentAtom.ResNum());
    double dGBrad = 0.0;
//    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
//    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
//    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
//                &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
//                &iElement, &iHybridization, sDesc);
    if (iGBparm < 3 || iGBparm == 6 || iGBparm == 8) {
      // Bondi or modified Bondi radii
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.2;

          // Make the modifications that hydrogen radii
          // depend upon the atoms they are bonded to.  
          // iGBparm=1 corresponds to Amber 6, JACS 122:2489 (2000);
          // iGBparm=2 adds the update of Biopolymers 56: 275 (2001)   
          if (currentAtom.Nbonds() > 0) {

	    // For multiply bonded Hydrogen atoms use the first bond for
            // determining modified GB radii.  WAT contains multiply bonded
	    // Hydrogen atoms so do not emit a warning.
            Atom const& aAtomA = top[currentAtom.Bond(0)];
            if (iGBparm == 1 || iGBparm == 2) {
              switch (aAtomA.Element()) {
                case Atom::CARBON: dGBrad = 1.3; break;      // Carbon
                case Atom::OXYGEN: dGBrad = 0.8; break;      // Oxygen
                case Atom::SULFUR: dGBrad = 0.8; break;      // Sulfur
                case Atom::NITROGEN:
                  if (iGBparm == 2) {
                    dGBrad = 1.3;
                  }
                  break;      // Nitrogen, mbondi
                case Atom::HYDROGEN:        // Special case: water hydrogen
                  if (aAtomA.Type() == "HW" || aAtomA.Type() == "hw") {
                                dGBrad = 0.8;
		  }
                  break;
                default : break;
              }
            }
	    else if (iGBparm == 6 || iGBparm == 8) {

              // Try Alexey's scheme
              if (currentAtom.Element() == Atom::NITROGEN) {
                dGBrad = 1.3;
                if (iGBparm == 8) {

                  // Update residue as appropriate
                  //if (saPAtom->iResidueIndex != iResidueIndex) {
                  //  iResidueIndex = saPAtom->iResidueIndex;
                  //  cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                  //                iResidueIndex - 1)->sName;
                  //  if (strlen(cPTemp) > 3) {
                  //    cPTemp += (strlen(cPTemp) - 3);
		  //  }
                  //}

		  // Adjust Arg HH and HE
                  if (currentRes.Name() == "ARG" && (currentAtom.Name().Match("HH*") || currentAtom.Name() == "HE"))
                  {
                  //if (!strcmp(cPTemp, "ARG") &&
		  //    !(strncmp(sAtomName(saPAtom->aAtom), "HH", 2) &&
                  //      strcmp(sAtomName(saPAtom->aAtom), "HE"))) {
                    dGBrad = 1.17;
                  }
                }
              }
            }
          }
          else {
            mprintf("Warning: Unbonded Hydrogen atom %s in %s.\n"
                    " Cannot determine the requested GB radius for this atom.\n"
                    " Writing the unmodified Bondi GB radius.\n",
                    *(currentAtom.Name()),
                    top.TruncResNameNum(currentAtom.ResNum()).c_str());
          }
          break;
        case Atom::CARBON:

	  // Use the mass of the carbon atom. We are testing for
          // carbons here. C1 == CH, C2 == CH2, C3 == CH3. UA carbons
          // have a larger radius (2.2), so we want to make sure that
          // the C1, C2, and C3 atom types _really_ correspond to UA
          // UA carbons. C1 atoms should have a mass of 12.01 + 1.01,
          // C2 should be 12.01 + 2.02, and C3 should be 12.01 + 3.03.
          // This mneumonic will not work for 13C named "C1". This is
          // a (hopefully) temporary kludge. FIXME not sure this first line is right
          //if (strncmp(sType, "C1", 2) && strncmp(sType, "C2", 2) && strncmp(sType, "C3", 2)) {
	  if (currentAtom.Type() != "C1" && currentAtom.Type() != "C2" && currentAtom.Type() != "C3") {
            dGBrad = 1.7;
	  }
          else if (currentAtom.Type() == "C1" && currentAtom.Mass() < 13.0) {
            dGBrad = 1.7;
	  }
	  else if (currentAtom.Type() == "C2" && currentAtom.Mass() < 14.0) {
	    dGBrad = 1.7;
	  }
	  else if (currentAtom.Type() == "C3" && currentAtom.Mass() < 15.0) {
	    dGBrad = 1.7;
	  }
          else {
            dGBrad = 2.2;
	  }
          break;
        case Atom::NITROGEN:
          dGBrad = 1.55;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.5;
          if (iGBparm == 8) {

            // Update residue as appropriate
            //if (saPAtom->iResidueIndex != iResidueIndex) {
            //  iResidueIndex = saPAtom->iResidueIndex;
            //  cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, iResidueIndex - 1)->sName;
            //  if (strlen(cPTemp) > 3) {
            //    cPTemp += (strlen(cPTemp) - 3);
	    //  }
            //}
	    
            // Adjust Asp OD and Glu OE, and terminal OXT
            if ( (currentRes.Name() == "ASP" && currentAtom.Name().Match("OD*")) ||
                 (currentRes.Name() == "AS4" && currentAtom.Name().Match("OD*")) ||
                 (currentRes.Name() == "GLU" && currentAtom.Name().Match("OE*")) ||
                 (currentRes.Name() == "GL4" && currentAtom.Name().Match("OE*")) ||
                 currentAtom.Name() == "OXT")
            {
            //if (!(strcmp(cPTemp, "ASP") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
            //      !(strcmp(cPTemp, "AS4") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
            //      !(strcmp(cPTemp, "GLU") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
            //      !(strcmp(cPTemp, "GL4") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
            //      (!strcmp(sAtomName(saPAtom->aAtom), "OXT") ||
            //      (i + 1 < iVarArrayElementCount(uUnit->vaAtoms) &&
            //       !strcmp(sAtomName(PVAI(uUnit->vaAtoms, SAVEATOMt, i + 1)->aAtom),
            //      "OXT")))) {
	      dGBrad = 1.4;
	    }
          }
          break;
        case Atom::FLUORINE:
          dGBrad = 1.5;
          break;
        case Atom::SILICON:
          dGBrad = 2.1;
          break;
        case Atom::PHOSPHORUS:
          dGBrad = 1.85;
          break;
        case Atom::SULFUR:
          dGBrad = 1.8;
          break;
        case Atom::CHLORINE:
          dGBrad = 1.7;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (iGBparm == 3) {

      // Radii from Huo & Kollman
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.15;
          break;
        case Atom::CARBON:
          dGBrad = 1.85;
          break;
        case Atom::NITROGEN:
          dGBrad = 1.64;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.53;
          break;
        case Atom::FLUORINE:
          dGBrad = 1.53;
          break;
        case Atom::PHOSPHORUS:
          dGBrad = 2.02;
          break;
        case Atom::SULFUR:
          dGBrad = 2.00;
          break;
        case Atom::CHLORINE:
          dGBrad = 1.97;
          break;
        case Atom::BROMINE:
          dGBrad = 2.03;
          break;
        case Atom::IODINE:
          dGBrad = 2.10;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (iGBparm == 7) {

      // Parse radii
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.00;
          break;
        case Atom::CARBON:
          dGBrad = 1.70;
          break;
        case Atom::NITROGEN:
          dGBrad = 1.50;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.40;
          break;
        case Atom::SULFUR:
          dGBrad = 1.85;
          break;
        default:
          dGBrad = 1.50;
          break;
          // Radii from J. Phys. Chem. 1994, 98, 1978-1988
      }
    }
    top.SetAtom(iat).SetGBradius( dGBrad );
  } // END loop over atoms
}

/** Assign GB radii and screening parameters based on the given radius set. */
int Cpptraj::Parm::Assign_GB_Radii(Topology& top, GB_RadiiType radiiSet)
{
  if (radiiSet == UNKNOWN_GB) {
    mprinterr("Error: Unknown GB radii set.\n");
    return 1;
  }
  mprintf("\tUsing GB radii set: %s\n", GB_RadiiTypeStr_[radiiSet]);
  top.SetGBradiiSet( std::string(GB_RadiiTypeStr_[radiiSet]) );
  assign_gb_radii( top, radiiSet );

  return 0;
}
