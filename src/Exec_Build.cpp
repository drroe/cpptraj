#include "Exec_Build.h"
#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h" // For casting DataSet_Parameters to ParameterSet
#include "Parm/AssignParams.h"
#include "Parm/GB_Params.h"
#include "Structure/Builder.h"
#include "Structure/Creator.h"
#include "Structure/Disulfide.h"
#include "Structure/HisProt.h"
#include "Structure/PdbCleaner.h"
#include "Structure/ResStatArray.h"
#include "Structure/SugarBuilder.h"
#include "Structure/Sugar.h"
#include "StructureCheck.h"
#include <set> // For warning about missing residue templates

/** CONSTRUCTOR */
Exec_Build::Exec_Build() :
  Exec(GENERAL),
  debug_(0),
  check_box_natom_(5000), // TODO make user specifiable
  check_structure_(true),
  addNonTemplateBonds_(false),
  sugarBuilder_(0)
{}

/** DESTRUCTOR */
Exec_Build::~Exec_Build() {
  if (sugarBuilder_ != 0)
    delete sugarBuilder_;
}

/** Map atoms in residue to template. */
std::vector<int> Exec_Build::MapAtomsToTemplate(Topology const& topIn,
                                                int rnum,
                                                DataSet_Coords* resTemplate,
                                                Cpptraj::Structure::Creator const& creator,
                                                std::vector<NameType>& sourceAtomNames,
                                                int& nTgtAtomsMissing)
{
  nTgtAtomsMissing = 0;
  std::vector<int> mapOut(resTemplate->Top().Natom(), -1);
  mapOut.reserve( resTemplate->Top().Natom() );
  Residue const& resIn = topIn.Res(rnum);
  // For each atom in topIn, find a template atom
  for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++)
  {
    NameType const& tgtName = topIn[itgt].Name();
    //mprintf("DEBUG: Search for atom %s\n", *tgtName);
    bool found = false;
    // Check if this atom has an alias.
    NameType alias;
    bool has_alias = creator.GetAlias( alias, tgtName );
//    if (creator.GetAlias( alias, tgtName )) {
//      mprintf("DEBUG: Atom %s alias is %s\n", *tgtName, *alias);
//    }
    // See if tgtName matches a reference (template) atom name.
    for (int iref = 0; iref != resTemplate->Top().Natom(); iref++)
    {
      NameType const& refName = resTemplate->Top()[iref].Name();
      if (refName == tgtName) {
        sourceAtomNames.push_back( tgtName );
        mapOut[iref] = itgt;
        found = true;
        break;
      }
    }
    if (!found && has_alias) {
      // See if alias matches a reference (template) atom name.
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        NameType const& refName = resTemplate->Top()[iref].Name();
        if (refName == alias) {
          sourceAtomNames.push_back( alias );
          mapOut[iref] = itgt;
          found = true;
          break;
        }
      } // END search template for alias
    } // END do alias search
    if (!found) {
      mprintf("Warning: Input atom %s was not mapped to a template atom.\n",
              topIn.TruncAtomResNameOnumId( itgt ).c_str());
      nTgtAtomsMissing++;
    }
  }
/*
  for (int iref = 0; iref != resTemplate->Top().Natom(); iref++)
  {
    // Find this atom name in topIn
    NameType const& refName = resTemplate->Top()[iref].Name();
    int iat = -1;
    for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++) {
      if ( refName == topIn[itgt].Name() ) {
        iat = itgt;
        break;
      }
    }
    
    mapOut.push_back( iat );
  }*/
  return mapOut;
}

/** Search in array of atom bonding pairs for given bonding pair. */
bool Exec_Build::hasBondingPair(IParray const& bpairs, Ipair const& bpair) {
  for (IParray::const_iterator it = bpairs.begin(); it != bpairs.end(); ++it)
    if (*it == bpair) return true;
  return false;
}

/** \return True if target residue is in array of residue connections. */
bool Exec_Build::resIsConnected(Iarray const& resConnections, int tgtRes) {
  for (Iarray::const_iterator it = resConnections.begin(); it != resConnections.end(); ++it)
    if (*it == tgtRes) return true;
  return false;
}

/** Use given templates to construct a final molecule. */
int Exec_Build::FillAtomsWithTemplates(Topology& topOut, Frame& frameOut,
                                       Topology const& topIn, Frame const& frameIn,
                                       Cpptraj::Structure::Creator const& creator)
const
{
  // Array of head/tail connect atoms for each residue
  Iarray resHeadAtoms;
  Iarray resTailAtoms;
  std::vector<Cpptraj::Structure::TerminalType> ResTermTypes;
  resHeadAtoms.reserve( topIn.Nres() );
  resTailAtoms.reserve( topIn.Nres() );
  ResTermTypes.reserve( topIn.Nres() );
  // Array of templates for each residue
  std::vector<DataSet_Coords*> ResTemplates;
  ResTemplates.reserve( topIn.Nres() );
  //typedef std::vector<Cpptraj::Structure::TerminalType> TermTypeArray;
  //TermTypeArray ResTypes;
  //ResTypes.reserve( topIn.Nres() );
  // Initial loop to try to match residues to templates
  int newNatom = 0;
  unsigned int n_no_template_found = 0;
  std::set<NameType> missing_templates;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& currentRes = topIn.Res(ires);
    if (debug_ > 0)
      mprintf("DEBUG: ---------- Processing Residue %s ---------- \n", topIn.TruncResNameNum(ires).c_str());
    int pres = ires - 1;
    int nres = ires + 1;
    // Determine if this is a terminal residue
    Cpptraj::Structure::TerminalType resTermType;
    if (currentRes.IsTerminal()) {
      //mprintf("DEBUG: End terminal due to TERMINAL status.\n");
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else if (ires == 0 && topIn.Nres() > 1) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
    } else if (pres > -1 && topIn.Res(pres).IsTerminal()) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
    } else if (nres < topIn.Nres() && (topIn.Res(nres).ChainID() != currentRes.ChainID())// ||
//                                        (topIn.Res(nres).OriginalResNum() == currentRes.OriginalResNum() &&
//                                         topIn.Res(nres).Icode()          != currentRes.Icode()
//                                        )
//                                      )
              )
    {
      //mprintf("DEBUG: End terminal due to chain ID.\n");
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else if (nres == topIn.Nres()) {
      //mprintf("DEBUG: End terminal due to last residue.\n");
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else {
      resTermType = Cpptraj::Structure::NON_TERMINAL;
    }
    if (debug_ > 0)
      mprintf("DEBUG: Residue type: %s terminal (IsTerminal=%i)\n", Cpptraj::Structure::terminalStr(resTermType), (int)currentRes.IsTerminal());
    ResTermTypes.push_back( resTermType );
    // Identify a template based on the residue name.
    DataSet_Coords* resTemplate = creator.IdTemplateFromResname(currentRes.Name(), resTermType);
    if (resTemplate == 0) {
      // Residue has no template.
      n_no_template_found++;
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
      missing_templates.insert( currentRes.Name() );
      newNatom += currentRes.NumAtoms();
      // Head and tail atoms are blank
      resHeadAtoms.push_back( -1 );
      resTailAtoms.push_back( -1 );
    } else {
      // Residue has a template.
      if (debug_ > 0)
        mprintf("\tTemplate %s being used for residue %s\n",
                resTemplate->legend(), topIn.TruncResNameOnumId(ires).c_str());
      // Save the head and tail atoms
      AssociatedData* ad = resTemplate->GetAssociatedData(AssociatedData::CONNECT);
      if (ad == 0) {
        mprintf("Warning: Unit '%s' does not have CONNECT data.\n", resTemplate->legend());
        resHeadAtoms.push_back( -1 );
        resTailAtoms.push_back( -1 );
      } else {
        AssociatedData_Connect const& CONN = static_cast<AssociatedData_Connect const&>( *ad );
        if (CONN.NconnectAtoms() < 2) {
          mprinterr("Error: Not enough connect atoms in unit '%s'\n", resTemplate->legend());
          return 1;
        }
        if (CONN.Connect()[0] > -1)
          resHeadAtoms.push_back( CONN.Connect()[0] + newNatom );
        else
          resHeadAtoms.push_back( -1 );
        if (CONN.Connect()[1] > -1)
          resTailAtoms.push_back( CONN.Connect()[1] + newNatom );
        else
          resTailAtoms.push_back( -1 );
      }
      // Update # of atoms
      newNatom += resTemplate->Top().Natom();
    }
    ResTemplates.push_back( resTemplate );
  }
  mprintf("\tFinal structure should have %i atoms.\n", newNatom);
  if (n_no_template_found > 0) {
    mprintf("Warning: No template was found for %u residues.\n", n_no_template_found);
    mprintf("Warning: Residue names:");
    for (std::set<NameType>::const_iterator rit = missing_templates.begin();
                                            rit != missing_templates.end(); ++rit)
      mprintf(" %s", *(*rit));
    mprintf("\n");
    mprintf("Warning: This may indicate that you have not loaded a library file\n"
            "Warning:   or have not loaded the correct force field.\n");
  }
  frameOut.SetupFrame( newNatom );
  // Clear frame so that AddXYZ can be used
  frameOut.ClearAtoms();

  // -----------------------------------
  // hasPosition - for each atom in topOut, status on whether atom in frameOut needs building
  Cpptraj::Structure::Builder::Barray hasPosition;
  hasPosition.reserve( newNatom );

  // Hold atom offsets needed when building residues
  Iarray AtomOffsets;
  AtomOffsets.reserve( topIn.Nres() );

  // For existing inter-residue bonding, use residue # and atom name since
  // atom numbering may change if atoms are added from templates.
  // TODO make this a class var so disulfide/sugar prep can use it.
  typedef std::pair<int,NameType> ResAtPair;
  typedef std::vector<ResAtPair> ResAtArray;
  ResAtArray detectedInterResBonds;

  // Hold template atom names corressponding to source atom names.
  typedef std::vector<NameType> NameArray;
  NameArray SourceAtomNames;
  SourceAtomNames.reserve( topIn.Natom() );

  // Loop for setting up atoms in the topology from residues or residue templates.
  int nRefAtomsMissing = 0;
  int nAtomsMissingTypes = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    if (debug_ > 0)
      mprintf("\tAdding atoms for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
    int atomOffset = topOut.Natom();
    //mprintf("DEBUG: atom offset is %i\n", atomOffset);
    DataSet_Coords* resTemplate = ResTemplates[ires];
    IParray intraResBonds;
    if (resTemplate == 0) {
      // ----- No template. Just add the atoms. ------------
      Residue const& currentRes = topIn.Res(ires);
      AtomOffsets.push_back( -1 );
      for (int itgt = currentRes.FirstAtom(); itgt != currentRes.LastAtom(); ++itgt)
      {
        // Track intra-residue bonds
        Atom sourceAtom = topIn[itgt];
        if (!sourceAtom.HasType())
          nAtomsMissingTypes++;
        SourceAtomNames.push_back( sourceAtom.Name() );
        int at0 = itgt - currentRes.FirstAtom() + atomOffset;
        for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
          if ( topIn[*bat].ResNum() == ires ) {
            // Intra-residue
            int at1 = *bat - currentRes.FirstAtom() + atomOffset;
            if (at1 > at0) {
              //mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, itgt+1, *bat + 1);
              intraResBonds.push_back( Ipair(at0, at1) );
            }
          } else {
            // Inter-residue. Only record if bonding to a previous residue.
            if (topIn[*bat].ResNum() < ires) {
              detectedInterResBonds.push_back( ResAtPair(ires, sourceAtom.Name()) );
              detectedInterResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), topIn[*bat].Name()) );
            }
          }
        }
        sourceAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
        topOut.AddTopAtom( sourceAtom, currentRes );
        frameOut.AddVec3( Vec3(frameIn.XYZ(itgt)) );
        hasPosition.push_back( true );
      }
    } else {
      // ----- A template exists for this residue. ---------
      Residue currentRes = topIn.Res(ires);
      // Use template residue name.
      // To match LEaP behavior, if the template name is > 3 characters,
      // truncate to the last 3 characters.
      NameType const& templateResName = resTemplate->Top().Res(0).Name();
      if (templateResName.len() < 4)
        currentRes.SetName( templateResName );
      else
        currentRes.SetName( NameType( (*templateResName) + ( templateResName.len() - 3) ) );
      // Map source atoms to template atoms.
      int nTgtAtomsMissing = 0;
      std::vector<int> map = MapAtomsToTemplate( topIn, ires, resTemplate, creator, SourceAtomNames, nTgtAtomsMissing );
      if (debug_ > 1) {
        mprintf("\t  Atom map:\n");
        // DEBUG - print map
        for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
          mprintf("\t\t%6i %6s =>", iref+1, *(resTemplate->Top()[iref].Name()));
          if (map[iref] == -1)
            mprintf(" No match\n");
          else
            mprintf(" %6i %6s\n", map[iref]+1, *(topIn[map[iref]].Name()));
        }
      }
      // Map template atoms back to source atoms.
      std::vector<int> pdb(currentRes.NumAtoms(), -1);
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        if (map[iref] != -1)
          pdb[map[iref]-currentRes.FirstAtom()] = iref;
      }
      bool atomsNeedBuilding = false;
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        // Track intra-residue bonds from the template.
        Atom templateAtom = resTemplate->Top()[iref];
        int at0 = iref + atomOffset;
        for (Atom::bond_iterator bat = templateAtom.bondbegin(); bat != templateAtom.bondend(); ++bat) {
          int at1 = *bat + atomOffset;
          if (at1 > at0) {
            //mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, iref+1, *bat + 1);
            intraResBonds.push_back( Ipair(at0, at1) );
          }
        }
        templateAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
        //mprintf("DEBUG: Adding template atom %s (elt %s)\n", templateAtom.c_str(), templateAtom.ElementName());
        topOut.AddTopAtom( templateAtom, currentRes );
        if (map[iref] == -1) {
          // Template atom not in input structure.
          frameOut.AddVec3( Vec3(0.0) );
          hasPosition.push_back( false );
          nRefAtomsMissing++;
          atomsNeedBuilding = true;
        } else {
          // Template atom was in input structure.
          int itgt = map[iref];
          frameOut.AddVec3( Vec3(frameIn.XYZ(itgt)) );
          hasPosition.push_back( true );
          //pdb[itgt-currentRes.FirstAtom()] = iref;
          // Check source atoms for inter-residue connections.
          Atom const& sourceAtom = topIn[itgt];
          for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
            if ( topIn[*bat].ResNum() < ires ) {
              // Use template atom names. Use saved names in case source name had an alias.
              detectedInterResBonds.push_back( ResAtPair(ires, templateAtom.Name()) );
              detectedInterResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), SourceAtomNames[*bat]) );
            }
          }
        }
      }
      // See if any current atoms were not mapped to reference template
      //int nTgtAtomsMissing = 0;
      //for (int itgt = 0; itgt != currentRes.NumAtoms(); itgt++) {
      //  if (pdb[itgt] == -1) {
      //    mprintf("Warning: Input atom %s was not mapped to a template atom.\n",
      //            topIn.TruncAtomResNameOnumId(currentRes.FirstAtom()+itgt).c_str());
      //    nTgtAtomsMissing++;
      //  }
      //}
      if (nTgtAtomsMissing > 0)
        mprintf("\t%i source atoms not mapped to template.\n", nTgtAtomsMissing);
      // Save atom offset if atoms need to be built
      if (atomsNeedBuilding)
        AtomOffsets.push_back( atomOffset );
      else
        AtomOffsets.push_back( -1 );
    } // END template exists
    // Add intra-residue bonds
    for (IParray::const_iterator it = intraResBonds.begin(); it != intraResBonds.end(); ++it)
    {
      //mprintf("DEBUG: Intra-res bond: Res %s atom %s to res %s atom %s\n",
      //        topOut.TruncResNameOnumId(topOut[it->first].ResNum()).c_str(), *(topOut[it->first].Name()),
      //        topOut.TruncResNameOnumId(topOut[it->second].ResNum()).c_str(), *(topOut[it->second].Name()));
      topOut.AddBond(it->first, it->second);
    }
  } // END loop over source residues
  mprintf("\t%i template atoms missing in source.\n", nRefAtomsMissing);
  if (nAtomsMissingTypes > 0) {
    mprinterr("Error: %i atoms are missing types, either because they did not have\n"
              "Error:  one initially or they could not be matched to a template.\n"
              "Error: This can happen if a parameter file is missing or a force field\n"
              "Error:  file has not been loaded.\n"
              "Error: Build cannot proceed unless all atoms have a type:\n",
              nAtomsMissingTypes);
    std::set<NameType> atoms_missing_types;
    for (int ires = 0; ires != topOut.Nres(); ires++) {
      for (int at = topOut.Res(ires).FirstAtom(); at != topOut.Res(ires).LastAtom(); ++at) {
        if ( !topOut[at].HasType() )
          atoms_missing_types.insert( topOut[at].Name() );
      }
    }
    mprinterr("Error: Atoms missing types:");
    for (std::set<NameType>::const_iterator ait = atoms_missing_types.begin();
                                            ait != atoms_missing_types.end(); ++ait)
      mprinterr(" %s", *(*ait));
    mprinterr("\n");
    if (debug_ > 0) {
      for (int ires = 0; ires != topOut.Nres(); ires++) {
        std::string missingTypes;
        for (int at = topOut.Res(ires).FirstAtom(); at != topOut.Res(ires).LastAtom(); ++at)
          if ( !topOut[at].HasType() )
            missingTypes.append(" " + topOut[at].Name().Truncated() );
        if (!missingTypes.empty())
          mprinterr("Error:\t%s missing types for%s\n", topOut.TruncResNameNum(ires).c_str(), missingTypes.c_str());
      }
    }
    return 1;
  }

  // -----------------------------------
  // DEBUG - Print primary connection atoms
  if (debug_ > 0) {
    for (unsigned int idx = 0; idx != ResTemplates.size(); idx++) {
      if (ResTemplates[idx] != 0) {
        mprintf("DEBUG: Template %s (%s)", ResTemplates[idx]->legend(), Cpptraj::Structure::terminalStr(ResTermTypes[idx]));
        if (resHeadAtoms[idx] > -1) mprintf(" head %s", topOut.AtomMaskName(resHeadAtoms[idx]).c_str());
        if (resTailAtoms[idx] > -1) mprintf(" tail %s", topOut.AtomMaskName(resTailAtoms[idx]).c_str());
        mprintf("\n");
      }
    }
  }

  // Keep track of which residues are connected.
  ResConnectArray ResidueConnections( topOut.Nres() );

  // Try to connect HEAD atoms to previous residue TAIL atoms.
  std::vector<IParray> resBondingAtoms(topOut.Nres());
  for (int ires = 1; ires < topOut.Nres(); ires++) {
    int pres = ires - 1;
    if (resHeadAtoms[ires] != -1) {
      if (ResTermTypes[ires] == Cpptraj::Structure::BEG_TERMINAL) {
        if (debug_ > 0)
          mprintf("DEBUG: Res %s is begin terminal, ignoring head atom.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
      } else {
        if (resTailAtoms[pres] != -1) {
          if (ResTermTypes[pres] == Cpptraj::Structure::END_TERMINAL) {
            if (debug_ > 0)
              mprintf("DEBUG: Res %s is end terminal, ignoring tail atom.\n",
                      topOut.TruncResNameOnumId(pres).c_str());
          } else {
            if (debug_ > 0)
              mprintf("DEBUG: Connecting HEAD atom %s to tail atom %s\n",
                      topOut.AtomMaskName(resHeadAtoms[ires]).c_str(),
                      topOut.AtomMaskName(resTailAtoms[pres]).c_str());
            resBondingAtoms[ires].push_back( Ipair(resHeadAtoms[ires], resTailAtoms[pres]) );
            resBondingAtoms[pres].push_back( Ipair(resTailAtoms[pres], resHeadAtoms[ires]) );
            resHeadAtoms[ires] = -1;
            resTailAtoms[pres] = -1;
            ResidueConnections[ires].push_back( pres );
            ResidueConnections[pres].push_back( ires );
          }
        }
      }
    }
  }

  // Report unused HEAD/TAIL atoms
  for (int ires = 0; ires != topOut.Nres(); ires++) { // TODO should be topIn?
    if (resHeadAtoms[ires] != -1)
      mprintf("Warning: Unused head atom %s\n", topOut.AtomMaskName(resHeadAtoms[ires]).c_str());
    if (resTailAtoms[ires] != -1)
      mprintf("Warning: Unused tail atom %s\n", topOut.AtomMaskName(resTailAtoms[ires]).c_str());
  }

  // Check detected inter-residue bonds
  for (ResAtArray::const_iterator it = detectedInterResBonds.begin();
                                  it != detectedInterResBonds.end(); ++it)
  {
    ResAtPair const& ra0 = *it;
    ++it;
    ResAtPair const& ra1 = *it;
    if (debug_ > 1)
      mprintf("DEBUG: Inter-res bond: Res %i atom %s to res %i atom %s\n",
              ra0.first+1, *(ra0.second),
              ra1.first+1, *(ra1.second));
    int at0 = topOut.FindAtomInResidue(ra0.first, ra0.second);
    if (at0 < 0) {
      mprinterr("Error: Atom %s not found in residue %i\n", *(ra0.second), ra0.first);
      return 1;
    }
    int at1 = topOut.FindAtomInResidue(ra1.first, ra1.second);
    if (at1 < 0) {
      mprinterr("Error: Atom %s not found in residue %i\n", *(ra1.second), ra1.first);
      return 1;
    }
    // Save detected inter-residue bonding atoms if not already added via
    // template connect atoms. Convention is atom belonging to the current
    // residue is first.
    // NOTE: Only checking at0/at1 here, which should be fine.
    if (!hasBondingPair(resBondingAtoms[ra0.first], Ipair(at0, at1))) {
      // Check if we already have a connection from ra0 to ra1.
      if (resIsConnected(ResidueConnections[ra0.first], ra1.first)) {
        mprintf("Warning: Residue %s already connected to residue %s; ignoring\n"
                "Warning:   potential detected bond %s - %s\n",
                topOut.TruncResNameNum(ra0.first).c_str(),
                topOut.TruncResNameNum(ra1.first).c_str(),
                topOut.AtomMaskName(at0).c_str(),
                topOut.AtomMaskName(at1).c_str());
      } else {
        if (addNonTemplateBonds_) { // FIXME need to change this or move the sugar/disulfide stuff to after
          mprintf("\tAdding non-template bond %s - %s\n",
                  topOut.AtomMaskName(at0).c_str(),
                  topOut.AtomMaskName(at1).c_str());
          resBondingAtoms[ra0.first].push_back( Ipair(at0, at1) );
          resBondingAtoms[ra1.first].push_back( Ipair(at1, at0) );
          ResidueConnections[ra0.first].push_back( ra1.first );
          ResidueConnections[ra1.first].push_back( ra0.first );
        } else {
          mprintf("Warning: Detected non-template bond %s - %s; not adding it.\n",
                  topOut.AtomMaskName(at0).c_str(),
                  topOut.AtomMaskName(at1).c_str());
        }
      }
    } //else
      //mprintf("DEBUG: Detected bond %s - %s already present.\n",
      //        topOut.AtomMaskName(at0).c_str(),
      //        topOut.AtomMaskName(at1).c_str());
    //if (!hasBondingPair(resBondingAtoms[ra1.first], Ipair(at1, at0)))
    //  resBondingAtoms[ra1.first].push_back( Ipair(at1, at0) );
    //else
    //  mprintf("DEBUG: Detected bond %s - %s already present.\n",
    //          topOut.AtomMaskName(at1).c_str(),
    //          topOut.AtomMaskName(at0).c_str());
  }

  // For each inter-residue bonding atom pair, check that they match expected
  // head/tail atoms.
/*  for (std::vector<IParray>::const_iterator rit = resBondingAtoms.begin();
                                            rit != resBondingAtoms.end(); ++rit)
  {
    long int ires = rit - resBondingAtoms.begin();
    for (IParray::const_iterator atPair = rit->begin(); atPair != rit->end(); ++atPair)
    {
      long int jres = topOut[atPair->second].ResNum();
      if (atPair->first == resHeadAtoms[ires]) {
        // HEAD atom. Should connect to other residue TAIL atom.
        if (atPair->second != resTailAtoms[jres]) {
          mprintf("Warning: Detected inter-res bond %s - %s HEAD does not match TAIL.\n",
                  topOut.AtomMaskName(atPair->first).c_str(),
                  topOut.AtomMaskName(atPair->second).c_str());
        }
      } else if (atPair->first == resTailAtoms[ires]) {
        // TAIL atom. Should connect to other residue HEAD atom.
        if (atPair->second != resHeadAtoms[jres]) {
          mprintf("Warning: Detected inter-res bond %s - %s TAIL does not match HEAD.\n",
                  topOut.AtomMaskName(atPair->first).c_str(),
                  topOut.AtomMaskName(atPair->second).c_str());
        }
      } else {
        mprintf("Warning: Atom %s is not a HEAD or TAIL atom.\n",
                topOut.AtomMaskName(atPair->first).c_str());
      }
    } // END loop over inter-residue bond pairs for this residue
  } // END loop over residues*/

  // Mark used HEAD/TAIL atoms in existing inter-residue bonds.
/*  for (std::vector<IParray>::const_iterator rit = resBondingAtoms.begin();
                                            rit != resBondingAtoms.end(); ++rit)
  {
    long int ires = rit - resBondingAtoms.begin();
    for (IParray::const_iterator atPair = rit->begin(); atPair != rit->end(); ++atPair)
    {
      if (atPair->first == resHeadAtoms[ires])
        resHeadAtoms[ires] = -1;
      else if (atPair->first == resTailAtoms[ires])
        resTailAtoms[ires] = -1;
    }
  }*/

  // DEBUG print inter-residue bonding atoms
  if (debug_ > 0) {
    for (std::vector<IParray>::const_iterator rit = resBondingAtoms.begin();
                                              rit != resBondingAtoms.end(); ++rit)
    {
      mprintf("\tResidue %s bonding atoms.\n", topOut.TruncResNameNum(rit-resBondingAtoms.begin()).c_str());
      for (IParray::const_iterator it = rit->begin(); it != rit->end(); ++it)
        mprintf("\t\t%s - %s\n", topOut.AtomMaskName(it->first).c_str(), topOut.AtomMaskName(it->second).c_str());
    }
  }

  // -----------------------------------
  // Do some error checking
  if (hasPosition.size() != (unsigned int)newNatom) {
    mprinterr("Internal Error: hasPosition size %zu != newNatom size %i\n", hasPosition.size(), newNatom);
    return 1;
  }
  if (AtomOffsets.size() != (unsigned int)topOut.Nres()) {
    mprinterr("Internal Error: AtomOffsets size %zu != newNres size %i\n", AtomOffsets.size(), topOut.Nres());
    return 1;
  }
  if (SourceAtomNames.size() != (unsigned int)topIn.Natom()) {
    mprinterr("Internal Error: Source atom names length %zu != # input atoms %i\n", SourceAtomNames.size(), topIn.Natom());
    return 1;
  }

  // -----------------------------------
  // Build using internal coords if needed.
//  std::vector<bool> resIsBuilt; // TODO is this needed?
//  resIsBuilt.reserve( AtomOffsets.size() );
//  for (Iarray::const_iterator it = AtomOffsets.begin(); it != AtomOffsets.end(); ++it) {
//    if ( *it < 0 )
//      resIsBuilt.push_back( true );
//    else
//      resIsBuilt.push_back( false );
//  }

  bool buildFailed = false;
  for (Iarray::const_iterator it = AtomOffsets.begin(); it != AtomOffsets.end(); ++it)
  {
    long int ires = it-AtomOffsets.begin();
    if (*it > -1) {
      if (debug_ > 0)
        mprintf("DEBUG: ***** BUILD residue %li %s *****\n", ires + 1,
                topOut.TruncResNameOnumId(ires).c_str());
      // Residue has atom offset which indicates it needs something built.
      Cpptraj::Structure::Builder structureBuilder;// = new Cpptraj::Structure::Builder();
      structureBuilder.SetDebug( debug_ );
      if (creator.HasMainParmSet())
        structureBuilder.SetParameters( creator.MainParmSetPtr() );
      // Generate internals from the template, update indices to this topology.
      DataSet_Coords* resTemplate = ResTemplates[ires];
      Frame templateFrame = resTemplate->AllocateFrame();
      resTemplate->GetFrame( 0, templateFrame );
      if (structureBuilder.GenerateInternals(templateFrame, resTemplate->Top(),
                                             std::vector<bool>(resTemplate->Top().Natom(), true)))
      {
        mprinterr("Error: Generate internals for residue template failed.\n");
        return 1;
      }
      structureBuilder.UpdateIndicesWithOffset( *it );
      //mprintf("DEBUG: Residue type: %s terminal\n", Cpptraj::Structure::terminalStr(*termType));
      // Is this residue connected to an earlier residue?
      for (IParray::const_iterator resBonds = resBondingAtoms[ires].begin();
                                   resBonds != resBondingAtoms[ires].end(); ++resBonds)
      {
        if (resBonds->second < resBonds->first) {
          if (debug_ > 0)
            mprintf("\t\tResidue connection: %s - %s\n",
                    topOut.AtomMaskName(resBonds->first).c_str(),
                    topOut.AtomMaskName(resBonds->second).c_str());
          topOut.AddBond(resBonds->first, resBonds->second);
          // Generate internals around the link
          if (structureBuilder.GenerateInternalsAroundLink(resBonds->first, resBonds->second,
                                                            frameOut, topOut, hasPosition, Cpptraj::Structure::Builder::BUILD))
          {
            mprinterr("Error: Assign torsions around inter-residue link %s - %s failed.\n",
                      topOut.AtomMaskName(resBonds->first).c_str(),
                      topOut.AtomMaskName(resBonds->second).c_str());
            return 1;
          }
        }
      }
      // Update internal coords from known positions
      if (structureBuilder.UpdateICsFromFrame( frameOut, topOut, hasPosition )) {
        mprinterr("Error: Failed to update internals with values from existing positions.\n");
        return 1;
      }
      if (structureBuilder.BuildFromInternals(frameOut, topOut, hasPosition)) {
        mprinterr("Error: Building residue %s failed.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
        buildFailed = true;
      }
    } else {
      // All atoms present. Just connect
      // Is this residue connected to an earlier residue?
      for (IParray::const_iterator resBonds = resBondingAtoms[ires].begin();
                                   resBonds != resBondingAtoms[ires].end(); ++resBonds)
      {
        if (resBonds->second < resBonds->first) {
          if (debug_ > 0)
            mprintf("\t\tResidue connection: %s - %s\n",
                    topOut.AtomMaskName(resBonds->first).c_str(),
                    topOut.AtomMaskName(resBonds->second).c_str());
          topOut.AddBond(resBonds->first, resBonds->second);
        }
      }
    }
  } // END loop over atom offsets

  // DEBUG - Print new top/coords
  if (debug_ > 1) {
    for (int iat = 0; iat != topOut.Natom(); iat++)
    {
      Residue const& res = topOut.Res( topOut[iat].ResNum() );
      const double* XYZ = frameOut.XYZ(iat);
      mprintf("%6i %6s %6i %6s (%i) %g %g %g\n",
              iat+1, *(topOut[iat].Name()), res.OriginalResNum(), *(res.Name()),
              (int)hasPosition[iat], XYZ[0], XYZ[1], XYZ[2]);
    }
  }

  // Finalize topology - determine molecules, dont renumber residues, dont assign default bond params
  topOut.CommonSetup(true, false, false);
  topOut.Summary();

  if (buildFailed) return 1;
  return 0;
}

/** Given an original topology and bonded atom indices, find those atoms
  * in another topology and ensure they are bonded.
  */
int Exec_Build::transfer_bonds(Topology& topOut, Topology const& topIn,
                               std::vector<BondType> const& bondsIn)
const
{
  for (std::vector<BondType>::const_iterator bnd = bondsIn.begin();
                                             bnd != bondsIn.end(); ++bnd)
  {
    // Get the original atom name and residue number
    Atom const& original_A1 = topIn[bnd->A1()];
    Atom const& original_A2 = topIn[bnd->A2()];
    mprintf("DEBUG: Original bond atoms %i (%s) %i (%s)\n",
            bnd->A1()+1, topIn.AtomMaskName(bnd->A1()).c_str(),
            bnd->A2()+1, topIn.AtomMaskName(bnd->A2()).c_str());
    // Find the atoms in the new topology
    int a1 = topOut.FindAtomInResidue( original_A1.ResNum(), original_A1.Name() );
    if (a1 < 0) {
      mprinterr("Error: Could not find atom %i (%s) in new topology.\n",
                bnd->A1()+1, topIn.AtomMaskName(bnd->A1()).c_str());
      return 1;
    }
    int a2 = topOut.FindAtomInResidue( original_A2.ResNum(), original_A2.Name() );
    if (a2 < 0) {
      mprinterr("Error: Could not find atom %i (%s) in new topology.\n",
                bnd->A2()+1, topIn.AtomMaskName(bnd->A2()).c_str());
      return 1;
    }
    // Add the bond to the new topology
    topOut.AddBond( a1, a2 );
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Exec_Build::Help()
void Exec_Build::Help() const
{
  mprintf("\tname <output COORDS> crdset <COORDS set> [frame <#>]\n"
          "\t[title <title>] [gb <radii>] [verbose <#>]\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n"
          "%s"
          "%s",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_,
          Cpptraj::Structure::HisProt::keywords_,
          Cpptraj::Structure::Disulfide::keywords_
         );
  mprintf("    <radii> =");
  for (int ig = 0; ig != (int)Cpptraj::Parm::UNKNOWN_GB; ig++)
    mprintf(" %s", Cpptraj::Parm::GbTypeKey((Cpptraj::Parm::GB_RadiiType)ig));
  mprintf("\n");
  mprintf("  Build complete topology and parameters from given crdset.\n");
}

// Exec_Build::Execute()
Exec::RetType Exec_Build::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get input coords
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify input COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* inCrdPtr = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (inCrdPtr == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  return BuildStructure(inCrdPtr, State.DSL(), State.Debug(), argIn);
}

/** Standalone execute. For DataIO_LeapRC. */
Exec::RetType Exec_Build::BuildStructure(DataSet* inCrdPtr, DataSetList& DSL, int debugIn, ArgList& argIn)
{
  t_total_.Start();
  if (inCrdPtr == 0) {
    mprinterr("Internal Error: Exec_Build::BuildStructure(): Null input coordinates.\n");
    return CpptrajState::ERR;
  }
  if (inCrdPtr->Group() != DataSet::COORDINATES) {
    mprinterr("Error: Set '%s' is not coordinates, cannot use for building.\n", inCrdPtr->legend());
    return CpptrajState::ERR;
  }
  debug_ = debugIn;
  int verbose = argIn.getKeyInt("verbose", 0);
  std::string title = argIn.GetStringKey("title");
  if (addNonTemplateBonds_)
    mprintf("\tNon-template bonds will be added if detected.\n");
  else
    mprintf("\tOnly bonding according to residue templates.\n");

  // TODO make it so this can be const (cant bc GetFrame)
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)inCrdPtr) );
  // Get frame from input coords
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);
  // Get modifiable topology
  //Topology& topIn = *(coords.TopPtr());
  //Topology const& topIn = coords.Top();
  Topology topIn = coords.Top(); // FIXME do not work on the copy, work on the top itself

  std::string solventResName = argIn.GetStringKey("solventresname", "HOH");
  mprintf("\tSolvent residue name: %s\n", solventResName.c_str());

  // Do histidine detection before H atoms are removed
  t_hisDetect_.Start();
  if (!argIn.hasKey("nohisdetect")) {
    Cpptraj::Structure::HisProt hisProt;
    if (hisProt.InitHisProt( argIn, debug_ )) {
      mprinterr("Error: Could not initialize histidine detection.\n");
      return CpptrajState::ERR;
    }
    hisProt.HisProtInfo();
    if (hisProt.DetermineHisProt( topIn )) {
      mprinterr("Error: HIS protonation detection failed.\n");
      return CpptrajState::ERR;
    }
  }
  t_hisDetect_.Stop();

  // Clean up structure
  t_clean_.Start();
  Cpptraj::Structure::PdbCleaner pdbCleaner;
  pdbCleaner.SetDebug( debug_ );
  if (pdbCleaner.InitPdbCleaner( argIn, solventResName, std::vector<int>() )) {
    mprinterr("Error: Could not init PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  if (pdbCleaner.SetupPdbCleaner( topIn )) {
    mprinterr("Error: Could not set up PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  pdbCleaner.PdbCleanerInfo();
  if (pdbCleaner.ModifyCoords(topIn, frameIn)) {
    mprinterr("Error: Could not clean PDB.\n");
    return CpptrajState::ERR;
  }
  t_clean_.Stop();

  // Set up Output coords
  std::string outset = argIn.GetStringKey("name");
  if (outset.empty()) {
    mprinterr("Error: Must specify output COORDS set with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet* outCrdPtr = DSL.AddSet( DataSet::COORDS, outset );
  if (outCrdPtr == 0) {
    mprinterr("Error: Could not allocate output COORDS set with name '%s'\n", outset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& crdout = static_cast<DataSet_Coords&>( *((DataSet_Coords*)outCrdPtr) );
  mprintf("\tOutput COORDS set: %s\n", crdout.legend());

  // GB radii set
  Cpptraj::Parm::GB_RadiiType gbradii = Cpptraj::Parm::MBONDI; // Default
  std::string gbset = argIn.GetStringKey("gb");
  if (!gbset.empty()) {
    gbradii = Cpptraj::Parm::GbTypeFromKey( gbset );
    if (gbradii == Cpptraj::Parm::UNKNOWN_GB) {
      mprinterr("Error: Unknown GB radii set: %s\n", gbset.c_str());
      return CpptrajState::ERR;
    }
  }
  mprintf("\tGB radii set: %s\n", Cpptraj::Parm::GbTypeStr(gbradii).c_str());

  // Get templates and parameter sets.
  t_get_templates_.Start();
  Cpptraj::Structure::Creator creator( debug_ );
  if (creator.InitCreator(argIn, DSL, debug_)) {
    return CpptrajState::ERR;
  }
  if (!creator.HasTemplates()) {
    mprintf("Warning: No residue templates loaded.\n");
  }
  if (!creator.HasMainParmSet()) {
    mprinterr("Error: No parameter sets.\n");
    return CpptrajState::ERR;
  }
  t_get_templates_.Stop();
  // FIXME hide behind ifdef?
  creator.TimingInfo(t_get_templates_.Total(), 2);

  // All residues start unknown
  Cpptraj::Structure::ResStatArray resStat( topIn.Nres() );
  std::vector<BondType> LeapBonds;

  // Disulfide search
  t_disulfide_.Start();
  if (!argIn.hasKey("nodisulfides")) {
    Cpptraj::Structure::Disulfide disulfide;
    if (disulfide.InitDisulfide( argIn, Cpptraj::Structure::Disulfide::ADD_BONDS, debug_ )) {
      mprinterr("Error: Could not init disulfide search.\n");
      return CpptrajState::ERR;
    }
    if (disulfide.SearchForDisulfides( resStat, topIn, frameIn, LeapBonds ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot searching for disulfides.\n");
  }
  t_disulfide_.Stop();

  // Handle sugars.
  t_sugar_.Start();
  // TODO should be on a residue by residue basis in FillAtomsWithTemplates
  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");
  if (sugarBuilder_ != 0) delete sugarBuilder_;
  sugarBuilder_ = 0;
  if (prepare_sugars) {
    sugarBuilder_ = new Cpptraj::Structure::SugarBuilder(debug_);
    // Init options
    if (sugarBuilder_->InitOptions( argIn.hasKey("hasglycam"),
                                    argIn.getKeyDouble("rescut", 8.0),
                                    argIn.getKeyDouble("bondoffset", 0.2),
                                    argIn.GetStringKey("sugarmask"),
                                    argIn.GetStringKey("determinesugarsby", "geometry"),
                                    argIn.GetStringKey("resmapfile") ))
    {
      mprinterr("Error: Sugar options init failed.\n");
      return CpptrajState::ERR;
    }
    bool splitres = !argIn.hasKey("nosplitres");
    if (splitres)
      mprintf("\tWill split off recognized sugar functional groups into separate residues.\n");
    else
      mprintf("\tNot splitting recognized sugar functional groups into separate residues.\n");
    bool c1bondsearch = !argIn.hasKey("noc1search");
    if (c1bondsearch)
      mprintf("\tWill search for missing bonds to sugar anomeric atoms.\n");
    else
      mprintf("\tNot searching for missing bonds to sugar anomeric atoms.\n");
    // May need to modify sugar structure/topology, either by splitting
    // C1 hydroxyls of terminal sugars into ROH residues, and/or by
    // adding missing bonds to C1 atoms.
    // This is done before any identification takes place since we want
    // to identify based on the most up-to-date topology.
    if (sugarBuilder_->FixSugarsStructure(topIn, frameIn,
                                          c1bondsearch, splitres, solventResName,
                                          LeapBonds))
    {
      mprinterr("Error: Sugar structure modification failed.\n");
      return CpptrajState::ERR;
    }
    if (sugarBuilder_->PrepareSugars(true, resStat, topIn, frameIn, LeapBonds))
    {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
  }
  t_sugar_.Stop();

  // Fill in atoms with templates
  t_fill_.Start();
  Topology topOut;
  topOut.SetDebug( debug_ );
  // TODO better default
  if (title.empty())
    title.assign( topIn.c_str() );
  topOut.SetParmName( title, FileName() );
  Frame frameOut;
  if (FillAtomsWithTemplates(topOut, frameOut, topIn, frameIn, creator)) {
    mprinterr("Error: Could not fill in atoms using templates.\n");
    return CpptrajState::ERR;
  }
  t_fill_.Stop();

  // Add the disulfide/sugar bonds
  if (transfer_bonds( topOut, topIn, LeapBonds )) {
    mprinterr("Error: Adding disulfide/sugar bonds failed.\n");
    return CpptrajState::ERR;
  }

  // Assign parameters. This will create the bond/angle/dihedral/improper
  // arrays as well.
  t_assign_.Start();
  Exec::RetType ret = CpptrajState::OK;
  Cpptraj::Parm::AssignParams AP;
  AP.SetDebug( debug_ );
  AP.SetVerbose( verbose );
  if ( AP.AssignParameters( topOut, *(creator.MainParmSetPtr()) ) ) {
    mprinterr("Error: Could not assign parameters for '%s'.\n", topOut.c_str());
    ret = CpptrajState::ERR;
  }
  // Assign GB parameters
  if (Cpptraj::Parm::Assign_GB_Radii( topOut, gbradii )) {
    mprinterr("Error: Could not assign GB parameters for '%s'\n", topOut.c_str());
    ret = CpptrajState::ERR;
  }
  // Create empty arrays for the TREE, JOIN, and IROTAT arrays
  topOut.AllocTreeChainClassification( );
  topOut.AllocJoinArray();
  topOut.AllocRotateArray();
  t_assign_.Stop();

  // Update coords 
  if (crdout.CoordsSetup( topOut, CoordinateInfo() )) { // FIXME better coordinate info
    mprinterr("Error: Could not set up output COORDS.\n");
    return CpptrajState::ERR;
  }
  crdout.SetCRD(0, frameOut);

  // Structure check
  if (check_structure_) {
    t_check_.Start();
    StructureCheck check;
    if (check.SetOptions( true, // image 
                          true, // check bonds
                          true,  // save problems
                          debug_, // debug
                          "*", // mask 1
                          "", // mask 2
                          0.8, // nonbond cut
                          1.15, // bond long offset
                          0.5, // bond short offset
                          -1, // pairlist cut (-1 for heuristic)
                          true, // ring check
                          0, // 0 = default ring check short distance cut
                          0, // 0 = default ring check long distance cut
                          0  // 0 = default ring check angle cut
        ))
    {
      mprinterr("Error: Structure check options failed.\n");
      return CpptrajState::ERR;
    }
    // For larger structures, automatically add a box
    bool box_added = false;
    if (!frameOut.BoxCrd().HasBox() && topOut.Natom() > check_box_natom_) {
      mprintf("\tAdding unit cell for check only.\n");
      box_added = true;
      // Get radii
      std::vector<double> Radii;
      Radii.reserve( topIn.Natom() );
      for (int atnum = 0; atnum != topIn.Natom(); ++atnum) {
        Radii.push_back( topOut.GetVDWradius(atnum) );
        //Radii.push_back( topIn[atnum].ParseRadius() );
        //Radii.push_back( 0.5 );
      }

      if (frameOut.SetOrthoBoundingBox(Radii, 1.0)) {
        mprinterr("Error: Setting orthogonal bounding box failed.\n");
        return CpptrajState::ERR;
      }
      frameOut.BoxCrd().PrintInfo();
    }
    if (check.Setup( topOut, frameOut.BoxCrd() )) {
      mprinterr("Error: Structure check setup failed.\n");
      return CpptrajState::ERR;
    }
    check.Mask1().MaskInfo();
    if (check.ImageOpt().ImagingEnabled())
      mprintf("\tImaging on.\n");
    else
      mprintf("\timaging off.\n");
    // TODO make file a user option
    CpptrajFile check_output;
    check_output.OpenWrite("");
    int Ntotal_problems = check.CheckOverlaps( frameOut );
    check.WriteProblemsToFile( &check_output, 1, topOut );
    Ntotal_problems += check.CheckBonds( frameOut );
    check.WriteProblemsToFile( &check_output, 1, topOut );
    Ntotal_problems += check.CheckRings( frameOut );
    check.WriteProblemsToFile( &check_output, 1, topOut );
    mprintf("\t%i total problems detected.\n", Ntotal_problems);
    // If box was added for check only, remove it
    if (box_added)
      frameOut.ModifyBox().SetNoBox();
    t_check_.Stop();
  }
  t_total_.Stop();

  t_total_.WriteTiming(1, "Build timing:");
  t_hisDetect_.WriteTiming    (2, "Histidine detection :", t_total_.Total());
  t_clean_.WriteTiming        (2, "Structure clean     :", t_total_.Total());
  t_get_templates_.WriteTiming(2, "Get templates/parms :", t_total_.Total());
  
  t_disulfide_.WriteTiming    (2, "Disulfide detection :", t_total_.Total());
  t_sugar_.WriteTiming        (2, "Sugar preparation   :", t_total_.Total());
  t_fill_.WriteTiming         (2, "Fill missing atoms  :", t_total_.Total());
  t_assign_.WriteTiming       (2, "Param./Top. gen.    :", t_total_.Total());
  t_check_.WriteTiming        (2, "Structure check     :", t_total_.Total());

  return ret;
}
