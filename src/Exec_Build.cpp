#include "Exec_Build.h"
#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h" // For casting DataSet_Parameters to ParameterSet
#include "Parm/GB_Params.h"
#include "Structure/Builder.h"
#include "Structure/Creator.h"
#include "Structure/PdbCleaner.h"

/** Map atoms in residue to template. */
std::vector<int> Exec_Build::MapAtomsToTemplate(Topology const& topIn,
                                                int rnum,
                                                DataSet_Coords* resTemplate,
                                                Cpptraj::Structure::Creator const& creator)
{
  std::vector<int> mapOut;
  mapOut.reserve( resTemplate->Top().Natom() );
  Residue const& resIn = topIn.Res(rnum);
  // For each atom in topIn, find a template atom
  for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++)
  {
    mprintf("DEBUG: Search for atom %s\n", *(topIn[itgt].Name()));
    // Did this atom have an alias
    NameType alias;
    if (creator.GetAlias( alias, topIn[itgt].Name() )) {
      mprintf("DEBUG: Atom %s alias is %s\n",
              *(topIn[itgt].Name()),
              *alias);
    }
  }

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
  }
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
  resHeadAtoms.reserve( topIn.Nres() );
  resTailAtoms.reserve( topIn.Nres() );
  // Array of templates for each residue
  std::vector<DataSet_Coords*> ResTemplates;
  ResTemplates.reserve( topIn.Nres() );
  //typedef std::vector<Cpptraj::Structure::TerminalType> TermTypeArray;
  //TermTypeArray ResTypes;
  //ResTypes.reserve( topIn.Nres() );
  // Initial loop to try to match residues to templates
  int newNatom = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& currentRes = topIn.Res(ires);
    if (debug_ > 0)
      mprintf("DEBUG: ---------- Processing Residue %s ---------- \n", topIn.TruncResNameNum(ires).c_str());
    int pres = ires - 1;
    int nres = ires + 1;
    // Determine if this is a terminal residue
    Cpptraj::Structure::TerminalType resTermType;
    if (ires == 0 && topIn.Nres() > 1) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
    } else if (currentRes.IsTerminal()) {
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else if (pres > -1 && topIn.Res(pres).IsTerminal()) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
    } else if (nres < topIn.Nres() && (topIn.Res(nres).ChainId() != currentRes.ChainId() ||
                                       topIn.Res(nres).Icode()   != currentRes.Icode()))
    {
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else if (nres == topIn.Nres()) {
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else {
      resTermType = Cpptraj::Structure::NON_TERMINAL;
    }
    if (debug_ > 0)
      mprintf("DEBUG: Residue type: %s terminal\n", Cpptraj::Structure::terminalStr(resTermType));
    // Identify a template based on the residue name.
    DataSet_Coords* resTemplate = creator.IdTemplateFromResname(currentRes.Name(), resTermType);
    if (resTemplate == 0) {
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
      newNatom += currentRes.NumAtoms();
      // Head and tail atoms are blank
      resHeadAtoms.push_back( -1 );
      resTailAtoms.push_back( -1 );
    } else {
      if (debug_ > 0)
        mprintf("\tTemplate %s being used for residue %s\n",
                resTemplate->legend(), topIn.TruncResNameOnumId(ires).c_str());
      // Save the head and tail atoms
      AssociatedData* ad = resTemplate->GetAssociatedData(AssociatedData::CONNECT);
      if (ad == 0) {
        mprintf("Warning: Unit '%s' does not have CONNECT data.\n", resTemplate->legend());
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
  typedef std::pair<int,NameType> ResAtPair;
  typedef std::vector<ResAtPair> ResAtArray;
  ResAtArray detectedInterResBonds;

  // Loop for setting up atoms in the topology from residues or residue templates.
  int nRefAtomsMissing = 0;
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
      std::vector<int> map = MapAtomsToTemplate( topIn, ires, resTemplate, creator );
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
          pdb[itgt-currentRes.FirstAtom()] = iref;
          // Check source atoms for inter-residue connections.
          Atom const& sourceAtom = topIn[itgt];
          for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
            if ( topIn[*bat].ResNum() < ires ) {
              detectedInterResBonds.push_back( ResAtPair(ires, sourceAtom.Name()) );
              detectedInterResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), topIn[*bat].Name()) );
            }
          }
        }
      }
      // See if any current atoms were not mapped to reference
      int nTgtAtomsMissing = 0;
      for (int itgt = 0; itgt != currentRes.NumAtoms(); itgt++)
        if (pdb[itgt] == -1)
          nTgtAtomsMissing++;
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

  // -----------------------------------
  // DEBUG - Print primary connection atoms
  if (debug_ > 0) {
    for (unsigned int idx = 0; idx != ResTemplates.size(); idx++) {
      if (ResTemplates[idx] != 0) {
        mprintf("DEBUG: Template %s", ResTemplates[idx]->legend());
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
      if (resTailAtoms[pres] != -1) {
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
        mprintf("\tAdding non-template bond %s - %s\n",
                topOut.AtomMaskName(at0).c_str(),
                topOut.AtomMaskName(at1).c_str());
        resBondingAtoms[ra0.first].push_back( Ipair(at0, at1) );
        resBondingAtoms[ra1.first].push_back( Ipair(at1, at0) );
        ResidueConnections[ra0.first].push_back( ra1.first );
        ResidueConnections[ra1.first].push_back( ra0.first );
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

// -----------------------------------------------------------------------------
// Exec_Build::Help()
void Exec_Build::Help() const
{
  mprintf("\tname <output COORDS> crdset <COORDS set> [frame <#>]\n"
          "\t[title <title>]\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n"
          "  Build complete topology and parameters from given crdset.\n",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_);
}

// Exec_Build::Execute()
Exec::RetType Exec_Build::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  std::string title = argIn.GetStringKey("title");
/*  // Atom scan direction
  std::string atomscandir = argIn.GetStringKey("atomscandir");
  if (!atomscandir.empty()) {
    if (atomscandir == "f")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_FORWARDS);
    else if (atomscandir == "b")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_BACKWARDS);
    else {
      mprinterr("Error: Unrecognized keyword for 'atomscandir' : %s\n", atomscandir.c_str());
      return CpptrajState::ERR;
    }
  }*/
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
  Cpptraj::Structure::PdbCleaner pdbCleaner;
  pdbCleaner.SetDebug( debug_ );
  if (pdbCleaner.InitPdbCleaner( argIn, "HOH", std::vector<int>() )) {
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

  // Output coords
  std::string outset = argIn.GetStringKey("name");
  if (outset.empty()) {
    mprinterr("Error: Must specify output COORDS set with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet* outCrdPtr = State.DSL().AddSet( DataSet::COORDS, outset );
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

/*  Carray Templates;
  std::string lib = argIn.GetStringKey("lib");
  if (lib.empty()) {
    mprintf("\tNo template(s) specified with 'lib'; using any loaded templates.\n");
    DataSetList sets = State.DSL().SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      // Should only be a single residue FIXME need new set type
      DataSet_Coords const& ds = static_cast<DataSet_Coords const&>( *(*it) );
      if ( ds.Top().Nres() == 1 )
        Templates.push_back( (DataSet_Coords*)(*it) );
    }
  } else {
    while (!lib.empty()) {
      DataSetList sets = State.DSL().SelectGroupSets( lib, DataSet::COORDINATES ); // TODO specific set type for units?
      if (sets.empty()) {
        mprintf("Warning: No sets corresponding to '%s'\n", lib.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
        {
          // Should only be a single residue FIXME need new set type
          DataSet_Coords const& ds = static_cast<DataSet_Coords const&>( *(*it) );
          if ( ds.Top().Nres() == 1 )
            Templates.push_back( (DataSet_Coords*)(*it) );
        }
      }
      lib = argIn.GetStringKey("lib");
    }
  }
  if (Templates.empty())
    mprintf("Warning: No residue templates loaded.\n");
  else {
    mprintf("\t%zu residue templates found:", Templates.size());
    for (std::vector<DataSet_Coords*>::const_iterator it = Templates.begin(); it != Templates.end(); ++it)
      mprintf(" %s", (*it)->legend());
    mprintf("\n");
  }*/

  // Get templates and parameter sets.
  Cpptraj::Structure::Creator creator;
  if (creator.InitCreator(argIn, State.DSL(), debug_)) {
    return CpptrajState::ERR;
  }
  if (!creator.HasTemplates()) {
    mprintf("Warning: No residue templates loaded.\n");
  }
  if (!creator.HasMainParmSet()) {
    mprinterr("Error: No parameter sets.\n");
    return CpptrajState::ERR;
  }
/*
  typedef std::vector<DataSet_Parameters*> Parray;
  Parray ParamSets;
  std::string parmset = argIn.GetStringKey("parmset");
  if (parmset.empty()) {
    mprintf("\tNo parameter set(s) specified with 'parmset'; using any loaded sets.\n");
    // See if there are any parameter sets.
    DataSetList sets = State.DSL().GetSetsOfType( "*", DataSet::PARAMETERS );
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
      ParamSets.push_back( (DataSet_Parameters*)(*it) );
  } else {
    while (!parmset.empty()) {
      DataSetList sets = State.DSL().GetSetsOfType( parmset, DataSet::PARAMETERS );
      if (sets.empty()) {
        mprintf("Warning: No parameter sets corresponding to '%s'\n", parmset.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
          ParamSets.push_back( (DataSet_Parameters*)(*it) );
      }
      parmset = argIn.GetStringKey("parmset");
    }
  }
  if (ParamSets.empty()) {
    mprinterr("Error: No parameter sets.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tParameter sets:\n");
  for (Parray::const_iterator it = ParamSets.begin(); it != ParamSets.end(); ++it)
    mprintf("\t  %s\n", (*it)->legend());

  // Combine parameters if needed
  DataSet_Parameters* mainParmSet = 0;
  bool free_parmset_mem = false;
  if (ParamSets.size() == 1)
    mainParmSet = ParamSets.front();
  else {
    free_parmset_mem = true;
    mprintf("\tCombining parameter sets.\n");
    Parray::const_iterator it = ParamSets.begin();
    mainParmSet = new DataSet_Parameters( *(*it) );
    ++it;
    ParameterSet::UpdateCount UC;
    for (; it != ParamSets.end(); ++it)
      mainParmSet->UpdateParamSet( *(*it), UC, State.Debug(), State.Debug() ); // FIXME verbose
  }*/

  // Fill in atoms with templates
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


  // Assign parameters. This will create the bond/angle/dihedral/improper
  // arrays as well.
  Exec::RetType ret = CpptrajState::OK;
  if ( topOut.AssignParams( *(creator.MainParmSetPtr())  ) ) {
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

  //if (free_parmset_mem && mainParmSet != 0) delete mainParmSet;

  // Update coords 
  if (crdout.CoordsSetup( topOut, CoordinateInfo() )) { // FIXME better coordinate info
    mprinterr("Error: Could not set up output COORDS.\n");
    return CpptrajState::ERR;
  }
  crdout.SetCRD(0, frameOut);

  return ret;
}
