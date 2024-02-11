#include "Exec_Build.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "Structure/Builder.h"
#include "Structure/GenerateConnectivityArrays.h"
#include "Structure/Zmatrix.h"
#include "Parm/GB_Params.h"
#include "AssociatedData_ResId.h"

/** Try to identify residue template DataSet from the given residue
  * name (from e.g. the PDB/Mol2/etc file).
  */
DataSet_Coords* Exec_Build::IdTemplateFromName(Carray const& Templates,
                                               NameType const& rname,
                                               Cpptraj::Structure::TerminalType termType)
{
  DataSet_Coords* out = 0;
  if (termType != Cpptraj::Structure::NON_TERMINAL) {
    // Looking for a terminal residue. Need to get sets with AssociatedData_ResId
    for (Carray::const_iterator it = Templates.begin(); it != Templates.end(); ++it) {
      AssociatedData* ad = (*it)->GetAssociatedData( AssociatedData::RESID );
      if (ad != 0) {
        AssociatedData_ResId const& resid = static_cast<AssociatedData_ResId const&>( *ad );
        if (rname == resid.ResName() && termType == resid.TermType()) {
          out = *it;
          break;
        }
      }
    }
  }
  if (out == 0) {
    // Terminal residue not found or non-terminal residue.
    if (termType != Cpptraj::Structure::NON_TERMINAL)
      mprintf("Warning: No terminal residue found for '%s'\n", *rname);
    // Assume Coords set aspect is what we need
    for (Carray::const_iterator it = Templates.begin(); it != Templates.end(); ++it) {
      if ( rname == NameType( (*it)->Meta().Aspect() ) ) {
        out = *it;
        break;
      }
    }
  }

  return out;
}

/** Map atoms in residue to template. */
std::vector<int> Exec_Build::MapAtomsToTemplate(Topology const& topIn,
                                                int rnum,
                                                DataSet_Coords* resTemplate)
{
  std::vector<int> mapOut;
  mapOut.reserve( resTemplate->Top().Natom() );
  Residue const& resIn = topIn.Res(rnum);
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

/** Use given templates to construct a final molecule. */
int Exec_Build::FillAtomsWithTemplates(Topology& topOut, Frame& frameOut,
                                       Carray const& Templates,
                                       Topology const& topIn, Frame const& frameIn,
                                       ParameterSet const& mainParmSet)
{
  // Array of templates for each residue
  std::vector<DataSet_Coords*> ResTemplates;
  ResTemplates.reserve( topIn.Nres() );
  typedef std::vector<Cpptraj::Structure::TerminalType> TermTypeArray;
  TermTypeArray ResTypes;
  ResTypes.reserve( topIn.Nres() );
  // Initial loop to try to match residues to templates
  int newNatom = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& currentRes = topIn.Res(ires);
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
    } else if (nres < topIn.Nres() && topIn.Res(nres).ChainId() != currentRes.ChainId()) {
      resTermType = Cpptraj::Structure::END_TERMINAL;
    } else {
      resTermType = Cpptraj::Structure::NON_TERMINAL;
    }
    mprintf("DEBUG: Residue type: %s terminal\n", Cpptraj::Structure::terminalStr(resTermType));
    // Identify a template based on the residue name.
    DataSet_Coords* resTemplate = IdTemplateFromName(Templates, currentRes.Name(), resTermType);
    if (resTemplate == 0) {
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
      newNatom += currentRes.NumAtoms();
    } else {
      mprintf("\tTemplate %s being used for residue %s\n",
              resTemplate->legend(), topIn.TruncResNameOnumId(ires).c_str());
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
  Cpptraj::Structure::Zmatrix::Barray hasPosition;
  hasPosition.reserve( newNatom );

  // Hold atom offsets needed when building residues
  typedef std::vector<int> Iarray;
  Iarray AtomOffsets;
  AtomOffsets.reserve( topIn.Nres() );

  // For inter-residue bonding, use residue # and atom name since
  // atom numbering may change if atoms are added from templates.
  typedef std::pair<int,NameType> ResAtPair;
  typedef std::vector<ResAtPair> ResAtArray;
  ResAtArray interResBonds;

  // For holding bonded atom pairs
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> IParray;

  // Loop for setting up atoms in the topology from residues or residue templates.
  int nRefAtomsMissing = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    mprintf("\tAdding atoms for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
    int atomOffset = topOut.Natom();
    mprintf("DEBUG: atom offset is %i\n", atomOffset);
    Residue const& currentRes = topIn.Res(ires);
    DataSet_Coords* resTemplate = ResTemplates[ires];
    IParray intraResBonds;
    if (resTemplate == 0) {
      // ----- No template. Just add the atoms. ------------
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
            // Inter-residue. Only record if bonding to the next residue.
            if (topIn[*bat].ResNum() > ires) {
              interResBonds.push_back( ResAtPair(ires, sourceAtom.Name()) );
              interResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), topIn[*bat].Name()) );
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
      // Map source atoms to template atoms.
      std::vector<int> map = MapAtomsToTemplate( topIn, ires, resTemplate );
      mprintf("\t  Atom map:\n");
      // DEBUG - print map
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        mprintf("\t\t%6i %6s =>", iref+1, *(resTemplate->Top()[iref].Name()));
        if (map[iref] == -1)
          mprintf(" No match\n");
        else
          mprintf(" %6i %6s\n", map[iref]+1, *(topIn[map[iref]].Name()));
      }
      // Map template atoms back to source atoms.
      std::vector<int> pdb(currentRes.NumAtoms(), -1);
      bool atomsNeedBuilding = false;
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        // Track intra-residue bonds
        Atom templateAtom = resTemplate->Top()[iref];
        int at0 = iref + atomOffset;
        for (Atom::bond_iterator bat = templateAtom.bondbegin(); bat != templateAtom.bondend(); ++bat) {
          //if ( resTemplate->Top()[*bat].ResNum() == ires ) {
            int at1 = *bat + atomOffset;
            if (at1 > at0) {
              //mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, iref+1, *bat + 1);
              intraResBonds.push_back( Ipair(at0, at1) );
            }
          //}
        }
        // TODO check connect atoms for inter-residue connections
        templateAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
        topOut.AddTopAtom( templateAtom, currentRes );
        if (map[iref] == -1) {
          frameOut.AddVec3( Vec3(0.0) );
          hasPosition.push_back( false );
          nRefAtomsMissing++;
          atomsNeedBuilding = true;
        } else {
          int itgt = map[iref];
          frameOut.AddVec3( Vec3(frameIn.XYZ(itgt)) );
          hasPosition.push_back( true );
          pdb[itgt-currentRes.FirstAtom()] = iref;
          // Check source atoms for inter-residue connections
          Atom const& sourceAtom = topIn[itgt];
          for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
            if ( topIn[*bat].ResNum() > ires ) {
              interResBonds.push_back( ResAtPair(ires, sourceAtom.Name()) );
              interResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), topIn[*bat].Name()) );
            }
          }
        }
      }
      // See if any current atoms were not mapped to reference
      int nTgtAtomsMissing = 0;
      for (int itgt = 0; itgt != currentRes.NumAtoms(); itgt++)
        if (pdb[itgt] == -1)
          nTgtAtomsMissing++;
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
  // Check inter-residue bonds
  std::vector<IParray> resBondingAtoms(topOut.Nres());
  for (ResAtArray::const_iterator it = interResBonds.begin(); it != interResBonds.end(); ++it)
  {
    ResAtPair const& ra0 = *it;
    ++it;
    ResAtPair const& ra1 = *it;
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
    // Save bonding atoms
    resBondingAtoms[ra0.first].push_back( Ipair(at0, at1) );
    resBondingAtoms[ra1.first].push_back( Ipair(at1, at0) );
  }
  // DEBUG print residue bonding atoms
  for (std::vector<IParray>::const_iterator rit = resBondingAtoms.begin();
                                            rit != resBondingAtoms.end(); ++rit)
  {
    mprintf("\tResidue %s bonding atoms.\n", topOut.TruncResNameNum(rit-resBondingAtoms.begin()).c_str());
    for (IParray::const_iterator it = rit->begin(); it != rit->end(); ++it)
      mprintf("\t\t%s - %s\n", topOut.AtomMaskName(it->first).c_str(), topOut.AtomMaskName(it->second).c_str());
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
      // Residue has atom offset which indicates it needs something built.
      Cpptraj::Structure::Builder structureBuilder;// = new Cpptraj::Structure::Builder();
      structureBuilder.SetDebug( 1 ); // DEBUG FIXME
      structureBuilder.SetParameters( &mainParmSet );
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

      mprintf("DEBUG: ***** BUILD residue %li %s *****\n", ires + 1,
              topOut.TruncResNameOnumId(ires).c_str());
      //mprintf("DEBUG: Residue type: %s terminal\n", Cpptraj::Structure::terminalStr(*termType));
      // Is this residue connected to an earlier residue?
      for (IParray::const_iterator resBonds = resBondingAtoms[ires].begin();
                                   resBonds != resBondingAtoms[ires].end(); ++resBonds)
      {
        if (resBonds->second < resBonds->first) {
          mprintf("\t\tResidue connection: %s - %s\n",
                  topOut.AtomMaskName(resBonds->first).c_str(),
                  topOut.AtomMaskName(resBonds->second).c_str());
          topOut.AddBond(resBonds->first, resBonds->second);
          // Generate internals around the link
          if (structureBuilder.GenerateInternalsAroundLink(resBonds->first, resBonds->second,
                                                            frameOut, topOut, hasPosition))
          {
            mprinterr("Error: Assign torsions around inter-residue link %s - %s failed.\n",
                      topOut.AtomMaskName(resBonds->first).c_str(),
                      topOut.AtomMaskName(resBonds->second).c_str());
            return 1;
          }
        }
      }
      // Update internal coords from known positions
      if (structureBuilder.UpdateICsFromFrame( frameOut, ires, topOut, hasPosition )) {
        mprinterr("Error: Failed to update Zmatrix with values from existing positions.\n");
        return 1;
      }
      // Convert to Zmatrix and assign missing atom positions
      Cpptraj::Structure::Zmatrix tmpz;
      tmpz.SetDebug( 1 ); // DEBUG
      if (structureBuilder.GetZmatrixFromInternals(tmpz, topOut)) {
        mprinterr("Error: Could not get Zmatrix from internals.\n");
        return 1;
      }
      if (tmpz.SetToFrame( frameOut, hasPosition )) {
        mprinterr("Error: Building residue %s failed.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
        buildFailed = true;
      }// else
       // resIsBuilt[ires] = true;
      //delete structureBuilder;
    }
  } // END loop over atom offsets

  // DEBUG - Print new top/coords
  for (int iat = 0; iat != topOut.Natom(); iat++)
  {
    Residue const& res = topOut.Res( topOut[iat].ResNum() );
    const double* XYZ = frameOut.XYZ(iat);
    mprintf("%6i %6s %6i %6s (%i) %g %g %g\n",
            iat+1, *(topOut[iat].Name()), res.OriginalResNum(), *(res.Name()),
            (int)hasPosition[iat], XYZ[0], XYZ[1], XYZ[2]);
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
          "\t[parmset <param set> ...] [lib <lib set> ...]\n"
          "\t[atomscandir {f|b}]\n"
         );
}

// Exec_Build::Execute()
Exec::RetType Exec_Build::Execute(CpptrajState& State, ArgList& argIn)
{
  // Atom scan direction
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
  }
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
  Topology const& topIn = coords.Top();

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

  // Get residue templates.
  Carray Templates;
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
  }

  // Get parameter sets.
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
  }

  // Fill in atoms with templates
  Topology topOut;
  Frame frameOut;
  if (FillAtomsWithTemplates(topOut, frameOut, Templates, topIn, frameIn, *mainParmSet)) {
    mprinterr("Error: Could not fill in atoms using templates.\n");
    return CpptrajState::ERR;
  }


  // Assign parameters. This will create the bond/angle/dihedral/improper
  // arrays as well.
  Exec::RetType ret = CpptrajState::OK;
  if ( topOut.AssignParams( *mainParmSet  ) ) {
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

  if (free_parmset_mem && mainParmSet != 0) delete mainParmSet;

  // Update coords 
  if (crdout.CoordsSetup( topOut, CoordinateInfo() )) { // FIXME better coordinate info
    mprinterr("Error: Could not set up output COORDS.\n");
    return CpptrajState::ERR;
  }
  crdout.SetCRD(0, frameOut);

  return ret;
}
