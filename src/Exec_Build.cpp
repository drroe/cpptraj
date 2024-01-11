#include "Exec_Build.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "Structure/GenerateConnectivityArrays.h"
#include "Structure/Zmatrix.h"
#include "Parm/GB_Params.h"

DataSet_Coords* Exec_Build::IdTemplateFromName(Carray const& Templates,
                                               NameType const& rname)
{
  DataSet_Coords* out = 0;
  // Assume Coords set aspect is what we need
  for (Carray::const_iterator it = Templates.begin(); it != Templates.end(); ++it) {
    if ( rname == NameType( (*it)->Meta().Aspect() ) ) {
      out = *it;
      break;
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
                                       Topology const& topIn, Frame const& frameIn)
{
  typedef std::vector<int> Iarray;
  typedef std::vector<Iarray> IIarray;
  // Determine which residues are terminal
  IIarray molTermResidues( topIn.Nmol() );
  // Track residue connections
  IIarray resConnections( topIn.Nres() );
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& currentRes = topIn.Res(ires);
    int n_connections = 0;
    for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); at++)
    {
      //if (topIn[at].Element() != Atom::HYDROGEN && topIn[at].Nbonds() > 1) {
        for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
        {
          if (topIn[*bat].ResNum() != topIn[at].ResNum()) {
            n_connections++;
            resConnections[ ires ].push_back( topIn[*bat].ResNum() );
          }
        }
      //}
    }
    mprintf("DEBUG: Res %s has %i connections.\n", topIn.TruncResNameNum(ires).c_str(), n_connections);
    if (n_connections == 1) {
      int molnum = topIn[currentRes.FirstAtom()].MolNum();
      molTermResidues[ molnum ].push_back( ires );
    }
  }
  // DEBUG Print residue connections
  mprintf("DEBUG: Residue connections:\n");
  for (IIarray::const_iterator rit = resConnections.begin(); rit != resConnections.end(); ++rit)
  {
    mprintf("DEBUG:\t\t%s to", topIn.TruncResNameNum(rit-resConnections.begin()).c_str());
    for (Iarray::const_iterator it = rit->begin(); it != rit->end(); ++it)
      mprintf(" %s", topIn.TruncResNameNum( *it ).c_str());
    mprintf("\n");
  }
  // Set residue terminal status
  enum TermType { BEG_TERMINAL = 0, REGULAR, END_TERMINAL };
  typedef std::vector<TermType> Tarray;
  Tarray residueTerminalStatus( topIn.Nres(), REGULAR );
  mprintf("\tTerminal residues:\n");
  for (IIarray::const_iterator mit = molTermResidues.begin(); mit != molTermResidues.end(); ++mit)
  {
    mprintf("\t\tMol %li:", mit - molTermResidues.begin() + 1);
    for (Iarray::const_iterator it = mit->begin(); it != mit->end(); ++it) {
      mprintf(" %s", topIn.TruncResNameNum( *it ).c_str());
      if (it == mit->begin())
        residueTerminalStatus[*it] = BEG_TERMINAL;
      else
        residueTerminalStatus[*it] = END_TERMINAL;
      mprintf("(%i)", (int)residueTerminalStatus[*it]);
    }
    mprintf("\n");
  }

  // Array of templates for each residue
  std::vector<DataSet_Coords*> ResTemplates;
  ResTemplates.reserve( topIn.Nres() );
  // Initial loop to try to match residues to templates
  int newNatom = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    // Identify a template based on the residue name.
    Residue const& currentRes = topIn.Res(ires);
    DataSet_Coords* resTemplate = IdTemplateFromName(Templates, currentRes.Name());
    if (resTemplate == 0) {
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameNum(ires).c_str());
      newNatom += currentRes.NumAtoms();
    } else {
      mprintf("\tTemplate %s being used for residue %s\n",
              resTemplate->legend(), topIn.TruncResNameNum(ires).c_str());
      newNatom += resTemplate->Top().Natom();
    }
    ResTemplates.push_back( resTemplate );
  }
  mprintf("\tFinal structure should have %i atoms.\n", newNatom);
  frameOut.SetupFrame( newNatom );
  // Clear frame so that AddXYZ can be used
  frameOut.ClearAtoms();

  // hasPosition - for each atom in topOut, status on whether atom in frameOut needs building
  Cpptraj::Structure::Zmatrix::Barray hasPosition;
  hasPosition.reserve( newNatom );

  // Hold Z-matrices for residues that have missing atoms
  typedef std::vector<Cpptraj::Structure::Zmatrix*> Zarray;
  Zarray ResZmatrices;
  ResZmatrices.reserve( topIn.Nres() );

  // For inter-residue bonding, use residue # and atom name since
  // atom numbering may change if atoms are added from templates.
  typedef std::pair<int,NameType> ResAtPair;
  typedef std::vector<ResAtPair> ResAtArray;
  ResAtArray interResBonds;

  // Loop for setting up atoms in the topology from residues or residue templates.
  int nRefAtomsMissing = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    mprintf("\tAdding atoms for residue %s\n", topIn.TruncResNameNum(ires).c_str());
    int atomOffset = topOut.Natom();
    Residue const& currentRes = topIn.Res(ires);
    DataSet_Coords* resTemplate = ResTemplates[ires];
    // For holding bonded atom pairs
    typedef std::pair<int,int> Ipair;
    typedef std::vector<Ipair> IParray;
    IParray intraResBonds;
    if (resTemplate == 0) {
      // No template. Just add the atoms.
      for (int itgt = currentRes.FirstAtom(); itgt != currentRes.LastAtom(); ++itgt)
      {
        // Track intra-residue bonds
        Atom sourceAtom = topIn[itgt];
        int at0 = itgt - currentRes.FirstAtom() + atomOffset;
        for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
          if ( topIn[*bat].ResNum() == ires ) {
            int at1 = *bat - currentRes.FirstAtom() + atomOffset;
            if (at1 > at0) {
              mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, itgt+1, *bat + 1);
              intraResBonds.push_back( Ipair(at0, at1) );
            }
          } else {
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
      // A template exists for this residue.
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
              mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, iref+1, *bat + 1);
              intraResBonds.push_back( Ipair(at0, at1) );
            }
          //}
        }
        // TODO connect atoms for inter-residue connections
        templateAtom.ClearBonds();
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
      // Save zmatrix if atoms need to be built
      if (atomsNeedBuilding) {
        Frame templateFrame = resTemplate->AllocateFrame();
        resTemplate->GetFrame( 0, templateFrame );
        Cpptraj::Structure::Zmatrix* zmatrix = new Cpptraj::Structure::Zmatrix();
        if (zmatrix->SetFromFrame( templateFrame, resTemplate->Top(), 0 )) {
          mprinterr("Error: Could not set up residue template zmatrix.\n");
          return 1;
        }
        zmatrix->print( resTemplate->TopPtr() );
        zmatrix->OffsetIcIndices( atomOffset );
        ResZmatrices.push_back( zmatrix );
        zmatrix->print( &topOut );
      } else
        ResZmatrices.push_back( 0 );
    } // END template exists
    // Add intra-residue bonds
    for (IParray::const_iterator it = intraResBonds.begin(); it != intraResBonds.end(); ++it)
    {
      mprintf("DEBUG: Intra-res bond: Res %s atom %s to res %s atom %s\n",
              topOut.TruncResNameNum(topOut[it->first].ResNum()).c_str(), *(topOut[it->first].Name()),
              topOut.TruncResNameNum(topOut[it->second].ResNum()).c_str(), *(topOut[it->second].Name()));
      topOut.AddBond(it->first, it->second);
    }
  } // END loop over source residues
  // Add inter-residue bonds
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
    topOut.AddBond(at0, at1);
  }
  mprintf("\t%i template atoms missing in source.\n", nRefAtomsMissing);

  // Build using internal coords if needed.
  bool buildFailed = false;
  for (Zarray::const_iterator it = ResZmatrices.begin(); it != ResZmatrices.end(); ++it)
  {
    Cpptraj::Structure::Zmatrix* zmatrix = *it;
    if (zmatrix != 0) {
      // Create zmatrix with previous residue if needed
      
      // Update zmatrix seeds
      if (zmatrix->AutoSetSeedsWithPositions( frameOut, topOut, it - ResZmatrices.begin(), hasPosition )) {
        mprinterr("Error: Could not set up seed atoms for Zmatrix.\n");
        buildFailed = true;
      } else {
        mprintf("DEBUG: Zmatrix for building residue %li %s\n", it - ResZmatrices.begin() + 1,
                topOut.TruncResNameNum(it - ResZmatrices.begin()).c_str());
        zmatrix->print(&topOut);
        zmatrix->SetDebug( 1 ); // DEBUG
        if (zmatrix->SetToFrame( frameOut, hasPosition )) {
          mprinterr("Error: Building residue %s failed.\n",
                    topOut.TruncResNameNum(it - ResZmatrices.begin()).c_str());
          buildFailed = true;
        }
      }
    }
  }

  // DEBUG - Print new top/coords
  for (int iat = 0; iat != topOut.Natom(); iat++)
  {
    Residue const& res = topOut.Res( topOut[iat].ResNum() );
    const double* XYZ = frameOut.XYZ(iat);
    mprintf("%6i %6s %6i %6s (%i) %g %g %g\n",
            iat+1, *(topOut[iat].Name()), res.OriginalResNum(), *(res.Name()),
            (int)hasPosition[iat], XYZ[0], XYZ[1], XYZ[2]);
  }

  // Clean up zmatrices
  for (Zarray::iterator it = ResZmatrices.begin(); it != ResZmatrices.end(); ++it)
    if (*it != 0) delete *it;

  // Finalize topology
  topOut.CommonSetup(); // TODO dont assign default bond parameters here
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

  // Fill in atoms with templates
  Topology topOut;
  Frame frameOut;
  if (FillAtomsWithTemplates(topOut, frameOut, Templates, topIn, frameIn)) {
    mprinterr("Error: Could not fill in atoms using templates.\n");
    return CpptrajState::ERR;
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
