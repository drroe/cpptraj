#include "Exec_Sequence.h"
#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h" // For casting DataSet_Parameters to ParameterSet
#include "Structure/Builder.h"
#include "Structure/Creator.h"

/** Generate and build the specified sequence. */
int Exec_Sequence::generate_sequence(DataSet_Coords* OUT,
                                     DataSetList const& DSL,
                                     Sarray const& main_sequence,
                                     Cpptraj::Structure::Creator const& creator)
const
{
  // First, get all units in order.
  typedef std::vector<DataSet_Coords*> Uarray;
  Uarray Units;
  Units.reserve( main_sequence.size() );

  int total_natom = 0;

  for (Sarray::const_iterator it = main_sequence.begin(); it != main_sequence.end(); ++it)
  {
    DataSet_Coords* unit = creator.IdTemplateFromName( *it );

    if (unit == 0) {
      mprinterr("Error: Unit '%s' not found.\n", it->c_str());
      return 1;
    }
    if (unit->Size() < 1) {
      mprinterr("Error: Unit '%s' is empty.\n", unit->legend());
      return 1;
    }
    if (unit->Size() > 1) {
      mprintf("Warning: Unit '%s' has more than 1 frame. Only using first frame.\n", unit->legend());
    }
    Units.push_back( unit );
    total_natom += unit->Top().Natom();
  } // END loop over sequence
  mprintf("\tFound %zu units.\n", Units.size());
  if (Units.empty()) {
    mprinterr("Error: No units.\n");
    return 1;
  }
//  for (unsigned int idx = 0; idx < Units.size(); idx++)
//    mprintf("\tUnit %s HEAD %i TAIL %i\n", Units[idx]->legend(), connectAt0[idx]+1, connectAt1[idx]+1);

  Topology topOut;
  topOut.SetDebug( debug_ );
  topOut.SetParmName( OUT->Meta().Name(), FileName() );
  Frame frameOut;
  mprintf("\tFinal structure should have %i atoms.\n", total_natom);
  frameOut.SetupFrame( total_natom );
  // Clear frame so that AddXYZ can be used
  frameOut.ClearAtoms();

  // hasPosition - for each atom in topOut, status on whether atom in frameOut needs building
  Cpptraj::Structure::Builder::Barray hasPosition;
  hasPosition.reserve( total_natom );

  // Hold atom offsets needed when building residues
  typedef std::vector<int> Iarray;
  Iarray AtomOffsets;
  AtomOffsets.reserve( Units.size() ); // FIXME assuming 1 residue per unit

  // For holding bonded atom pairs
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> IParray;

  // Loop for setting up atoms in the topology from units
  IParray interResBonds;
  int prevTailAtom = -1;
  for (unsigned int idx = 0; idx < Units.size(); idx++)
  {
    int atomOffset = topOut.Natom();
    int resOffset = topOut.Nres();
    DataSet_Coords* unit = Units[idx];
    // Needs to have connect associated data
    AssociatedData* ad = unit->GetAssociatedData(AssociatedData::CONNECT);
    if (ad == 0) {
      mprinterr("Error: Unit '%s' does not have CONNECT data.\n", unit->legend());
      return 1;
    }
    AssociatedData_Connect const& CONN = static_cast<AssociatedData_Connect const&>( *ad );
    if (CONN.NconnectAtoms() < 2) {
      mprinterr("Error: Not enough connect atoms in unit '%s'\n", unit->legend());
      return 1;
    }
    int headAtom = CONN.Connect()[0] + atomOffset;
    int tailAtom = CONN.Connect()[1] + atomOffset;
    mprintf("\tAdding atoms for unit %s (head %i tail %i)\n", unit->legend(), headAtom+1, tailAtom+1);
    //mprintf("DEBUG: atom offset is %i\n", atomOffset);
    // Add the unit atoms. Only the first unit has known position.
    bool atomPosKnown = (idx == 0);
    Frame unitFrm = unit->AllocateFrame();
    unit->GetFrame(0, unitFrm);
    IParray intraResBonds;
    for (int itgt = 0; itgt < unit->Top().Natom(); itgt++)
    {
      Atom sourceAtom = unit->Top()[itgt];
      Residue currentRes = unit->Top().Res( sourceAtom.ResNum() );
      currentRes.SetOriginalNum( currentRes.OriginalResNum() + resOffset );
      // Save the intra-residue bonds
      int at0 = itgt + atomOffset;
      for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
        int at1 = *bat + atomOffset;
        if (at1 > at0) {
          if (debug_ > 1)
            mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, itgt+1, *bat + 1);
          intraResBonds.push_back( Ipair(at0, at1) );
        }
      }
      sourceAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
      topOut.AddTopAtom( sourceAtom, currentRes );
      frameOut.AddVec3( Vec3(unitFrm.XYZ(itgt)) );
      hasPosition.push_back( atomPosKnown );
    }
    // Add intra-residue bonds
    for (IParray::const_iterator it = intraResBonds.begin(); it != intraResBonds.end(); ++it)
    {
      //mprintf("DEBUG: Intra-res bond: Res %s atom %s to res %s atom %s\n",
      //        topOut.TruncResNameOnumId(topOut[it->first].ResNum()).c_str(), *(topOut[it->first].Name()),
      //        topOut.TruncResNameOnumId(topOut[it->second].ResNum()).c_str(), *(topOut[it->second].Name()));
      topOut.AddBond(it->first, it->second);
    }
    // Connect HEAD atom of this residue to TAIL of previous residue
    if (idx > 0) {
      if (prevTailAtom < 0 || headAtom < 0) {
        mprinterr("Error: Could not find connect atoms for previous residue (%i) and/or this residue (%i)\n",
                  prevTailAtom+1, headAtom+1);
        return 1;
      }
      if (debug_ > 1)
        mprintf("Will add bond between %i and %i\n", prevTailAtom+1, headAtom+1);
      // To preserve compat. with LEaP, make first atom the head atom.
      interResBonds.push_back( Ipair(headAtom, prevTailAtom) );
    } else
      // Placeholder
      interResBonds.push_back( Ipair(-1, -1) );
    prevTailAtom = tailAtom;

    // All units after the first need building
    if (atomPosKnown)
      AtomOffsets.push_back( -1 );
    else
      AtomOffsets.push_back( atomOffset );
  } // END loop over units

  // Build
  bool buildFailed = false;
  for (Iarray::const_iterator it = AtomOffsets.begin(); it != AtomOffsets.end(); ++it)
  {
    long int ires = it-AtomOffsets.begin();
    if (*it > -1) {
      // Residue has atom offset which indicates it needs something built.
      Cpptraj::Structure::Builder structureBuilder;// = new Cpptraj::Structure::Builder();
      structureBuilder.SetDebug( debug_ );
      if (creator.HasMainParmSet())
        structureBuilder.SetParameters( creator.MainParmSetPtr() );
      // Generate internals from the template, update indices to this topology.
      DataSet_Coords* unit = Units[ires];
      Frame unitFrm = unit->AllocateFrame();
      unit->GetFrame(0, unitFrm);
      if (structureBuilder.GenerateInternals(unitFrm, unit->Top(),
                                              std::vector<bool>(unit->Top().Natom(), true)))
      {
        mprinterr("Error: Generate internals for unit failed.\n");
        return 1;
      }
      structureBuilder.UpdateIndicesWithOffset( *it );

      mprintf("DEBUG: ***** BUILD unit %li %s *****\n", ires + 1,
              topOut.TruncResNameOnumId(ires).c_str());

      // Connect unit
      mprintf("DEBUG: Linking atoms %s and %s\n",
              topOut.AtomMaskName(interResBonds[ires].first).c_str(),
              topOut.AtomMaskName(interResBonds[ires].second).c_str());
      topOut.AddBond( interResBonds[ires].first, interResBonds[ires].second );
      // Generate internals around the link
      if (structureBuilder.GenerateInternalsAroundLink(interResBonds[ires].first, interResBonds[ires].second,
                                                       frameOut, topOut, hasPosition,
                                                       Cpptraj::Structure::Builder::SEQUENCE))
      {
        mprinterr("Error: Assign torsions around inter-unit link %s - %s failed.\n",
                  topOut.AtomMaskName(interResBonds[ires].first).c_str(),
                  topOut.AtomMaskName(interResBonds[ires].second).c_str());
        return 1;
      }
      // Update internal coords from known positions
      // NOTE: By defintion, there are no known positions.
      //if (structureBuilder.UpdateICsFromFrame( frameOut, topOut, hasPosition )) {
      //  mprinterr("Error: Failed to update Zmatrix with values from existing positions.\n");
      //  return 1;
      //}
      if (structureBuilder.BuildSequenceFromInternals(frameOut, topOut, hasPosition,
                                                      interResBonds[ires].first,
                                                      interResBonds[ires].second))
      {
        mprinterr("Error: Building residue %s failed.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
        buildFailed = true;
      }
      //delete structureBuilder;
    }
  }

  // Finalize topology - determine molecules, dont renumber residues, assign default bond params
  topOut.CommonSetup(true, false, true);
  topOut.Summary();

  OUT->CoordsSetup(topOut, frameOut.CoordsInfo());
  OUT->AddFrame( frameOut );

  if (buildFailed) return 1;
/*
  Topology combinedTop;
  combinedTop.SetDebug( debug_ );
  combinedTop.AppendTop( Units.front()->Top() );
  //combinedTop.SetParmBox( Units // TODO
  combinedTop.Brief("Sequence topology:");

  Frame CombinedFrame = Units.front()->AllocateFrame();
  Units.front()->GetFrame(0, CombinedFrame);

  using namespace Cpptraj::Structure;
  Builder builder;
  builder.SetDebug( debug_ );
  for (unsigned int idx = 1; idx < Units.size(); idx++) {
    mprintf("\tConnect %s atom %i to %s atom %i\n",
            Units[idx-1]->legend(), connectAt1[idx-1]+1,
            Units[idx]->legend(),   connectAt0[idx]  +1);
    Frame mol1frm = Units[idx]->AllocateFrame();
    Units[idx]->GetFrame(0, mol1frm);
    int bondat0 = connectAt1[idx-1];
    int bondat1 = connectAt0[idx];
    if (bondat0 < 0 || bondat1 < 0) {
      mprinterr("Error: Invalid connect atom(s) between %s atom %i to %s atom %i\n",
                Units[idx-1]->legend(), bondat0+1, Units[idx]->legend(), bondat1+1);
      return 1;
    }
    if (builder.Combine( combinedTop, CombinedFrame, Units[idx]->Top(), mol1frm,
                         connectAt1[idx-1], connectAt0[idx] )) {
      mprinterr("Error: Sequence combine between units %u %s and %u %s failed.\n",
                idx, Units[idx-1]->legend(), idx+1, Units[idx]->legend());
      return 1;
    }
  }

  // Generate angles and dihedrals
  //if (combinedTop.Nbonds() > 0) {
  //  if (Cpptraj::Structure::GenerateBondAngleTorsionArrays(combinedTop)) {
  //    mprinterr("Error: Angle generation failed.\n");
  //    return 1;
  //  }
  //}

  OUT->CoordsSetup(combinedTop, CombinedFrame.CoordsInfo());
  OUT->AddFrame( CombinedFrame );*/

  return 0;
}

// Exec_Sequence::Help()
void Exec_Sequence::Help() const
{
  mprintf("\tname <output set name> <unit0> <unit1> ...\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n"
          "  Create a molecule from a sequence of units.\n",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_);
}


// Exec_Sequence::Execute()
Exec::RetType Exec_Sequence::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();

  // Args
  Cpptraj::Structure::Creator creator(debug_);
  if (creator.InitCreator(argIn, State.DSL(), debug_)) {
    return CpptrajState::ERR;
  }

/*
  Sarray LibSetNames;
  std::string libsetname = argIn.GetStringKey("libset");
  while (!libsetname.empty()) {
    LibSetNames.push_back( libsetname );
    libsetname = argIn.GetStringKey("libset");
  } */
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: No output set name specified with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* OUT = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, MetaData(dsname));
  if (OUT == 0) {
    mprinterr("Error: Could not create output COORDS set named '%s'\n", dsname.c_str());
    return CpptrajState::ERR;
  }

  // Get the actual sequence from remaining args.
  Sarray main_sequence;
  ArgList remaining = argIn.RemainingArgs();
  std::string unit = remaining.GetStringNext();
  while (!unit.empty()) {
    main_sequence.push_back( unit );
    unit = remaining.GetStringNext();
  }
  if (main_sequence.empty()) {
    mprinterr("Error: No units specified.\n");
    return CpptrajState::ERR;
  }

  // Info
/*  if (!LibSetNames.empty()) {
    mprintf("\tLibrary set names:");
    for (Sarray::const_iterator it = LibSetNames.begin(); it != LibSetNames.end(); ++it)
      mprintf(" %s", it->c_str());
    mprintf("\n");
  }*/
  mprintf("\tMain sequence:");
  for (Sarray::const_iterator it = main_sequence.begin(); it != main_sequence.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");
  mprintf("\tOutput set name : %s\n", OUT->legend());

  // Execute
  if (generate_sequence(OUT, State.DSL(), main_sequence, creator)) {
    mprinterr("Error: Could not generate sequence.\n");
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
