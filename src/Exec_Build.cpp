#include "Exec_Build.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "Structure/GenerateAngles.h"
#include "Structure/Zmatrix.h"

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
  std::vector<Vec3> XYZ; // FIXME should use frameOut
  Cpptraj::Structure::Zmatrix::Barray hasPosition;
  int nAtomsMissing = 0;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    // Identify a template based on the residue name.
    DataSet_Coords* resTemplate = IdTemplateFromName(Templates, topIn.Res(ires).Name());
    if (resTemplate == 0) {
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameNum(ires).c_str());
    } else {
      mprintf("\tTemplate %s being used for residue %s\n",
              resTemplate->legend(), topIn.TruncResNameNum(ires).c_str());
      // Map atoms to template atoms
      std::vector<int> map = MapAtomsToTemplate( topIn, ires, resTemplate );
      mprintf("\t  Atom map:\n");
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        mprintf("\t\t%6i %6s =>", iref+1, *(resTemplate->Top()[iref].Name()));
        if (map[iref] == -1)
          mprintf(" No match\n");
        else
          mprintf(" %6i %6s\n", map[iref]+1, *(topIn[map[iref]].Name()));
      }
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        topOut.AddTopAtom( resTemplate->Top()[iref], topIn.Res(ires) );
        if (map[iref] == -1) {
          XYZ.push_back( Vec3(0.0) );
          hasPosition.push_back( false );
          nAtomsMissing++;
        } else {
          XYZ.push_back( Vec3(frameIn.XYZ(map[iref])) );
          hasPosition.push_back( true );
        }
      }
      // DEBUG
      Frame templateFrame = resTemplate->AllocateFrame();
      resTemplate->GetFrame( 0, templateFrame );
      Cpptraj::Structure::Zmatrix zmatrix;
      if (zmatrix.SetFromFrame( templateFrame, resTemplate->Top(), 0 )) {
        mprinterr("Error: Could not set up residue template zmatrix.\n");
        return 1;
      }
      zmatrix.print();
      // If no atoms missing just fill in the residue
/*      if (nAtomsMissing == 0) {
        for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
          topOut.AddTopAtom( resTemplate->Top()[iref], topIn.Res(ires) );
          XYZ.push_back( Vec3(frameIn.XYZ(map[iref])) );
        }
      } else {
        mprintf("\tTrying to fill in missing atoms.\n");
        // one or more atoms missing. Try to use Zmatrix to fill it in
        
         // DEBUG*/

    }
  }
  mprintf("\t%i atoms missing.\n", nAtomsMissing);
  for (int iat = 0; iat != topOut.Natom(); iat++)
  {
    Residue const& res = topOut.Res( topOut[iat].ResNum() );
    mprintf("%6i %6s %6i %6s (%i) %g %g %g\n",
            iat+1, *(topOut[iat].Name()), res.OriginalResNum(), *(res.Name()),
            (int)hasPosition[iat], XYZ[iat][0], XYZ[iat][1], XYZ[iat][2]);
  }

  return 0;
}

// Exec_Build::Help()
void Exec_Build::Help() const
{
  mprintf("\tcrdset <COORDS set> [frame <#>] [parmset <param set> ...]\n");
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
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
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
  Topology& topIn = *(coords.TopPtr());

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

  // Generate angles/dihedrals
  if (Cpptraj::Structure::GenerateAngles(topIn)) {
    mprinterr("Error: Could not generate angles/dihedrals for '%s'\n", topIn.c_str());
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

  // Update parameters
  Exec::RetType ret = CpptrajState::OK;
  if ( topIn.UpdateParams( *mainParmSet  ) ) {
    mprinterr("Error: Could not update parameters for '%s'.\n", topIn.c_str());
    ret = CpptrajState::ERR;
  }

  if (free_parmset_mem && mainParmSet != 0) delete mainParmSet;

  return ret;
}
