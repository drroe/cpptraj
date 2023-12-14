#include "Exec_Build.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "Structure/GenerateAngles.h"

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
  typedef std::vector<DataSet_Coords*> Carray;
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
