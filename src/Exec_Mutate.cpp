#include "Exec_Mutate.h"
#include "CpptrajStdio.h"
#include "Structure/Creator.h"

// Exec_Mutate::Help()
void Exec_Mutate::Help() const
{
  mprintf("\tcrdset <COORDS set> resmask <mask>\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_);
}

// Exec_Mutate::Execute()
Exec::RetType Exec_Mutate::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  Cpptraj::Structure::Creator creator( debug_ );
  if (creator.InitCreator(argIn, State.DSL(), debug_)) {
    return CpptrajState::ERR;
  }
  if (!creator.HasTemplates()) {
    mprinterr("Error: No residue templates loaded.\n");
  }

  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS dataset name with 'crdset'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: No COORDS set with name %s found.\n", setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());

  std::string resmask = argIn.GetStringKey("resmask");
  if (resmask.empty()) {
    mprinterr("Error: Specify mask of residues to mutate with 'resmask'\n");
    return CpptrajState::ERR;
  }

  std::string templateName = argIn.GetStringKey("to");
  if (templateName.empty()) {
    mprinterr("Error: Specify template name to mutate to with 'to'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* UNIT = creator.IdTemplateFromName( templateName );
  if (UNIT == 0) {
    mprinterr("Error: Could not get template for '%s'\n", templateName.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tMutate residues selected by '%s' to '%s'\n", resmask.c_str(), UNIT->legend());

  return CpptrajState::OK;
}
