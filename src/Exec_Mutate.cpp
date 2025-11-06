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

  AtomMask mask;
  if (mask.SetMaskString( resmask )) {
    mprinterr("Error: Could not set mask '%s'\n", resmask.c_str());
    return CpptrajState::ERR;
  }
  if (CRD->Top().SetupIntegerMask( mask )) {
    mprinterr("Error: Could not setup mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }
  if (mask.None()) {
    mprinterr("Error: Nothing selected by mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }
  //mask.MaskInfo();
  std::vector<int> resnums = CRD->Top().ResnumsSelectedBy( mask );
  mprintf("\t%zu residues selected by '%s'\n", resnums.size(), mask.MaskString());

  AtomMask toRemove;
  toRemove.SetNatoms( CRD->Top().Natom() );
  for (std::vector<int>::const_iterator rnum = resnums.begin(); rnum != resnums.end(); ++rnum)
  {
    Residue const& currentRes = CRD->Top().Res( *rnum );
    int atomsToRemove = 0;
    for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); at++)
    {
      // Does this atom exist in the template?
      int idx = UNIT->Top().FindAtomInResidue(0, CRD->Top()[at].Name());
      if (idx < 0) {
        toRemove.AddSelectedAtom( at );
        atomsToRemove++;
      }
      if (idx > -1)
        mprintf("DEBUG: Found atom %s in template.\n", CRD->Top().AtomMaskName(at).c_str());
      else
        mprintf("DEBUG: Atom %s not in template.\n", CRD->Top().AtomMaskName(at).c_str());
    }
    mprintf("DEBUG: Removing %i atoms from residue %i\n", atomsToRemove, *rnum + 1 );
    if (currentRes.NumAtoms() - atomsToRemove < 1) {
      mprinterr("Error: Number of atoms to remove %i >= number of atoms in residue %s (%i)\n",
                atomsToRemove, CRD->Top().TruncResNameNum(*rnum).c_str(), currentRes.NumAtoms());
      return CpptrajState::ERR;
    }
  }
  toRemove.InvertMask();

  Topology* newTop = CRD->Top().modifyStateByMask( toRemove );
  if (newTop == 0) {
    mprinterr("Error: Could not remove atoms from '%s'\n", CRD->legend());
    return CpptrajState::ERR;
  }
  newTop->Summary();

  if (newTop != 0) delete newTop;

  return CpptrajState::OK;
}
