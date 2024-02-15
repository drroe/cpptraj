#include "Creator.h"
#include "GenerateConnectivityArrays.h" // For setting atom scan direction
#include "../ArgList.h"
#include "../AssociatedData_ResId.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h" // TODO new coords type
#include "../DataSet_Parameters.h"
#include "../DataSetList.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Creator::Creator() :
  mainParmSet_(0),
  debug_(0),
  free_parmset_mem_(false)
{}

/** DESTRUCTOR */
Creator::~Creator() {
  if (mainParmSet_ != 0 && free_parmset_mem_)
    delete mainParmSet_;
}

const char* Creator::parm_keywords_ = "parmset <parameter setname>";

const char* Creator::template_keywords_ = "lib <template setname>";

const char* Creator::other_keywords_ = "atomscandir {f|b}";

/** Initialize */
int Creator::InitCreator(ArgList& argIn, DataSetList const& DSL, int debugIn)
{
  debug_ = debugIn;

  // Atom scan direction
  std::string atomscandir = argIn.GetStringKey("atomscandir");
  if (!atomscandir.empty()) {
    if (atomscandir == "f")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_FORWARDS);
    else if (atomscandir == "b")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_BACKWARDS);
    else {
      mprinterr("Error: Unrecognized keyword for 'atomscandir' : %s\n", atomscandir.c_str());
      return 1;
    }
  }

  if (getTemplates(argIn, DSL)) return 1;
  if (getParameterSets(argIn, DSL)) return 1;

  return 0;
}

/** Try to identify residue template DataSet from the given residue
  * name (from e.g. the PDB/Mol2/etc file).
  */
DataSet_Coords* Creator::IdTemplateFromName(NameType const& rname,
                                            TerminalType termType)
const
{
  DataSet_Coords* out = 0;
  if (termType != Cpptraj::Structure::NON_TERMINAL) {
    // Looking for a terminal residue. Need to get sets with AssociatedData_ResId
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
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
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
      if ( rname == NameType( (*it)->Meta().Aspect() ) ) {
        out = *it;
        break;
      }
    }
  }

  return out;
}

/** Get templates */
int Creator::getTemplates(ArgList& argIn, DataSetList const& DSL) {
  // Clear existing templates
  Templates_.clear();
  std::string lib = argIn.GetStringKey("lib");
  if (lib.empty()) {
    mprintf("\tNo template(s) specified with 'lib'; using any loaded templates.\n");
    DataSetList sets = DSL.SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      // Should only be a single residue FIXME need new set type
      DataSet_Coords const& ds = static_cast<DataSet_Coords const&>( *(*it) );
      if ( ds.Top().Nres() == 1 )
        Templates_.push_back( (DataSet_Coords*)(*it) );
    }
  } else {
    while (!lib.empty()) {
      DataSetList sets = DSL.SelectGroupSets( lib, DataSet::COORDINATES ); // TODO specific set type for units?
      if (sets.empty()) {
        mprintf("Warning: No sets corresponding to '%s'\n", lib.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
        {
          // Should only be a single residue FIXME need new set type
          //DataSet_Coords const& ds = static_cast<DataSet_Coords const&>( *(*it) );
          //if ( ds.Top().Nres() == 1 )
            Templates_.push_back( (DataSet_Coords*)(*it) );
        }
      }
      lib = argIn.GetStringKey("lib");
    }
  }
  if (!Templates_.empty()) {
    mprintf("\t%zu residue templates found:", Templates_.size());
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it)
      mprintf(" %s", (*it)->legend());
    mprintf("\n");
  }

  return 0;
}

/** Get parameter sets. */
int Creator::getParameterSets(ArgList& argIn, DataSetList const& DSL) {
  // Clear any existing set
  if (mainParmSet_ != 0) {
    if (free_parmset_mem_) delete mainParmSet_;
  }
  mainParmSet_ = 0;
  free_parmset_mem_ = false;
  // Look for parmset args
  typedef std::vector<DataSet_Parameters*> Parray;
  Parray ParamSets;
  std::string parmset = argIn.GetStringKey("parmset");
  if (parmset.empty()) {
    mprintf("\tNo parameter set(s) specified with 'parmset'; using any loaded sets.\n");
    // See if there are any parameter sets.
    DataSetList sets = DSL.GetSetsOfType( "*", DataSet::PARAMETERS );
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
      ParamSets.push_back( (DataSet_Parameters*)(*it) );
  } else {
    while (!parmset.empty()) {
      DataSetList sets = DSL.GetSetsOfType( parmset, DataSet::PARAMETERS );
      if (sets.empty()) {
        mprintf("Warning: No parameter sets corresponding to '%s'\n", parmset.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
          ParamSets.push_back( (DataSet_Parameters*)(*it) );
      }
      parmset = argIn.GetStringKey("parmset");
    }
  }
  //if (ParamSets.empty()) {
  //  mprinterr("Error: No parameter sets.\n");
  //  return CpptrajState::ERR;
  //}
  if (!ParamSets.empty()) {
    mprintf("\tParameter sets:\n");
    for (Parray::const_iterator it = ParamSets.begin(); it != ParamSets.end(); ++it)
      mprintf("\t  %s\n", (*it)->legend());

    // Combine parameters if needed

    if (ParamSets.size() == 1)
      mainParmSet_ = ParamSets.front();
    else {
      free_parmset_mem_ = true;
      mprintf("\tCombining parameter sets.\n");
      Parray::const_iterator it = ParamSets.begin();
      mainParmSet_ = new DataSet_Parameters( *(*it) );
      ++it;
      ParameterSet::UpdateCount UC;
      for (; it != ParamSets.end(); ++it)
        mainParmSet_->UpdateParamSet( *(*it), UC, debug_, debug_ ); // FIXME verbose
    }
  }
  return 0;
}

