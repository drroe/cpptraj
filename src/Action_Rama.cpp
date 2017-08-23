#include "Action_Rama.h"
#include "CpptrajStdio.h"

// Action_Rama::Help()
void Action_Rama::Help() const {

}

const char* Action_Rama::TypeKeys_[] = {"alpha", "beta", "pp2", "left", 0};

// Action_Rama::Init()
Action::RetType Action_Rama::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  Nframe_ = 0;
  // Only phi and psi for now
  dihSearch_.SearchFor(MetaData::PHI);
  dihSearch_.SearchFor(MetaData::PSI);
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string temp = actionArgs.GetStringKey("sumout");
  if (temp.empty() && outfile_ != 0)
    temp = outfile_->DataFilename().Full() + ".sum";
  sumFile_ = init.DFL().AddDataFile( temp );
  DataFile* totalout = init.DFL().AddDataFile( actionArgs.GetStringKey("totalout"), actionArgs );
  //assignout_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("assignout"), "SS assignment");
  BB_N_ = actionArgs.GetStringKey("namen", "N");
  BB_H_ = actionArgs.GetStringKey("nameh", "H");
  BB_C_ = actionArgs.GetStringKey("namec", "C");
  BB_O_ = actionArgs.GetStringKey("nameo", "O");
  BB_CA_ = actionArgs.GetStringKey("nameca", "CA");
  // Type parameters
  std::string typearg = actionArgs.GetStringKey("type");
  while (!typearg.empty()) {
    // Comma-separated
    ArgList Args(typearg, ",");
    if (Args.Nargs() < 2) {
      mprinterr("Error: Expected 'type <type>,<args>', got '%s'\n", typearg.c_str());
      return Action::ERR;
    }
    Type currentType = NTYPES;
    for (int i = 0; i < (int)NTYPES; i++)
      if (Args[0].compare(TypeKeys_[i]) == 0) {
        currentType = (Type)i;
        break;
      }
    if (currentType == NTYPES) {
      mprinterr("Error: Unrecognized type '%s'\n", Args[0].c_str());
      return Action::ERR;
    }
    for (int i = 1; i != Args.Nargs(); i++)
    {
      // <arg>=<val>
      ArgList A0(Args[i], "=");
      if (A0.Nargs() != 2) {
        mprinterr("Error: Bad argument '%s'; expect <arg>=<value>\n", Args[i].c_str());
        return Action::ERR;
      }
      if (A0[0] == "phioff") phiOff_[currentType] = A0.getNextDouble(10.0);
      else if (A0[0] == "psioff") psiOff_[currentType] = A0.getNextDouble(10.0);
    }
    typearg = actionArgs.GetStringKey("type");
  } // END loop over type args
    
  
  // Get masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );
  // Set up the data set
  dsetname_ = actionArgs.GetStringNext();
  if (dsetname_.empty())
    dsetname_ = init.DSL().GenerateDefaultName("RAMA");
  // Type total v time
  for (int i = 0; i < (int)NTYPES; i++) {
    ds_[i] = init.DSL().AddSet( DataSet::INTEGER, MetaData(dsetname_, TypeKeys_[i]) );
    if (ds_[i] == 0) return Action::ERR;
    if (totalout != 0) totalout->AddDataSet( ds_[i] );
  }

  mprintf("    RAMACHANDRAN PLOT:");
  mprintf("\tUsing residues in mask '%s'\n", Mask_.MaskString());
  //mprintf("\t Atom Names: N=%s H=%s C=%s
  dihSearch_.PrintTypes();

  return Action::OK;
}

// Action_Rama::Setup()
Action::RetType Action_Rama::Setup(ActionSetup& setup)
{
  if (setup.Top().SetupCharMask( Mask_ )) return Action::ERR;
  if ( Mask_.None() ) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  // Loop over residues, see which ones are selected
  Range actualRange;
  for (Topology::res_iterator res = setup.Top().ResStart(); res != setup.Top().ResEnd(); ++res)
  {
    if ( Mask_.AtomsInCharMask( res->FirstAtom(), res->LastAtom() ) )
      actualRange.AddToRange( res - setup.Top().ResStart() );
  }
  mprintf("\t%i residues selected.\n", actualRange.Size());
  // Find dihedrals
  if (dihSearch_.FindDihedrals(setup.Top(), actualRange))
    return Action::SKIP;
  mprintf("\t%i dihedrals.\n", dihSearch_.Ndihedrals());
  if (dihSearch_.Ndihedrals() < 1) return Action::SKIP;
  // Loop over dihedrals
  
  return Action::OK;
}

// Action_Rama::DoAction()
Action::RetType Action_Rama::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::ERR;
}
