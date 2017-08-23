#include "Action_Rama.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
#include "Constants.h"

Action_Rama::Action_Rama() {
  Phi_[ALPHA]    =  -57.8;
  Psi_[ALPHA]    =  -47.0;
  Phi_[LEFT]     =   57.8;
  Psi_[LEFT]     =   47.0;
  Phi_[PP2]      =  -75.0;
  Psi_[PP2]      =  145.0;
  Phi_[HAIRPIN]  = -100.0;
  Psi_[HAIRPIN]  =  130.0;
  Phi_[EXTENDED] = -150.0;
  Psi_[EXTENDED] =  155.0;
  std::fill(phiOff_, phiOff_+NTYPES, 10.0);
  std::fill(psiOff_, psiOff_+NTYPES, 10.0);
}

// Action_Rama::Help()
void Action_Rama::Help() const {

}

const char* Action_Rama::TypeKeys_[] = {"alpha", "left", "pp2", "hairpin", "extended", 0};

// Action_Rama::Init()
Action::RetType Action_Rama::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  masterDSL_ = init.DslPtr();
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
  mprintf("\n"); // for PrintTypes
  for (int i = 0; i < (int)NTYPES; i++)
    mprintf("\t%8s : %8.3f +/- %8.3f  %8.3f +/- %8.3f\n",
            TypeKeys_[i], Phi_[i], phiOff_[i], Psi_[i], psiOff_[i]);

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
  for (Rarray::iterator r = residues_.begin(); r != residues_.end(); ++r)
    r->SetActive( false );
  int nActive = 0;
  // Loop over dihedrals
  DihedralSearch::mask_it dih = dihSearch_.begin();
  while (dih != dihSearch_.end()) {
    if ( (dih+1) == dihSearch_.end() ) break;
    // Relies on phi/psi pairs being together
    // Assume if consecutive res nums do not match we need to advance.
    if ( dih->ResNum() != (dih+1)->ResNum() )
      ++dih;
    else {
      bool hasPhi = (dih->Type() == MetaData::PHI || (dih+1)->Type() == MetaData::PHI);
      bool hasPsi = (dih->Type() == MetaData::PSI || (dih+1)->Type() == MetaData::PSI);
      if (hasPhi && hasPsi) {
/*
        ResMapType::iterator it = resMap_.lower_bound( dih->ResNum() );
        if (it == resMap_.end() || it->first != dih->ResNum())
        {
          // New residue
          DataSet* ds = masterDSL_->AddSet(DataSet::INTEGER,MetaData(dsetname_,dih->ResNum()+1));
          if (ds == 0) return Action::ERR;
          resMap_.insert(it, std::pair<int,Res>(dih->ResNum(), Res(ds, *dih, *(dih+1))));
        }
*/
        if ( dih->ResNum() >= (int)residues_.size() )
          residues_.resize( dih->ResNum() + 1 );
        Res& res = residues_[dih->ResNum()];
        res.SetActive( true );
        nActive++;
        if (res.Data() == 0) {
          res.SetData(masterDSL_->AddSet(DataSet::INTEGER,MetaData(dsetname_,dih->ResNum()+1)));
          if (res.Data() == 0) return Action::ERR;
          if (outfile_ != 0) outfile_->AddDataSet( res.Data() );
        }
        res.SetMasks( *dih, *(dih+1) );
      }
      dih += 2;
    }
  } // END loop over dihedrals
  mprintf("\t%i residues active.\n", nActive);
  if (nActive < 1) return Action::SKIP;

  return Action::OK;
}

// Action_Rama::DoAction()
Action::RetType Action_Rama::DoAction(int frameNum, ActionFrame& frm)
{
  for (Rarray::const_iterator res = residues_.begin(); res != residues_.end(); ++res)
  {
    if (res->IsActive()) {
      double phi = Torsion( frm.Frm().XYZ(res->Phi().A0()),
                            frm.Frm().XYZ(res->Phi().A1()),
                            frm.Frm().XYZ(res->Phi().A2()),
                            frm.Frm().XYZ(res->Phi().A3()) );
      double psi = Torsion( frm.Frm().XYZ(res->Psi().A0()),
                            frm.Frm().XYZ(res->Psi().A1()),
                            frm.Frm().XYZ(res->Psi().A2()),
                            frm.Frm().XYZ(res->Psi().A3()) );
      // TODO Radians
      phi *= Constants::RADDEG;
      psi *= Constants::RADDEG;
      // Determine Rama. region
      int currentType = -1;
      for (int i = 0; i < (int)NTYPES; i++) {
        if (phi > Phi_[i] - phiOff_[i] &&
            phi < Phi_[i] + phiOff_[i] &&
            psi > Psi_[i] - psiOff_[i] &&
            psi < Psi_[i] + psiOff_[i])
        {
          currentType = i;
          break;
        }
      }
      res->Data()->Add(frameNum, &currentType);
    }
  }
  return Action::OK;
}
