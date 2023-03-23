#include "Exec_MEAD.h"
#include "CpptrajStdio.h"
#include "MeadInterface.h"
#include "StringRoutines.h"

// Exec_MEAD::Help()
void Exec_MEAD::Help() const
{
  mprintf("\t[ogm <ngridpoints>,<gridspacing> ...]\n"
          "\t[crdset <COORDS set>] [radmode {gb|parse|vdw}]\n"
          "\t[name <output set name>] [out <file>]\n"
          "\t[verbose <#>]\n"
          "\t{ potential [epsin <epsilon in>] [epsext <epsilon out>]\n"
          "\t            [fpt <X>,<Y>,<Z> ...] |\n"
          "\t  solvate [epsin <epsilon in>] [epssol <solvent epsilon>]\n"
          "\t          [epsvac <vacuum epsilon>] [solrad <probe radius>]\n"
          "\t          [sterln <ion layer thickness>] [ionicstr <ionic strength>]\n"
          "\t          [temp <temperature>]\n"
          "\t}\n"
         );
}

/** Check MEAD is properly set up. */
int Exec_MEAD::CheckMead(Cpptraj::MeadInterface const& MEAD) {
  // Sanity checks
  if (!MEAD.HasFDM()) {
    mprinterr("Error: No MEAD grid allocated.\n");
    return 1;
  }
  if (!MEAD.HasAtoms()) {
    mprinterr("Error: No MEAD atoms allocated.\n");
    return 1;
  }
  return 0;
}

/** MEAD solvate. */
int Exec_MEAD::Solvate(Cpptraj::MeadInterface& MEAD, ArgList& argIn, DataSet* outset)
const
{
  if (CheckMead( MEAD )) return 1;

  double epsin = argIn.getKeyDouble("epsin", 1);
  double epssol = argIn.getKeyDouble("epssol", 80);
  double epsvac = argIn.getKeyDouble("epsvac", 1);
  double solrad = argIn.getKeyDouble("solrad", 1.4);
  double sterln = argIn.getKeyDouble("sterln", 2.0);
  double ionicstr = argIn.getKeyDouble("ionicstr", 0.0);
  double temperature = argIn.getKeyDouble("temp", 300.0);

  double Esolv = 0;
  int err = MEAD.Solvate(Esolv, epsin, epssol, epsvac, solrad, sterln, ionicstr, temperature);

  if (err != 0) return 1;
  outset->Add(0, &Esolv);
  return 0;
}

/** MEAD potential. */
int Exec_MEAD::Potential(Cpptraj::MeadInterface& MEAD, ArgList& argIn, DataSet_Vector_Scalar& outset) const {
  if (CheckMead( MEAD )) return 1;

  double epsin = argIn.getKeyDouble("epsin", 1);
  double epsext = argIn.getKeyDouble("epsext", 80);

  std::vector<Vec3> fieldPoints;
  std::string fptstr = argIn.GetStringKey("fpt");
  while (!fptstr.empty()) {
    ArgList fptarg(fptstr, ",");
    if (fptarg.Nargs() != 3) {
      mprinterr("Error: Malformed fpt arg '%s', expected <X>,<Y>,<Z>\n", fptstr.c_str());
      return 1;
    }
    for (int ii = 0; ii != 3; ii++) {
      if (!validDouble(fptarg[ii])) {
        mprinterr("Error: '%s' is not a valid double.\n", fptarg[ii].c_str());
        return 1;
      }
    }
    fieldPoints.push_back( Vec3( convertToDouble(fptarg[0]),
                                 convertToDouble(fptarg[1]),
                                 convertToDouble(fptarg[2]) ) );
    fptstr = argIn.GetStringKey("fpt");
  }
  if (fieldPoints.empty()) {
    mprintf("Warning: No field points specified.\n");
  } else {
    if (MEAD.Potential(outset, epsin, epsext, fieldPoints)) {
      mprinterr("Error: Could not process MEAD field points.\n");
      return 1;
    }
  }

  return 0;
}

// Exec_MEAD::Execute()
Exec::RetType Exec_MEAD::Execute(CpptrajState& State, ArgList& argIn)
{
  Cpptraj::MeadInterface MEAD;
  int verbose = argIn.getKeyInt("verbose", 0);
  MEAD.MeadVerbosity( verbose );

  std::string ogmstr = argIn.GetStringKey("ogm");
  while (!ogmstr.empty()) {
    // Format: N,spacing TODO centering
    ArgList ogmarg( ogmstr, "," );
    if (ogmarg.Nargs() < 2) {
      mprinterr("Error: Malformed ogm key; expected <N>,<spacing>\n");
      return CpptrajState::ERR;
    }
    int ngridpts = convertToInteger( ogmarg[0] );
    double spacing = convertToDouble( ogmarg[1] );
    mprintf("\tAdding grid of %i points, spacing %g\n", ngridpts, spacing);
    if (MEAD.AddGrid(ngridpts, spacing, Vec3(0,0,0))) {
      mprinterr("Error: Adding MEAD grid.\n");
      return CpptrajState::ERR;
    }
    ogmstr = argIn.GetStringKey("ogm");
  }

  Cpptraj::MeadInterface::Radii_Mode radiiMode;
  std::string radmode = argIn.GetStringKey("radii");
  if (radmode == "gb")
    radiiMode = Cpptraj::MeadInterface::GB;
  else if (radmode == "parse")
    radiiMode = Cpptraj::MeadInterface::PARSE;
  else if (radmode == "vdw")
    radiiMode = Cpptraj::MeadInterface::VDW;
  else
    radiiMode = Cpptraj::MeadInterface::GB;

  switch (radiiMode) {
    case Cpptraj::MeadInterface::GB : mprintf("\tUsing GB radii.\n"); break;
    case Cpptraj::MeadInterface::PARSE : mprintf("\tUsing PARSE radii.\n"); break;
    case Cpptraj::MeadInterface::VDW : mprintf("\tUsing VDW radii.\n"); break;
  }

  std::string outSetName = argIn.GetStringKey("name");

  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: No COORDS set specified.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return CpptrajState::ERR;
  }
  // TODO multiple frames
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' has no frames.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  Frame frameIn = CRD->AllocateFrame();
  CRD->GetFrame(0, frameIn);
  if (MEAD.SetupAtoms( CRD->Top(), frameIn, radiiMode )) {
    mprinterr("Error: Setting up frame/topology failed.\n");
    return CpptrajState::ERR; 
  }

  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );

  MEAD.Print();

  int err = 0;
  if (argIn.hasKey("potential")) {
    // Allocate output set
    if (outSetName.empty())
      outSetName = State.DSL().GenerateDefaultName("POTENTIAL");
    DataSet_Vector_Scalar* outset = (DataSet_Vector_Scalar*)State.DSL().AddSet( DataSet::VECTOR_SCALAR, outSetName );
    if (outset == 0) {
      mprinterr("Error: Could not allocate output set '%s' for potential\n", outSetName.c_str());
      return CpptrajState::ERR;
    }
    if (outfile != 0)
      outfile->AddDataSet( (DataSet*)outset );
    err = Potential( MEAD, argIn, *outset );
  } else if (argIn.hasKey("solvate")) {
    // Allocate output set
    if (outSetName.empty())
      outSetName = State.DSL().GenerateDefaultName("SOLVATE");
    DataSet* outset = State.DSL().AddSet( DataSet::DOUBLE, outSetName );
    if (outset == 0) {
      mprinterr("Error: Could not allocate output set '%s' for solvate\n", outSetName.c_str());
      return CpptrajState::ERR;
    }
    if (outfile != 0)
      outfile->AddDataSet( outset );
    err = Solvate( MEAD, argIn, outset );
  } else {
    mprinterr("Error: No MEAD calculation keywords given.\n");
    err = 1;
  }

  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;    
}
