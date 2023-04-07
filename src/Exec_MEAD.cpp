#include "Exec_MEAD.h"
#include "CpptrajStdio.h"
#include "MeadInterface.h"
#include "MultiFlexResults.h"
#include "StringRoutines.h"
#include "Structure/TitrationData.h"
#include "Structure/TitratableSite.h"
#include "Mead/MeadGrid.h"
#include "Mead/MeadOpts.h"

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
          "\t          [temp <temperature>] [rxngrid <input grid>]\n"
          "\t}\n"
         );
}

/** Check MEAD is properly set up. */
int Exec_MEAD::CheckMead(Cpptraj::MeadInterface const& MEAD) {
  // Sanity checks
  //if (!MEAD.HasFDM()) {
  //  mprinterr("Error: No MEAD grid allocated.\n");
  //  return 1;
  //}
  if (!MEAD.HasAtoms()) {
    mprinterr("Error: No MEAD atoms allocated.\n");
    return 1;
  }
  return 0;
}

/** MEAD solvate. */
int Exec_MEAD::Solvate(Cpptraj::MeadInterface& MEAD, Cpptraj::Mead::MeadGrid const& ogm, ArgList& argIn, DataSet* outset, DataSet_3D* rgrid)
const
{
  using namespace Cpptraj::Mead;
  //if (CheckMead( MEAD )) return 1; FIXME add check

  MeadOpts Opts;
  Opts.SetEpsIn(argIn.getKeyDouble("epsin", 1));
  Opts.SetEpsExt(argIn.getKeyDouble("epssol", 80));
  Opts.SetEpsVac(argIn.getKeyDouble("epsvac", 1));
  Opts.SetSolRad(argIn.getKeyDouble("solrad", 1.4));
  Opts.SetSterLn(argIn.getKeyDouble("sterln", 2.0));
  Opts.SetIonicStr(argIn.getKeyDouble("ionicstr", 0.0));
  Opts.SetTemperature(argIn.getKeyDouble("temp", 300.0));

  double Esolv = 0;
  int err = MEAD.Solvate(Esolv, Opts, ogm, rgrid);

  if (err != 0) return 1;
  outset->Add(0, &Esolv);
  return 0;
}

/** MEAD potential. */
int Exec_MEAD::Potential(Cpptraj::MeadInterface& MEAD, Cpptraj::Mead::MeadGrid const& ogm, ArgList& argIn, DataSet_Vector_Scalar& outset) const {
  //if (CheckMead( MEAD )) return 1; FIXME add check
  using namespace Cpptraj::Mead;
 
  MeadOpts Opts; 
  Opts.SetEpsIn(argIn.getKeyDouble("epsin", 1));
  Opts.SetEpsExt(argIn.getKeyDouble("epsext", 80));

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
    if (MEAD.Potential(outset, Opts, ogm, fieldPoints)) {
      mprinterr("Error: Could not process MEAD field points.\n");
      return 1;
    }
  }

  return 0;
}

/** Run multiflex. */
int Exec_MEAD::MultiFlex(Cpptraj::MeadInterface& MEAD, Cpptraj::Mead::MeadGrid const& ogm,
                         Cpptraj::Mead::MeadGrid const& mgm, 
                         ArgList& argIn, Topology const& topIn, Frame const& frameIn,
                         int iradiimode, Cpptraj::MultiFlexResults const& results)
const
{
  using namespace Cpptraj::Structure;

  Cpptraj::MeadInterface::Radii_Mode radiiMode = (Cpptraj::MeadInterface::Radii_Mode)iradiimode;

  using namespace Cpptraj::Mead;
  MeadOpts Opts;
  std::string sitesFileName = argIn.GetStringKey("sites");
  std::string sitesDirName = argIn.GetStringKey("sitesdir");
  Opts.SetEpsIn(argIn.getKeyDouble("epsin", 1));
  Opts.SetEpsExt(argIn.getKeyDouble("epssol", 80));
  Opts.SetSolRad(argIn.getKeyDouble("solrad", 1.4));
  Opts.SetSterLn(argIn.getKeyDouble("sterln", 2.0));
  Opts.SetIonicStr(argIn.getKeyDouble("ionicstr", 0.0));

  mprintf("\tSites file : %s\n", sitesFileName.c_str());
  mprintf("\tSites dir  : %s\n", sitesDirName.c_str());

  if (sitesFileName.empty()) {
    mprinterr("Error: No sites file provided.\n");
    return 1;
  }

  TitrationData titrationData;

  if (titrationData.LoadTitrationData( sitesFileName, sitesDirName )) {
    mprinterr("Error: Could not load titration data.\n");
    return 1;
  }

  if (MEAD.MultiFlex(results, Opts, ogm, mgm, topIn, frameIn, titrationData, radiiMode)) {
    mprinterr("Error: Multiflex failed.\n");
    return 1;
  } 

  return 0;
}

/** Add level to a Mead grid. */
int Exec_MEAD::addGridLevel(Cpptraj::Mead::MeadGrid& ogm, std::string const& ogmstr) {
  using namespace Cpptraj::Mead;
  // Format: N,spacing[,centering]
  ArgList ogmarg( ogmstr, "," );
  if (ogmarg.Nargs() < 2) {
    mprinterr("Error: Malformed ogm key; expected <N>,<spacing>\n");
    return 1;
  }
  int ngridpts = convertToInteger( ogmarg[0] );
  double spacing = convertToDouble( ogmarg[1] );
  if ( ogmarg.Nargs() == 2) {
    // No centering, default to origin
    mprintf("\tAdding grid of %i points, spacing %g, center on origin.\n", ngridpts, spacing);
    if (ogm.AddGrid(ngridpts, spacing, Vec3(0,0,0))) {
      mprinterr("Error: Adding MEAD grid.\n");
      return 1;
    }
  } else if ( ogmarg.Nargs() == 3) {
    MeadGrid::Center_Mode gc;
    // Center using string
    if (ogmarg[2] == "o")
      gc = MeadGrid::C_ON_ORIGIN;
    else if (ogmarg[2] == "i")
      gc = MeadGrid::C_ON_CENT_OF_INTR;
    else if (ogmarg[2] == "g")
      gc = MeadGrid::C_ON_GEOM_CENT;
    else {
      mprinterr("Error: Expected grid centering argument to be 'o', 'i', or 'g', got '%s'\n",
                ogmarg[2].c_str());
      return 1;
    }
    mprintf("\tAdding grid of %i points, spacing %g, center on %s\n", ngridpts, spacing,
            MeadGrid::Center_ModeStr(gc));
    if (ogm.AddGrid(ngridpts, spacing, gc)) {
      mprinterr("Error: Adding MEAD grid.\n");
      return 1;
    }
  } else {
    mprinterr("Error: Invalid centering argument for grid.\n");
    return 1;
  }
  return 0;
}

// Exec_MEAD::Execute()
Exec::RetType Exec_MEAD::Execute(CpptrajState& State, ArgList& argIn)
{
  using namespace Cpptraj;
  MeadInterface MEAD;
  int verbose = argIn.getKeyInt("verbose", 0);
  MEAD.MeadVerbosity( verbose );

  using namespace Cpptraj::Mead;
  MeadGrid ogm, mgm;

  std::string ogmstr = argIn.GetStringKey("ogm");
  while (!ogmstr.empty()) {
    if (addGridLevel(ogm, ogmstr)) return CpptrajState::ERR;
    ogmstr = argIn.GetStringKey("ogm");
  }
  ogm.Print();
  if (!ogm.IsSetup()) {
    mprinterr("Error: No grid set up (ogm).\n");
    return CpptrajState::ERR;
  }

  std::string mgmstr = argIn.GetStringKey("mgm");
  while (!mgmstr.empty()) {
    if (addGridLevel(mgm, mgmstr)) return CpptrajState::ERR;
    mgmstr = argIn.GetStringKey("mgm");
  }
  if (!mgm.IsSetup())
    mgm = ogm;
  mgm.Print();

  MeadInterface::Radii_Mode radiiMode;
  std::string radmode = argIn.GetStringKey("radii");
  if (radmode == "gb")
    radiiMode = MeadInterface::GB;
  else if (radmode == "parse")
    radiiMode = MeadInterface::PARSE;
  else if (radmode == "vdw")
    radiiMode = MeadInterface::VDW;
  else
    radiiMode = MeadInterface::GB;

  switch (radiiMode) {
    case MeadInterface::GB : mprintf("\tUsing GB radii.\n"); break;
    case MeadInterface::PARSE : mprintf("\tUsing PARSE radii.\n"); break;
    case MeadInterface::VDW : mprintf("\tUsing VDW radii.\n"); break;
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
    err = Potential( MEAD, ogm, argIn, *outset );
  } else if (argIn.hasKey("solvate")) {
    // Allocate output set(s)
    if (outSetName.empty())
      outSetName = State.DSL().GenerateDefaultName("SOLVATE");
    DataSet* outset = State.DSL().AddSet( DataSet::DOUBLE, outSetName );
    if (outset == 0) {
      mprinterr("Error: Could not allocate output set '%s' for solvate\n", outSetName.c_str());
      return CpptrajState::ERR;
    }
    if (outfile != 0)
      outfile->AddDataSet( outset );
    std::string rxngrid = argIn.GetStringKey("rxngrid");
    DataSet_3D* rgrid = 0;
    if (!rxngrid.empty()) {
      DataSet* ds = State.DSL().FindSetOfGroup( rxngrid, DataSet::GRID_3D );
      if (ds == 0) {
        mprinterr("Error: No grid set with name '%s' found.\n", rxngrid.c_str());
        return CpptrajState::ERR;
      }
      rgrid = (DataSet_3D*)ds;
    }
    err = Solvate( MEAD, ogm, argIn, outset, rgrid );
  } else if (argIn.hasKey("multiflex")) {
    /// Allocate output sets
    if (outSetName.empty())
      outSetName = State.DSL().GenerateDefaultName("MULTIFLEX");
    MultiFlexResults results;
    if (results.Allocate(State.DSL(), outSetName)) return CpptrajState::ERR;
    DataFile* ssiout = State.DFL().AddDataFile( argIn.GetStringKey("ssiout"), argIn );
    results.AddSetsToFile( outfile, ssiout );
    if (results.CreateOutputFiles(State.DFL(),
                                  argIn.GetStringKey("pkint"),
                                  argIn.GetStringKey("summ"),
                                  argIn.GetStringKey("gfile")))
    {
      mprinterr("Error: Could not create MEAD output files for multiflex.\n");
      return CpptrajState::ERR;
    }
    err = MultiFlex( MEAD, ogm, mgm, argIn, CRD->Top(), frameIn, (int)radiiMode, results );
  } else {
    mprinterr("Error: No MEAD calculation keywords given.\n");
    err = 1;
  }

  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;    
}
