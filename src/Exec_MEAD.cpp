#include "Exec_MEAD.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"
#include "Structure/SiteData.h"
#include "Structure/TitratableSite.h"
#include "Mead/MeadGrid.h"
#include "Mead/MeadOpts.h"
#include "Mead/MultiFlexResults.h"
#include "Mead/MeadCalc.h"
#include "Mead/MeadCalc_Solvate.h"
#include "Mead/MeadCalc_Potential.h"
#include "Mead/MeadCalc_Multiflex.h"
#include <cmath> //sqrt

using namespace Cpptraj::Mead;

Exec_MEAD::Exec_MEAD() :
  Exec(GENERAL),
  MEAD_(0)
{}

Exec_MEAD::~Exec_MEAD() {
  if (MEAD_ != 0) delete MEAD_;
}

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
int Exec_MEAD::CheckMead(MeadGrid const& ogm) const {
  // Sanity checks
  if (!ogm.IsSetup()) {
    mprinterr("Error: No MEAD grid set up (ogm).\n");
    return 1;
  }
  if (MEAD_ == 0 || !MEAD_->HasAtoms()) {
    mprinterr("Error: No MEAD atoms allocated.\n");
    return 1;
  }
  return 0;
}

/** MEAD solvate. */
int Exec_MEAD::Solvate(CpptrajState& State, ArgList& argIn, MeadGrid const& ogm,
                       std::string const& dsname, DataFile* outfile)
const
{
  // Allocate output set(s)
  std::string outSetName = dsname;
  if (outSetName.empty())
    outSetName = State.DSL().GenerateDefaultName("SOLVATE");
  DataSet* outset = State.DSL().AddSet( DataSet::DOUBLE, outSetName );
  if (outset == 0) {
    mprinterr("Error: Could not allocate output set '%s' for solvate\n", outSetName.c_str());
    return 1;
  }
  if (outfile != 0)
    outfile->AddDataSet( outset );
  std::string rxngrid = argIn.GetStringKey("rxngrid");
  DataSet_3D* rgrid = 0;
  if (!rxngrid.empty()) {
    DataSet* ds = State.DSL().FindSetOfGroup( rxngrid, DataSet::GRID_3D );
    if (ds == 0) {
      mprinterr("Error: No grid set with name '%s' found.\n", rxngrid.c_str());
      return 1;
    }
    rgrid = (DataSet_3D*)ds;
  }

  if (CheckMead( ogm )) return 1;

  MeadOpts Opts;
  Opts.SetEpsIn(argIn.getKeyDouble("epsin", 1));
  Opts.SetEpsExt(argIn.getKeyDouble("epssol", 80));
  Opts.SetEpsVac(argIn.getKeyDouble("epsvac", 1));
  Opts.SetSolRad(argIn.getKeyDouble("solrad", 1.4));
  Opts.SetSterLn(argIn.getKeyDouble("sterln", 2.0));
  Opts.SetIonicStr(argIn.getKeyDouble("ionicstr", 0.0));
  Opts.SetTemperature(argIn.getKeyDouble("temp", 300.0));

  double Esolv = 0;
  int err = ((MeadCalc_Solvate*)MEAD_)->Solvate(Esolv, Opts, ogm, rgrid);

  if (err != 0) return 1;
  outset->Add(0, &Esolv);
  return 0;
}

/** MEAD potential. */
int Exec_MEAD::Potential(CpptrajState& State, ArgList& argIn,
                         MeadGrid const& ogm,
                         std::string const& dsname, DataFile* outfile)
const
{
  // Allocate output set
  std::string outSetName = dsname;
  if (outSetName.empty())
    outSetName = State.DSL().GenerateDefaultName("POTENTIAL");
  DataSet_Vector_Scalar* outset = (DataSet_Vector_Scalar*)State.DSL().AddSet( DataSet::VECTOR_SCALAR, outSetName );
  if (outset == 0) {
    mprinterr("Error: Could not allocate output set '%s' for potential\n", outSetName.c_str());
    return 1;
  }
  if (outfile != 0)
      outfile->AddDataSet( (DataSet*)outset );
  if (CheckMead( ogm )) return 1;
 
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
    if ( ((MeadCalc_Potential*)MEAD_)->Potential(*outset, Opts, ogm, fieldPoints)) {
      mprinterr("Error: Could not process MEAD field points.\n");
      return 1;
    }
  }

  return 0;
}

/** Run multiflex. */
int Exec_MEAD::MultiFlex(CpptrajState& State, ArgList& argIn, 
                         MeadGrid const& ogm, MeadGrid const& mgm, 
                         Topology const& topIn, Frame const& frameIn,
                         std::string const& dsname, DataFile* outfile
                         )
const
{
  // Allocate output sets
  std::string outSetName = dsname;
  if (outSetName.empty())
    outSetName = State.DSL().GenerateDefaultName("MULTIFLEX");
  MultiFlexResults results;
  if (results.CreateSets(State.DSL(), outSetName)) return 1;
  DataFile* ssiout = State.DFL().AddDataFile( argIn.GetStringKey("ssiout"), argIn );
  results.AddSetsToFile( outfile, ssiout );
  if (results.CreateOutputFiles(State.DFL(),
                                argIn.GetStringKey("pkint"),
                                argIn.GetStringKey("summ"),
                                argIn.GetStringKey("gfile")))
  {
    mprinterr("Error: Could not create MEAD output files for multiflex.\n");
    return 1;
  }

  if (CheckMead( ogm )) return 1;
  using namespace Cpptraj::Structure;

  MeadOpts Opts;
  std::string sitesFileName = argIn.GetStringKey("sites");
  std::string sitesDirName = argIn.GetStringKey("sitesdir");
  int siteIdx = argIn.getKeyInt("site", 0) - 1; // User args begin at 1
  Opts.SetEpsIn(argIn.getKeyDouble("epsin", 1));
  Opts.SetEpsExt(argIn.getKeyDouble("epssol", 80));
  Opts.SetSolRad(argIn.getKeyDouble("solrad", 1.4));
  Opts.SetSterLn(argIn.getKeyDouble("sterln", 2.0));
  Opts.SetIonicStr(argIn.getKeyDouble("ionicstr", 0.0));

  mprintf("\tSites file : %s\n", sitesFileName.c_str());
  mprintf("\tSites dir  : %s\n", sitesDirName.c_str());
  if (siteIdx == -1)
    mprintf("\tCalculating all sites.\n");
  else
    mprintf("\tOnly calculating site %i\n", siteIdx + 1);

  //if (sitesFileName.empty()) {
  //  mprinterr("Error: No sites file provided.\n");
  //  return 1;
  //}

  SiteData titrationData;

  if (sitesFileName.empty()) {
    if (titrationData.LoadSiteDirectory( sitesDirName )) {
      mprinterr("Error: Could not load titration sites data from directory.\n");
      return 1;
    }
    if (titrationData.SetupSitesFromTop( topIn )) {
      mprinterr("Error: Could not set up titratable sites from topology '%s'\n", topIn.c_str());
      return 1;
    }
  } else {
    if (titrationData.LoadMeadSiteData( sitesFileName, sitesDirName, topIn )) {
      mprinterr("Error: Could not load MEAD titration sites data.\n");
      return 1;
    }
  }
  if (titrationData.NoSites()) {
    mprinterr("Error: No sites to calculate titration for.\n");
    return 1;
  }

  if (((MeadCalc_Multiflex*)MEAD_)->MultiFlex(results, Opts, ogm, mgm, topIn, frameIn, titrationData, siteIdx)) {
    mprinterr("Error: Multiflex failed.\n");
    return 1;
  }
  ((MeadCalc_Multiflex*)MEAD_)->TotalTime().WriteTiming();

  return 0;
}

/** Add level to a Mead grid. */
int Exec_MEAD::addGridLevel(MeadGrid& ogm, std::string const& ogmstr) {
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

/** Set up MEAD grid based on furthest atom-atom distance. */
int Exec_MEAD::setup_grid_from_coords(MeadGrid& ogm, Frame const& frameIn) {
  double maxdist2 = 0;
  int maxat1 = -1;
  int maxat2 = -1;
  for (int at1 = 0; at1 < frameIn.Natom(); at1++) {
    const double* xyz1 = frameIn.XYZ(at1);
    for (int at2 = at1 + 1; at2 < frameIn.Natom(); at2++) {
      const double* xyz2 = frameIn.XYZ(at2);
      double dist2 = DIST2_NoImage(xyz1, xyz2);
      if (dist2 > maxdist2) {
        maxdist2 = dist2;
        maxat1 = at1;
        maxat2 = at2;
      }
    }
  }
  double max = sqrt(maxdist2);
  mprintf("\tMax distance is %g Ang. between atoms %i and %i.\n", max, maxat1+1, maxat2+1);
  max = int(max + 5);
  int ival = (int)max;
  ival = ival % 2;
  if (ival == 0)
    max = max + 1;

  mprintf("\tMAX= %g\n", max);
  ogm.AddGrid(max, 8.0, MeadGrid::C_ON_GEOM_CENT);
  ogm.AddGrid(max, 2.0, MeadGrid::C_ON_CENT_OF_INTR);
  ogm.AddGrid(max, 0.5, MeadGrid::C_ON_CENT_OF_INTR);

  return 0;
}
 


// Exec_MEAD::Execute()
Exec::RetType Exec_MEAD::Execute(CpptrajState& State, ArgList& argIn)
{
  using namespace Cpptraj;
  if (MEAD_ != 0) delete MEAD_;
  enum McalcType { POTENTIAL = 0, SOLVATE, MULTIFLEX };
  McalcType calcType;
  if (argIn.hasKey("potential")) {
    calcType = POTENTIAL;
    MEAD_ = new MeadCalc_Potential();
  } else if (argIn.hasKey("solvate")) {
    calcType = SOLVATE;
    MEAD_ = new MeadCalc_Solvate();
  } else if (argIn.hasKey("multiflex")) {
    calcType = MULTIFLEX;
    MEAD_ = new MeadCalc_Multiflex();
  } else {
    mprinterr("Error: No valid MEAD keyword specified (potential, solvate, multiflex).\n");
    return CpptrajState::ERR;
  }

  int verbose = argIn.getKeyInt("verbose", 0);
  MEAD_->MeadVerbosity( verbose );

  MeadGrid ogm, mgm;
  std::string ogmstr = argIn.GetStringKey("ogm");
  while (!ogmstr.empty()) {
    if (addGridLevel(ogm, ogmstr)) return CpptrajState::ERR;
    ogmstr = argIn.GetStringKey("ogm");
  }
  
  std::string mgmstr = argIn.GetStringKey("mgm");
  while (!mgmstr.empty()) {
    if (addGridLevel(mgm, mgmstr)) return CpptrajState::ERR;
    mgmstr = argIn.GetStringKey("mgm");
  }
  
  MeadCalc::Radii_Mode radiiMode;
  std::string radmode = argIn.GetStringKey("radii");
  if (radmode == "gb")
    radiiMode = MeadCalc::GB;
  else if (radmode == "parse")
    radiiMode = MeadCalc::PARSE;
  else if (radmode == "vdw")
    radiiMode = MeadCalc::VDW;
  else
    radiiMode = MeadCalc::GB;

  switch (radiiMode) {
    case MeadCalc::GB : mprintf("\tUsing GB radii.\n"); break;
    case MeadCalc::PARSE : mprintf("\tUsing PARSE radii.\n"); break;
    case MeadCalc::VDW : mprintf("\tUsing VDW radii.\n"); break;
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
  if (MEAD_->SetupAtoms( CRD->Top(), frameIn, radiiMode )) {
    mprinterr("Error: Setting up frame/topology failed.\n");
    return CpptrajState::ERR; 
  }

  // Set up grid if needed
  if (!ogm.IsSetup()) {
    mprintf("\tPerforming automatic grid setup.\n");
    if (setup_grid_from_coords(ogm, frameIn)) {
      mprinterr("Error: Automatic grid setup failed.\n");
      return CpptrajState::ERR;
    }
  }
  ogm.Print();

  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );

  MEAD_->Print();

  int err = 0;
  if (calcType == POTENTIAL) {
    err = Potential( State, argIn, ogm, outSetName, outfile );
  } else if (calcType == SOLVATE) {
    err = Solvate( State, argIn, ogm, outSetName, outfile );
  } else if (calcType == MULTIFLEX) {
    if (!mgm.IsSetup())
      mgm = ogm;
    mgm.Print();
    err = MultiFlex( State, argIn, ogm, mgm, CRD->Top(), frameIn, outSetName, outfile );
  } else {
    mprinterr("Error: No MEAD calculation keywords given.\n");
    err = 1;
  }

  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;    
}
