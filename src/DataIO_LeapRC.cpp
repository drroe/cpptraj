#include "DataIO_LeapRC.h"
#include "AssociatedData_ResId.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataIO_AmberFF.h"
#include "DataIO_AmberFrcmod.h"
#include "DataIO_AmberLib.h"
#include "DataIO_AmberPrep.h"
#include "DataIO_Coords.h"
#include "DataSet_NameMap.h"
#include "DataSet_Parameters.h"
#include "Exec_Build.h"
#include "StringRoutines.h" // ToLower
#include <cstdlib> //getenv

/// CONSTRUCTOR
DataIO_LeapRC::DataIO_LeapRC()
{

}

/** Track already loaded parm files. */
DataIO_LeapRC::Sarray DataIO_LeapRC::paramFiles_ = Sarray();

/** Track already loaded lib/prep files. */
DataIO_LeapRC::Sarray DataIO_LeapRC::libFiles_ = Sarray();

// DataIO_LeapRC::ID_DataFormat()
bool DataIO_LeapRC::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  bool isLeaprc = false;
  // Scan the first 5 non-blank non-comment lines
  int nLinesScanned = 0;
  while (nLinesScanned < 5 && !isLeaprc) {
    const char* ptr = infile.NextLine();
    if (ptr == 0) break;
    if (ptr[0] == '\0' || ptr[0] == '#') continue;
    nLinesScanned++;
    // LEaP commands are case-insensitive
    ArgList line(ToLower(std::string(ptr)), " \t");
    if (line.Nargs() > 0) {
      if (line[0] == "logfile")
        isLeaprc = true;
      else if (line[0] == "source")
        isLeaprc = true;
      else if (line[0] == "addatomtypes")
        isLeaprc = true;
      // TODO more commands
    }
  }
  infile.CloseFile();
  return isLeaprc;
}

// DataIO_LeapRC::ReadHelp()
void DataIO_LeapRC::ReadHelp()
{

}

// DataIO_LeapRC::processReadArgs()
int DataIO_LeapRC::processReadArgs(ArgList& argIn)
{

  return 0;
}

/** Check if file already loaded. */
bool DataIO_LeapRC::check_already_loaded(Sarray const& files, std::string const& filename) {
  for (Sarray::const_iterator it = files.begin(); it != files.end(); ++it)
    if (*it == filename) return true;
  return false;
}

/** First look for filename, then look for AMBERHOME/dir/filename. */
std::string DataIO_LeapRC::find_path(std::string const& filename,
                                     std::string const& dir)
const
{
  if (File::Exists( filename )) return filename;
  if (amberhome_.empty()) {
    mprinterr("Error: '%s' not found.\n", filename.c_str());
    return filename;
  }
  std::string amberpath = amberhome_ + dir + filename;
  if (File::Exists( amberpath )) return amberpath;
  mprinterr("Error: '%s' not found.\n", amberpath.c_str());
  return amberpath;
}

/** LEaP loadAmberParams command. */
int DataIO_LeapRC::LoadAmberParams(std::string const& filename, DataSetList& dsl,
                                   std::string const& dsname,
                                   NHarrayType const& atomHybridizations)
const
{
  DataSet* paramSet = 0;
  // TODO detect this better
  ArgList fargs( filename, "." );
  if (fargs.hasKey("frcmod")) {
    mprintf("\tLoading force field modifications from '%s'\n", filename.c_str());
    DataIO_AmberFrcmod infile;
    infile.SetDebug( debug_ );
    if (infile.ReadData( find_path(filename, "parm/"), dsl, dsname)) {
      mprinterr("Error: Could not load force field modifications from '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() != 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Expected only 1 parameter set added, got %u\n", infile.Nadded());
      return 1;
    }
    paramSet = infile.added_back();
  } else {
    if (check_already_loaded(paramFiles_, filename)) {
      mprintf("Warning: Force field %s has already been loaded, skipping.\n", filename.c_str());
      return 0;
    } else {
      mprintf("\tLoading force field from '%s'\n", filename.c_str());
      DataIO_AmberFF infile;
      if (infile.ReadData( find_path(filename, "parm/"), dsl, dsname)) {
        mprinterr("Error: Could not load force field from '%s'\n", filename.c_str());
        return 1;
      }
      paramFiles_.push_back( filename );
      if (infile.Nadded() != 1) {
        mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Expected only 1 parameter set added, got %u\n", infile.Nadded());
        return 1;
      }
      paramSet = infile.added_back();
    }
  }
  if (paramSet == 0) {
    mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Parameter set is null.\n");
    return 1;
  }
  // Update hybridizations for parameter atom types
  //for (DataIO::set_iterator ds = paramDSL.begin(); ds != paramDSL.end(); ++ds)
  //{
    if ( paramSet->Type() == DataSet::PARAMETERS ) {
      DataSet_Parameters& param = static_cast<DataSet_Parameters&>( *paramSet );
      mprintf("\tUpdating atom hybridizations in set %s\n", param.legend());
      for (ParmHolder<AtomType>::iterator it = param.AT().begin();
                                          it != param.AT().end(); ++it)
      {
        NHarrayType::const_iterator ah = atomHybridizations.find( it->first[0] );
        if (ah == atomHybridizations.end())
          mprintf("Warning: No hybridization set for atom type '%s'\n", *(it->first[0]));
        else
          it->second.SetHybridization( ah->second );
      }
    } else {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Set %s is not parameter set.\n",
                paramSet->legend());
    }
  //}
  return 0;
}

/** LEaP loadOff command. */
int DataIO_LeapRC::LoadOFF(std::string const& filename, DataSetList& dsl, std::string const& dsname, DSarray& units) const {
  if (check_already_loaded(libFiles_, filename)) {
    mprintf("Warning: Library %s has already been loaded, skipping.\n", filename.c_str());
  } else {
    DataIO_AmberLib infile;
    infile.SetDebug( debug_ );
    // Allow lib to overwrite e.g. something from previous prep
    ArgList tmpArgs("allowoverwrite");
    infile.processReadArgs(tmpArgs);
    if (infile.ReadData( find_path(filename, "lib/"), dsl, dsname)) {
      mprinterr("Error: Could not load library file '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() < 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadOFF(): No unit sets added.\n");
      return 1;
    }
    for (DataIO::set_iterator it = infile.added_begin(); it != infile.added_end(); ++it)
      units.push_back( *it );
    libFiles_.push_back( filename );
  }
  return 0;
}

/** LEaP loadAmberPrep command. */
int DataIO_LeapRC::LoadAmberPrep(std::string const& filename, DataSetList& dsl, std::string const& dsname, DSarray& units) const {
  if (check_already_loaded(libFiles_, filename)) {
    mprintf("Warning: Prep file %s has already been loaded, skipping.\n", filename.c_str());
  } else {
    DataIO_AmberPrep infile;
    infile.SetDebug( debug_ );
    if (infile.ReadData( find_path(filename, "prep/"), dsl, dsname)) {
      mprinterr("Error: Could not load prep file '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() < 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberPrep(): No unit sets added.\n");
      return 1;
    }
    for (DataIO::set_iterator it = infile.added_begin(); it != infile.added_end(); ++it)
      units.push_back( *it );
    libFiles_.push_back( filename );
  }
  return 0;
}

/** LEaP addAtomTypes command. */
int DataIO_LeapRC::AddAtomTypes(NHarrayType& atomHybridizations, BufferedLine& infile)
const
{
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0) mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // Some entries (like LP and EP) are not required to have elements.
          // Set the hybridization index to 1 or 2.
          int hidx;
          if (aline.Nargs() == 3)
            hidx = 2;
          else if (aline.Nargs() == 2)
            hidx = 1;
          else {
            mprinterr("Error: Malformed addAtomTypes entry %s\n", tmp.c_str());
            return 1;
          }
          AtomType::HybridizationType ht;
          if (aline[hidx] == "sp3")
            ht = AtomType::SP3;
          else if (aline[hidx] == "sp2")
            ht = AtomType::SP2;
          else if (aline[hidx] == "sp")
            ht = AtomType::SP;
          else {
            mprintf("Warning: Unknown hybridization in addAtomTypes entry %s\n", tmp.c_str());
            ht = AtomType::UNKNOWN_HYBRIDIZATION;
          }
          NameType atype(aline[0]);
          NHarrayType::iterator it = atomHybridizations.lower_bound( atype );
          if (it == atomHybridizations.end() || it->first != atype) {
            it = atomHybridizations.insert( it, NHpairType(atype, ht) );
          } else {
            mprintf("Warning: Duplicate entry for '%s' in addAtomTypes.", *atype);
            if (it->second != ht) {
              mprintf(" Overwriting.\n");
              mprintf("Warning: Line is %s\n", tmp.c_str());
              it->second = ht;
            } else
              mprintf("\n");
          }
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addAtomTypes command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
      //mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addAtomTypes command.\n");
    return 1;
  }
  mprintf("\tRead %zu atom hybridizations.\n", atomHybridizations.size());
  return 0;
}

/** Find unit among loaded units. */
/*DataSet* DataIO_LeapRC::findUnit(std::string const& unitName) const {
  for (DSarray::const_iterator it = units_.begin(); it != units_.end(); ++it)
  {
    if ( (*it)->Meta().Aspect() == unitName ) {
      mprintf("DEBUG: findUnit: Aspect match: %s %s\n", unitName.c_str(), (*it)->legend());
      return *it;
    }
  }
  for (DSarray::const_iterator it = units_.begin(); it != units_.end(); ++it)
  {
    if ( (*it)->Meta().Name() == unitName ) { 
      mprintf("DEBUG: findUnit: Name match: %s %s\n", unitName.c_str(), (*it)->legend());
      return *it;
    }
  }
  return 0;
}*/

/** Associate PDB residue name map with COORDS set. */ // TODO call every time a new unit is loaded?
void DataIO_LeapRC::addPdbResMapToUnit(DataSet_Coords* unit, PdbResMapType const& prm)
const
{
  // Does this unit already have associated ResId?
  AssociatedData* ad = unit->GetAssociatedData( AssociatedData::RESID );
  if (ad == 0) {
    AssociatedData_ResId resid( prm.pdbName_, prm.termType_ );
    unit->AssociateData( &resid );
    //if (debug_ > 0) {
      mprintf("DEBUG: Found unit %s", unit->legend());
      resid.Ainfo();
      mprintf("\n");
    //}
  } else {
    // TODO check if map is different
    mprintf("DEBUG: Unit %s already has resmap. ", unit->legend());
    ad->Ainfo();
    mprintf("\n");
  }
}

/** LEaP addPdbResMap command. */
int DataIO_LeapRC::AddPdbResMap(BufferedLine& infile, PdbResMapArray& pdbResMapIn)
const
{
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0) mprintf("DEBUG: addPdbResMap: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // 3 tokens: terminal type (0=beg 1=end), PDB name, unit name
          if (aline.Nargs() < 2 || aline.Nargs() > 3) {
            mprinterr("Error: Malformed entry in addPdbResMap: %s\n", tmp.c_str());
            return 1;
          }
          Cpptraj::Structure::TerminalType termType = Cpptraj::Structure::NON_TERMINAL;
          int pdbidx = 0;
          int unitidx = 1;
          if (aline.Nargs() == 3) {
            if (aline[0] == "0")
              termType = Cpptraj::Structure::BEG_TERMINAL;
            else if (aline[0] == "1")
              termType = Cpptraj::Structure::END_TERMINAL;
            else
              mprintf("Warning: Unrecognized terminal type in addPdbResMap: %s\n", aline[0].c_str());
            pdbidx = 1;
            unitidx = 2;
          }
          //if (termType != Cpptraj::Structure::NON_TERMINAL) {
            PdbResMapType prm;
            prm.termType_ = termType;
            prm.pdbName_ = aline[pdbidx];
            prm.unitName_ = aline[unitidx];
            pdbResMapIn.push_back( prm );
            // Find among loaded units
            //DataSet* unitSet = findUnit( aline[unitidx] );
            //if (unitSet == 0) {
            //  mprintf("Warning: Unit '%s' was not found among loaded units.\n",
            //          aline[unitidx].c_str());
            //} else {
            //  //DataSet_Coords& crd = static_cast<DataSet_Coords&>( *ds );
            //  AssociatedData_ResId resid( aline[pdbidx], termType );
            //  unitSet->AssociateData( &resid );
            //  if (debug_ > 0) {
            //    mprintf("DEBUG: Found unit %s", unitSet->legend());
            //    resid.Ainfo();
            //    mprintf("\n");
            //  }
            //}
          //}
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addPdbResMap command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
    //mprintf("DEBUG: END OF LINE: addPdbResMap: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addPdbResMap command.\n");
    return 1;
  }
  return 0;
}
/** LEaP addPdbAtomMap command. */
int DataIO_LeapRC::AddPdbAtomMap(std::string const& dsname, DataSetList& DSL, BufferedLine& infile)
const
{
  MetaData meta(dsname, "atommap");
  DataSet* ds = DSL.CheckForSet(meta);
  if (ds == 0) {
    ds = DSL.AddSet(DataSet::NAMEMAP, meta);
    if (ds == 0) return 1;
  }
  DataSet_NameMap& namemap = static_cast<DataSet_NameMap&>( *ds );
  mprintf("DEBUG: Name map set: %s\n", namemap.legend());
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0)
            mprintf("DEBUG: addPdbAtomMap: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // 2 tokens: Old name, new name
          if (aline.Nargs() != 2) {
            mprinterr("Error: Malformed entry in addPdbAtomMap: %s\n", tmp.c_str());
            return 1;
          }
          if (debug_ > 0) mprintf("DEBUG: old= %s  new= %s\n", aline[0].c_str(), aline[1].c_str());
          namemap.AddNameMap( aline[0], aline[1] );
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addPdbAtomMap command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
    //mprintf("DEBUG: END OF LINE: addPdbResMap: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addPdbAtomMap command.\n");
    return 1;
  }
  return 0;
}

/** Load mol2 as a COORDS set */
int DataIO_LeapRC::LoadMol2(ArgList const& argIn, DataSetList& dsl) const {
  ArgList args( argIn.ArgLineStr(), " =" );
  //mprintf("DEBUG: LoadMol2\n");
  //args.PrintList();
  // Should be at least 3 args: NAME loadmol2 FILE
  if (args.Nargs() < 3) {
    mprinterr("Error: Expected at least <NAME> = loadmol2 <FILE>, got: %s\n", argIn.ArgLine());
    return 1;
  }
  DataIO_Coords coordsIn;
  coordsIn.SetDebug( debug_ );
  if (coordsIn.ReadData( args[2], dsl, args[0] )) {
    mprinterr("Error: Could not load structure from '%s' into '%s'\n",
              args[2].c_str(), args[0].c_str());
    return 1;
  }
  mprintf("\tLoaded file '%s' into '%s'\n",
          args[2].c_str(), args[0].c_str());
  return 0;
}

/** Load PDB and build it. */
int DataIO_LeapRC::LoadPDB(ArgList const& argIn, DataSetList& dsl) const {
  ArgList args( argIn.ArgLineStr(), " =" );
  // Should be at least 3 args: NAME loadpdb FILE
  if (args.Nargs() < 3) {
    mprinterr("Error: Expected at least <NAME> = loadpdb <FILE>, got: %s\n", argIn.ArgLine());
    return 1;
  }
  DataIO_Coords coordsIn;
  coordsIn.SetDebug( debug_ );
  DataSetList tmpdsl;
  if (coordsIn.ReadData( args[2], tmpdsl, args[0] )) {
    mprinterr("Error: Could not load structure from '%s' into '%s'\n",
              args[2].c_str(), args[0].c_str());
    return 1;
  }
  if (coordsIn.Nadded() != 1) {
    mprinterr("Internal Error: DataIO_LeapRC::LoadPDB(): Expected 1 COORDS set loaded from PDB, got %u\n", coordsIn.Nadded());
    return 1;
  }
  mprintf("\tLoaded file '%s' into '%s'\n", args[2].c_str(), args[0].c_str());

  // Assume we want the remaining args to be passed to the builder.
  args.MarkArg(0);
  args.MarkArg(1);
  args.MarkArg(2);
  ArgList tmparg = args.RemainingArgs();

  tmparg.AddArg("name " + args[0]);
  Exec_Build build;
  Exec::RetType ret = build.BuildStructure( coordsIn.added_back(), dsl, debug_, tmparg );
  if (ret == CpptrajState::ERR) {
    mprinterr("Error: Build of '%s' failed.\n", args[2].c_str());
    return 1;
  }
  return 0;
}

/** Save specified unit to a topology and restart file. */
int DataIO_LeapRC::SaveAmberParm(std::string const& unitName, ArgList& line, DataSetList const& dsl)
const
{
  // saveamberparm <unit> <topfile> <restartfile>
  if (unitName.empty()) {
    mprinterr("Error: 'saveamberparm' : Unit name missing.\n");
    return 1;
  }
  std::string topName = line.GetStringNext();
  if (topName.empty()) {
    mprinterr("Error: 'saveamberparm' : Topology name missing.\n");
    return 1;
  }
  std::string crdName = line.GetStringNext();
  if (crdName.empty()) {
    mprinterr("Error: 'saveamberparm' : Coordinates name missing.\n");
    return 1;
  }
  mprintf("\tSaving unit '%s' to topology '%s' and coordinates '%s'\n",
          unitName.c_str(), topName.c_str(), crdName.c_str());
  // Get the unit
  return 0;
}

/// Move sets from paramDSL to dsl
/*static inline int addSetsToList(DataSetList& dsl, DataSetList& paramDSL)
{
  // Add data sets to the main data set list
  for (DataSetList::const_iterator ds = paramDSL.begin(); ds != paramDSL.end(); ++ds) {
    DataSet* dtmp = dsl.CheckForSet( (*ds)->Meta() );
    if (dtmp != 0) {
      mprinterr("Error: Set '%s' already exists.\n", (*ds)->legend());
      return 1;
    }
    dsl.AddSet( *ds );
  }
  paramDSL.SetHasCopies( true );
  return 0;
}*/

/** Read (source) a leaprc (input) file. */
int DataIO_LeapRC::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  atomHybridizations_.clear();
  units_.clear();
  // First, need to determine where the Amber FF files are
  const char* env = getenv("AMBERHOME");
  if (env != 0)
    amberhome_ = std::string(env) + "/dat/leap/";
  else {
    mprintf("Warning: AMBERHOME is not set. Determining FF file location based on leaprc file.\n");
    // Try to guess based on where the leaprc file is
    FileName leapcmddir( fname.DirPrefix_NoSlash() );
    if (leapcmddir.Base() == "oldff") {
      FileName leapcmddir2( leapcmddir.DirPrefix_NoSlash() );
      amberhome_ = leapcmddir2.DirPrefix();
    } else
      amberhome_ = leapcmddir.DirPrefix();
  }
  if (!amberhome_.empty())
    mprintf("\tForce field files located in '%s'\n", amberhome_.c_str());

  if (Source(fname, dsl, dsname)) {
    mprinterr("Error: Could not read LEaP input from '%s'\n", fname.full());
    return 1;
  }

  return 0;
}

/** Execute leap source command */
int DataIO_LeapRC::Source(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead( find_path(fname.Full(), "cmd/") )) {
    mprinterr("Error: Could not open leaprc file '%s'\n", fname.full());
    return 1;
  }
  mprintf("\tReading LEaP input from '%s'\n", fname.base());
  enum LeapCmdType { LOADAMBERPARAMS = 0, LOADOFF, LOADAMBERPREP, ADDATOMTYPES,
                     ADDPDBRESMAP, ADDPDBATOMMAP, LOADMOL2, LOADPDB, SOURCE,
                     QUIT, SAVEAMBERPARM, UNKNOWN_CMD };
  //DataSetList paramDSL;
  //DataSetList unitDSL;
  //NHarrayType atomHybridizations;
  //PdbResMapArray pdbResMap;
  int err = 0;
  const char* ptr = infile.Line();
  // FIXME need to convert to all lowercase for matching commands; leap allows
  //       mixed case
  while (ptr != 0) {
    if (ptr[0] != '\0' && ptr[0] != '#') {
      // Note if this line contains an equals sign.
      bool has_equals = false;
      for (const char* p = ptr; *p != '\0'; ++p) {
        if (*p == '=') {
          has_equals = true;
          break;
        }
      }

      ArgList line( ptr, " =\t" );
      mprintf("\tLEAP> %s\n", ptr);

      LeapCmdType leapcmd = UNKNOWN_CMD;
      // Look through all args, lowercase, for recognized commands.
      int pos = -1;
      for (int arg = 0; arg != line.Nargs(); arg++) {
        std::string argStr = ToLower( line[arg] );
        if      (argStr == "loadamberparams") { pos = arg; leapcmd = LOADAMBERPARAMS; break; }
        else if (argStr == "loadoff"        ) { pos = arg; leapcmd = LOADOFF; break; }
        else if (argStr == "loadamberprep"  ) { pos = arg; leapcmd = LOADAMBERPREP; break; }
        else if (argStr == "addatomtypes"   ) { pos = arg; leapcmd = ADDATOMTYPES; break; }
        else if (argStr == "addpdbresmap"   ) { pos = arg; leapcmd = ADDPDBRESMAP; break; }
        else if (argStr == "addpdbatommap"  ) { pos = arg; leapcmd = ADDPDBATOMMAP; break; }
        else if (argStr == "loadmol2"       ) { pos = arg; leapcmd = LOADMOL2; break; }
        else if (argStr == "loadpdb"        ) { pos = arg; leapcmd = LOADPDB; break; }
        else if (argStr == "source"         ) { pos = arg; leapcmd = SOURCE; break; }
        else if (argStr == "quit"           ) { pos = arg; leapcmd = QUIT; break; }
        else if (argStr == "saveamberparm"  ) { pos = arg; leapcmd = SAVEAMBERPARM; break; }
      }

      err = 0;
      if (leapcmd == LOADAMBERPARAMS)
        err = LoadAmberParams(line.GetStringKey(line[pos]), dsl, dsname, atomHybridizations_ );
      else if (leapcmd == SAVEAMBERPARM)
        err = SaveAmberParm(line.GetStringKey(line[pos]), line, dsl);
      else if (leapcmd == LOADOFF)
        err = LoadOFF( line.GetStringKey(line[pos]), dsl, dsname, units_ );
      else if (leapcmd == LOADAMBERPREP)
        err = LoadAmberPrep( line.GetStringKey(line[pos]), dsl, dsname, units_ );
      else if (leapcmd == ADDATOMTYPES)
        err = AddAtomTypes(atomHybridizations_, infile);
      else if (leapcmd == ADDPDBRESMAP) {
        err = AddPdbResMap(infile, pdbResMap_);
        if (err == 0) {
          // Update units with pdb residue map info
          //for (DataSetList::const_iterator ds = unitDSL.begin(); ds != paramDSL.end(); ++ds)
          //{
          //  if ( (*ds)->Group() == DataSet::COORDINATES ) {
          for (PdbResMapArray::const_iterator it = pdbResMap_.begin();
                                              it != pdbResMap_.end(); ++it)
          {
            // Find the unit in unit DSL
            DataSet* ds = dsl.CheckForSet( MetaData(dsname, it->unitName_) );
            if (ds != 0 && ds->Group() == DataSet::COORDINATES) {
              addPdbResMapToUnit( static_cast<DataSet_Coords*>(ds), *it );
            }
            /*if (ds == 0) {
              mprintf("Warning: Unit '%s' was not found among loaded units.\n", it->unitName_.c_str());
            } else {
              DataSet_Coords& crd = static_cast<DataSet_Coords&>( *ds );
              AssociatedData_ResId resid( it->pdbName_, it->termType_ );
              crd.AssociateData( &resid );
              if (debug_ > 0) {
                mprintf("DEBUG: Found unit %s", crd.legend());
                resid.Ainfo();
                mprintf("\n");
              }
            }*/
          }
        }
      } else if (leapcmd == ADDPDBATOMMAP)
        err = AddPdbAtomMap(dsname, dsl, infile);
      else if (leapcmd == LOADMOL2)
        err = LoadMol2(line, dsl);
      else if (leapcmd == LOADPDB)
        err = LoadPDB(line, dsl);
      else if (leapcmd == SOURCE) {
        std::string fname1 = line.GetStringKey("source");
        if (fname1.empty()) {
          mprinterr("Error: No filename given for 'source'\n");
          return 1;
        } else if (fname1 == fname.Full()) {
          mprinterr("Error: File '%s' attempting to source itself '%s'\n", 
                    fname.full(), fname1.c_str());
          return 1;
        }
        err = Source(fname1, dsl, dsname);
      } else if (leapcmd == QUIT) {
        // Do not read any more.
        mprintf("\tEncountered 'quit' in leaprc file, not reading any more.\n");
        break;
      } else {
        // Unrecognized so far. See if this is a unit alias (interpret as 'alias = unit')
        if (has_equals) {
          if (line.Nargs() == 2) {
            mprintf("DEBUG: %s = %s\n", line[0].c_str(), line[1].c_str());
            // Find the unit to make a copy of
            DataSet* ds0 = dsl.CheckForSet( MetaData(dsname, line[1]) );
            if (ds0 == 0) {
              mprinterr("Error: Could not find unit '%s' to copy to '%s'\n", line[1].c_str(), line[0].c_str());
              return 1;
            }
            DataSet_Coords& crd0 = static_cast<DataSet_Coords&>( *ds0 );
            // Allocate copy
            DataSet* ds1 = dsl.AddSet( DataSet::COORDS, MetaData(dsname, line[0]) );
            if (ds1 == 0) {
              mprinterr("Error: Could not allocate unit '%s' for '%s'\n", line[0].c_str(), line[1].c_str());
              return 1;
            }
            DataSet_Coords& crd1 = static_cast<DataSet_Coords&>( *ds1 );
            if (crd1.CoordsSetup( crd0.Top(), crd0.CoordsInfo() )) {
              mprinterr("Error: Could not set up unit '%s' for '%s'\n", line[0].c_str(), line[1].c_str());
              return 1;
            }
            crd1.Allocate( DataSet::SizeArray(1, 1) );
            // Copy
            Frame tmpFrm = crd0.AllocateFrame();
            crd0.GetFrame(0, tmpFrm);
            crd1.SetCRD(0, tmpFrm );
            // Copy associated data
            //crd1.CopyAssociatedDataFrom( crd0 );
            mprintf("DEBUG: Created unit set %s\n", crd1.legend());
            // See if there is a PDB residue name map
            for (PdbResMapArray::const_iterator it = pdbResMap_.begin();
                                                it != pdbResMap_.end(); ++it)
            {
              if (it->unitName_ == line[0]) {
                addPdbResMapToUnit( &crd1, *it );
                break;
              }
            }
          }
        } else {
          mprintf("Warning: Skipping unhandled LEaP command line: %s\n", ptr);
        }
      }
    }
    if (err != 0) break;
    ptr = infile.Line();
  }
  infile.CloseFile();



  // Add data sets to the main data set list
  //if (addSetsToList(dsl, paramDSL)) return err+1;

  //if (addSetsToList(dsl, unitDSL)) return err+1;

  return err;
}

// DataIO_LeapRC::WriteHelp()
void DataIO_LeapRC::WriteHelp()
{

}

// DataIO_LeapRC::processWriteArgs()
int DataIO_LeapRC::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_LeapRC::WriteData()
int DataIO_LeapRC::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
