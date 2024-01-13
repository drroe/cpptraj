#include "DataIO_LeapRC.h"
#include "AssociatedData_ResId.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataIO_AmberFF.h"
#include "DataIO_AmberFrcmod.h"
#include "DataIO_AmberLib.h"
#include "DataIO_AmberPrep.h"
#include "DataSet_Parameters.h"
#include <cstdlib> //getenv

/// CONSTRUCTOR
DataIO_LeapRC::DataIO_LeapRC()
{

}

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
    ArgList line(ptr, " \t");
    if (line.Nargs() > 0) {
      if (line[0] == "logFile" || line[0] == "logfile")
        isLeaprc = true;
      else if (line[0] == "source")
        isLeaprc = true;
      else if (line[0] == "addAtomTypes" || line[0] == "addatomtypes")
        isLeaprc = true;
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

/** LEaP loadAmberParams command. */
int DataIO_LeapRC::LoadAmberParams(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  // TODO detect this better
  ArgList fargs( filename, "." );
  if (fargs.hasKey("frcmod")) {
    mprintf("\tLoading force field modifications from '%s'\n", filename.c_str());
    DataIO_AmberFrcmod infile;
    if (infile.ReadData(amberhome_ + "parm/" + filename, dsl, dsname)) {
      mprinterr("Error: Could not load force field modifications from '%s'\n", filename.c_str());
      return 1;
    }
  } else {
    mprintf("\tLoading force field from '%s'\n", filename.c_str());
    DataIO_AmberFF infile;
    if (infile.ReadData(amberhome_ + "parm/" + filename, dsl, dsname)) {
      mprinterr("Error: Could not load force field from '%s'\n", filename.c_str());
      return 1;
    }
  }
  return 0;
}

/** LEaP loadOff command. */
int DataIO_LeapRC::LoadOFF(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  DataIO_AmberLib infile;
  if (infile.ReadData(amberhome_ + "lib/" + filename, dsl, dsname)) {
    mprinterr("Error: Could not load library file '%s'\n", filename.c_str());
    return 1;
  }
  return 0;
}

/** LEaP loadAmberPrep command. */
int DataIO_LeapRC::LoadAmberPrep(std::string const& filename, DataSetList& dsl, std::string const& dsname) const {
  DataIO_AmberPrep infile;
  if (infile.ReadData(amberhome_ + "lib/" + filename, dsl, dsname)) {
    mprinterr("Error: Could not load prep file '%s'\n", filename.c_str());
    return 1;
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
      if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
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

/** LEaP addPdbResMap command. */
int DataIO_LeapRC::AddPdbResMap(PdbResMapArray& pdbResMap, BufferedLine& infile)
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
      if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          mprintf("DEBUG: addPdbResMap: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // 3 tokens: terminal type (0=beg 1=end), PDB name, unit name
          if (aline.Nargs() != 3) {
            mprinterr("Error: Malformed entry in addPdbResMap: %s\n", tmp.c_str());
            return 1;
          }
          Cpptraj::Structure::TerminalType termType = Cpptraj::Structure::NON_TERMINAL;
          if (aline[0] == "0")
            termType = Cpptraj::Structure::BEG_TERMINAL;
          else if (aline[0] == "1")
            termType = Cpptraj::Structure::END_TERMINAL;
          else
            mprintf("Warning: Unrecognized terminal type in addPdbResMap: %s\n", aline[0].c_str());
          if (termType != Cpptraj::Structure::NON_TERMINAL) {
            PdbResMapType prm;
            prm.termType_ = termType;
            prm.pdbName_ = aline[1];
            prm.unitName_ = aline[2];
            pdbResMap.push_back( prm );
          }
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
      //mprintf("DEBUG: addPdbResMap: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addPdbResMap command.\n");
    return 1;
  }
  return 0;
}

/// Move sets from paramDSL to dsl
static inline int addSetsToList(DataSetList& dsl, DataSetList& paramDSL)
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
}

// DataIO_LeapRC::ReadData()
int DataIO_LeapRC::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
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
  mprintf("\tForce field files located in '%s'\n", amberhome_.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) {
    mprinterr("Error: Could not open leaprc file '%s'\n", fname.full());
    return 1;
  }
  DataSetList paramDSL;
  DataSetList unitDSL;
  NHarrayType atomHybridizations;
  PdbResMapArray pdbResMap;
  int err = 0;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (ptr[0] != '\0' && ptr[0] != '#') {
      ArgList line( ptr, " \t" );
      if (line.Contains("loadAmberParams"))
        err = LoadAmberParams( line.GetStringKey("loadAmberParams"), paramDSL, dsname );
      else if (line.Contains("loadamberparams"))
        err = LoadAmberParams( line.GetStringKey("loadamberparams"), paramDSL, dsname );
      else if (line.Contains("loadOff"))
        err = LoadOFF( line.GetStringKey("loadOff"), unitDSL, dsname );
      else if (line.Contains("loadoff"))
        err = LoadOFF( line.GetStringKey("loadoff"), unitDSL, dsname );
      else if (line.Contains("loadAmberPrep"))
        err = LoadAmberPrep( line.GetStringKey("loadAmberPrep"), unitDSL, dsname );
      else if (line.Contains("loadamberprep"))
        err = LoadAmberPrep( line.GetStringKey("loadamberprep"), unitDSL, dsname );
      else if (line.Contains("addAtomTypes") || line.Contains("addatomtypes"))
        err = AddAtomTypes(atomHybridizations, infile);
      else if (line.Contains("addPdbResMap") || line.Contains("addpdbresmap"))
        err = AddPdbResMap(pdbResMap, infile);
      else {
        // Does this line contain an equals sign?
        bool has_equals = false;
        for (const char* p = ptr; *p != '\0'; ++p) {
          if (*p == '=') {
            has_equals = true;
            break;
          }
        }
        // See if this is a unit alias (interpret as 'alias = unit')
        if (has_equals) {
          ArgList equals(ptr, " =\t");
          if (equals.Nargs() == 2) {
            mprintf("DEBUG: %s = %s\n", equals[0].c_str(), equals[1].c_str());
            // Find the unit to make a copy of
            DataSet* ds0 = unitDSL.CheckForSet( MetaData(dsname, equals[1]) );
            if (ds0 == 0) {
              mprinterr("Error: Could not find unit '%s' to copy to '%s'\n", equals[1].c_str(), equals[0].c_str());
              return 1;
            }
            DataSet_Coords& crd0 = static_cast<DataSet_Coords&>( *ds0 );
            // Allocate copy
            DataSet* ds1 = unitDSL.AddSet( DataSet::COORDS, MetaData(dsname, equals[0]) );
            if (ds1 == 0) {
              mprinterr("Error: Could not allocate unit '%s' for '%s'\n", equals[0].c_str(), equals[1].c_str());
              return 1;
            }
            DataSet_Coords& crd1 = static_cast<DataSet_Coords&>( *ds1 );
            if (crd1.CoordsSetup( crd0.Top(), crd0.CoordsInfo() )) {
              mprinterr("Error: Could not set up unit '%s' for '%s'\n", equals[0].c_str(), equals[1].c_str());
              return 1;
            }
            crd1.Allocate( DataSet::SizeArray(1, 1) );
            // Copy
            Frame tmpFrm = crd0.AllocateFrame();
            crd0.GetFrame(0, tmpFrm);
            crd1.SetCRD(0, tmpFrm );
            // Copy associated data
            crd1.CopyAssociatedDataFrom( crd0 );
            mprintf("DEBUG: Created unit set %s\n", crd1.legend());
          }
        }
      }
    }
    if (err != 0) break;
    ptr = infile.Line();
  }
  infile.CloseFile();

  // Update hybridizations for parameter atom types
  for (DataSetList::const_iterator ds = paramDSL.begin(); ds != paramDSL.end(); ++ds)
  {
    if ( (*ds)->Type() == DataSet::PARAMETERS ) {
      DataSet_Parameters& param = static_cast<DataSet_Parameters&>( *(*ds) );
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
    }
  }
  // Update units with pdb residue map info
  //for (DataSetList::const_iterator ds = unitDSL.begin(); ds != paramDSL.end(); ++ds)
  //{
  //  if ( (*ds)->Group() == DataSet::COORDINATES ) {
  for (PdbResMapArray::const_iterator it = pdbResMap.begin();
                                      it != pdbResMap.end(); ++it)
  {
    // Find the unit in unit DSL
    DataSet* ds = unitDSL.CheckForSet( MetaData(dsname, it->unitName_) );
    if (ds == 0) {
      mprintf("Warning: Unit '%s' was not found among loaded units.\n", it->unitName_.c_str());
    } else {
      DataSet_Coords& crd = static_cast<DataSet_Coords&>( *ds );
      AssociatedData_ResId resid( it->pdbName_, it->termType_ );
      crd.AssociateData( &resid );
      mprintf("DEBUG: Found unit %s", crd.legend());
      resid.Ainfo();
      mprintf("\n");
    }
  }

  // Add data sets to the main data set list
  if (addSetsToList(dsl, paramDSL)) return err+1;

  if (addSetsToList(dsl, unitDSL)) return err+1;

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
