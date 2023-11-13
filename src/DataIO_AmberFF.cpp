#include "DataIO_AmberFF.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_Parameters.h"
#include <cstdio> // sscanf

/// CONSTRUCTOR
DataIO_AmberFF::DataIO_AmberFF()
{

}

// DataIO_AmberFF::ID_DataFormat()
bool DataIO_AmberFF::ID_DataFormat(CpptrajFile& infile)
{
  //if (infile.OpenFile()) return false;
  //std::string line = infile.GetLine(); // Title
  //infile.CloseFile();
  //bool isLib = (line == "!!index array str");
  return false;
}

// DataIO_AmberFF::ReadHelp()
void DataIO_AmberFF::ReadHelp()
{
  mprintf("\tnbset <nonbond set name> : Nonbonded set name to use when multiple nonbond parameter sets.\n");
}

// DataIO_AmberFF::processReadArgs()
int DataIO_AmberFF::processReadArgs(ArgList& argIn)
{
  nbsetname_ = argIn.GetStringKey("nbset");
  return 0;
}

/** Read symbols delimited by - and space. */
int DataIO_AmberFF::read_symbols(const char* ptrIn, std::vector<std::string>& symbols, int nsymbols)
{
  int isymbol = 0;
  bool char_has_been_read = false;
  for (const char* ptr = ptrIn; *ptr != '\0'; ++ptr)
  {
    if (*ptr == '-') {
      isymbol++;
      char_has_been_read = false;
    } else if (*ptr == ' ' && isymbol + 1 == nsymbols && char_has_been_read) {
      return (ptr - ptrIn);
    } else {
      symbols[isymbol] += *ptr;
      if (*ptr != ' ') char_has_been_read = true;
    }
  }
  return -1;
}

/// Hold a set of nonbonded parameters
class NonbondSet {
  public:
    NonbondSet(std::string const& n) : name_(n) {}

    std::string name_;          ///< Name of set parameters
    ParmHolder<LJparmType> LJ_; ///< Hold LJ 6-12 parameters
};

// DataIO_AmberFF::ReadData()
int DataIO_AmberFF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  static const int MAXSYMLEN = 16;
  // Allocate data set
  MetaData md( dsname );
  DataSet* ds = dsl.CheckForSet( md );
  if (ds != 0) {
    if (ds->Type() != DataSet::PARAMETERS) {
      mprinterr("Error: Set '%s' does not have parameters, cannot append.\n", ds->legend());
      return 1;
    }
    mprintf("\tAdding to existing set %s\n", ds->legend());
  } else {
    ds = dsl.AddSet( DataSet::PARAMETERS, md );
    if (ds == 0) return 1;
  }
  DataSet_Parameters& prm = static_cast<DataSet_Parameters&>( *ds ); 

  typedef std::vector<NonbondSet> NbSetArrayType;
  NbSetArrayType NBsets;

  // Read title
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read anything from Amber FF file %s\n", fname.full());
    return 1;
  }
  std::string title(ptr);
  mprintf("\tTitle: %s\n", title.c_str());
  // Read file
  bool ljedit = false;
  enum SectionType { ATYPE = 0, HYDROPHILIC, BOND, ANGLE, DIHEDRAL, IMPROPER, 
                     LJ1012, NB_EQUIV, NONBOND, UNKNOWN };
  SectionType section = ATYPE;
  ptr = infile.Line();
  while (ptr != 0) {
    // Advance to first non-space char
    while (*ptr == ' ' && *ptr != '\0') ++ptr;
    mprintf("DEBUG: First char: %c (%i)\n", *ptr, (int)*ptr);
    if (*ptr == '\0') {
      // Section Change
      if (section != UNKNOWN) {
        if (section == NONBOND) {
          // Do a lookahead to see if there are multiple NB sets.
          // It will either be another set, END, or LJEDIT
          ptr = infile.Line();
          while (*ptr == ' ' && *ptr != '\0') ++ptr;
          std::string nbline(ptr);
          if (nbline == "END") {
            mprintf("END\n");
            section = UNKNOWN;
          } else if (nbline == "LJEDIT") {
            ljedit = true;
            section = UNKNOWN;
          } // Otherwise assume another nonbond section
        } else {
          mprintf("SECTION %i change to %i\n", (int)section, (int)section + 1);
          section = (SectionType)((int)section + 1);
        }
      }
      // Special cases
      if (section == HYDROPHILIC) {
        // Special case: hydrophilic atom types. One line.
        // FORMAT(20(A2,2X))
        ptr = infile.Line();
        mprintf("DEBUG: Hydrophilic: %s\n", ptr);
        // Take advantage of the fact that we expect whitespace-delimiters
        ArgList hsymbols( ptr, " " );
        for (int iarg = 0; iarg != hsymbols.Nargs(); ++iarg)
          prm.AddHydrophilicAtomType( NameType( hsymbols[iarg] ) );
        section = (SectionType)((int)section + 1);
        // Look ahead. Older parm files have the hydrophilic section delimited
        // with a newline.
        ptr = infile.Line();
        while (*ptr == ' ' && *ptr != '\0') ++ptr;
        if (*ptr != '\0') continue;
      } else if (section == NONBOND) {
        // Special case: first read the line.
        // LABEL , KINDNB
        // FORMAT(A4,6X,A2)
        if (NBsets.empty())
          ptr = infile.Line();
        char nb_label[MAXSYMLEN], nb_kind[MAXSYMLEN];
        sscanf(ptr, "%s %s", nb_label, nb_kind);
        mprintf("DEBUG: NB label= %s  NB kind = %s\n", nb_label, nb_kind);
        if (nb_kind[0] != 'R' || nb_kind[1] != 'E') {
          mprinterr("Error: Nonbond parameters are not of type 'RE' (Rmin, Epsilon).\n");
          return 1;
        }
        NBsets.push_back( NonbondSet( std::string(nb_label) ) );
      }
    } else if (section == ATYPE) {
      // Input for atom symbols and masses
      // Format (A2,2X,F10.2x,f10.2)
      mprintf("DEBUG: Atype: %s\n", ptr);
      char kndsym[MAXSYMLEN];
      double amass;
      double atpol;
      int nscan = sscanf(ptr, "%s %lf %lf", kndsym, &amass, &atpol);
      ParameterHolders::RetType ret;
      if (nscan == 3) {
        ret = prm.AT().AddParm( TypeNameHolder(kndsym),
                                AtomType(amass, atpol),
                                true );
      } else if (nscan == 2) {
        // Only mass
        ret = prm.AT().AddParm( TypeNameHolder(kndsym),
                                AtomType(amass),
                                true );
      } else {
        mprinterr("Error: Expected atom type, mass, polarizability, got only %i columns.\n", nscan);
        return 1;
      }
      if (ret == ParameterHolders::UPDATED)
        mprintf("Warning: Redefining atom type %s\n", kndsym);
    } else if (section == BOND) {
      // Bond parameters
      // IBT , JBT , RK , REQ
      // FORMAT(A2,1X,A2,2F10.2)
      mprintf("DEBUG: Bond: %s\n", ptr);
      std::vector<std::string> symbols(2);
      int pos = read_symbols(ptr, symbols, 2);
      if (pos < 0) {
        mprinterr("Error: Could not read symbols for bond from %s\n", ptr);
        return 1;
      }
      mprintf("DEBUG: %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), ptr+pos);
      double RK, REQ;
      int nscan = sscanf(ptr+pos, "%lf %lf", &RK, &REQ);
      if (nscan != 2) {
        mprinterr("Error: Expected RK, REQ, got only %i elements\n", nscan);
        return 1;
      }
      TypeNameHolder types(2);
      types.AddName( symbols[0] );
      types.AddName( symbols[1] );
      prm.BP().AddParm(types, BondParmType(RK, REQ), false);
    } else if (section == ANGLE) {
      // Angle parameters
      // ITT , JTT , KTT , TK , TEQ
      // FORMAT(A2,1X,A2,1X,A2,2F10.2)
      mprintf("DEBUG: Angle: %s\n", ptr);
      std::vector<std::string> symbols(3);
      int pos = read_symbols(ptr, symbols, 3);
      if (pos < 0) {
        mprinterr("Error: Could not read symbols for angle from %s\n", ptr);
        return 1;
      }
      mprintf("DEBUG: %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), ptr+pos);
      double TK, TEQ;
      int nscan = sscanf(ptr+pos, "%lf %lf", &TK, &TEQ);
      if (nscan != 2) {
        mprinterr("Error: Expected TK, TEQ, got only %i elements\n", nscan);
        return 1;
      }
      TypeNameHolder types(3);
      types.AddName( symbols[0] );
      types.AddName( symbols[1] );
      types.AddName( symbols[2] );
      prm.AP().AddParm(types, AngleParmType(TK, TEQ), false);
    } else if (section == DIHEDRAL) {
      // Dihedral parameters
      // IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
      // FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
      // If IPT .eq. 'X ' .and. LPT .eq. 'X ' then any dihedrals in the
      // system involving the atoms "JPT" and and "KPT" are assigned 
      // the same parameters.  This is called the general dihedral type
      // and is of the form "X "-"JPT"-"KPT"-"X ".
      // IDIVF is the factor by which the torsional barrier is divided.
      // Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
      // details. Basically, the actual torsional potential is
      //   (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
      mprintf("DEBUG: Dihedral: %s\n", ptr);
      std::vector<std::string> symbols(4);
      int pos = read_symbols(ptr, symbols, 4);
      if (pos < 0) {
        mprinterr("Error: Could not read symbols for dihedral from %s\n", ptr);
        return 1;
      }
      mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
      int IDIVF;
      double PK, PHASE, PN;
      int nscan = sscanf(ptr+pos, "%i %lf %lf %lf", &IDIVF, &PK, &PHASE, &PN);
      if (nscan != 4) {
        mprinterr("Error: Expected IDIVF, PK, PHASE, PN, got only %i elements\n", nscan);
        return 1;
      }
      TypeNameHolder types(4);
      types.AddName( symbols[0] );
      types.AddName( symbols[1] );
      types.AddName( symbols[2] );
      types.AddName( symbols[3] );
      prm.DP().AddParm(types, DihedralParmType(PK / (double)IDIVF, PN, PHASE), false);
    } else if (section == IMPROPER) {
      // Improper parameters
      // IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
      // FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
      mprintf("DEBUG: Improper: %s\n", ptr);
      std::vector<std::string> symbols(4);
      int pos = read_symbols(ptr, symbols, 4);
      if (pos < 0) {
        mprinterr("Error: Could not read symbols for improper from %s\n", ptr);
        return 1;
      }
      mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
      double PK, PHASE, PN;
      int nscan = sscanf(ptr+pos, "%lf %lf %lf", &PK, &PHASE, &PN);
      if (nscan != 3) {
        mprinterr("Error: Expected PK, PHASE, PN, got only %i elements\n", nscan);
        return 1;
      }
      TypeNameHolder types(4);
      types.AddName( symbols[0] );
      types.AddName( symbols[1] );
      types.AddName( symbols[2] );
      types.AddName( symbols[3] );
      prm.IP().AddParm(types, DihedralParmType(PK, PN, PHASE), false);
    } else if (section == LJ1012) {
      // Lennard-Jones 10-12 hydrogen bonding term.
      // According to the docs, the format should be:
      // KT1 , KT2 , A , B , ASOLN , BSOLN , HCUT , IC
      // FORMAT(2X,A2,2X,A2,2x,5F10.2,I2)
      // In practive (in e.g. parm91.dat), the actual format appears to be:
      // KT1, KT2, A, B, HCUT
      char KT1[MAXSYMLEN], KT2[MAXSYMLEN];
      double A, B, HCUT;
      //double ASOLN, BSOLN, HCUT;
      //int IC;
      mprintf("DEBUG: LJ 10-12: %s\n", ptr);
      //int nscan = sscanf(ptr, "%s %s %lf %lf %lf %lf %lf %i\n",
      int nscan = sscanf(ptr, "%s %s %lf %lf %lf", KT1, KT2, &A, &B, &HCUT);
      if (nscan < 5) {
        mprinterr("Error: Expected at least KT1, KT2, A, B, HCUT, got only %i elements\n", nscan);
        return 1;
      }
      TypeNameHolder types(2);
      types.AddName( KT1 );
      types.AddName( KT2 );
      prm.HB().AddParm(types, HB_ParmType(A, B, HCUT), false);
    } else if (section == NB_EQUIV) {
      mprintf("DEBUG: Nonbond equiv: %s\n", ptr);
    } else if (section == NONBOND) {
      // ***** ONLY IF KINDNB .EQ. 'RE' *****
      // LTYNB , R , EDEP
      mprintf("DEBUG: Nonbond: %s\n", ptr);
      char LTYNB[MAXSYMLEN];
      double R, EDEP;
      double R14, E14;
      int nscan = sscanf(ptr, "%s %lf %lf %lf %lf", LTYNB, &R, &EDEP, &R14, &E14);
      
      if (nscan == 3 || nscan == 5) {
        // symbol, rmin, epsilon
        NBsets.back().LJ_.AddParm( TypeNameHolder(LTYNB), LJparmType(R, EDEP), false );
        ParmHolder<AtomType>::iterator it = prm.AT().GetParam( TypeNameHolder(LTYNB) );
        if (it == prm.AT().end()) {
          mprintf("Warning: Nonbond parameters defined for previously undefined type '%s'."
                          " Skipping.\n", LTYNB);
        } else {
          it->second.SetLJ().SetRadius( R );
          it->second.SetLJ().SetDepth( EDEP );
        }
      }
    } else if (ljedit) {
      mprintf("DEBUG: LJedit: %s\n", ptr);
    }
      
    ptr = infile.Line();
  } // END loop over file.
  // Deal with nonbond and equivalence
  if (!NBsets.empty()) {
    int nbsetidx = 0;
    if (NBsets.size() > 1 || !nbsetname_.empty()) {
      // Choose from multiple nbsets
      if (nbsetname_.empty()) {
        mprinterr("Error: Parm set in file '%s' contains %zu nonbonded parameter sets\n"
                  "Error:  but no set has been specified.\n", fname.full(), NBsets.size());
        mprinterr("Error: Need to specify one of");
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
          mprinterr(" %s", it->name_.c_str());
        mprinterr("\n");
      } else {
        nbsetidx = -1;
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it) {
          if (it->name_ == nbsetname_) {
            nbsetidx = it - NBsets.begin();
            break;
          }
        }
        if (nbsetidx < 0) {
          mprinterr("Error: Nonbonded set '%s' not found.\n", nbsetname_.c_str());
          mprinterr("Error: Need to specify one of");
          for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
            mprinterr(" %s", it->name_.c_str());
          mprinterr("\n");
          return 1;
        }
      }
    }
    mprintf("\tUsing nonbonded parm set: %s\n", NBsets[nbsetidx].name_.c_str());
  }

  prm.Debug(); // TODO debug level
  infile.CloseFile();
  return 0;
}

// DataIO_AmberFF::WriteHelp()
void DataIO_AmberFF::WriteHelp()
{

}

// DataIO_AmberFF::processWriteArgs()
int DataIO_AmberFF::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberFF::WriteData()
int DataIO_AmberFF::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
