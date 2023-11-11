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

}

// DataIO_AmberFF::processReadArgs()
int DataIO_AmberFF::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberFF::ReadData()
int DataIO_AmberFF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
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
  enum SectionType { ATYPE = 0, HYDROPHILIC, UNKNOWN };
  SectionType section = ATYPE;
  ptr = infile.Line();
  while (ptr != 0) {
    if (*ptr == '\0') {
      // Section Change
      if (section != UNKNOWN) {
        section = (SectionType)((int)section + 1);
      }
      // Special case; hydrophilic atom types. One line.
      // FORMAT(20(A2,2X))
      if (section == HYDROPHILIC) {
        ptr = infile.Line();
        mprintf("DEBUG: Hydrophilic: %s\n", ptr);
        // Take advantage of the fact that we expect whitespace-delimiters
        ArgList hsymbols( ptr, " " );
        for (int iarg = 0; iarg != hsymbols.Nargs(); ++iarg)
          prm.AddHydrophilicAtomType( NameType( hsymbols[iarg] ) );
        /*
        char htype[5];
        htype[4] = '\0';
        std::string hline(ptr);
        for (unsigned int idx = 0; idx < hline.size(); idx+=4) {
          htype[0] = hline[idx  ];
          htype[1] = hline[idx+1];
          htype[2] = hline[idx+2];
          htype[3] = hline[idx+3];
          //mprintf("DEBUG:\t%s\n", htype);
          prm.AddHydrophilicAtomType( NameType(htype) );
        }*/
        section = (SectionType)((int)section + 1);
      }
    } else if (section == ATYPE) {
      // Input for atom symbols and masses
      // Format (A2,2X,F10.2x,f10.2)
      mprintf("DEBUG: Atype: %s\n", ptr);
      char kndsym[16];
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
        mprinterr("Error: Expected atom type, mass, polarizability, got only %i\n", nscan);
        return 1;
      }
      if (ret == ParameterHolders::UPDATED)
        mprintf("Warning: Redefining atom type %s\n", kndsym);
    }
    ptr = infile.Line();
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
