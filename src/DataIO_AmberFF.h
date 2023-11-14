#ifndef INC_DATAIO_AMBERFF_H
#define INC_DATAIO_AMBERFF_H
#include "DataIO.h"
/// <Enter description of DataIO_AmberFF here>
class DataIO_AmberFF : public DataIO {
  public:
    DataIO_AmberFF();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberFF(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    static int read_symbols(const char*, std::vector<std::string>&, int);
    int writeParameterSet(ParameterSet const&) const;

    std::string nbsetname_; ///< Nonbonded parameter set name to use
};
#endif
