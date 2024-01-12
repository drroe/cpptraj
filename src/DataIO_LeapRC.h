#ifndef INC_DATAIO_LEAPRC_H
#define INC_DATAIO_LEAPRC_H
#include "DataIO.h"
class BufferedLine;
/// Read parameters and units from a leap rc file 
class DataIO_LeapRC : public DataIO {
  public:
    DataIO_LeapRC();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_LeapRC(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    typedef std::pair<NameType, AtomType::HybridizationType> NHpairType;
    typedef std::vector<NHpairType> NHarrayType;

    int LoadAmberParams(std::string const&, DataSetList&, std::string const&) const;
    int LoadOFF(std::string const&, DataSetList&, std::string const&) const;
    int LoadAmberPrep(std::string const&, DataSetList&, std::string const&) const;
    int AddAtomTypes(NHarrayType&, BufferedLine&) const;

    std::string amberhome_;
};
#endif
