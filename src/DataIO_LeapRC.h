#ifndef INC_DATAIO_LEAPRC_H
#define INC_DATAIO_LEAPRC_H
#include "DataIO.h"
#include "Structure/StructureEnum.h"
#include <map>
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
    typedef std::map<NameType, AtomType::HybridizationType> NHarrayType;
    typedef std::vector<std::string> Sarray;

    struct PdbResMapType {
      std::string unitName_;
      NameType pdbName_;
      Cpptraj::Structure::TerminalType termType_;
    };
    typedef std::vector<PdbResMapType> PdbResMapArray;

    int LoadAmberParams(std::string const&, DataSetList&, std::string const&) const;
    int LoadOFF(std::string const&, DataSetList&, std::string const&) const;
    int LoadAmberPrep(std::string const&, DataSetList&, std::string const&) const;
    int AddAtomTypes(NHarrayType&, BufferedLine&) const;
    int AddPdbResMap(PdbResMapArray&, BufferedLine&) const;
    int AddPdbAtomMap(std::string const&, DataSetList&, BufferedLine&) const;
    int LoadMol2(ArgList const&, DataSetList&) const;
    /// Used to check if a parm/lib file was already loaded.
    static inline bool check_already_loaded(Sarray const&, std::string const&);
    /// \return either file or Amberhome/dir/file
    std::string find_path(std::string const&, std::string const&) const;

    std::string amberhome_;
    static Sarray paramFiles_; ///< Track amber FF param files loaded from leaprc files
    static Sarray libFiles_;   ///< Track amber library/prep files loaded from leaprc files
};
#endif
