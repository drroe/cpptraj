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
    typedef std::vector<DataSet*> DSarray;

    struct PdbResMapType {
      std::string unitName_;
      NameType pdbName_;
      Cpptraj::Structure::TerminalType termType_;
    };
    typedef std::vector<PdbResMapType> PdbResMapArray;

    int LoadAmberParams(std::string const&, DataSetList&, std::string const&, NHarrayType const&) const;
    int LoadOFF(std::string const&, DataSetList&, std::string const&, DSarray&) const;
    int LoadAmberPrep(std::string const&, DataSetList&, std::string const&, DSarray&) const;
    int AddAtomTypes(NHarrayType&, BufferedLine&) const;
    int AddPdbResMap(BufferedLine&, PdbResMapArray&) const;
    int AddPdbAtomMap(std::string const&, DataSetList&, BufferedLine&) const;
    int LoadMol2(ArgList const&, DataSetList&) const;
    int LoadPDB(ArgList const&, DataSetList&) const;
    int SaveAmberParm(std::string const&, ArgList&, DataSetList const& dsl) const;
    int Source(FileName const&, DataSetList&, std::string const&);
    /// Add PDB residue map to COORDS unit
    void addPdbResMapToUnit(DataSet_Coords*, PdbResMapType const&, bool) const;
    /// Add PDB residue map to COORDS unit, no update
    void addPdbResMapToUnit(DataSet_Coords*, PdbResMapType const&) const;
    /// Add PDB residue map to COORDS unit, allow update
    void updatePdbResMapToUnit(DataSet_Coords*, PdbResMapType const&) const;
    /// \return Previously loaded unit set with given name
    //DataSet* findUnit(std::string const&) const;
    /// Used to check if a parm/lib file was already loaded.
    static inline bool check_already_loaded(Sarray const&, std::string const&);
    /// \return either file or Amberhome/dir/file
    std::string find_path(std::string const&, std::string const&) const;

    std::string amberhome_;
    NHarrayType atomHybridizations_; ///< Store hybridizations for atom types
    PdbResMapArray pdbResMap_; ///< Hold PDB residue name map
    DSarray units_;            ///< Hold COORDS sets which have been added as units
    static Sarray paramFiles_; ///< Track amber FF param files loaded from leaprc files
    static Sarray libFiles_;   ///< Track amber library/prep files loaded from leaprc files
};
#endif
