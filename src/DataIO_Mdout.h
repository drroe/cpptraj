#ifndef INC_DATAIO_MDOUT_H
#define INC_DATAIO_MDOUT_H
#include "DataIO.h"
#include <map>
/// Read energies from Amber MDOUT files.
class DataIO_Mdout : public DataIO {
  public:
    DataIO_Mdout();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Mdout(); }
    static void ReadHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&)   { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    typedef std::vector<std::string> Sarray;
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> DSarray;
    typedef std::map<std::string, unsigned int> NameIdxMap;
    typedef std::pair<std::string, unsigned int> NameIdxPair;
    /// These are common energy field types.
    enum FieldType { ETOT = 0, EPTOT, GMAX,        BOND,
                     ANGLE,    DIHED, VDWAALS,     EEL,        EGB,    EPB, ECAVITY, EDISPER,
                     VDW14,    EEL14, RESTRAINT,   EAMBER,     DENSITY,
                     RMS,      EKTOT, ESURF,       EAMD_BOOST, VOLUME, TEMP,
                     PRESS,    DVDL,  N_FIELDTYPES };
    /// These are the aspects for each FieldType
    static const char* fieldTypeStr_[];
    /// This array holds the aspect for each current field set.
    Sarray SetAspects_;
    /// This array holds the offset for each current field set (if applicable).
    Darray SetOffsets_;
    /// This array holds the DataSets for each current field.
    DSarray EneSets_;
    /// Map field names to indices into energy sets.
    NameIdxMap termIdxMap_;

    static FieldType getEindex(Sarray const&);
};
#endif
