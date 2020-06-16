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

    /// \return Index corresponding to given energy term name
    unsigned int getTermIdx(std::string const&);
    /// Add value to the specified set.
    int AddData(unsigned int, double, std::string const&, DataSetList&);
    /// Parse amber energy terms from given line
    int GetAmberEterms(const char*, std::string const&, DataSetList&);

    /// This array holds the aspect for each current field set.
    Sarray SetAspects_;
    /// This array holds the offset for each current field set (if applicable).
    Darray SetOffsets_;
    /// This array holds the DataSets for each current field.
    DSarray EneSets_;
    /// Map field names to indices into energy sets.
    NameIdxMap termIdxMap_;
    /// Time step for this file (MD only)
    double dt_;
    /// Initial time for this file (MD only)
    double t0_;
    /// Imin value for this file; determines MD, minimization, or post-process
    int imin_;
    /// Minimization step (imin==1 or 5 only)
    int minStep_;
    /// MD step (imin==0 only)
    int nstep_;
    /// Write frequency in steps (imin==1 or 0 only)
    int ntpr_;
    /// False unless we are actually reading data
    bool reachedNstep_;
};
#endif
