#ifndef INC_HB_HBDATA_H
#define INC_HB_HBDATA_H
#include "Bridge.h"
#include "Hbond.h"
#include <vector>
#include <set>
#include <map>
#include <string>
class ArgList;
class CpptrajFile;
class DataFile;
class DataFileList;
class DataSet;
class DataSet_2D;
class DataSet_integer;
class DataSetList;
class Topology;
namespace Cpptraj {
namespace HB {
/// Hold hydrogen bond calculation data
class HbData {
  public:
    /// CONSTRUCTOR
    HbData();
    /// Process data-related args
    int ProcessArgs(ArgList&, DataFileList&);
    /// Initialize hydrogen bond data 
    int InitHbData(DataSetList*, std::string const&);
    /// Set pointer to current Topology
    void SetCurrentParm(Topology const*);
    /// Output HBond data
    void PrintHbData();

    /// Add a solute-solute hydrogen bond
    void AddUU(double, double, int, int, int, int, int);
    /// Finish calc for a Frame and increment total # frames
    void IncrementNframes();
  private:
    enum MatrixNormType { NORM_NONE = 0, NORM_FRAMES, NORM_RESMAX };

    typedef std::vector<int> Iarray;
    typedef std::pair<int,int> Hpair;
    typedef std::map<Hpair,Hbond> UUmapType;
    typedef std::map<int,Hbond> UVmapType;
//    typedef std::map< int,std::set<int> > RmapType;
    typedef std::map< std::set<int>,Bridge > BmapType;
    typedef std::vector<Hbond> Harray;
    typedef std::map<int,int> IdxMapType;

    class bridgeSorter;

    /// \return legend for hydrogen bond series set
    static inline std::string CreateHBlegend(Topology const&, int, int, int);
    /// \return DataSet Index for solute-solute hbond
    int UU_Set_Idx(int, int) const;
    /// \return Solute-solute hydrogen bond time series with legend set
    DataSet_integer* UUset(int, int, int);
    /// \return String containing estimated memory usage
    std::string MemoryUsage(size_t, size_t, size_t) const;
    /// Finish hbond time series
    void FinishSeries(DataSet_integer*, unsigned int);
    /// Ensure all hbond time series have same # frames
    void UpdateSeries();
    /// Write summary by parts file header
    static void summary_Parts_header(CpptrajFile*, unsigned int);
    /// Print parts summary for given hbond
    void summary_Parts(CpptrajFile*, Hbond const&) const;
    // ---------------------------------
    DataSetList* masterDSL_;  ///< Pointer to the master DataSetList
    Topology const* CurrentParm_; ///< Pointer to current Topology

    UUmapType UU_Map_;        ///< Map solute donorH/acceptor pair to UU hbond
    UVmapType UV_Map_;        ///< Map solute donorH or solute acceptor to UV hbond
//    RmapType solvent2solute_; ///< Map solvent res # to solute residues it is bound to each frame
    BmapType BridgeMap_; ///< Map residues involved in bridging to # frames bridge present
    IdxMapType DidxMap_; ///< Map solute hydrogen donor atom # to index (series only)
    IdxMapType AidxMap_; ///< Map solute acceptor atom # to index (series only).

    DataSet* NumHbonds_;     ///< Hold # UU hbonds per frame.
    DataSet* NumSolvent_;    ///< Hold # UV hbonds per frame.
    DataSet* NumBridge_;     ///< Hold # solute-solvent bridges per frame.
    DataSet* BridgeID_;      ///< Hold info on each bridge per frame.
    DataSet_2D* UU_matrix_byRes_; ///< Record # hbonds between each residue pair.
    DataFile* nhbout_;       ///< File to write # hbonds vs time to
    DataFile* UUseriesout_;  ///< File to write UU time series to.
    DataFile* UVseriesout_;  ///< File to write UV time series to.
    DataFile* Bseriesout_;   ///< File to write bridge time series to.
    DataFile* uuResMatrixFile_; ///< UU hbond matrix file
    CpptrajFile* avgout_;    ///< File to write UU averages to.
    CpptrajFile* solvout_;   ///< File to write UV averages to.
    CpptrajFile* bridgeout_; ///< File to write bridge totals to.

    Iarray splitFrames_;         ///< For calculating hydrogen bonds by parts
    std::string hbsetname_;   ///< Hydrogen bond data set name
    unsigned int Nframes_;    ///< Total # of frames going into the calculation
    MatrixNormType UUmatByRes_norm_;
    int nuuhb_;               ///< Number of UU hydrogen bonds for the current frame.
    int nuvhb_;               ///< Number of UV hydrogen bonds for the current frame.
    int nbridge_;             ///< Number of UV bridges for the current frame.
    bool series_;             ///< If true, track hbond time series TODO
    bool Bseries_;            ///< If true, track bridge time series
    bool calcSolvent_;        ///< If true track solvent hbonds TODO
    bool seriesUpdated_;      ///< If false time series need to be updated to have same # frames
    bool useAtomNum_;         ///< If true include atom numbers in labels/legends
    bool bridgeByAtom_;       ///< If true determine bridging by atom.
    bool do_uuResMatrix_;     ///< If true calculate UU matrix
};

}
}
#endif
