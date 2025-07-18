#ifndef INC_HB_HBDATA_H
#define INC_HB_HBDATA_H
#include "Bridge.h"
#include "Hbond.h"
#include "../Parallel.h" // FIXME Comm should be in a namespace
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
class DataSet_MatrixDbl;
class DataSetList;
class Topology;
namespace Cpptraj {
namespace HB {
/// Hold hydrogen bond calculation data
class HbData {
    typedef std::vector<int> Iarray;
  public:
    /// CONSTRUCTOR
    HbData();
    /// DESTRUCTOR
    ~HbData();
    /// Set debug level
    void SetDebug(int);
    /// Process data-related args
    int ProcessArgs(ArgList&, DataFileList&, bool);
    /// Initialize hydrogen bond data 
    int InitHbData(DataSetList*, std::string const&);
    /// Set pointer to current Topology
    int SetCurrentParm(Topology const*, Iarray const&, Iarray const&, Iarray const&, Iarray const&);
    /// Output HBond data
    void PrintHbData();

    /// Print options to stdout
    void PrintHbDataOpts() const;
    /// \return True if solvent hydrogen bonds are being calculated
    bool CalcSolvent() const { return calcSolvent_; }
    /// \return True if ignoring hbonds between atoms in same molecule
    bool NoIntramol() const { return noIntramol_; }
    /// \return True if saving hydrogen bond time series data
    bool Series() const { return series_; }
    /// \return True if calculating interaction matrix
    bool InteractionMatrix() const { return (UU_matrix_byRes_ != 0); }
    /// \return Debug level
    int Debug() const { return debug_; }
    /// \return Pointer to current topology
    Topology const* CurrentParmPtr() const { return CurrentParm_; }
    /// \return String containing estimated memory usage
    std::string MemoryUsage(size_t, size_t, size_t) const;

    /// Add a solute-solute hydrogen bond
    void AddUU(double, double, int, int, int, int, int);
    /// Add a solute-solvent hydrogen bond, note potential bridges
    void AddUV(double, double, int, int, int, int, bool, int);
    /// Finish calc for a Frame and increment total # frames, do bridge calc if needed
    void IncrementNframes(int, int);
#   ifdef MPI
    /// Set the across traj communicator
    void SetTrajComm(Parallel::Comm const&);
    /// Sync data to the master process
    int SyncToMaster();
#   endif
  private:
    enum MatrixNormType { NORM_NONE = 0, NORM_FRAMES, NORM_RESMAX };

    static const int ID_SOLVENT_;
    static const int ID_ION_;

    typedef std::pair<int,int> Hpair;
    typedef std::map<Hpair,Hbond> UUmapType;
    typedef std::map<int,Hbond> UVmapType;
    typedef std::map< int,std::set<int> > RmapType;
    typedef std::map< std::set<int>,Bridge > BmapType;
    typedef std::vector<Hbond> Harray;
    typedef std::map<int,int> IdxMapType;

    class bridgeSorter;

    /// Set up solute/solute interaction matrix
    int SetupInteractionMatrix(Iarray const&, Iarray const&, Iarray const&);
    /// \return legend for hydrogen bond series set
    static inline std::string CreateHBlegend(Topology const&, int, int, int);
    /// \return DataSet Index for solute-solute hbond
    int UU_Set_Idx(int, int) const;
    /// \return Solute-solute hydrogen bond time series with legend set
    DataSet_integer* UUset(int, int, int);
    /// \return legend for bridge based on indices
    static inline std::string CreateBridgeLegend(std::string const&, std::set<int> const&);
    /// Calculate Bridges for the current frame
    void BridgeCalc(int, int);
    /// Finish hbond time series
    void FinishSeries(DataSet_integer*, unsigned int);
    /// Ensure all hbond time series have same # frames
    void UpdateSeries();
    /// Write summary by parts file header
    static void summary_Parts_header(CpptrajFile*, unsigned int);
    /// Print parts summary for given hbond
    void summary_Parts(CpptrajFile*, Hbond const&) const;
#   ifdef MPI
    /// \return Array containing # hbonds on each rank
    static std::vector<int> GetRankNhbonds(int,Parallel::Comm const&);
    /// Flatten given Hbond into arrays
    static void HbondToArray(std::vector<double>&, std::vector<int>&, Hbond const&);
#   endif
    // ---------------------------------
    DataSetList* masterDSL_;  ///< Pointer to the master DataSetList
    Topology const* CurrentParm_; ///< Pointer to current Topology

    UUmapType UU_Map_;        ///< Map solute donorH/acceptor pair to UU hbond
    UVmapType UV_Map_;        ///< Map solute donorH or solute acceptor to UV hbond
    RmapType solvent2solute_; ///< Map solvent res # to solute residues it is bound to each frame
    BmapType BridgeMap_; ///< Map residues involved in bridging to # frames bridge present
    IdxMapType DidxMap_; ///< Map solute hydrogen donor atom # to index (series only)
    IdxMapType AidxMap_; ///< Map solute acceptor atom # to index (series only).

    DataSet* NumHbonds_;               ///< Hold # UU hbonds per frame.
    DataSet* NumSolvent_;              ///< Hold # UV hbonds per frame.
    DataSet* NumBridge_;               ///< Hold # solute-solvent bridges per frame.
    DataSet* BridgeID_;                ///< Hold info on each bridge per frame.
    DataSet_2D* UU_matrix_byRes_;      ///< Record # hbonds between each residue pair.
    DataSet_MatrixDbl* UU_norm_byRes_; ///< For normalizing the max possible # hbonds by residue
    DataFile* nhbout_;                 ///< File to write # hbonds vs time to
    DataFile* UUseriesout_;            ///< File to write UU time series to.
    DataFile* UVseriesout_;            ///< File to write UV time series to.
    DataFile* Bseriesout_;             ///< File to write bridge time series to.
    DataFile* uuResMatrixFile_;        ///< UU hbond matrix file
    CpptrajFile* avgout_;              ///< File to write UU averages to.
    CpptrajFile* solvout_;             ///< File to write UV averages to.
    CpptrajFile* bridgeout_;           ///< File to write bridge totals to.

    Iarray splitFrames_;         ///< For calculating hydrogen bonds by parts
    std::string hbsetname_;   ///< Hydrogen bond data set name
    unsigned int Nframes_;    ///< Total # of frames going into the calculation
    MatrixNormType UUmatByRes_norm_;
    int debug_;               ///< Debug level
    int nuuhb_;               ///< Number of UU hydrogen bonds for the current frame.
    int nuvhb_;               ///< Number of UV hydrogen bonds for the current frame.
//    int nbridge_;             ///< Number of UV bridges for the current frame.
    bool series_;             ///< If true, track hbond time series TODO
    bool Bseries_;            ///< If true, track bridge time series
    bool calcSolvent_;        ///< If true track solvent hbonds TODO
    bool seriesUpdated_;      ///< If false time series need to be updated to have same # frames
    bool useAtomNum_;         ///< If true include atom numbers in labels/legends
    bool bridgeByAtom_;       ///< If true determine bridging by atom.
    bool do_uuResMatrix_;     ///< If true calculate UU matrix
    bool noIntramol_;         ///< If true ignore intramolecular hydrogen bonds/bridges.
#   ifdef MPI
    Parallel::Comm trajComm_; ///< Across-trajectory communicator
#   endif
};
}
}
#endif
