#ifndef INC_HB_HBDATA_H
#define INC_HB_HBDATA_H
#include <vector>
#include <map>
#include <string>
class DataFile;
class DataSet_2D;
class DataSet_integer;
class DataSetList;
class Topology;
namespace Cpptraj {
namespace HB {
class Hbond;
/// Hold hydrogen bond calculation data
class HbData {
  public:
    /// CONSTRUCTOR
    HbData();
    /// Initialize hydrogen bond data 
    void InitHbData(DataSetList*, std::string const&, DataFile*);
    /// Set pointer to current Topology
    void SetCurrentParm(Topology const*);

    /// Add a solute-solute hydrogen bond
    void AddUU(double, double, int, int, int, int, int);
  private:
    typedef std::vector<int> Iarray;
    typedef std::pair<int,int> Hpair;
    typedef std::map<Hpair,Hbond> UUmapType;
    typedef std::map<int,int> IdxMapType;

    /// \return legend for hydrogen bond series set
    static inline std::string CreateHBlegend(Topology const&, int, int, int);
    /// \return DataSet Index for solute-solute hbond
    int UU_Set_Idx(int, int) const;
    /// \return Solute-solute hydrogen bond time series with legend set
    DataSet_integer* UUset(int, int, int);

    // ---------------------------------
    DataSetList* masterDSL_;  ///< Pointer to the master DataSetList
    Topology const* CurrentParm_; ///< Pointer to current Topology

    UUmapType UU_Map_;        ///< Map solute donorH/acceptor pair to UU hbond
    IdxMapType DidxMap_; ///< Map solute hydrogen donor atom # to index (series only)
    IdxMapType AidxMap_; ///< Map solute acceptor atom # to index (series only).

    DataFile* UUseriesout_;  ///< File to write UU time series to.
    DataSet_2D* UU_matrix_byRes_; ///< Record # hbonds between each residue pair.

    Iarray splitFrames_;         ///< For calculating hydrogen bonds by parts
    std::string hbsetname_;   ///< Hydrogen bond data set name
    bool series_;             ///< If true, track hbond time series TODO
};

}
}
#endif
