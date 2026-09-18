#ifndef INC_ANALYSIS_MDANCE_H
#define INC_ANALYSIS_MDANCE_H
#include "Analysis.h"
#include "ExtendedSimilarity.h" // Has the same metrics
#include "TrajectoryFile.h"
#ifdef HAS_EIGEN
#include "Mdance/types.h"
#endif 
class DataSet_Coords;
/// Interface to MDANCE
/** MDANCE: Molecular Dynamics Analysis with N-aray Clustering Ensembles.
  */ 
class Analysis_MDANCE : public Analysis {
  public:
    Analysis_MDANCE();
    ~Analysis_MDANCE();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_MDANCE(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static const char* kinitKeys_[];
    static const char* kinitStr_[];
    static const Cpptraj::Mdance::MD::KinitType kinitTypes_[];

    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> ClusterArray;

    void getClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                            TrajectoryFile::TrajFormatType&) const;
    void writeClusterTraj(ClusterArray const&) const;
    void writeCenterTraj() const;


    int debug_;
    DataSet_Coords* coords_;
    int kClusters_;
    int percentage_;
    int vthresh_;
    ExtendedSimilarity::MetricType metric_;
    Cpptraj::Mdance::MD::KinitType kinit_;
    AtomMask mask_;

    DataSet* cnumvtime_;                        ///< Hold cluster number for each frame
    DataSet_Coords* centers_;                   ///< Hold cluster centers
    std::string clusterfile_;                   ///< Cluster trajectory base filename.
    std::string centerfile_;
    TrajectoryFile::TrajFormatType clusterfmt_; ///< Cluster trajectory format.
    TrajectoryFile::TrajFormatType centerfmt_;
};
#endif
