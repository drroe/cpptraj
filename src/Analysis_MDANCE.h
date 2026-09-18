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
#   ifdef HAS_EIGEN
    static const Cpptraj::Mdance::MD::KinitType kinitTypes_[];
#   endif
    typedef std::vector<int> Iarray;
    /// Hold frames/center for MDANCE clustering
    class MdCluster {
      public:
        MdCluster() {}
        Iarray const& Frames() const { return frames_; }
        Frame const& Ctr() const { return ctr_; }
        unsigned int size() const { return frames_.size(); }
        void push_back(int i) { frames_.push_back( i ); }
        void SetCtr( Frame const& f ) { ctr_ = f; }

        bool operator<(MdCluster const& rhs) const {
          if (frames_.size() > rhs.frames_.size())
            return true; // TODO absolute cluster number
          else
            return false;
        }
      private:
        Iarray frames_;
        Frame ctr_;
    };
    typedef std::vector<MdCluster> ClusterArray;

    void getClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                            TrajectoryFile::TrajFormatType&) const;
    void writeClusterTraj(ClusterArray const&) const;
    void writeCenterTraj(ClusterArray const&) const;
    void writeSummary(CpptrajFile&, ClusterArray const&, unsigned int) const;

    int debug_;
    DataSet_Coords* coords_;
    int kClusters_;
    int percentage_;
    int vthresh_;
    ExtendedSimilarity::MetricType metric_;
#   ifdef HAS_EIGEN
    Cpptraj::Mdance::MD::KinitType kinit_;
#   endif
    AtomMask mask_;

    DataSet* cnumvtime_;                        ///< Hold cluster number for each frame
    DataSet_Coords* centers_;                   ///< Hold cluster centers
    std::string clusterfile_;                   ///< Cluster trajectory base filename.
    std::string centerfile_;
    TrajectoryFile::TrajFormatType clusterfmt_; ///< Cluster trajectory format.
    TrajectoryFile::TrajFormatType centerfmt_;
    CpptrajFile* outfile_;
    bool sort_; ///< If true, sort by cluster population (default)
};
#endif
