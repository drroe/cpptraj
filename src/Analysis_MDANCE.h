#ifndef INC_ANALYSIS_MDANCE_H
#define INC_ANALYSIS_MDANCE_H
#include "Analysis.h"
#include "ExtendedSimilarity.h" // Has the same metrics 
class DataSet_Coords;
/// Interface to MDANCE
/** MDANCE: Molecular Dynamics Analysis with N-aray Clustering Ensembles.
  */ 
class Analysis_MDANCE : public Analysis {
  public:
    Analysis_MDANCE();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_MDANCE(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int debug_;
    DataSet_Coords* coords_;
    int kClusters_;
    ExtendedSimilarity::MetricType metric_;
};
#endif
