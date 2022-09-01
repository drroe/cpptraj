#ifndef INC_ANALYSIS_TWOPARTICLEDIFFUSION_H
#define INC_ANALYSIS_TWOPARTICLEDIFFUSION_H
#include "Analysis.h"
#include "ImageOption.h"
/// <Enter description of Analysis_TwoParticleDiffusion here>
class Analysis_TwoParticleDiffusion : public Analysis {
  public:
    Analysis_TwoParticleDiffusion();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TwoParticleDiffusion(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    ImageOption imageOpt_;
    DataSet_Coords* coords_;
    AtomMask mask_;     ///< Selected atoms for calculation
    DataSet* outDrr_;   ///< Output Drr matrix
    DataSet* outDtt_;   ///< Output Dtt matrix
    double rmax_;       ///< R dimension max
    double rstep_;      ///< R dimension step
    int maxlag_;        ///< Maximum lag to calculate for
};
#endif
