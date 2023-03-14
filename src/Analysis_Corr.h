#ifndef INC_ANALYSIS_CORR_H
#define INC_ANALYSIS_CORR_H
#include "Analysis.h"
class DataSet_1D;
// Class: Analysis_Corr
/// Calculate autocorrelation or correlation of dataset(s).
class Analysis_Corr : public Analysis {
  public:
    Analysis_Corr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Corr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int cosine_dihedral_autocorr(DataSet_1D const&);

    DataSet *D1_;
    DataSet *D2_;
    DataSet* Ct_;
    DataSet* Coeff_;
    int lagmax_;
    bool usefft_;
    bool calc_covar_;
    bool cosphi_;
};
#endif
