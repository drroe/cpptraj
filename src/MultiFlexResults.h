#ifndef INC_MULTIFLEXRESULTS_H
#define INC_MULTIFLEXRESULTS_H
#include <string>
class DataSetList;
class DataSet;
namespace Cpptraj {
/// Hold results from MultiFlex calc within MeadInterface class
class MultiFlexResults {
  public:
    MultiFlexResults();
    int Allocate(DataSetList&, std::string const&);
    void AddSiteResult(int, std::string const&, int, double, double, double) const;
  private:
    DataSet* ssi_matrix_;    ///< Site-site interaction matrix
    DataSet* pkInt_;         ///< Intrinsic pKa for each site
    DataSet* delta_pK_self_; ///< Self contribution to the intrinsic pKa for each site
    DataSet* delta_pK_back_; ///< Background contribution to the intrinsic pKa for each site
    DataSet* siteNames_;     ///< Hold name of each site.
};
}
#endif
