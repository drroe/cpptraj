#ifndef INC_MULTIFLEXRESULTS_H
#define INC_MULTIFLEXRESULTS_H
#include <string>
#include <vector>
class DataSetList;
class DataSet;
class DataFile;
namespace Cpptraj {
/// Hold results from MultiFlex calc within MeadInterface class
class MultiFlexResults {
  public:
    MultiFlexResults();
    /// Allocate the data sets
    int Allocate(DataSetList&, std::string const&);
    /// Put 1d sets into an output file
    void AddSetsToFile(DataFile*, DataFile*) const;
    /// Allocate space in each set
    void AllocateSets(unsigned int) const;
    /// Add site results
    void AddSiteResult(int, std::string const&, int, double, double, double) const;
    /// Add site-site matrix 
    void AddSiteSiteMatrix(std::vector< std::vector<double> > const&) const;
  private:
    DataSet* ssi_matrix_;    ///< Site-site interaction matrix
    DataSet* pkInt_;         ///< Intrinsic pKa for each site
    DataSet* delta_pK_self_; ///< Self contribution to the intrinsic pKa for each site
    DataSet* delta_pK_back_; ///< Background contribution to the intrinsic pKa for each site
    DataSet* siteNames_;     ///< Hold name of each site.
};
}
#endif
