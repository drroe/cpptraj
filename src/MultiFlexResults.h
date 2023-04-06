#ifndef INC_MULTIFLEXRESULTS_H
#define INC_MULTIFLEXRESULTS_H
#include <string>
#include <vector>
class DataSetList;
class DataSet;
class DataFile;
class DataFileList;
class CpptrajFile;
namespace Cpptraj {
/// Hold results from MultiFlex calc within MeadInterface class
class MultiFlexResults {
  public:
    MultiFlexResults();
    /// Allocate the data sets
    int Allocate(DataSetList&, std::string const&);
    /// Create text output files
    int CreateOutputFiles(DataFileList&, std::string const&, std::string const&, std::string const&);
    /// Put 1d sets into an output file
    void AddSetsToFile(DataFile*, DataFile*) const;
    /// Allocate space in each set
    void AllocateSets(unsigned int) const;
    /// Add site results
    void AddSiteResult(int, std::string const&, int, double, double, double) const;
    /// Add site-site matrix 
    void AddSiteSiteMatrix(std::vector< std::vector<double> > const&) const;

    CpptrajFile* PkIntFile() const { return pkintfile_; }
    CpptrajFile* SummFile() const { return summfile_; }
    CpptrajFile* Gfile() const { return gfile_; }

    DataSet* PkIntSet() const { return pkInt_; }
    DataSet* Delta_pK_SelfSet() const { return delta_pK_self_; }
    DataSet* Delta_pK_BackSet() const { return delta_pK_back_; }
    DataSet* SiteNamesSet() const { return siteNames_; }
  private:
    DataSet* ssi_matrix_;    ///< Site-site interaction matrix
    DataSet* pkInt_;         ///< Intrinsic pKa for each site
    DataSet* delta_pK_self_; ///< Self contribution to the intrinsic pKa for each site
    DataSet* delta_pK_back_; ///< Background contribution to the intrinsic pKa for each site
    DataSet* siteNames_;     ///< Hold name of each site.
    CpptrajFile* pkintfile_; ///< File to write pkint values to in MEAD format.
    CpptrajFile* summfile_;  ///< File to write self and background contributions to the intrinsic pK in MEAD format.
    CpptrajFile* gfile_;     ///< File to write site-site interactions in units of charge squared per length in MEAD format.
};
}
#endif
