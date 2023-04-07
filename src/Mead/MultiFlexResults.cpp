#include "MultiFlexResults.h"
#include "../DataSetList.h"
#include "../DataFileList.h"
#include "../StringRoutines.h"
#include "../CpptrajStdio.h"
#include "../DataSet_2D.h"
#include "../DataFile.h"

using namespace Cpptraj::Mead;
// MultiFlexResults class
MultiFlexResults::MultiFlexResults() :
  ssi_matrix_(0),
  pkInt_(0),
  delta_pK_self_(0),
  delta_pK_back_(0),
  siteNames_(0),
  pkintfile_(0),
  summfile_(0),
  gfile_(0)
{}

static inline int ds_alloc_err(const char* desc) {
  mprinterr("Error: Could not allocate %s\n");
  return 1;
}

int MultiFlexResults::Allocate(DataSetList& dsl, std::string const& dsname) {
  if (ssi_matrix_ != 0) {
    mprinterr("Internal Error: MultiFlexResults::Allocate: Already allocated.\n");
    return 1;
  }
  ssi_matrix_ = dsl.AddSet( DataSet::MATRIX_DBL, MetaData(dsname, "ssi") );
  if (ssi_matrix_ == 0) return ds_alloc_err("site-site interaction matrix");
  ssi_matrix_->SetupFormat().SetFormatType( TextFormat::SCIENTIFIC );
  pkInt_ = dsl.AddSet( DataSet::DOUBLE, MetaData(dsname, "pkint") );
  if (pkInt_ == 0) return ds_alloc_err("intrinsic pKa array");
  delta_pK_self_ = dsl.AddSet( DataSet::DOUBLE, MetaData(dsname, "dpkself") );
  if (delta_pK_self_ == 0) return ds_alloc_err("intrinsic pKa self contribution array");
  delta_pK_back_ = dsl.AddSet( DataSet::DOUBLE, MetaData(dsname, "dpkback") );
  if (delta_pK_back_ == 0) return ds_alloc_err("intrinsic pKa background contribution array");
  siteNames_ = dsl.AddSet( DataSet::STRING, MetaData(dsname, "sitenames" ) );
  if (siteNames_ == 0) return ds_alloc_err("site names array");

  return 0;
}

int MultiFlexResults::CreateOutputFiles(DataFileList& dfl, std::string const& pkintname, std::string const& summname, std::string const& gname)
{
  pkintfile_ = dfl.AddCpptrajFile(pkintname, "intrinsic pKa values", DataFileList::TEXT, true);
  summfile_ = dfl.AddCpptrajFile(summname, "intrinsic pKa summary", DataFileList::TEXT, true);
  gfile_ = dfl.AddCpptrajFile(gname, "site-site interactions", DataFileList::TEXT, true);
  if (pkintfile_ == 0 || summfile_ == 0 || gfile_ == 0) {
    mprinterr("Error: Could not create multiflex output files.\n");
    return 1;
  }
  return 0;
}

void MultiFlexResults::AddSetsToFile(DataFile* outfile, DataFile* ssiout)
const
{
  if (outfile != 0) {
    outfile->AddDataSet( pkInt_ );
    outfile->AddDataSet( delta_pK_self_ );
    outfile->AddDataSet( delta_pK_back_ );
    outfile->AddDataSet( siteNames_ );
  }
  if (ssiout != 0) {
    ssiout->AddDataSet(ssi_matrix_);
  }
}

void MultiFlexResults::AllocateSets(unsigned int sizeIn)
const
{
  DataSet::SizeArray dsize(1, sizeIn);
  pkInt_->Allocate( dsize );
  delta_pK_self_->Allocate( dsize );
  delta_pK_back_->Allocate( dsize );
  siteNames_->Allocate( dsize );

  // Allocate matrix
  DataSet_2D& mat = static_cast<DataSet_2D&>( *ssi_matrix_ );
  mat.AllocateTriangle( sizeIn );
}

void MultiFlexResults::AddSiteResult(int idx,
                                     std::string const& sitenameIn,
                                     int originalResNumIn,
                                     double pkInt_in,
                                     double delta_pK_self_in,
                                     double delta_pK_back_in)
const
{
  std::string sitename = sitenameIn + "-" + integerToString(originalResNumIn);
  siteNames_->Add(idx, sitename.c_str());
  pkInt_->Add(idx, &pkInt_in);
  delta_pK_self_->Add(idx, &delta_pK_self_in);
  delta_pK_back_->Add(idx, &delta_pK_back_in);
} 

void MultiFlexResults::AddSiteSiteMatrix(std::vector< std::vector<double> > const& ssi)
const
{
  DataSet_2D& mat = static_cast<DataSet_2D&>( *ssi_matrix_ );

  for (unsigned int i = 0; i < ssi.size(); i++) {
    for (unsigned int j = i + 1; j < ssi.size(); j++) {
      mat.SetElement(i, j, ssi[i][j]);
    }
  }
}
