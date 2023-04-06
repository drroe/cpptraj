#include "MultiFlexResults.h"
#include "DataSetList.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;
// MultiFlexResults class
MultiFlexResults::MultiFlexResults() :
  ssi_matrix_(0),
  pkInt_(0),
  delta_pK_self_(0),
  delta_pK_back_(0),
  siteNames_(0)
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

