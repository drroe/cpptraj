#include "HbData.h"
#include "Hbond.h"
#include "../DataFile.h"
#include "../DataSetList.h"
#include "../DataSet_integer.h"
#include "../DataSet_2D.h"

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
HbData::HbData() :
  masterDSL_(0),
  CurrentParm_(0),
  UUseriesout_(0),
  UU_matrix_byRes_(0),
  series_(false)
{}

/** Set pointer to the master DataSetList, set HB data set name. */
void HbData::InitHbData(DataSetList* dslPtr, std::string const& setNameIn, DataFile* UUseriesoutIn) {
  masterDSL_ = dslPtr;
  if (setNameIn.empty())
    hbsetname_ = masterDSL_->GenerateDefaultName("HB");
  else
    hbsetname_ = setNameIn;
  if (series_) {
    UUseriesout_ = UUseriesoutIn;
    masterDSL_->SetDataSetsPending(true);
  }
}

/** Set pointer to current Topology */
void HbData::SetCurrentParm( Topology const* topPtr ) {
  CurrentParm_ = topPtr;
}

/** Create legend for hydrogen bond based on given atoms. */
std::string HbData::CreateHBlegend(Topology const& topIn, int a_atom, int h_atom, int d_atom)
{
  if (a_atom == -1)
    return (topIn.TruncResAtomName(h_atom) + "-V");
  else if (d_atom == -1)
    return (topIn.TruncResAtomName(a_atom) + "-V");
  else
    return (topIn.TruncResAtomName(a_atom) + "-" +
            topIn.TruncResAtomName(d_atom) + "-" +
            topIn[h_atom].Name().Truncated());
}

/** Determine solute-solute hbond index for backwards compatibility:
  *   hbidx = (donorIndex * #acceptors) + acceptorIndex
  */
int HbData::UU_Set_Idx(int a_atom, int h_atom) const {
  IdxMapType::const_iterator it = DidxMap_.find( h_atom );
  int didx = it->second;
  it = AidxMap_.find( a_atom );
  int aidx = it->second;
  return (didx * AidxMap_.size()) + aidx;
}

/** \return solute-solute hydrogen bond time series set with legend set. */
DataSet_integer* HbData::UUset(int a_atom, int h_atom, int d_atom) {
  DataSet_integer* ds = (DataSet_integer*)
    masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,"solutehb",UU_Set_Idx(a_atom,h_atom)));
  if (UUseriesout_ != 0) UUseriesout_->AddDataSet( ds );
  ds->SetLegend( CreateHBlegend(*CurrentParm_, a_atom, h_atom, d_atom) );
  return ds;
}

/** Add or update a solute-solute hydrogen bond with given angle/distance. */
void HbData::AddUU(double dist, double angle, int fnum, int a_atom, int h_atom, int d_atom, int onum)
{
  // Index UU hydrogen bonds by DonorH-Acceptor
  Hpair hbidx(h_atom, a_atom);
  UUmapType::iterator it = UU_Map_.lower_bound( hbidx );
  if (it == UU_Map_.end() || it->first != hbidx)
  {
//      mprintf("DBG1: NEW hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
    DataSet_integer* ds = 0;
    if (series_) {
      ds = UUset(a_atom, h_atom, d_atom);
    }
    it = UU_Map_.insert(it, std::pair<Hpair,Hbond>(hbidx,Hbond(ds, a_atom, h_atom, d_atom, splitFrames_)));
  } else {
//      mprintf("DBG1: OLD hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
  }
  it->second.Update(dist, angle, fnum, splitFrames_, onum);
  if (UU_matrix_byRes_ != 0) {
    int a_res = (*CurrentParm_)[a_atom].ResNum();
    int d_res = (*CurrentParm_)[d_atom].ResNum();
    UU_matrix_byRes_->UpdateElement(a_res, d_res, 1.0);
  }
}

