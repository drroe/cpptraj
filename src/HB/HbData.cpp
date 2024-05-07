#include "HbData.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataFile.h"
#include "../DataFileList.h"
#include "../DataSetList.h"
#include "../DataSet_integer.h"
#include "../DataSet_2D.h"
#include "../StringRoutines.h"
#include <algorithm> //sort

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
HbData::HbData() :
  masterDSL_(0),
  CurrentParm_(0),
  UU_matrix_byRes_(0),
  UUseriesout_(0),
  UVseriesout_(0),
  Bseriesout_(0),
  avgout_(0),
  solvout_(0),
  bridgeout_(0),
  Nframes_(0),
  UUmatByRes_norm_(NORM_FRAMES),
  series_(false),
  Bseries_(false),
  calcSolvent_(false),
  seriesUpdated_(false),
  useAtomNum_(false),
  bridgeByAtom_(false)
{}

/** Process data-related args */
int HbData::ProcessArgs(ArgList& actionArgs, DataFileList& DFL) {
  series_ = actionArgs.hasKey("series");
  if (series_) {
    UUseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("uuseries"), actionArgs);
    UVseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("uvseries"), actionArgs);
  }
  Bseries_ = actionArgs.hasKey("bseries");
  if (Bseries_) {
    Bseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("bseriesfile"), actionArgs);
  }
  bool do_uuResMatrix = actionArgs.hasKey("uuresmatrix");
  DataFile* uuResMatrixFile = 0; // FIXME class var ?
  if (do_uuResMatrix) {
    uuResMatrixFile = DFL.AddDataFile(actionArgs.GetStringKey("uuresmatrixout"),
                                             ArgList("nosquare2d"),
                                             actionArgs);
    std::string uuResMatrixNorm = actionArgs.GetStringKey("uuresmatrixnorm");
    if (!uuResMatrixNorm.empty()) {
      if (uuResMatrixNorm == "none")
        UUmatByRes_norm_ = NORM_NONE;
      else if (uuResMatrixNorm == "frames")
        UUmatByRes_norm_ = NORM_FRAMES;
      //else if (uuResMatrixNorm == "resmax") // DISABLE for now, needs more testing
      //  UUmatByRes_norm_ = NORM_RESMAX;
      else {
        mprinterr("Error: Invalid keyword for 'uuresmatrixnorm'.\n");
        return 1;
      }
    }
    
  }
  std::string avgname = actionArgs.GetStringKey("avgout");
  std::string solvname = actionArgs.GetStringKey("solvout");
  if (solvname.empty()) solvname = avgname;
  std::string bridgename = actionArgs.GetStringKey("bridgeout");
  if (bridgename.empty()) bridgename = solvname;
  avgout_ = DFL.AddCpptrajFile(avgname, "Avg. solute-solute HBonds");
  if (calcSolvent_) {
    solvout_ = DFL.AddCpptrajFile(solvname,"Avg. solute-solvent HBonds");
    bridgeout_ = DFL.AddCpptrajFile(bridgename,"Solvent bridging info");
  }

  useAtomNum_ = actionArgs.hasKey("printatomnum");
  bridgeByAtom_ = actionArgs.hasKey("bridgebyatom");

  // Hbond split analysis
  std::string splitarg = actionArgs.GetStringKey("splitframe");
  if (!splitarg.empty()) {
    ArgList splits( splitarg, "," );
    if (splits.Nargs() < 1) {
      mprinterr("Error: Invalid argument for 'splitframe': %s\n", splitarg.c_str());
      return 1;
    }
    int sf = splits.getNextInteger(-1); // User frame #s start at 1
    while (sf > 0) {
      // Since user frame #s start at 1, subtract 1 for actual frame index.
      sf--;
      // Check that frame arg is valid
      if (!splitFrames_.empty()) {
        if (sf <= splitFrames_.back()) {
          mprinterr("Error: 'splitframe' #s must be in increasing order.\n");
          return 1;
        }
      }
      splitFrames_.push_back( sf );
      sf = splits.getNextInteger(-1);
    }
    if ((int)splitFrames_.size() < splits.Nargs()) {
      mprinterr("Error: Invalid split frame arguments.\n");
      splits.CheckForMoreArgs();
      return 1;
    }
  }

  return 0;
}

/** Set pointer to the master DataSetList, set HB data set name. */
void HbData::InitHbData(DataSetList* dslPtr, std::string const& setNameIn) {
  masterDSL_ = dslPtr;
  if (setNameIn.empty())
    hbsetname_ = masterDSL_->GenerateDefaultName("HB");
  else
    hbsetname_ = setNameIn;
  if (series_ || Bseries_) {
    masterDSL_->SetDataSetsPending(true);
  }
  seriesUpdated_ = false;
  Nframes_ = 0;
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

/** Estimate the memory usage of the hbond command. */
std::string HbData::MemoryUsage(size_t n_uu_pairs, size_t n_uv_pairs, size_t nFrames) const
{
  static const size_t sizeHbond = sizeof(Hbond);
  // NOTE: Assuming an overhead of 32 bytes per map element.
  static const size_t sizeElt = 32;
  static const size_t sizeUUmapElt = sizeElt + sizeof(Hpair) + sizeHbond;
  static const size_t sizeUVmapElt = sizeElt + sizeof(int) + sizeHbond;
  static const size_t sizeBRmapElt = sizeElt + sizeof(std::set<int>) + sizeof(Bridge);
  // Solute-solute pairs
  size_t memTotal = sizeof(UUmapType) + (n_uu_pairs * sizeUUmapElt);
  // Solute-solvent pairs
  memTotal += sizeof(UVmapType) + (n_uv_pairs * sizeUVmapElt);
  // Time series TODO bridge series
  if (series_ && nFrames > 0) {
    size_t seriesSet = (nFrames * sizeof(int)) + sizeof(DataSet_integer);
    memTotal += (seriesSet * (n_uu_pairs + n_uv_pairs));
  }
  // Solute-solvent bridges
  // Cannot really estimate bridging, so always base it on BridgeMap_
  memTotal += sizeof(BmapType);
  for (BmapType::const_iterator it = BridgeMap_.begin(); it != BridgeMap_.end(); ++it)
    memTotal += (sizeBRmapElt + it->first.size()*sizeof(int));
  // Matrices
  if (UU_matrix_byRes_ != 0)
    memTotal += UU_matrix_byRes_->MemUsageInBytes();
 
  return ByteString( memTotal, BYTE_DECIMAL );
}

/** Ensure DataSet time series has at least N frames, fill out if necessary. */
void HbData::FinishSeries(DataSet_integer* data, unsigned int N) {
  static const int ZERO = 0;
  if (data != 0 && N > 0) {
    if ( data->Size() < N ) {
      data->Add( N-1, &ZERO );
#     ifdef MPI
      data->SetNeedsSync( false );
#     endif
    }
  }
}

/** Ensure all time series data is up-to-date with Nframes.
  * Should only be called once.
  */
void HbData::UpdateSeries() {
  if (seriesUpdated_) return;
  if (series_ && Nframes_ > 0) {
    for (UUmapType::iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb)
      FinishSeries(hb->second.Data(), Nframes_);
    for (UVmapType::iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb)
      FinishSeries(hb->second.Data(), Nframes_);
  }
  if (Bseries_ && Nframes_ > 0) {
    for (BmapType::iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
      FinishSeries( b->second.Data(), Nframes_ );
  }
  seriesUpdated_ = true;
}

/** Print header for summary by parts. */
void HbData::summary_Parts_header(CpptrajFile* avgout, unsigned int nParts)
{
  if (nParts < 1) return;
  for (unsigned int idx = 0; idx != nParts; idx++) {
    std::string spart(integerToString(idx+1));
    std::string frames( "Frames"  + spart);
    std::string frac(   "Frac"    + spart);
    std::string avgdist("AvgDist" + spart);
    std::string avgang( "AvgAng"  + spart);
    avgout->Printf(" %8s %12s %12s %12s", frames.c_str(), frac.c_str(), avgdist.c_str(), avgang.c_str());
  }
}

/** Print summary by parts for given hbond. */
void HbData::summary_Parts(CpptrajFile* avgout, Hbond const& hb) const {
  for (unsigned int idx = 0; idx != hb.Nparts(); idx++)
    avgout->Printf(" %8i %12.4f %12.4f %12.4f",
                   hb.PartFrames(idx), hb.PartFrac(idx, Nframes_),
                   hb.PartDist(idx).mean(), hb.PartAngle(idx).mean()*Constants::RADDEG);
}

/** Used to associate Bridge with solute atoms/residues and sort Bridges by frames. */
class HbData::bridgeSorter {
  public:
    /// CONSTRUCTOR - list of solute atoms/residues, Bridge info
    bridgeSorter(std::set<int> const& uIdx, Bridge const& bridge) :
      uIdx_(uIdx), bridge_(bridge) {}
    /// Used to sort by bridge # frames
    bool operator<(bridgeSorter const& rhs) const {
      if (bridge_.Frames() == rhs.bridge_.Frames())
        return (uIdx_ < rhs.uIdx_);
      else
        return (bridge_.Frames() > rhs.bridge_.Frames());
    }
    /// \return List of solute atom/residue #s
    std::set<int> const& Uidx() const { return uIdx_; }
    /// \return Bridging info
    Bridge const& Binfo()       const { return bridge_; }
  private:
    std::set<int> uIdx_; ///< Hold solute atom/residue #s
    Bridge bridge_;      ///< Hold bridging info.
};

/** Print average occupancies over all frames for all detected Hbonds. */
void HbData::PrintHbData() {
  Harray HbondList; // For sorting
  std::string Aname, Hname, Dname;

  // Final memory usage
  mprintf("    HBOND: Actual memory usage is %s\n",
          MemoryUsage(UU_Map_.size(), UV_Map_.size(), Nframes_).c_str());
  mprintf("\t%zu solute-solute hydrogen bonds.\n", UU_Map_.size());
  if (calcSolvent_) {
   mprintf("\t%zu solute-solvent hydrogen bonds.\n", UV_Map_.size());
   mprintf("\t%zu unique solute-solvent bridging interactions.\n", BridgeMap_.size());
  }
# ifdef TIMER
  t_uu_.WriteTiming(      2,"Solute-Solute   :",t_action_.Total());
  if (calcSolvent_) {
    t_uv_.WriteTiming(    2,"Solute-Solvent  :",t_uv_.Total());
    t_bridge_.WriteTiming(2,"Bridging waters :",t_action_.Total());
  }
  t_action_.WriteTiming(1,"Total:");
# endif
  // Ensure all series have been updated for all frames.
  UpdateSeries();
  // Matrix normalization
  if (UU_matrix_byRes_ != 0) {
    if (UUmatByRes_norm_ == NORM_FRAMES) {
      double norm = 1.0 / ((double)Nframes_);
      for (unsigned int r = 0; r != UU_matrix_byRes_->Nrows(); r++)
        for (unsigned int c = 0; c != UU_matrix_byRes_->Ncols(); c++)
          UU_matrix_byRes_->SetElement(c, r, UU_matrix_byRes_->GetElement(c, r) * norm);
    }
  }

  if (CurrentParm_ == 0) return;
  // Calculate necessary column width for strings based on how many residues.
  // ResName+'_'+ResNum+'@'+AtomName | NUM = 4+1+R+1+4 = R+10
  int NUM = DigitWidth( CurrentParm_->Nres() ) + 10;
  // If useAtomNum_ +'_'+AtomNum += 1+A
  if (useAtomNum_) NUM += ( DigitWidth( CurrentParm_->Natom() ) + 1 );

  // Solute Hbonds 
  if (avgout_ != 0) { 
    // Place all detected Hbonds in a list and sort.
    for (UUmapType::const_iterator it = UU_Map_.begin(); it != UU_Map_.end(); ++it) {
      HbondList.push_back( it->second );
      // Calculate average distance and angle for this hbond.
      HbondList.back().CalcAvg();
    }
    UU_Map_.clear();
    // Sort and Print 
    sort( HbondList.begin(), HbondList.end() );
    avgout_->Printf("%-*s %*s %*s %8s %12s %12s %12s", NUM, "#Acceptor", 
                    NUM, "DonorH", NUM, "Donor", "Frames", "Frac", "AvgDist", "AvgAng");
    if (!splitFrames_.empty())
      summary_Parts_header(avgout_, splitFrames_.size()+1);
    avgout_->Printf("\n");
    for (Harray::const_iterator hbond = HbondList.begin(); hbond != HbondList.end(); ++hbond ) 
    {
      double avg = ((double)hbond->Frames()) / ((double) Nframes_);
      Aname = CurrentParm_->TruncResAtomName(hbond->A());
      Hname = CurrentParm_->TruncResAtomName(hbond->H());
      Dname = CurrentParm_->TruncResAtomName(hbond->D());
      if (useAtomNum_) {
        Aname.append("_" + integerToString(hbond->A()+1));
        Hname.append("_" + integerToString(hbond->H()+1));
        Dname.append("_" + integerToString(hbond->D()+1));
      }
      avgout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     hbond->Frames(), avg, hbond->Dist(), hbond->Angle());
      if (!splitFrames_.empty())
        summary_Parts(avgout_, *hbond);
      avgout_->Printf("\n");
    }
  }

  // Solute-solvent Hbonds 
  if (solvout_ != 0 && calcSolvent_) {
    HbondList.clear();
    for (UVmapType::const_iterator it = UV_Map_.begin(); it != UV_Map_.end(); ++it) {
      HbondList.push_back( it->second );
      // Calculate average distance and angle for this hbond.
      HbondList.back().CalcAvg();
    }
    UV_Map_.clear();
    sort( HbondList.begin(), HbondList.end() );
    // Calc averages and print
    solvout_->Printf("#Solute-Solvent Hbonds:\n");
    solvout_->Printf("%-*s %*s %*s %8s %12s %12s %12s", NUM, "#Acceptor", 
                     NUM, "DonorH", NUM, "Donor", "Count", "Frac", "AvgDist", "AvgAng");
    if (!splitFrames_.empty())
      summary_Parts_header(solvout_, splitFrames_.size()+1);
    solvout_->Printf("\n");
    for (Harray::const_iterator hbond = HbondList.begin(); hbond != HbondList.end(); ++hbond )
    {
      // Average has slightly diff meaning since for any given frame multiple
      // solvent can bond to the same solute.
      double avg = ((double)hbond->Frames()) / ((double) Nframes_);
      if (hbond->A()==-1) // Solvent acceptor
        Aname = "SolventAcc";
      else {
        Aname = CurrentParm_->TruncResAtomName(hbond->A());
        if (useAtomNum_) Aname.append("_" + integerToString(hbond->A()+1));
      }
      if (hbond->D()==-1) { // Solvent donor
        Dname = "SolventDnr";
        Hname = "SolventH";
      } else {
        Dname = CurrentParm_->TruncResAtomName(hbond->D());
        Hname = CurrentParm_->TruncResAtomName(hbond->H());
        if (useAtomNum_) {
          Dname.append("_" + integerToString(hbond->D()+1));
          Hname.append("_" + integerToString(hbond->H()+1));
        }
      }
      solvout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     hbond->Frames(), avg, hbond->Dist(), hbond->Angle());
      if (!splitFrames_.empty())
        summary_Parts(solvout_, *hbond);
      solvout_->Printf("\n");
    }
    HbondList.clear();
  }

  // BRIDGING INFO
  if (bridgeout_ != 0 && calcSolvent_) {
    if (bridgeByAtom_)
      bridgeout_->Printf("#Bridging Solute Atoms:\n");
    else
      bridgeout_->Printf("#Bridging Solute Residues:\n");
    // Place bridging values in a vector for sorting
    typedef std::vector<bridgeSorter> Bvec;
    Bvec bridgevector;
    bridgevector.reserve( BridgeMap_.size() );
    for (BmapType::const_iterator it = BridgeMap_.begin();
                                  it != BridgeMap_.end(); ++it)
      bridgevector.push_back( bridgeSorter(it->first, it->second) );
    std::sort( bridgevector.begin(), bridgevector.end() );
    for (Bvec::const_iterator bv = bridgevector.begin(); bv != bridgevector.end(); ++bv)
    {
      if (bridgeByAtom_) {
        bridgeout_->Printf("Bridge Atm");
        for (std::set<int>::const_iterator atm = bv->Uidx().begin();
                                           atm != bv->Uidx().end(); ++atm)
          bridgeout_->Printf(" %s", CurrentParm_->TruncAtomNameNum(*atm).c_str());
      } else {
        bridgeout_->Printf("Bridge Res");
        for (std::set<int>::const_iterator res = bv->Uidx().begin();
                                           res != bv->Uidx().end(); ++res)
          bridgeout_->Printf(" %i:%s", *res+1, CurrentParm_->Res( *res ).Name().Formatted(4).c_str());
      }
      bridgeout_->Printf(", %i frames.", bv->Binfo().Frames());
      if (!splitFrames_.empty()) {
        bridgeout_->Printf(" Parts:");
        for (unsigned int idx = 0; idx != bv->Binfo().Nparts(); idx++)
          bridgeout_->Printf(" %i", bv->Binfo().PartFrames(idx));
      }
      bridgeout_->Printf("\n");
    } 
  }
}

