#include "HbData.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataFile.h"
#include "../DataFileList.h"
#include "../DataSetList.h"
#include "../DataSet_integer.h"
#include "../DataSet_2D.h"
#include "../DataSet_MatrixDbl.h"
#include "../StringRoutines.h"
#include <algorithm> //sort

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
HbData::HbData() :
  masterDSL_(0),
  CurrentParm_(0),
  NumHbonds_(0),
  NumSolvent_(0),
  NumBridge_(0),
  BridgeID_(0),
  UU_matrix_byRes_(0),
  UU_norm_byRes_(0),
  nhbout_(0),
  UUseriesout_(0),
  UVseriesout_(0),
  Bseriesout_(0),
  uuResMatrixFile_(0),
  avgout_(0),
  solvout_(0),
  bridgeout_(0),
  Nframes_(0),
  UUmatByRes_norm_(NORM_FRAMES),
  debug_(0),
  nuuhb_(0),
  nuvhb_(0),
//  nbridge_(0),
  series_(false),
  Bseries_(false),
  calcSolvent_(false),
  seriesUpdated_(false),
  useAtomNum_(false),
  bridgeByAtom_(false),
  do_uuResMatrix_(false),
  noIntramol_(false)
{}

/** DESTRUCTOR */
HbData::~HbData() {
  if (UU_norm_byRes_ != 0) delete UU_norm_byRes_;
}

const int HbData::ID_SOLVENT_ = -1;
const int HbData::ID_ION_ = -2;

/** Set debug level */
void HbData::SetDebug(int debugIn) {
  debug_ = debugIn;
}

/** Process data-related args */
int HbData::ProcessArgs(ArgList& actionArgs, DataFileList& DFL) {
  nhbout_ = DFL.AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  series_ = actionArgs.hasKey("series");
  if (series_) {
    UUseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("uuseries"), actionArgs);
    UVseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("uvseries"), actionArgs);
  }
  Bseries_ = actionArgs.hasKey("bseries");
  if (Bseries_) {
    Bseriesout_ = DFL.AddDataFile(actionArgs.GetStringKey("bseriesfile"), actionArgs);
  }
  do_uuResMatrix_ = actionArgs.hasKey("uuresmatrix");
  uuResMatrixFile_ = 0;
  if (do_uuResMatrix_) {
    uuResMatrixFile_ = DFL.AddDataFile(actionArgs.GetStringKey("uuresmatrixout"),
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

  calcSolvent_ = actionArgs.hasKey("solvent");

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
  noIntramol_ = actionArgs.hasKey("nointramol");
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

/** Print options to stdout. Should only be called after ProcessArgs()
  * and InitHbData().
  */
void HbData::PrintHbDataOpts() const {
  if (calcSolvent_)
    mprintf("\tWill look for solvent-solute hydrogen bonds.\n");
  else
    mprintf("\tOnly looking for solute-solute hydrogen bonds.\n");
  if (noIntramol_)
    mprintf("\tOnly looking for intermolecular hydrogen bonds.\n");
  if (nhbout_ != 0) 
    mprintf("\tWriting # Hbond v time results to %s\n", nhbout_->DataFilename().full());
  if (avgout_ != 0)
    mprintf("\tWriting Hbond avgs to %s\n",avgout_->Filename().full());
  if (!splitFrames_.empty()) {
    mprintf("\tWill split analysis at frames:");
    for (Iarray::const_iterator it = splitFrames_.begin(); it != splitFrames_.end(); ++it)
      mprintf(" %i", *it + 1);
    mprintf("\n");
  }
  if (useAtomNum_)
    mprintf("\tAtom numbers will be written to output.\n");
  if (series_) {
    mprintf("\tTime series data for each hbond will be saved for analysis.\n");
    if (UUseriesout_ != 0) mprintf("\tWriting solute-solute time series to %s\n",
                                   UUseriesout_->DataFilename().full());
    if (UVseriesout_ != 0) mprintf("\tWriting solute-solvent time series to %s\n",
                                   UVseriesout_->DataFilename().full());
  }
  if (Bseries_) {
    mprintf("\tTime series data for each bridge will be saved for analysis.\n");
    if (Bseriesout_ != 0)
      mprintf("\tWriting bridge time series to '%s'\n", Bseriesout_->DataFilename().full());
  }
  if (UU_matrix_byRes_ != 0) {
    mprintf("\tCalculating solute-solute residue matrix: %s\n", UU_matrix_byRes_->legend());
    if (uuResMatrixFile_ != 0)
      mprintf("\tWriting solute-solute residue matrix to '%s'\n", uuResMatrixFile_->DataFilename().full());
    if (UUmatByRes_norm_ == NORM_NONE)
      mprintf("\tNot normalizing solute-solute residue matrix.\n");
    else if (UUmatByRes_norm_ == NORM_FRAMES)
      mprintf("\tNormalizing solute-solute residue matrix by frames.\n");
    else if (UUmatByRes_norm_ == NORM_RESMAX)
      mprintf("\tNormalizing solute-solute residue matrix by max possible h. bonds between residues.\n");
  }
}

/** Set pointer to the master DataSetList, set HB data set name. */
int HbData::InitHbData(DataSetList* dslPtr, std::string const& setNameIn) {
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
  nuuhb_ = 0;
  nuvhb_ = 0;
//  nbridge_ = 0;
  solvent2solute_.clear();

  NumHbonds_ = masterDSL_->AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UU"));
  if (NumHbonds_ == 0) return 1;
  if (nhbout_ != 0) nhbout_->AddDataSet( NumHbonds_ );
  if (calcSolvent_) {
    NumSolvent_ = masterDSL_->AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UV"));
    if (NumSolvent_ == 0) return 1;
    if (nhbout_ != 0) nhbout_->AddDataSet( NumSolvent_ );
    NumBridge_ = masterDSL_->AddSet(DataSet::INTEGER, MetaData(hbsetname_, "Bridge"));
    if (NumBridge_ == 0) return 1;
    if (nhbout_ != 0) nhbout_->AddDataSet( NumBridge_ );
    BridgeID_ = masterDSL_->AddSet(DataSet::STRING, MetaData(hbsetname_, "ID"));
    if (BridgeID_ == 0) return 1;
    if (nhbout_ != 0) nhbout_->AddDataSet( BridgeID_ );
  }
  if (do_uuResMatrix_) {
    UU_matrix_byRes_ = (DataSet_2D*)
                       masterDSL_->AddSet(DataSet::MATRIX_DBL, MetaData(hbsetname_, "UUresmat"));
    if (UU_matrix_byRes_ == 0) return 1;
    UU_matrix_byRes_->ModifyDim(Dimension::X).SetLabel("Res");
    UU_matrix_byRes_->ModifyDim(Dimension::Y).SetLabel("Res");
    if (uuResMatrixFile_ != 0)
      uuResMatrixFile_->AddDataSet( UU_matrix_byRes_ );
    UU_norm_byRes_ = new DataSet_MatrixDbl();
  }

  return 0;
}

/** Do any Topology-related setup. Set pointer to current Topology.
  * \param topPtr Pointer to current topology.
  * \param both_indices Indices of heavy atoms that can be hbond donor and acceptor
  * \param donor_indices Indices of heavy atoms that can only be hbond donor
  * \param donor_h_indices Indices of donor hydrogens
  * \param acceptor_indices Indices of heavy atoms that can only be hbond acceptor.
  */
int HbData::SetCurrentParm( Topology const* topPtr,
                             Iarray const& both_indices,
                             Iarray const& donor_indices,
                             Iarray const& donor_h_indices,
                             Iarray const& acceptor_indices )
{
  CurrentParm_ = topPtr;

  // For backwards compat. store donor/acceptor indices for determining data set index.
  if (series_) {
    // Donor hydrogen indices
    //std::sort(at_array.begin(), at_array.end()); // Should already be sorted
    int didx = 0;
    for (Iarray::const_iterator at = donor_h_indices.begin();
                                at != donor_h_indices.end(); ++at)
      DidxMap_[*at] = didx++;
    // Acceptor indices
    Iarray at_array;
    at_array.reserve( acceptor_indices.size() + both_indices.size() );
    at_array = both_indices;
    for (Iarray::const_iterator it = acceptor_indices.begin(); it != acceptor_indices.end(); ++it)
      at_array.push_back( *it );
    std::sort(at_array.begin(), at_array.end());
    int aidx = 0;
    for (Iarray::const_iterator at = at_array.begin(); at != at_array.end(); ++at)
      AidxMap_[*at] = aidx++;
  }
  // Set up interaction matrix if needed
  if (UU_matrix_byRes_ != 0) {
    if (SetupInteractionMatrix(both_indices, donor_indices, acceptor_indices)) {
      mprinterr("Error: Could not set up hydrogen bond interaction matrix.\n");
      return 1;
    }
  }
  return 0;
}

/** Set up hbond UU matrix */
int HbData::SetupInteractionMatrix(Iarray const& both_indices,
                                   Iarray const& donorOnly_indices,
                                   Iarray const& acceptorOnly_indices)
{
  if (UU_matrix_byRes_ == 0) return 0;

  unsigned int bothEnd = both_indices.size();
  Iarray Both;
  Both.reserve( both_indices.size() + donorOnly_indices.size() );
  Both = both_indices;
  for (Iarray::const_iterator it = donorOnly_indices.begin(); it != donorOnly_indices.end(); ++it)
    Both.push_back( *it );

  // Sanity check
  if (CurrentParm_ == 0) {
    mprinterr("Internal Error: HbData::SetupInteractionMatrix(): CurrentParm_ is null.\n");
    return 1;
  }

  // Find highest solute index
  int highest_U_idx = 0;
  for (int rnum = 0; rnum < CurrentParm_->Nres(); rnum++)
  {
    // Find molecule number
    int at0 = CurrentParm_->Res(rnum).FirstAtom();
    int mnum = (*CurrentParm_)[at0].MolNum();
    if (!CurrentParm_->Mol(mnum).IsSolvent())
      highest_U_idx = rnum;
  }
  mprintf("\tHighest solute residue # = %i\n", highest_U_idx+1);
  if (UU_matrix_byRes_->AllocateHalf(highest_U_idx+1)) {
    mprinterr("Error: Allocating solute-solute hbond matrix failed.\n");
    return 1;
  }
  // Set up normalization matrix if necesary
  if (UUmatByRes_norm_ == NORM_RESMAX) {
    UU_norm_byRes_->AllocateHalf( UU_matrix_byRes_->Nrows() );
    // Outer loop over solute donors
    for (unsigned int sidx0 = 0; sidx0 < Both.size(); sidx0++)
    {
      int atidx0 = Both[sidx0];
      int mol0 = (*CurrentParm_)[atidx0].MolNum();
      int rnum0 = (*CurrentParm_)[atidx0].ResNum();
      // Loop over solute sites that can be both donor and acceptor
      for (unsigned int sidx1 = sidx0 + 1; sidx1 < bothEnd; sidx1++)
      {
        int atidx1 = Both[sidx1];
        if (!noIntramol_ || mol0 != (*CurrentParm_)[atidx1].MolNum()) {
          int rnum1 = (*CurrentParm_)[atidx1].ResNum();
          // The max # of hbonds between sites depends on # hydrogens.
          // However, given H1-X-H2 Y, you tend to have either Y..H1-X or
          // Y..H2-X, not both, so do based on sites only.
          //UU_norm_byRes_.UpdateElement(rnum0, rnum1, Site0.n_hydrogens() + Site1.n_hydrogens());
          UU_norm_byRes_->UpdateElement(rnum0, rnum1, 2.0);
        }
      } // END loop over solute sites that are D/A
      // Loop over solute acceptor-only
      for (Iarray::const_iterator a_atom = acceptorOnly_indices.begin();
                                  a_atom != acceptorOnly_indices.end(); ++a_atom)
      {
        if (!noIntramol_ || mol0 != (*CurrentParm_)[*a_atom].MolNum()) {
          int rnum1 = (*CurrentParm_)[*a_atom].ResNum();
          //UU_norm_byRes_.UpdateElement(rnum0, rnum1, Site0.n_hydrogens());
          UU_norm_byRes_->UpdateElement(rnum0, rnum1, 1.0);
        }
      } // END loop over acceptor atoms
    } // END loop over all D/A sites
    mprintf("DEBUG: Residue normalization matrix:\n");
    for (unsigned int rnum0 = 0; rnum0 < UU_norm_byRes_->Nrows(); rnum0++) {
      for (unsigned int rnum1 = rnum0; rnum1 < UU_norm_byRes_->Ncols(); rnum1++) {
        mprintf("\t%6i %6i %f\n", rnum0+1, rnum1+1, UU_norm_byRes_->GetElement(rnum1, rnum0));
      }
    }
  }
  return 0;
}

/** Create legend for hydrogen bond based on given atoms. */
std::string HbData::CreateHBlegend(Topology const& topIn, int a_atom, int h_atom, int d_atom)
{
  if (a_atom == ID_SOLVENT_)
    return (topIn.TruncResAtomName(h_atom) + "-V");
  else if (a_atom == ID_ION_)
    return (topIn.TruncResAtomName(h_atom) + "-I");
  else if (d_atom == ID_SOLVENT_)
    return (topIn.TruncResAtomName(a_atom) + "-V");
  else if (d_atom == ID_ION_)
    return (topIn.TruncResAtomName(a_atom) + "-I");
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
  nuuhb_++;
}

/** Add or update solute-solvent hydrogen bond, track bridging water. */
void HbData::AddUV(double dist, double angle, int fnum,
                   int a_atom, int h_atom, int d_atom, bool udonor, int onum)
{
  int hb_id;
  std::string aspect;
  if (d_atom == h_atom) {
    hb_id = ID_ION_;
    aspect = "ionhb";
  } else {
    hb_id = ID_SOLVENT_;
    aspect = "solventhb";
  }
  // TODO return if not calcSolvent_?
  int hbidx, solventres, soluteres;
  // TODO: Option to use solvent mol num?
  if (udonor) {
    // Index U-H .. V hydrogen bonds by solute H atom.
    hbidx = h_atom;
    if (bridgeByAtom_)
      soluteres = h_atom;
    else
      soluteres = (*CurrentParm_)[d_atom].ResNum();
    solventres = (*CurrentParm_)[a_atom].ResNum();
  } else {
    // Index U .. H-V hydrogen bonds by solute A atom.
    hbidx = a_atom;
    if (bridgeByAtom_)
      soluteres = a_atom;
    else
      soluteres = (*CurrentParm_)[a_atom].ResNum();
    solventres = (*CurrentParm_)[d_atom].ResNum();
  }
  solvent2solute_[solventres].insert( soluteres );
  UVmapType::iterator it = UV_Map_.lower_bound( hbidx );
  if (it == UV_Map_.end() || it->first != hbidx)
  {
//      mprintf("DBG1: NEW hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
    DataSet_integer* ds = 0;
    if (series_) {
      ds = (DataSet_integer*)
           masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_, aspect, hbidx));
      if (UVseriesout_ != 0) UVseriesout_->AddDataSet( ds );
    }
    Hbond hb;
    if (udonor) { // Do not care about which solvent acceptor
      if (ds != 0) ds->SetLegend( CreateHBlegend(*CurrentParm_, hb_id, h_atom, d_atom) );
      hb = Hbond(ds, hb_id, h_atom, d_atom, splitFrames_);
    } else {           // Do not care about which solvent donor
      if (ds != 0) ds->SetLegend( CreateHBlegend(*CurrentParm_, a_atom, hb_id, hb_id) );
      hb = Hbond(ds, a_atom, hb_id, hb_id, splitFrames_);
    }
    it = UV_Map_.insert(it, std::pair<int,Hbond>(hbidx,hb));
  } else {
//      mprintf("DBG1: OLD hbond : %8i .. %8i - %8i\n", a_atom+1,h_atom+1,d_atom+1);
  }
  it->second.Update(dist, angle, fnum, splitFrames_, onum);
  nuvhb_++;
}

/** Create legend for bridge based on given indices. */
std::string HbData::CreateBridgeLegend(std::string const& prefix, std::set<int> const& indices)
{
  std::string blegend(prefix);
  for (std::set<int>::const_iterator brs = indices.begin(); brs != indices.end(); ++brs)
    blegend.append("_" + integerToString(*brs + 1));
  return blegend;
}

/** Calculate bridging water. */
void HbData::BridgeCalc(int frameNum, int trajoutNum) {
    int numHB = 0;
    std::string bridgeID;
    for (RmapType::const_iterator bridge = solvent2solute_.begin();
                                  bridge != solvent2solute_.end(); ++bridge)
    {
      // bridge->first is solvent residue number.
      // bridge->second is a set of solute residue numbers the solvent
      // residue is bound to.
      // If solvent molecule is bound to 2 or more different solute residues,
      // it is bridging. 
      if ( bridge->second.size() > 1) {
        bool isBridge = true;
        if (noIntramol_) {
          // If all residues belong to the same molecule and 'nointramol',
          // do not consider this bridge.
          int firstmol = -1;
          unsigned int nequal = 1;
          for (std::set<int>::const_iterator res = bridge->second.begin();
                                             res != bridge->second.end(); ++res)
          {
            int currentMol;
            if (bridgeByAtom_)
              currentMol = (*CurrentParm_)[*res].MolNum();
            else
              currentMol = (*CurrentParm_)[CurrentParm_->Res(*res).FirstAtom()].MolNum();
            if ( firstmol == -1 )
              firstmol = currentMol;
            else if (currentMol == firstmol)
              ++nequal;
          }
          isBridge = (nequal < bridge->second.size());
        }
        if (isBridge) {
          // numHB is used to track the number of bridges
          ++numHB;
          // Bridging Solvent residue number
          bridgeID.append(integerToString( bridge->first+1 ) + "(");
          // Loop over solute residues this solvent is bound to.
          for (std::set<int>::const_iterator res = bridge->second.begin();
                                             res != bridge->second.end(); ++res)
            // Solute residue number being bridged
            bridgeID.append( integerToString( *res+1 ) + "+" );
          bridgeID.append("),");
          // Find bridge in map based on this combo of residues (bridge->second)
          BmapType::iterator b_it = BridgeMap_.lower_bound( bridge->second );
          if (b_it == BridgeMap_.end() || b_it->first != bridge->second) {
            // New Bridge
            DataSet_integer* bds = 0; 
            if (Bseries_) {
              bds = (DataSet_integer*)
                masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,CreateBridgeLegend("bridge",bridge->second),BridgeMap_.size()));
              // Create a legend from the indices.
              bds->SetLegend( CreateBridgeLegend( "B", bridge->second ) );
              if (Bseriesout_ != 0) Bseriesout_->AddDataSet( bds );
            }
            b_it = BridgeMap_.insert( b_it, std::pair<std::set<int>,Bridge>(bridge->second, Bridge(bds, splitFrames_)) );
          }
          // Increment bridge #frames
          b_it->second.Update(frameNum, splitFrames_, trajoutNum);
        }
      }
    } // END LOOP OVER solvent2solute_
    if (bridgeID.empty())
      bridgeID.assign("None");
    NumBridge_->Add(frameNum, &numHB);
    BridgeID_->Add(frameNum, bridgeID.c_str());
//#   ifdef TIMER
//    t_bridge_.Stop();
//#   endif
}

/** Finish hbond calc for a Frame. */
void HbData::IncrementNframes(int frameNum, int trajoutNum) {
  if (NumHbonds_ != 0)
    NumHbonds_->Add( frameNum, &nuuhb_ );
  if (calcSolvent_) {
    if (NumSolvent_ != 0) NumSolvent_->Add( frameNum, &nuvhb_ );
    BridgeCalc(frameNum, trajoutNum);
//    if (NumBridge_ != 0) NumBridge_->Add( Nframes_, &nbridge_ );
  }
  nuuhb_ = 0;
  nuvhb_ = 0;
//  nbridge_ = 0;
  solvent2solute_.clear();
  Nframes_++;
}

/** Estimate the memory usage of the hbond command. */
std::string HbData::MemoryUsage(size_t n_uu_pairs, size_t n_uv_pairs, size_t nFramesIn) const
{
  size_t nFrames;
  if (nFramesIn < 1)
    nFrames = masterDSL_->MaxFrames();
  else
    nFrames = nFramesIn;
  mprintf("DEBUG: nuu %zu nuv %zu nframes %zu\n", n_uu_pairs, n_uv_pairs, nFrames);
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
  if (UU_matrix_byRes_ != 0) {
    memTotal += UU_matrix_byRes_->MemUsageInBytes();
    memTotal += UU_norm_byRes_->MemUsageInBytes();
  }
 
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
      if (hbond->A() < 0) {
        // Solvent/ion acceptor
        if (hbond->A() == ID_SOLVENT_)
          Aname = "SolventAcc";
        else
          Aname = "IonAcc";
      } else {
        Aname = CurrentParm_->TruncResAtomName(hbond->A());
        if (useAtomNum_) Aname.append("_" + integerToString(hbond->A()+1));
      }
      if (hbond->D() < 0) {
        // Solvent/ion donor
        if (hbond->D() == ID_SOLVENT_) {
          Dname = "SolventDnr";
          Hname = "SolventH";
        } else {
          Dname = "Ion";
          Hname = "Ion";
        }
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

