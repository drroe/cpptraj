#include "HbParallel.h"
#include "Hbond.h"

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
HbParallel::HbParallel() {}

#ifdef MPI
/** Determine how many hydrogen bonds are on each rank. */
std::vector<int> HbParallel::GetRankNhbonds( int num_hb, Parallel::Comm const& commIn )
{
  std::vector<int> nhb_on_rank;
  if (commIn.Master())
    nhb_on_rank.resize( commIn.Size() );
  commIn.GatherMaster( &num_hb, 1, MPI_INT, &nhb_on_rank[0] );
  return nhb_on_rank;
}

/** Add Hbond class to flat arrays. */
void HbParallel::HbondToArray(std::vector<double>& Dvals, std::vector<int>& Ivals, Hbond const& hb)
{
  Dvals.push_back( hb.Dist() );
  Dvals.push_back( hb.Angle() );
  for (unsigned int idx = 0; idx != hb.Nparts(); idx++) {
    Dvals.push_back( hb.PartDist(idx).nData() );
    Dvals.push_back( hb.PartDist(idx).mean() );
    Dvals.push_back( hb.PartDist(idx).M2() );
    Dvals.push_back( hb.PartAngle(idx).nData() );
    Dvals.push_back( hb.PartAngle(idx).mean() );
    Dvals.push_back( hb.PartAngle(idx).M2() );
  }
  Ivals.push_back( hb.A() );
  Ivals.push_back( hb.H() );
  Ivals.push_back( hb.D() );
  Ivals.push_back( hb.Frames() );
}

/** PARALLEL NOTES:
  * The following tags are used for MPI send/receive:
  *   1300  : Array containing hbond double info on rank.
  *   1301  : Array containing hbond integer info on rank.
  *   1302  : Number of bridges to expect from rank.
  *   1303  : Array containing bridge integer info on rank.
  *   1304+X: Array of hbond X series info from rank.
  */
int HbParallel::SyncToMaster(int& Nframes_, Iarray const& splitFrames_) const {
  // Make sure all time series are updated at this point.
  UpdateSeries();
  // TODO consolidate # frames / offset calc code with Action_NAstruct
  // Get total number of frames
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &Nframes_, 1, MPI_INT, &rank_frames[0] );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      Nframes_ += rank_frames[rank];
  }
  // Convert rank frames to offsets.
  std::vector<int> rank_offsets( trajComm_.Size(), 0 );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      rank_offsets[rank] = rank_offsets[rank-1] + rank_frames[rank-1];
  }

  // Need to send hbond data from all ranks to master.
  std::vector<double> Dvals;           // Hold dist_ and angle_ for each hbond, (as well as n_/mean_/M2_ for dist/angle each part)
  std::vector<int> Ivals;              // Hold A_, H_, D_, and frames_ for each hbond
  unsigned int dvalsPerHbond;
  unsigned int nParts;
  if (!splitFrames_.empty()) {
    nParts = splitFrames_.size() + 1;
    dvalsPerHbond = 2 + (nParts * 6);
  } else {
    nParts = 0;
    dvalsPerHbond = 2;
  }
  // Need to know how many hbonds on each process.
  std::vector<int> nuu_on_rank = GetRankNhbonds( UU_Map_.size(), trajComm_ );
  std::vector<int> nuv_on_rank = GetRankNhbonds( UV_Map_.size(), trajComm_ );
  if (trajComm_.Master()) {
    // MASTER RANK
    for (int rank = 1; rank < trajComm_.Size(); rank++)
    {
      int n_total_on_rank = nuu_on_rank[rank] + nuv_on_rank[rank];
      if (n_total_on_rank > 0)
      {
        std::vector<DataSet_integer*> Svals;
        if (series_) Svals.reserve( n_total_on_rank );
        Dvals.resize( dvalsPerHbond * n_total_on_rank );
        Ivals.resize( 4             * n_total_on_rank );
        trajComm_.Recv( &(Dvals[0]), Dvals.size(), MPI_DOUBLE, rank, 1300 );
        trajComm_.Recv( &(Ivals[0]), Ivals.size(), MPI_INT,    rank, 1301 );
        // Loop over all received hydrogen bonds
        // Dvals = dist, angle
        // Avals = A, H, D, frames
        const int* IV = &Ivals[0];
        const double* DV = &Dvals[0];
        // UU Hbonds
        for (int in = 0; in != nuu_on_rank[rank]; in++, IV += 4, DV += dvalsPerHbond)
        {
          Hpair hbidx(IV[1], IV[0]);
          UUmapType::iterator it = UU_Map_.lower_bound( hbidx );
          DataSet_integer* ds = 0;
          if (it == UU_Map_.end() || it->first != hbidx)
          {
            // Hbond on rank that has not been found on master
            if (series_)
              ds = UUset(IV[0], IV[1], IV[2]);
            it = UU_Map_.insert(it, std::pair<Hpair,Hbond>(hbidx,Hbond(DV[0],DV[1],ds,IV[0],IV[1],IV[2],IV[3])));
            it->second.SetupParts(nParts, DV+2);
          } else {
            // Hbond on rank and master. Update on master.
            it->second.Combine(DV[0], DV[1], IV[3]);
            it->second.CombineParts(nParts, DV+2);
            ds = it->second.Data();
          }
          Svals.push_back( ds );
        }
        // UV Hbonds
        for (int in = 0; in != nuv_on_rank[rank]; in++, IV += 4, DV += dvalsPerHbond)
        {
          int hbidx;
          if (IV[1] != -1)
            // Index U-H .. V hydrogen bonds by solute H atom.
            hbidx = IV[1];
          else
            // Index U .. H-V hydrogen bonds by solute A atom.
            hbidx = IV[0];
          DataSet_integer* ds = 0;
          UVmapType::iterator it = UV_Map_.lower_bound( hbidx );
          if (it == UV_Map_.end() || it->first != hbidx)
          {
            // Hbond on rank that has not been found on master
            if (series_) {
              ds = (DataSet_integer*)
                   masterDSL_->AddSet(DataSet::INTEGER, MetaData(hbsetname_,"solventhb",hbidx));
              ds->SetLegend( CreateHBlegend(*CurrentParm_, IV[0], IV[1], IV[2]) );
              if (UVseriesout_ != 0) UVseriesout_->AddDataSet( ds );
            }
            it = UV_Map_.insert(it, std::pair<int,Hbond>(hbidx,Hbond(DV[0],DV[1],ds,IV[0],IV[1],IV[2],IV[3])));
            it->second.SetupParts(nParts, DV+2);
          } else {
            // Hbond on rank and master. Update on master.
            it->second.Combine(DV[0], DV[1], IV[3]);
            it->second.CombineParts(nParts, DV+2);
            ds = it->second.Data();
          }
          Svals.push_back( ds );
        }
        // Update all time series
        if (series_) {
          for (int in = 0; in != n_total_on_rank; in++) {
            DataSet_integer* ds = Svals[in]; 
            //ds->Resize( Nframes_ );
            //int* d_beg = ds->Ptr() + rank_offsets[ rank ];
            //rprintf("Resizing hbond series data to %i, starting frame %i, # frames %i from rank %i (%i)\n",
            //        Nframes_, rank_offsets[rank], rank_frames[rank], rank, 1304 + in);
            ds->Recv(Nframes_, rank_offsets[ rank ], rank_frames[ rank ],
                     rank, 1304 + in, trajComm_);
            //trajComm_.Recv( d_beg, rank_frames[ rank ], MPI_INT, rank, 1304 + in );
            ds->SetNeedsSync( false );
          }
        }
      }
    } // END master loop over ranks
    // At this point we have all hbond sets from all ranks. Mark all HB sets
    // smaller than Nframes_ as synced and ensure the time series has been
    // updated to reflect overall # frames.
    if (series_) {
      for (UUmapType::iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb)
        FinishSeries( hb->second.Data(), Nframes_ );
      for (UVmapType::iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb)
        FinishSeries( hb->second.Data(), Nframes_ );
    }
  } else {
    // NON-MASTER RANK
    if (!UU_Map_.empty()) {
      unsigned int ntotal = UU_Map_.size() + UV_Map_.size();
      Dvals.reserve( ntotal * dvalsPerHbond );
      Ivals.reserve( ntotal * 4 );
      // Store UU bonds in flat arrays.
      for (UUmapType::const_iterator it = UU_Map_.begin(); it != UU_Map_.end(); ++it) {
        HbondToArray(Dvals, Ivals, it->second);
      }
      // Store UV bonds in flat arrays
      for (UVmapType::const_iterator it = UV_Map_.begin(); it != UV_Map_.end(); ++it) {
        HbondToArray(Dvals, Ivals, it->second);
      }
      trajComm_.Send( &(Dvals[0]), Dvals.size(), MPI_DOUBLE, 0, 1300 );
      trajComm_.Send( &(Ivals[0]), Ivals.size(), MPI_INT,    0, 1301 );
      // Send series data to master
      if (series_) {
        int in = 0; // For tag
        for (UUmapType::const_iterator hb = UU_Map_.begin(); hb != UU_Map_.end(); ++hb, in++) {
          //rprintf("Sending %zu frames to master (%i).\n", hb->second.Data()->Size(), 1304+in);
          hb->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          hb->second.Data()->SetNeedsSync( false );
        }
        for (UVmapType::const_iterator hb = UV_Map_.begin(); hb != UV_Map_.end(); ++hb, in++) {
          //rprintf("Sending %zu frames to master (%i).\n", hb->second.Data()->Size(), 1304+in);
          hb->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          hb->second.Data()->SetNeedsSync( false );
        }
      }
    }
  } // END COMMUNICATING HBOND DATA TO MASTER

  if (calcSolvent_) {
    // Sync bridging data
    // iArray will contain for each bridge: Nres, res1, ..., resN, Frames[, Npart1, ..., NpartN]
    std::vector<int> iArray;
    int iSize;
    if (trajComm_.Master()) {
      // MASTER RANK
      for (int rank = 1; rank < trajComm_.Size(); rank++)
      {
        std::vector<DataSet_integer*> Svals;
        // Receive size of iArray
        trajComm_.Recv( &iSize,           1, MPI_INT, rank, 1302 );
        //mprintf("DEBUG: Receiving %i bridges from rank %i\n", iSize, rank);
        iArray.resize( iSize );
        trajComm_.Recv( &(iArray[0]), iSize, MPI_INT, rank, 1303 );
        unsigned int idx = 0;
        while (idx < iArray.size()) {
          std::set<int> residues;
          unsigned int i2 = idx + 1;
          for (int ir = 0; ir != iArray[idx]; ir++, i2++)
            residues.insert( iArray[i2] );
          BmapType::iterator b_it = BridgeMap_.lower_bound( residues );
          DataSet_integer* bds = 0;
          if (b_it == BridgeMap_.end() || b_it->first != residues ) {
            // Bridge not found on master. Create new Bridge.
            if (Bseries_) {
              bds = (DataSet_integer*)
                masterDSL_->AddSet(DataSet::INTEGER,MetaData(hbsetname_,CreateBridgeLegend("bridge",residues),BridgeMap_.size()));
                // Create a legend from the indices.
                bds->SetLegend( CreateBridgeLegend( "B", residues ) );
            }
            b_it = BridgeMap_.insert( b_it, std::pair<std::set<int>,Bridge>(residues, Bridge(bds, iArray[i2])) );
            b_it->second.SetupParts(nParts, &iArray[0] + i2 + 1);
            //if (bds != 0) mprintf("DEBUG: '%s' was not on master.\n", bds->legend());
            if (Bseriesout_ != 0) Bseriesout_->AddDataSet( bds );
          } else {
            // Bridge on master and rank. Increment bridge #frames.
            bds = b_it->second.Data();
            b_it->second.Combine( iArray[i2] );
            b_it->second.CombineParts(nParts, &iArray[0] + i2 + 1);
            //if (bds != 0) mprintf("DEBUG: '%s' was already on master.\n", bds->legend());
          }
          Svals.push_back( bds );
          idx = i2 + 1 + nParts;
        }
        // Update all time series
        if (Bseries_) {
          for (unsigned int in = 0; in != Svals.size(); in++) {
            DataSet_integer* ds = Svals[in]; 
            //ds->Resize( Nframes_ );
            //int* d_beg = ds->Ptr() + rank_offsets[ rank ];
            //rprintf("Receiving %i frames of bridge series data for %s, starting frame %i, # frames %i from rank %i (%i)\n",
            //        Nframes_, ds->legend(), rank_offsets[rank], rank_frames[rank], rank, 1304 + in);
            ds->Recv(Nframes_, rank_offsets[ rank ], rank_frames[ rank ],
                     rank, 1304 + in, trajComm_);
            //trajComm_.Recv( d_beg, rank_frames[ rank ], MPI_INT, rank, 1304 + in );
            ds->SetNeedsSync( false );
          }
        }
      } // END LOOP OVER MASTER RANKS
      // At this point we have all bridges from all ranks. Mark all bridge sets
      // smaller than Nframes_ as synced and ensure the time series has been
      // updated to reflect overall # frames.
      if (Bseries_) {
        for (BmapType::iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
          FinishSeries( b->second.Data(), Nframes_ );
      }
    } else {
       // NON-MASTER
       // Construct bridge info array.
       for (BmapType::const_iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
       {
         iArray.push_back( b->first.size() ); // # of bridging res
         for ( std::set<int>::const_iterator r = b->first.begin(); r != b->first.end(); ++r)
           iArray.push_back( *r ); // Bridging res
         iArray.push_back( b->second.Frames() ); // # frames
         if (nParts > 0) {
           for (unsigned int part = 0; part != nParts; part++)
             iArray.push_back( b->second.PartFrames(part) );
         }
      }
      // Since the size of each bridge can be different (i.e. differing #s of
      // residues may be bridged), first send size of the transport array.
      iSize = (int)iArray.size();
      trajComm_.Send( &iSize,           1, MPI_INT, 0, 1302 );
      trajComm_.Send( &(iArray[0]), iSize, MPI_INT, 0, 1303 );
      // Send series data to master
      if (Bseries_) {
        int in = 0; // For tag
        for (BmapType::const_iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b, in++) {
          //rprintf("Sending %zu frames of %s to master (%i).\n", b->second.Data()->Size(), b->second.Data()->legend(), 1304+in);
          b->second.Data()->Send( 0, 1304 + in, trajComm_ );
          //trajComm_.Send( hb->second.Data()->Ptr(), hb->second.Data()->Size(), MPI_INT, 0, 1304 + in );
          b->second.Data()->SetNeedsSync( false );
        }
      }
    }
  } // END COMMUNICATING BRIDGE DATA TO MASTER
  return 0;
}
#endif /* MPI */

