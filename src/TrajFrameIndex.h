#ifndef INC_TRAJFRAMEINDEX_H
#define INC_TRAJFRAMEINDEX_H
/// Class used to find absolute position inside a group of trajectories.
template <class T> class TrajFrameIndex {
  public:
    /// CONSTRUCTOR
    TrajFrameIndex() : currentTrajNum_(-1), trajHasChanged_(false) {}
    /// \return Internal index for current traj given a global index 
    int FindIndex(int, T const&);
    /// \return Total # of readable frames in all trajectories
    int MaxFrames(T const&) const;
    /// \return size in bytes used
    unsigned int DataSize() const { return sizeof(int) + sizeof(bool); }
    /// \return true if a different traj was accessed after last call to FindIndex()
    bool TrajHasChanged() const { return trajHasChanged_; }
    /// \return # of last accessed trajectory
    int CurrentTrajNum() const { return currentTrajNum_; }
  private:
    int currentTrajNum_;  ///< # of currently open input trajectory.
    bool trajHasChanged_; ///< True if new traj opened after last call to FindIndex()
};

/** \return Internal index within current trajectory that corresponds to global index */
template <class T> int TrajFrameIndex<T>::FindIndex(int idx, T const& trajArray)
{
  // Determine which trajectory has the desired index.
  int globalOffset = 0; // Internal offset for converting global index to traj index.
  int currentMax = 0;   // Index after which we are in next trajectory.
  int desiredTrajNum = 0;
  for (; desiredTrajNum < (int)trajArray.size(); ++desiredTrajNum) {
    currentMax += trajArray[desiredTrajNum]->Traj().Counter().TotalReadFrames();
    if (idx < currentMax) break;
    globalOffset += trajArray[desiredTrajNum]->Traj().Counter().TotalReadFrames();
  }
  if (desiredTrajNum == (int)trajArray.size()) return -1;
  trajHasChanged_ = (desiredTrajNum != currentTrajNum_);
  currentTrajNum_ = desiredTrajNum;
  // Convert desired index into trajectory internal index.
  int idx_in_traj = idx - globalOffset;
  return trajArray[currentTrajNum_]->Traj().Counter().IdxToFrame(idx_in_traj);
  //int internalIdx = ((idx - globalOffset) * Offsets_[currentTrajNum_]) +
  //                  Starts_[currentTrajNum_];
  //return internalIdx;
}

/** \return total number of readable frames in all trajectories */
template <class T> int TrajFrameIndex<T>::MaxFrames(T const& trajArray) const {
  int maxFrames = 0;
  for (unsigned int idx = 0; idx != trajArray.size(); idx++)
    maxFrames += trajArray[idx]->Traj().Counter().TotalReadFrames();
  return maxFrames;
}
#endif
