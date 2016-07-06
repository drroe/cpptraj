#include "DataSet_Cmatrix.h"
#include "CpptrajStdio.h"

void DataSet_Cmatrix::PrintElements() const {
  // NOTE: Matrix is square upper triangle, Nrows == Ncols
  for (unsigned int row = 0; row != Nrows(); row++)
    for (unsigned int col = row + 1; col != Nrows(); col++)
      mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
}

/** Set up sieving info as necessary and set up cluster matrix based on actual
  * number of frames to be clustered.
  */
int DataSet_Cmatrix::SetupWithSieve(ClusterDist* CdistIn, size_t sizeIn, int sieveIn, int iseed)
{
  if (CdistIn == 0) {
    mprinterr("Internal Error: DataSet_Cmatrix::SetupWithSieve called with empty ClusterDist.\n");
    return 1;
  }
  metricDescription_.assign( CdistIn->Description() );
  if (sievedFrames_.SetSieve( sieveIn, sizeIn, iseed )) return 1;
  if (AllocateCmatrix( sievedFrames_.ActualNframes() )) return 1;
  if (SetCdist(CdistIn)) return 1;
  if (sievedFrames_.Type() != ClusterSieve::NONE)
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %i sieved frames.\n",
            sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes());
  else
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  return 0;
}

#ifdef MPI
/** Set up sieving info in parallel and set up cluster matrix based on
  * actual number of frames to be clustered. Not valid for no sieve.
  */
int DataSet_Cmatrix::SetupWithParallelSieve(ClusterDist* CdistIn, size_t sizeIn,
                                            int sieveIn, int iseed,
                                            Parallel::Comm const& commIn)
{
  if (sieveIn == 1) {
    mprinterr("Internal Error: SetupWithParallelSieve called with sieve == 1\n");
    return 1;
  }
  if (CdistIn == 0) {
    mprinterr("Internal Error: SetupWithParallelSieve called with empty ClusterDist.\n");
    return 1;
  }
  metricDescription_.assign( CdistIn->Description() );
  ClusterSieve::SieveErr err = sievedFrames_.SetParallelSieve( sieveIn, sizeIn, iseed, commIn );
  switch (err) {
    case ClusterSieve::OK: err = ClusterSieve::OK; break;
    case ClusterSieve::NO_FRAMES: mprinterr("Error: No input frames.\n"); break;
    case ClusterSieve::SIEVE_LARGER_THAN_THREADS:
      mprinterr("Error: Number of threads %i is greater than sieve value %i.\n",
                commIn.Size(), sieveIn);
      break;
    case ClusterSieve::NO_SIEVE: mprinterr("Error: No sieve not valid in parallel.\n"); break;
    case ClusterSieve::OTHER: mprinterr("Error: Could not set up sieve.\n"); break;
  }
  if (err != ClusterSieve::OK) return 1;
  if (AllocateCmatrix( sievedFrames_.ActualNframes() )) return 1;
  if (SetCdist(CdistIn)) return 1;
  if (sievedFrames_.Type() != ClusterSieve::NONE)
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %i sieved frames.\n",
            sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes());
  else
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  return 0;
}
#endif

/** Set up sieve info from an array that contains 'T' if the frame was sieved
  * out and 'F' otherwise.
  */
int DataSet_Cmatrix::SetSieveFromArray(std::vector<char> const& sieveStatus, int sieveIn)
{
  if (sieveStatus.empty()) return 1;
  // Setup sieve class
  if (sievedFrames_.SetSieve( sieveIn, sieveStatus )) {
    mprinterr("Error: Could not set sieve from cluster matrix file.\n");
    return 1;
  }
  mprintf("\tSet up %s: %u original frames, %u actual frames, %u elements, sieve=%i\n",
          legend(), sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes(), Nelements(), sieveIn);
  return 0;
}
