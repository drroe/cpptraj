#include "wrapper_WaterShell.cuh"
#include "kWaterShell.cuh"
#include "cuda_box.cuh"
#include "../CpptrajStdio.h"
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#include "../HipDefinitions.h"
#endif

static inline int calc_nblocks(int ntotal, int nthreadsPerBlock)
{
  int nblocks = ntotal / nthreadsPerBlock;
  if ( (ntotal % nthreadsPerBlock) != 0 )
    nblocks++;
  return nblocks;
}

/** Report any cuda errors. */
static inline int Cuda_check(cudaError_t err, const char* desc) {
  //cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    mprintf("Warning: CUDA Error %s: %s\n", desc, cudaGetErrorString(err));
    mprinterr("Error: CUDA Error %s: %s\n", desc, cudaGetErrorString(err));
    //return 1;
  }
  return 0;
}

/** Calculate distances between pairs of atoms and bin them into a 1D histogram. */
int Cpptraj_GPU_WaterShell(int& nlower, int& nupper,
                           CpptrajGpu::FpType lowerCut2, CpptrajGpu::FpType upperCut2,
                           const CpptrajGpu::FpType* xyz1, int N1,
                           const CpptrajGpu::FpType* xyz2, int N2,
                           ImageOption::Type imageType,
                           CpptrajGpu::HostBox<CpptrajGpu::FpType> const& box)
{
  int* device_counts;
  Cuda_check(cudaMalloc(((void**)(&device_counts)), 2 * sizeof(int)), "Allocating watershell bins");
  Cuda_check(cudaMemset( device_counts, 0, 2*sizeof(int) ), "Setting watershell bins to 0");

  CpptrajGpu::FpType* device_xyz1;
  Cuda_check(cudaMalloc(((void**)(&device_xyz1)), N1 * 3 * sizeof(CpptrajGpu::FpType)), "Allocating xyz1");
  Cuda_check(cudaMemcpy(device_xyz1, xyz1, N1 * 3 * sizeof(CpptrajGpu::FpType), cudaMemcpyHostToDevice), "Copying xyz1");

  CpptrajGpu::FpType* device_xyz2;
  Cuda_check(cudaMalloc(((void**)(&device_xyz2)), N2 * 3 * sizeof(CpptrajGpu::FpType)), "Allocating xyz2");
  Cuda_check(cudaMemcpy(device_xyz2, xyz2, N2 * 3 * sizeof(CpptrajGpu::FpType), cudaMemcpyHostToDevice), "Copying xyz2");

  cuda_box<CpptrajGpu::FpType> gpuBox;
  if ( gpuBox.Setup( imageType, box.BoxLengths(), box.Ucell(), box.Frac() ) ) {
    Cuda_check( gpuBox.LastErr(), gpuBox.LastErrDesc() );
    return 1;
  }

  // Determine number of blocks
  unsigned int BLOCKDIM = CpptrajGpu::MaxBlockDim_2D();

  dim3 threadsPerBlock(BLOCKDIM, BLOCKDIM);
  dim3 numBlocks(calc_nblocks(N1, threadsPerBlock.x), calc_nblocks(N2, threadsPerBlock.y));
  mprintf("#Atoms = %i, %i; Threads per block = %i, %i;  #Blocks = %i, %i\n",
          N1, N2, threadsPerBlock.x, threadsPerBlock.y, numBlocks.x, numBlocks.y);

  // Launch kernel
  // Must have Non-overlapping coords
  mprintf("DEBUG: before launch: nlower= %i nupper= %i\n", nlower, nupper);
  switch (imageType) {
      case ImageOption::NONORTHO:
        kWaterShell_nonOrtho<<<numBlocks, threadsPerBlock>>>(
          device_xyz1, N1, device_xyz2, N2, gpuBox.FracDev(), gpuBox.UcellDev(), lowerCut2, upperCut2, device_counts);
      break;
      case ImageOption::ORTHO:
        kWaterShell_Ortho<<<numBlocks, threadsPerBlock>>>(
          device_xyz1, N1, device_xyz2, N2, gpuBox.BoxDev(), lowerCut2, upperCut2, device_counts);
      break;
      case ImageOption::NO_IMAGE:
        kWaterShell_NoImage<<<numBlocks, threadsPerBlock>>>(
          device_xyz1, N1, device_xyz2, N2, lowerCut2, upperCut2, device_counts);
      break;
      //default:
      //  mprinterr("Internal Error: kernel_rdf: Unhandled image type.\n");
      //  return 1;
  }
  mprintf("DEBUG: after launch: nlower= %i nupper= %i\n", nlower, nupper);

  // Error check
  Cuda_check(cudaGetLastError(), "watershell kernel launch");

  // Copy the result back
  int local_counts[2];
  cudaMemcpy(local_counts, device_counts, 2*sizeof(int), cudaMemcpyDeviceToHost);
  nlower = local_counts[0];
  nupper = local_counts[1];

  // Free device memory
  cudaFree(device_counts);
  cudaFree(device_xyz1);
  cudaFree(device_xyz2);

  return 0;
}
