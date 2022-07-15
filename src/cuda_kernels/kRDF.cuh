#ifndef INC_KRDF_CUH
#define INC_KRDF_CUH
#include "../Gpu.h"
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#endif
// ----- Device RDF Kernels ----------------------
// RDF no imaging
__global__ void kBinDistances_nonOverlap_NoImage(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                               CpptrajGpu::FpType, CpptrajGpu::FpType);
// RDF ortho imaging
__global__ void kBinDistances_nonOverlap_Ortho(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                               const CpptrajGpu::FpType*, CpptrajGpu::FpType, CpptrajGpu::FpType);
// RDF nonortho imaging
__global__ void kBinDistances_nonOverlap_nonOrtho(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                                  const CpptrajGpu::FpType*, const CpptrajGpu::FpType*, CpptrajGpu::FpType, CpptrajGpu::FpType);
#endif
