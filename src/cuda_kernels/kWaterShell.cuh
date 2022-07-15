#ifndef INC_KWATERSHELL_CUH
#define INC_KWATERSHELL_CUH
#include "../Gpu.h"
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#endif
// ----- Device WaterShell Kernels ---------------
// WaterShell no imaging
__global__ void kWaterShell_NoImage(const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                    CpptrajGpu::FpType, CpptrajGpu::FpType, int*);
// WaterShell ortho imaging
__global__ void kWaterShell_Ortho(const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                  const CpptrajGpu::FpType*,
                                  CpptrajGpu::FpType, CpptrajGpu::FpType, int*);
// WaterShell nonortho imaging
__global__ void kWaterShell_nonOrtho(const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                     const CpptrajGpu::FpType*, const CpptrajGpu::FpType*,
                                     CpptrajGpu::FpType, CpptrajGpu::FpType, int*);
#endif
