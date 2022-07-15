#ifndef INC_WRAPPER_WATERSHELL_CUH
#define INC_WRAPPER_WATERSHELL_CUH
#include "../Gpu.h"
#include "../ImageOption.h"
int Cpptraj_GPU_WaterShell(int&, int&, CpptrajGpu::FpType, CpptrajGpu::FpType,
                    const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                    ImageOption::Type, const CpptrajGpu::FpType*,
                    const CpptrajGpu::FpType*, const CpptrajGpu::FpType*);
#endif
