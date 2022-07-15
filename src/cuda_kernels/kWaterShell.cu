#include "kWaterShell.cuh"
#include "NonOrtho_dist2.cuh"
#include "ortho_dist2.cuh"
#include <cstdio> // DEBUG

// -----------------------------------------------------------------------------
/** Calculate # waters in 1st and 2nd solvation shells based on distance cutoffs. */
__global__ void kWaterShell_NoImage(const CpptrajGpu::FpType* xyz1, int N1, const CpptrajGpu::FpType* xyz2, int N2,
                                    CpptrajGpu::FpType lowerCut2, CpptrajGpu::FpType upperCut2, int* VatomShell)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    CpptrajGpu::FpType a1x = xyz1[idx1  ];
    CpptrajGpu::FpType a1y = xyz1[idx1+1];
    CpptrajGpu::FpType a1z = xyz1[idx1+2];

    int idx2 = a2 * 3;
    CpptrajGpu::FpType x = a1x - xyz2[idx2  ];
    CpptrajGpu::FpType y = a1y - xyz2[idx2+1];
    CpptrajGpu::FpType z = a1z - xyz2[idx2+2];

    CpptrajGpu::FpType dist2 = (x*x) + (y*y) + (z*z);
    if (dist2 < upperCut2) {
      VatomShell[a1] = 2;
      if (dist2 < lowerCut2) {
        printf("a1= %i  a2= %i  dist= %f lower.\n", a1+1, a2+1, sqrt(dist2));
        VatomShell[a1] = 1;
      }
    }
  }
}

/** Calculate # waters in 1st and 2nd solvation shells based on distance cutoffs. */
__global__ void kWaterShell_Ortho(const CpptrajGpu::FpType* xyz1, int N1, const CpptrajGpu::FpType* xyz2, int N2,
                                  const CpptrajGpu::FpType* box,
                                  CpptrajGpu::FpType lowerCut2, CpptrajGpu::FpType upperCut2, int* VatomShell)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    CpptrajGpu::FpType a1x = xyz1[idx1  ];
    CpptrajGpu::FpType a1y = xyz1[idx1+1];
    CpptrajGpu::FpType a1z = xyz1[idx1+2];

    int idx2 = a2 * 3;
    CpptrajGpu::FpType a2x = xyz2[idx2  ];
    CpptrajGpu::FpType a2y = xyz2[idx2+1];
    CpptrajGpu::FpType a2z = xyz2[idx2+2];

    CpptrajGpu::FpType dist2 = ortho_dist2<CpptrajGpu::FpType>(a1x, a1y, a1z, a2x, a2y, a2z, box);
    if (dist2 < upperCut2) {
      VatomShell[a1] = 2;
      if (dist2 < lowerCut2)
        VatomShell[a1] = 1;
    }
  }
}

/** Calculate # waters in 1st and 2nd solvation shells based on distance cutoffs. */
__global__ void kWaterShell_nonOrtho(const CpptrajGpu::FpType* xyz1, int N1, const CpptrajGpu::FpType* xyz2, int N2,
                                     const CpptrajGpu::FpType* frac, const CpptrajGpu::FpType* ucell,
                                     CpptrajGpu::FpType lowerCut2, CpptrajGpu::FpType upperCut2, int* VatomShell)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    CpptrajGpu::FpType a1x = xyz1[idx1  ];
    CpptrajGpu::FpType a1y = xyz1[idx1+1];
    CpptrajGpu::FpType a1z = xyz1[idx1+2];
    CpptrajGpu::FpType f1x = frac[0]*a1x + frac[1]*a1y + frac[2]*a1z;
    CpptrajGpu::FpType f1y = frac[3]*a1x + frac[4]*a1y + frac[5]*a1z;
    CpptrajGpu::FpType f1z = frac[6]*a1x + frac[7]*a1y + frac[8]*a1z;

    int idx2 = a2 * 3;
    CpptrajGpu::FpType a2x = xyz2[idx2  ];
    CpptrajGpu::FpType a2y = xyz2[idx2+1];
    CpptrajGpu::FpType a2z = xyz2[idx2+2];
    CpptrajGpu::FpType f2x = frac[0]*a2x + frac[1]*a2y + frac[2]*a2z;
    CpptrajGpu::FpType f2y = frac[3]*a2x + frac[4]*a2y + frac[5]*a2z;
    CpptrajGpu::FpType f2z = frac[6]*a2x + frac[7]*a2y + frac[8]*a2z;

    CpptrajGpu::FpType dist2 =  NonOrtho_dist2<CpptrajGpu::FpType>(f2x, f2y, f2z, f1x ,f1y, f1z, ucell);
    if (dist2 < upperCut2) {
      VatomShell[a1] = 2;
      if (dist2 < lowerCut2)
        VatomShell[a1] = 1;
    }
  }
}
