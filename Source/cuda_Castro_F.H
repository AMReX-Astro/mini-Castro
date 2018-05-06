
#include <AMReX_BLFort.H>
#include <AMReX_Device.H>

#ifndef _cuda_Castro_F_
#define _cuda_Castro_F_

#ifdef AMREX_USE_CUDA
extern "C" {

__device__ void ca_ctoprim(const int* lo, const int* hi,
                                    const amrex::Real* u, const int* u_lo, const int* u_hi,
                                    const amrex::Real* q, const int* q_lo, const int* q_hi,
                                    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);


__global__ static void cuda_ca_ctoprim(const int* lo, const int* hi,
                                    const amrex::Real* u, const int* u_lo, const int* u_hi,
                                    const amrex::Real* q, const int* q_lo, const int* q_hi,
                                    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi)
{
   int blo[3];
   int bhi[3];
   for (int k = lo[2] + blockIdx.z * blockDim.z + threadIdx.z; k <= hi[2]; k += blockDim.z * gridDim.z) {
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo[1] + blockIdx.y * blockDim.y + threadIdx.y; j <= hi[1]; j += blockDim.y * gridDim.y) {
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo[0] + blockIdx.x * blockDim.x + threadIdx.x; i <= hi[0]; i += blockDim.x * gridDim.x) {
         blo[0] = i;
         bhi[0] = i;
         ca_ctoprim(blo, bhi, u, u_lo, u_hi, q, q_lo, q_hi, qaux, qa_lo, qa_hi);
       }
     }
   }
}

__device__ void ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state));


__global__ static void cuda_ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state))
{
   int blo[3];
   int bhi[3];
   for (int k = lo[2] + blockIdx.z * blockDim.z + threadIdx.z; k <= hi[2]; k += blockDim.z * gridDim.z) {
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo[1] + blockIdx.y * blockDim.y + threadIdx.y; j <= hi[1]; j += blockDim.y * gridDim.y) {
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo[0] + blockIdx.x * blockDim.x + threadIdx.x; i <= hi[0]; i += blockDim.x * gridDim.x) {
         blo[0] = i;
         bhi[0] = i;
         ca_compute_temp
    (blo, bhi, BL_FORT_FAB_VAL_3D(state));
       }
     }
   }
}

__device__ void ca_divu
    (const int* lo, const int* hi,
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(div));


__global__ static void cuda_ca_divu
    (const int* lo, const int* hi,
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(div))
{
   int blo[3];
   int bhi[3];
   for (int k = lo[2] + blockIdx.z * blockDim.z + threadIdx.z; k <= hi[2]; k += blockDim.z * gridDim.z) {
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo[1] + blockIdx.y * blockDim.y + threadIdx.y; j <= hi[1]; j += blockDim.y * gridDim.y) {
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo[0] + blockIdx.x * blockDim.x + threadIdx.x; i <= hi[0]; i += blockDim.x * gridDim.x) {
         blo[0] = i;
         bhi[0] = i;
         ca_divu
    (blo, bhi, dx, BL_FORT_FAB_VAL_3D(q), BL_FORT_FAB_VAL_3D(div));
       }
     }
   }
}


}
#endif

#endif
