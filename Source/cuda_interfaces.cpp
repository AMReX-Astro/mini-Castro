#include <Castro.H>
#include <cuda_Castro_F.H>
#include <AMReX_BLFort.H>
#include <AMReX_Device.H>

__global__ void cuda_ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state))
{

   int blo[3];
   int bhi[3];
   get_loop_bounds(blo, bhi, lo, hi);
   ca_compute_temp
    (blo, bhi, BL_FORT_FAB_VAL_3D(state));
}
