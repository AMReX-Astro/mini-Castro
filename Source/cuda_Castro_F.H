#ifdef AMREX_USE_CUDA
extern "C" {
__device__ void ca_ctoprim(const int* lo, const int* hi,
                                    const amrex::Real* u, const int* u_lo, const int* u_hi,
                                    const amrex::Real* q, const int* q_lo, const int* q_hi,
                                    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

__global__ void cuda_ca_ctoprim(const int* lo, const int* hi,
                                    const amrex::Real* u, const int* u_lo, const int* u_hi,
                                    const amrex::Real* q, const int* q_lo, const int* q_hi,
                                    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

__device__ void ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state));

__global__ void cuda_ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state));

}
#endif
