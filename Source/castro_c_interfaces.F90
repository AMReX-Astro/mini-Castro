module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
  use mempool_module, only: bl_allocate, bl_deallocate
#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, &
                       cudaStreamSynchronize, cudaDeviceSynchronize, dim3, cuda_stream_kind
    use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams
#endif

contains

  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi,idx) &
                                     bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_enforce_consistent_e
#endif

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

#ifdef CUDA

    attributes(device) :: state

    integer, managed, pointer :: lo_d(:), hi_d(:)
    integer, managed, pointer :: s_lo_d(:), s_hi_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(s_lo_d, 1, 3)
    call bl_allocate(s_hi_d, 1, 3)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_enforce_consistent_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(s_lo_d)
    call bl_deallocate(s_hi_d)

#else

    call enforce_consistent_e(lo, hi, state, s_lo, s_hi)

#endif

  end subroutine ca_enforce_consistent_e



  subroutine ca_compute_temp(lo, hi, state, s_lo, s_hi, idx) &
                             bind(C, name="ca_compute_temp")

    use castro_util_module, only: compute_temp
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_compute_temp
#endif

    implicit none

    integer,  intent(in   ) :: lo(3),hi(3)
    integer,  intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ) :: idx

#ifdef CUDA

    attributes(device) :: state

    integer, managed, pointer :: lo_d(:), hi_d(:)
    integer, managed, pointer :: s_lo_d(:), s_hi_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)

    call bl_allocate(s_lo_d, 1, 3)
    call bl_allocate(s_hi_d, 1, 3)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_compute_temp<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(s_lo_d)
    call bl_deallocate(s_hi_d)

#else

    call compute_temp(lo, hi, state, s_lo, s_hi)

#endif

  end subroutine ca_compute_temp


  subroutine ca_reset_internal_e(lo, hi, u, u_lo, u_hi, verbose, idx) &
                                 bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_reset_internal_e
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer, intent(in)     :: idx

#ifdef CUDA

    attributes(device) :: u

    integer, managed, pointer :: lo_d(:), hi_d(:)
    integer, managed, pointer :: u_lo_d(:), u_hi_d(:)
    integer, managed, pointer :: verbose_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(u_lo_d, 1, 3)
    call bl_allocate(u_hi_d, 1, 3)
    call bl_allocate(verbose_d, 1, 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(verbose_d, verbose, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_reset_internal_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, u, u_lo_d, u_hi_d, verbose_d(1))

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(u_lo_d)
    call bl_deallocate(u_hi_d)
    call bl_deallocate(verbose_d)

#else

    call reset_internal_e(lo, hi, u, u_lo, u_hi, verbose)

#endif

  end subroutine ca_reset_internal_e


  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi, idx) &
                                  bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_normalize_species
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer,  intent(in   ) :: idx

#ifdef CUDA

    attributes(device) :: u

    integer, managed, pointer :: lo_d(:), hi_d(:)
    integer, managed, pointer :: u_lo_d(:), u_hi_d(:)

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(u_lo_d, 1, 3)
    call bl_allocate(u_hi_d, 1, 3)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_normalize_species<<<numBlocks, numThreads, 0, stream>>>(u, u_lo_d, u_hi_d, lo_d, hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(u_lo_d)
    call bl_deallocate(u_hi_d)

#else

    call normalize_species(u, u_lo, u_hi, lo, hi)

#endif

  end subroutine ca_normalize_species


  subroutine ca_enforce_minimum_density(uin, uin_lo, uin_hi, &
                                        uout, uout_lo, uout_hi, &
                                        vol, vol_lo, vol_hi, &
                                        lo, hi, frac_change, verbose, idx) &
                                        bind(C, name="ca_enforce_minimum_density")

    use advection_util_module, only: enforce_minimum_density
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_enforce_minimum_density
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose
    integer,  intent(in   ) ::  uin_lo(3),  uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) ::  vol_lo(3),  vol_hi(3)

    real(rt), intent(in   ) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in   ) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt), intent(inout) :: frac_change
    integer,  intent(in   ) :: idx

#ifdef CUDA

    attributes(device) :: uin, uout, vol

    integer,  managed, pointer :: lo_d(:), hi_d(:)
    integer,  managed, pointer :: uin_lo_d(:), uin_hi_d(:)
    integer,  managed, pointer :: uout_lo_d(:), uout_hi_d(:)
    integer,  managed, pointer :: vol_lo_d(:), vol_hi_d(:)
    real(rt), managed, pointer :: frac_change_d(:)
    integer,  managed, pointer :: verbose_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    real(rt) :: local_frac_change

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(uin_lo_d, 1, 3)
    call bl_allocate(uin_hi_d, 1, 3)
    call bl_allocate(uout_lo_d, 1, 3)
    call bl_allocate(uout_hi_d, 1, 3)
    call bl_allocate(vol_lo_d, 1, 3)
    call bl_allocate(vol_hi_d, 1, 3)
    call bl_allocate(frac_change_d, 1, 1)
    call bl_allocate(verbose_d, 1, 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uin_lo_d, uin_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(uin_hi_d, uin_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uout_lo_d, uout_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(uout_hi_d, uout_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(vol_lo_d, vol_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(vol_hi_d, vol_hi, 3, cudaMemcpyHostToDevice, stream)

    local_frac_change = frac_change

    cuda_result = cudaMemcpyAsync(frac_change_d, local_frac_change, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(verbose_d, verbose, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_enforce_minimum_density<<<numBlocks, numThreads, 0, stream>>>(uin, uin_lo_d, uin_hi_d, &
                                                                            uout, uout_lo_d, uout_hi_d, &
                                                                            vol, vol_lo_d, vol_hi_d, &
                                                                            lo_d, hi_d, frac_change_d(1), verbose_d(1))

    cuda_result = cudaMemcpyAsync(local_frac_change, frac_change_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize(stream)

    frac_change = min(frac_change, local_frac_change)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(uin_lo_d)
    call bl_deallocate(uin_hi_d)
    call bl_deallocate(uout_lo_d)
    call bl_deallocate(uout_hi_d)
    call bl_deallocate(vol_lo_d)
    call bl_deallocate(vol_hi_d)
    call bl_deallocate(frac_change_d)
    call bl_deallocate(verbose_d)

#else

    call enforce_minimum_density(uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 vol, vol_lo, vol_hi, &
                                 lo, hi, frac_change, verbose)

#endif

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi, idx) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    integer,  intent(in) :: idx

    call check_initial_species(lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species


  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi, idx) bind(C, name = "ca_ctoprim")

    use advection_util_module, only: ctoprim
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_ctoprim
#endif

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: uin_lo(3), uin_hi(3)
    integer,  intent(in) :: q_lo(3), q_hi(3)
    integer,  intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    integer,  intent(in   ) :: idx

#ifdef CUDA

    attributes(device) :: uin, q, qaux

    integer, managed, pointer :: lo_d(:), hi_d(:)
    integer, managed, pointer :: uin_lo_d(:), uin_hi_d(:)
    integer, managed, pointer :: q_lo_d(:), q_hi_d(:)
    integer, managed, pointer :: qa_lo_d(:), qa_hi_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(uin_lo_d, 1, 3)
    call bl_allocate(uin_hi_d, 1, 3)
    call bl_allocate(q_lo_d, 1, 3)
    call bl_allocate(q_hi_d, 1, 3)
    call bl_allocate(qa_lo_d, 1, 3)
    call bl_allocate(qa_hi_d, 1, 3)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uin_lo_d, uin_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(uin_hi_d, uin_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(q_lo_d, q_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(q_hi_d, q_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(qa_lo_d, qa_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(qa_hi_d, qa_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_ctoprim<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, &
                                                            uin, uin_lo_d, uin_hi_d, &
                                                            q,     q_lo_d,   q_hi_d, &
                                                            qaux, qa_lo_d,  qa_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(uin_lo_d)
    call bl_deallocate(uin_hi_d)
    call bl_deallocate(q_lo_d)
    call bl_deallocate(q_hi_d)
    call bl_deallocate(qa_lo_d)
    call bl_deallocate(qa_hi_d)

#else

    call ctoprim(lo, hi, &
                 uin, uin_lo, uin_hi, &
                 q,     q_lo,   q_hi, &
                 qaux, qa_lo,  qa_hi)

#endif

  end subroutine ca_ctoprim


  subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(c, name="ca_dervel")

    use castro_util_module, only: dervel
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_dervel
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3), nv
    integer,  intent(in   ) :: d_lo(3), d_hi(3), nc
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(3,2,nc)
    real(rt), intent(in   ) :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ) :: level, grid_no

#ifdef CUDA

    attributes(device) :: vel, dat

    integer,  managed, pointer :: lo_d(:), hi_d(:)
    integer,  managed, pointer :: nv_d(:), nc_d(:)
    integer,  managed, pointer :: v_lo_d(:), v_hi_d(:)
    integer,  managed, pointer :: d_lo_d(:), d_hi_d(:)
    integer,  managed, pointer :: domlo_d(:), domhi_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)
    real(rt), managed, pointer :: delta_d(:), xlo_d(:), time_d(:), dt_d(:)
    integer,  managed, pointer :: level_d(:), grid_no_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    ! Note that this stream calculation is not ideal because there are
    ! potentially multiple tiles per box.

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(nv_d, 1, 1)
    call bl_allocate(nc_d, 1, 1)
    call bl_allocate(v_lo_d, 1, 3)
    call bl_allocate(v_hi_d, 1, 3)
    call bl_allocate(d_lo_d, 1, 3)
    call bl_allocate(d_hi_d, 1, 3)
    call bl_allocate(domlo_d, 1, 3)
    call bl_allocate(domhi_d, 1, 3)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, nc)
    call bl_allocate(delta_d, 1, 3)
    call bl_allocate(xlo_d, 1, 3)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(dt_d, 1, 1)
    call bl_allocate(level_d, 1, 1)
    call bl_allocate(grid_no_d, 1, 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(nv_d, nv, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(v_lo_d, v_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(v_hi_d, v_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(d_lo_d, d_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(d_hi_d, d_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(domlo_d, domlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(domhi_d, domhi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(delta_d, delta, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_dervel<<<numBlocks, numThreads, 0, stream>>>(vel,v_lo_d,v_hi_d,nv_d(1), &
                                                           dat,d_lo_d,d_hi_d,nc_d(1), &
                                                           lo_d,hi_d,domlo_d,domhi_d, &
                                                           delta_d,xlo_d,time_d(1),dt_d(1), &
                                                           bc_d,level_d(1),grid_no_d(1))

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(nv_d)
    call bl_deallocate(nc_d)
    call bl_deallocate(v_lo_d)
    call bl_deallocate(v_hi_d)
    call bl_deallocate(d_lo_d)
    call bl_deallocate(d_hi_d)
    call bl_deallocate(domlo_d)
    call bl_deallocate(domhi_d)
    call bl_deallocate(bc_d)
    call bl_deallocate(delta_d)
    call bl_deallocate(xlo_d)
    call bl_deallocate(time_d)
    call bl_deallocate(dt_d)
    call bl_deallocate(level_d)
    call bl_deallocate(grid_no_d)

#else

    call dervel(vel,v_lo,v_hi,nv, &
                dat,d_lo,d_hi,nc,lo,hi,domlo, &
                domhi,delta,xlo,time,dt,bc,level,grid_no)

#endif

  end subroutine ca_dervel



  subroutine ca_derpres(p,p_lo,p_hi,np, &
                        u,u_lo,u_hi,nc,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(c, name="ca_derpres")

    use castro_util_module, only: derpres
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_derpres
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3), np
    integer,  intent(in   ) :: u_lo(3), u_hi(3), nc
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),np)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc)
    real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
    integer,  intent(in   ) :: bc(3,2,nc), level, grid_no

#ifdef CUDA

    attributes(device) :: p, u

    integer,  managed, pointer :: lo_d(:), hi_d(:)
    integer,  managed, pointer :: np_d(:), nc_d(:)
    integer,  managed, pointer :: p_lo_d(:), p_hi_d(:)
    integer,  managed, pointer :: u_lo_d(:), u_hi_d(:)
    integer,  managed, pointer :: domlo_d(:), domhi_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)
    real(rt), managed, pointer :: dx_d(:), xlo_d(:), time_d(:), dt_d(:)
    integer,  managed, pointer :: level_d(:), grid_no_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(np_d, 1, 1)
    call bl_allocate(nc_d, 1, 1)
    call bl_allocate(p_lo_d, 1, 3)
    call bl_allocate(p_hi_d, 1, 3)
    call bl_allocate(u_lo_d, 1, 3)
    call bl_allocate(u_hi_d, 1, 3)
    call bl_allocate(domlo_d, 1, 3)
    call bl_allocate(domhi_d, 1, 3)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, nc)
    call bl_allocate(dx_d, 1, 3)
    call bl_allocate(xlo_d, 1, 3)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(dt_d, 1, 1)
    call bl_allocate(level_d, 1, 1)
    call bl_allocate(grid_no_d, 1, 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(np_d, np, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(p_lo_d, p_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(p_hi_d, p_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(domlo_d, domlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(domhi_d, domhi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_derpres<<<numBlocks, numThreads, 0, stream>>>(p,p_lo_d,p_hi_d,np_d(1), &
                                                            u,u_lo_d,u_hi_d,nc_d(1), &
                                                            lo_d,hi_d,domlo_d,domhi_d, &
                                                            dx_d,xlo_d,time_d(1),dt_d(1), &
                                                            bc_d,level_d(1),grid_no_d(1))

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(np_d)
    call bl_deallocate(nc_d)
    call bl_deallocate(p_lo_d)
    call bl_deallocate(p_hi_d)
    call bl_deallocate(u_lo_d)
    call bl_deallocate(u_hi_d)
    call bl_deallocate(domlo_d)
    call bl_deallocate(domhi_d)
    call bl_deallocate(bc_d)
    call bl_deallocate(dx_d)
    call bl_deallocate(xlo_d)
    call bl_deallocate(time_d)
    call bl_deallocate(dt_d)
    call bl_deallocate(level_d)
    call bl_deallocate(grid_no_d)

#else

    call derpres(p,p_lo,p_hi,np, &
                 u,u_lo,u_hi,nc,lo,hi,domlo, &
                 domhi,dx,xlo,time,dt,bc,level,grid_no)

#endif

  end subroutine ca_derpres

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt,idx) bind(C, name="ca_estdt")

    use timestep_module, only: estdt
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_estdt
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt
    integer,  intent(in   ) :: idx

#ifdef CUDA

    attributes(device) :: u

    integer,  managed, pointer :: lo_d(:), hi_d(:)
    integer,  managed, pointer :: u_lo_d(:), u_hi_d(:)
    real(rt), managed, pointer :: dx_d(:)
    real(rt), managed, pointer :: dt_loc_d(:)

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    real(rt) :: dt_loc

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(u_lo_d, 1, 3)
    call bl_allocate(u_hi_d, 1, 3)
    call bl_allocate(dx_d, 1, 3)
    call bl_allocate(dt_loc_d, 1, 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)

    dt_loc = dt

    cuda_result = cudaMemcpyAsync(dt_loc_d, dt, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_estdt<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, u, u_lo_d, u_hi_d, dx_d, dt_loc_d(1))

    cuda_result = cudaMemcpyAsync(dt_loc, dt_loc_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize(stream)

    dt = min(dt, dt_loc)

    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(u_lo_d)
    call bl_deallocate(u_hi_d)
    call bl_deallocate(dx_d)
    call bl_deallocate(dt_loc_d)

#else

    call estdt(lo, hi, u, u_lo, u_hi, dx, dt)

#endif

  end subroutine ca_estdt



  subroutine ca_mol_single_stage(time, &
                                 lo, hi, domlo, domhi, &
                                 uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 q, q_lo, q_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 update, updt_lo, updt_hi, &
                                 dx, dt, &
                                 flux1, f1_lo, f1_hi, &
                                 flux2, f2_lo, f2_hi, &
                                 flux3, f3_lo, f3_hi, &
                                 area1, a1_lo, a1_hi, &
                                 area2, a2_lo, a2_hi, &
                                 area3, a3_lo, a3_hi, &
                                 vol, vol_lo, vol_hi, &
                                 courno, verbose, idx) bind(C, name="ca_mol_single_stage")

    use mol_module, only: prepare_for_fluxes, prepare_profile, construct_flux
    use advection_util_module, only: construct_hydro_update
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_prepare_for_fluxes, cuda_construct_flux, cuda_construct_hydro_update, cuda_prepare_profile
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose, idx
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: updt_lo(3), updt_hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)
    integer,  intent(in   ) :: a1_lo(3), a1_hi(3)
    integer,  intent(in   ) :: a2_lo(3), a2_hi(3)
    integer,  intent(in   ) :: a3_lo(3), a3_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1), f1_lo(2):f1_hi(2), f1_lo(3):f1_hi(3), NVAR)
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1), f2_lo(2):f2_hi(2), f2_lo(3):f2_hi(3), NVAR)
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1), f3_lo(2):f3_hi(2), f3_lo(3):f3_hi(3), NVAR)
    real(rt), intent(in   ) :: area1(a1_lo(1):a1_hi(1), a1_lo(2):a1_hi(2), a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: area2(a2_lo(1):a2_hi(1), a2_lo(2):a2_hi(2), a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: area3(a3_lo(1):a3_hi(1), a3_lo(2):a3_hi(2), a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in   ) :: dx(3), dt, time
    real(rt), intent(inout) :: courno

    integer :: ngf
    integer :: It_lo(3), It_hi(3)
    integer :: st_lo(3), st_hi(3)
    integer :: g_lo(3), g_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer :: k_lo(3), k_hi(3)
    integer :: idir

    real(rt), pointer :: q1(:,:,:,:)
    real(rt), pointer :: q2(:,:,:,:)
    real(rt), pointer :: q3(:,:,:,:)

    real(rt), pointer :: flatn(:,:,:)
    real(rt), pointer :: div(:,:,:)

    real(rt), pointer :: qm(:,:,:,:,:)
    real(rt), pointer :: qp(:,:,:,:,:)

    real(rt), pointer :: sxm(:,:,:,:), sym(:,:,:,:), szm(:,:,:,:)
    real(rt), pointer :: sxp(:,:,:,:), syp(:,:,:,:), szp(:,:,:,:)

#ifdef CUDA
    attributes(device) :: uin, uout, q, qaux, update, flux1, flux2, flux3, area1, area2, area3, vol
    attributes(managed) :: q1, q2, q3, flatn, div, sxm, sxp, sym, syp, szm, szp, qm, qp

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    real(rt), managed, pointer :: time_d(:)
    integer,  managed, pointer :: lo_d(:), hi_d(:), domlo_d(:), domhi_d(:)
    integer,  managed, pointer :: uin_lo_d(:), uin_hi_d(:)
    integer,  managed, pointer :: uout_lo_d(:), uout_hi_d(:)
    integer,  managed, pointer :: q_lo_d(:), q_hi_d(:)
    integer,  managed, pointer :: qa_lo_d(:), qa_hi_d(:)
    integer,  managed, pointer :: updt_lo_d(:), updt_hi_d(:)
    real(rt), managed, pointer :: dx_d(:), dt_d(:)
    integer,  managed, pointer :: f1_lo_d(:), f1_hi_d(:)
    integer,  managed, pointer :: f2_lo_d(:), f2_hi_d(:)
    integer,  managed, pointer :: f3_lo_d(:), f3_hi_d(:)
    integer,  managed, pointer :: a1_lo_d(:), a1_hi_d(:)
    integer,  managed, pointer :: a2_lo_d(:), a2_hi_d(:)
    integer,  managed, pointer :: a3_lo_d(:), a3_hi_d(:)
    integer,  managed, pointer :: vol_lo_d(:), vol_hi_d(:)
    integer,  managed, pointer :: k_lo_d(:), k_hi_d(:)
    integer,  managed, pointer :: st_lo_d(:), st_hi_d(:)
    integer,  managed, pointer :: It_lo_d(:), It_hi_d(:)
    integer,  managed, pointer :: g_lo_d(:), g_hi_d(:)
    real(rt), managed, pointer :: courno_d(:)
    integer,  managed, pointer :: verbose_d(:)
    integer,  managed, pointer :: idir_d(:)

    real(rt) :: courno_loc
#endif

    ngf = 1

    It_lo = lo - 1
    It_hi = hi + 1

    st_lo = lo - 2
    st_hi = hi + 2

    gd_lo = lo
    gd_hi = hi + 1

    g_lo = lo - ngf
    g_hi = hi + ngf

    call bl_allocate(q1, f1_lo(1), f1_hi(1), f1_lo(2), f1_hi(2), f1_lo(3), f1_hi(3), 1, NGDNV)
    call bl_allocate(q2, f2_lo(1), f2_hi(1), f2_lo(2), f2_hi(2), f2_lo(3), f2_hi(3), 1, NGDNV)
    call bl_allocate(q3, f3_lo(1), f3_hi(1), f3_lo(2), f3_hi(2), f3_lo(3), f3_hi(3), 1, NGDNV)

    call bl_allocate(flatn, q_lo(1), q_hi(1), q_lo(2), q_hi(2), q_lo(3), q_hi(3))

    call bl_allocate(div, g_lo(1), g_hi(1), g_lo(2), g_hi(2), g_lo(3), g_hi(3))

    call bl_allocate(qm, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, NQ, 1, 3)
    call bl_allocate(qp, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, NQ, 1, 3)

    call bl_allocate(sxm, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)
    call bl_allocate(sxp, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)
    call bl_allocate(sym, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)
    call bl_allocate(syp, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)
    call bl_allocate(szm, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)
    call bl_allocate(szp, st_lo(1), st_hi(1), st_lo(2), st_hi(2), st_lo(3), st_hi(3), 1, NQ)

#ifdef CUDA

    call bl_allocate(time_d, 1, 1)
    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(domlo_d, 1, 3)
    call bl_allocate(domhi_d, 1, 3)
    call bl_allocate(uin_lo_d, 1, 3)
    call bl_allocate(uin_hi_d, 1, 3)
    call bl_allocate(uout_lo_d, 1, 3)
    call bl_allocate(uout_hi_d, 1, 3)
    call bl_allocate(q_lo_d, 1, 3)
    call bl_allocate(q_hi_d, 1, 3)
    call bl_allocate(qa_lo_d, 1, 3)
    call bl_allocate(qa_hi_d, 1, 3)
    call bl_allocate(updt_lo_d, 1, 3)
    call bl_allocate(updt_hi_d, 1, 3)
    call bl_allocate(dx_d, 1, 3)
    call bl_allocate(dt_d, 1, 1)
    call bl_allocate(f1_lo_d, 1, 3)
    call bl_allocate(f1_hi_d, 1, 3)
    call bl_allocate(f2_lo_d, 1, 3)
    call bl_allocate(f2_hi_d, 1, 3)
    call bl_allocate(f3_lo_d, 1, 3)
    call bl_allocate(f3_hi_d, 1, 3)
    call bl_allocate(a1_lo_d, 1, 3)
    call bl_allocate(a1_hi_d, 1, 3)
    call bl_allocate(a2_lo_d, 1, 3)
    call bl_allocate(a2_hi_d, 1, 3)
    call bl_allocate(a3_lo_d, 1, 3)
    call bl_allocate(a3_hi_d, 1, 3)
    call bl_allocate(vol_lo_d, 1, 3)
    call bl_allocate(vol_hi_d, 1, 3)
    call bl_allocate(k_lo_d, 1, 3)
    call bl_allocate(k_hi_d, 1, 3)
    call bl_allocate(st_lo_d, 1, 3)
    call bl_allocate(st_hi_d, 1, 3)
    call bl_allocate(It_lo_d, 1, 3)
    call bl_allocate(It_hi_d, 1, 3)
    call bl_allocate(g_lo_d, 1, 3)
    call bl_allocate(g_hi_d, 1, 3)
    call bl_allocate(courno_d, 1, 1)
    call bl_allocate(verbose_d, 1, 1)
    call bl_allocate(idir_d, 1, 1)

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(domlo_d, domlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(domhi_d, domhi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uin_lo_d, uin_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(uin_hi_d, uin_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uout_lo_d, uout_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(uout_hi_d, uout_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(q_lo_d, q_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(q_hi_d, q_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(qa_lo_d, qa_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(qa_hi_d, qa_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(updt_lo_d, updt_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(updt_hi_d, updt_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(f1_lo_d, f1_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(f1_hi_d, f1_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(f2_lo_d, f2_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(f2_hi_d, f2_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(f3_lo_d, f3_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(f3_hi_d, f3_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(a1_lo_d, a1_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(a1_hi_d, a1_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(a2_lo_d, a2_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(a2_hi_d, a2_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(a3_lo_d, a3_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(a3_hi_d, a3_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(vol_lo_d, vol_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(vol_hi_d, vol_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(st_lo_d, st_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(st_hi_d, st_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(It_lo_d, It_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(It_hi_d, It_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(g_lo_d, g_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(g_hi_d, g_hi, 3, cudaMemcpyHostToDevice, stream)

    courno_loc = courno

    cuda_result = cudaMemcpyAsync(courno_d, courno_loc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(verbose_d, verbose, 1, cudaMemcpyHostToDevice, stream)

    ! Construct edge states as inputs to flux construction

    k_lo = g_lo
    k_hi = g_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call cuda_prepare_for_fluxes<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, dt_d(1), dx_d, courno_d(1), &
                                                                       q, flatn, q_lo_d, q_hi_d, &
                                                                       div, g_lo_d, g_hi_d, &
                                                                       qaux, qa_lo_d, qa_hi_d, &
                                                                       sxm, sxp, sym, syp, szm, szp, st_lo_d, st_hi_d, &
                                                                       qm, qp, It_lo_d, It_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    cuda_result = cudaMemcpyAsync(courno_loc, courno_d, 1, cudaMemcpyDeviceToHost, stream)

    call cuda_prepare_profile<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, &
                                                                    q, flatn, q_lo_d, q_hi_d, &
                                                                    sxm, sxp, sym, syp, szm, szp, st_lo_d, st_hi_d, &
                                                                    qm, qp, It_lo_d, It_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    ! Compute F^x

    idir = 1

    cuda_result = cudaMemcpyAsync(idir_d, idir, 1, cudaMemcpyHostToDevice, stream)

    k_lo = f1_lo
    k_hi = f1_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call cuda_construct_flux<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, domlo_d, domhi_d, dx_d, dt_d(1), idir_d(1), &
                                                                   div, g_lo_d, g_hi_d, &
                                                                   uin, uin_lo_d, uin_hi_d, &
                                                                   qm, qp, It_lo_d, It_hi_d, &
                                                                   flux1, q1, f1_lo_d, f1_hi_d, &
                                                                   area1, a1_lo_d, a1_hi_d, &
                                                                   qaux, qa_lo_d, qa_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    ! Compute F^y

    idir = 2

    cuda_result = cudaMemcpyAsync(idir_d, idir, 1, cudaMemcpyHostToDevice, stream)

    k_lo = f2_lo
    k_hi = f2_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call cuda_construct_flux<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, domlo_d, domhi_d, dx_d, dt_d(1), idir_d(1), &
                                                                   div, g_lo_d, g_hi_d, &
                                                                   uin, uin_lo_d, uin_hi_d, &
                                                                   qm, qp, It_lo_d, It_hi_d, &
                                                                   flux2, q2, f2_lo_d, f2_hi_d, &
                                                                   area2, a2_lo_d, a2_hi_d, &
                                                                   qaux, qa_lo_d, qa_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    ! Compute F^z

    idir = 3

    cuda_result = cudaMemcpyAsync(idir_d, idir, 1, cudaMemcpyHostToDevice, stream)

    k_lo = f3_lo
    k_hi = f3_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call cuda_construct_flux<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, domlo_d, domhi_d, dx_d, dt_d(1), idir_d(1), &
                                                                   div, g_lo_d, g_hi_d, &
                                                                   uin, uin_lo_d, uin_hi_d, &
                                                                   qm, qp, It_lo_d, It_hi_d, &
                                                                   flux3, q3, f3_lo_d, f3_hi_d, &
                                                                   area3, a3_lo_d, a3_hi_d, &
                                                                   qaux, qa_lo_d, qa_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    ! Create an update source term based on the flux divergence.

    k_lo = lo
    k_hi = hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call cuda_construct_hydro_update<<<numBlocks, numThreads, 0, stream>>>(k_lo_d, k_hi_d, dx_d, dt_d(1), &
                                                                           flux1, q1, f1_lo_d, f1_hi_d, &
                                                                           flux2, q2, f2_lo_d, f2_hi_d, &
                                                                           flux3, q3, f3_lo_d, f3_hi_d, &
                                                                           area1, a1_lo_d, a1_hi_d, &
                                                                           area2, a2_lo_d, a2_hi_d, &
                                                                           area3, a3_lo_d, a3_hi_d, &
                                                                           vol, vol_lo_d, vol_hi_d, &
                                                                           update, updt_lo_d, updt_hi_d)

    cuda_result = cudaStreamSynchronize(stream)

    courno = max(courno, courno_loc)

#else

    ! Construct edge states as inputs to flux construction

    k_lo = g_lo
    k_hi = g_hi
    call prepare_for_fluxes(k_lo, k_hi, dt, dx, courno, &
                            q, flatn, q_lo, q_hi, &
                            div, g_lo, g_hi, &
                            qaux, qa_lo, qa_hi, &
                            sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                            qm, qp, It_lo, It_hi)

    call prepare_profile(k_lo, k_hi, &
                         q, flatn, q_lo, q_hi, &
                         sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                         qm, qp, It_lo, It_hi)

    ! Compute F^x

    idir = 1
    k_lo = f1_lo
    k_hi = f1_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, g_lo, g_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, It_lo, It_hi, &
                        flux1, q1, f1_lo, f1_hi, &
                        area1, a1_lo, a1_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^y

    idir = 2
    k_lo = f2_lo
    k_hi = f2_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, g_lo, g_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, It_lo, It_hi, &
                        flux2, q2, f2_lo, f2_hi, &
                        area2, a2_lo, a2_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^z

    idir = 3
    k_lo = f3_lo
    k_hi = f3_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, g_lo, g_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, It_lo, It_hi, &
                        flux3, q3, f3_lo, f3_hi, &
                        area3, a3_lo, a3_hi, &
                        qaux, qa_lo, qa_hi)

    ! Create an update source term based on the flux divergence.

    k_lo = lo
    k_hi = hi
    call construct_hydro_update(k_lo, k_hi, dx, dt, &
                                flux1, q1, f1_lo, f1_hi, &
                                flux2, q2, f2_lo, f2_hi, &
                                flux3, q3, f3_lo, f3_hi, &
                                area1, a1_lo, a1_hi, &
                                area2, a2_lo, a2_hi, &
                                area3, a3_lo, a3_hi, &
                                vol, vol_lo, vol_hi, &
                                update, updt_lo, updt_hi)

#endif

    call bl_deallocate(q1)
    call bl_deallocate(q2)
    call bl_deallocate(q3)

    call bl_deallocate(flatn)
    call bl_deallocate(div)

    call bl_deallocate(qm)
    call bl_deallocate(qp)

    call bl_deallocate(sxm)
    call bl_deallocate(sxp)
    call bl_deallocate(sym)
    call bl_deallocate(syp)
    call bl_deallocate(szm)
    call bl_deallocate(szp)

#ifdef CUDA

    call bl_deallocate(time_d)
    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(domlo_d)
    call bl_deallocate(domhi_d)
    call bl_deallocate(uin_lo_d)
    call bl_deallocate(uin_hi_d)
    call bl_deallocate(uout_lo_d)
    call bl_deallocate(uout_hi_d)
    call bl_deallocate(q_lo_d)
    call bl_deallocate(q_hi_d)
    call bl_deallocate(qa_lo_d)
    call bl_deallocate(qa_hi_d)
    call bl_deallocate(updt_lo_d)
    call bl_deallocate(updt_hi_d)
    call bl_deallocate(dx_d)
    call bl_deallocate(dt_d)
    call bl_deallocate(f1_lo_d)
    call bl_deallocate(f1_hi_d)
    call bl_deallocate(f2_lo_d)
    call bl_deallocate(f2_hi_d)
    call bl_deallocate(f3_lo_d)
    call bl_deallocate(f3_hi_d)
    call bl_deallocate(a1_lo_d)
    call bl_deallocate(a1_hi_d)
    call bl_deallocate(a2_lo_d)
    call bl_deallocate(a2_hi_d)
    call bl_deallocate(a3_lo_d)
    call bl_deallocate(a3_hi_d)
    call bl_deallocate(vol_lo_d)
    call bl_deallocate(vol_hi_d)
    call bl_deallocate(k_lo_d)
    call bl_deallocate(k_hi_d)
    call bl_deallocate(st_lo_d)
    call bl_deallocate(st_hi_d)
    call bl_deallocate(It_lo_d)
    call bl_deallocate(It_hi_d)
    call bl_deallocate(g_lo_d)
    call bl_deallocate(g_hi_d)
    call bl_deallocate(courno_d)
    call bl_deallocate(verbose_d)
    call bl_deallocate(idir_d)

#endif

  end subroutine ca_mol_single_stage



  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx, &
                        vol,v_lo,v_hi,mass,idx) bind(c, name='ca_summass')

    use bl_constants_module, only: ZERO
    use amrex_fort_module, only: rt => amrex_real
    use castro_util_module, only: summass
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_summass
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), idx
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(inout) :: mass

#ifdef CUDA

    attributes(device) :: rho, vol

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    real(rt)         :: mass_loc
    real(rt), managed, pointer :: mass_d(:)
    integer,  managed, pointer :: lo_d(:), hi_d(:)
    integer,  managed, pointer :: r_lo_d(:), r_hi_d(:)
    integer,  managed, pointer :: v_lo_d(:), v_hi_d(:)
    real(rt), managed, pointer :: dx_d(:)

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call bl_allocate(mass_d, 1, 1)
    call bl_allocate(lo_d, 1, 3)
    call bl_allocate(hi_d, 1, 3)
    call bl_allocate(r_lo_d, 1, 3)
    call bl_allocate(r_hi_d, 1, 3)
    call bl_allocate(v_lo_d, 1, 3)
    call bl_allocate(v_hi_d, 1, 3)
    call bl_allocate(dx_d, 1, 3)

    mass_loc = ZERO

    cuda_result = cudaMemcpyAsync(mass_d, mass_loc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(r_lo_d, r_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(r_hi_d, r_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(v_lo_d, v_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(v_hi_d, v_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_summass<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, rho, r_lo_d, r_hi_d, &
                                                            dx_d, vol, v_lo_d, v_hi_d, mass_d(1))

    cuda_result = cudaMemcpyAsync(mass_loc, mass_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize(stream)

    mass = mass + mass_loc

    call bl_deallocate(mass_d)
    call bl_deallocate(lo_d)
    call bl_deallocate(hi_d)
    call bl_deallocate(r_lo_d)
    call bl_deallocate(r_hi_d)
    call bl_deallocate(v_lo_d)
    call bl_deallocate(v_hi_d)
    call bl_deallocate(dx_d)

#else

    mass = ZERO

    call summass(lo,hi,rho,r_lo,r_hi,dx, &
                 vol,v_lo,v_hi,mass)

#endif

  end subroutine ca_summass



  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real
    use bc_fill_module, only: hypfill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer :: blo(3), bhi(3)
    integer :: adv_lo(3), adv_hi(3)

    adv_lo(1) = adv_l1
    adv_lo(2) = adv_l2
    adv_lo(3) = adv_l3
    adv_hi(1) = adv_h1
    adv_hi(2) = adv_h2
    adv_hi(3) = adv_h3

    blo = adv_lo
    bhi = adv_hi

    call hypfill(blo, bhi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real
    use bc_fill_module, only: denfill

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer :: blo(3), bhi(3)
    integer :: adv_lo(3), adv_hi(3)

    adv_lo(1) = adv_l1
    adv_lo(2) = adv_l2
    adv_lo(3) = adv_l3
    adv_hi(1) = adv_h1
    adv_hi(2) = adv_h2
    adv_hi(3) = adv_h3

    blo = adv_lo
    bhi = adv_hi

    call denfill(blo, bhi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_denfill

end module c_interface_modules
