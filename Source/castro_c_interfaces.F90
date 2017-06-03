module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR
  use amrex_fort_module, only: rt => amrex_real

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

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: s_lo_d(3), s_hi_d(3)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_enforce_consistent_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

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

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: s_lo_d(3), s_hi_d(3)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_compute_temp<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

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

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: u_lo_d(3), u_hi_d(3)
    integer, device :: verbose_d

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(verbose_d, verbose, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_reset_internal_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, u, u_lo_d, u_hi_d, verbose_d)

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

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: u_lo_d(3), u_hi_d(3)

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_normalize_species<<<numBlocks, numThreads, 0, stream>>>(u, u_lo_d, u_hi_d, lo_d, hi_d)

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

    integer,  device :: lo_d(3), hi_d(3)
    integer,  device :: uin_lo_d(3), uin_hi_d(3)
    integer,  device :: uout_lo_d(3), uout_hi_d(3)
    integer,  device :: vol_lo_d(3), vol_hi_d(3)
    real(rt), device :: frac_change_d
    integer,  device :: verbose_d

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    real(rt) :: local_frac_change

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

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
                                                                            lo_d, hi_d, frac_change_d, verbose_d)

    cuda_result = cudaMemcpyAsync(local_frac_change, frac_change_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize(stream)

    frac_change = min(frac_change, local_frac_change)

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

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    integer, intent(in)     :: idx

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

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: uin_lo_d(3), uin_hi_d(3)
    integer, device :: q_lo_d(3), q_hi_d(3)
    integer, device :: qa_lo_d(3), qa_hi_d(3)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

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

    integer,  device :: lo_d(3), hi_d(3)
    integer,  device :: nv_d, nc_d
    integer,  device :: v_lo_d(3), v_hi_d(3)
    integer,  device :: d_lo_d(3), d_hi_d(3)
    integer,  device :: domlo_d(3), domhi_d(3)
    integer,  device :: bc_d(3,2,nc)
    real(rt), device :: delta_d(3), xlo_d(3), time_d, dt_d
    integer,  device :: level_d, grid_no_d

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    ! Note that this stream calculation is not ideal because there are
    ! potentially multiple tiles per box.

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(nv_d, nv, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(v_lo_d, v_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(v_hi_d, v_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(d_lo_d, d_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(d_hi_d, d_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(delta_d, delta, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call dervel<<<numBlocks, numThreads, 0, stream>>>(vel,v_lo_d,v_hi_d,nv_d, &
                                                      dat,d_lo_d,d_hi_d,nc_d, &
                                                      lo_d,hi_d,domlo_d,domhi_d, &
                                                      delta_d,xlo_d,time_d,dt_d, &
                                                      bc_d,level_d,grid_no_d)

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

    integer,  device :: lo_d(3), hi_d(3)
    integer,  device :: np_d, nc_d
    integer,  device :: p_lo_d(3), p_hi_d(3)
    integer,  device :: u_lo_d(3), u_hi_d(3)
    integer,  device :: domlo_d(3), domhi_d(3)
    integer,  device :: bc_d(3,2,nc)
    real(rt), device :: dx_d(3), xlo_d(3), time_d, dt_d
    integer,  device :: level_d, grid_no_d

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(np_d, np, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(p_lo_d, p_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(p_hi_d, p_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call derpres<<<numBlocks, numThreads, 0, stream>>>(p,p_lo_d,p_hi_d,np_d, &
                                                       u,u_lo_d,u_hi_d,nc_d, &
                                                       lo_d,hi_d,domlo_d,domhi_d, &
                                                       dx_d,xlo_d,time_d,dt_d, &
                                                       bc_d,level_d,grid_no_d)

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

    integer,  device :: lo_d(3), hi_d(3)
    integer,  device :: u_lo_d(3), u_hi_d(3)
    real(rt), device :: dx_d(3)
    real(rt), device :: dt_loc_d

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    real(rt) :: dt_loc

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(u_lo_d, u_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(u_hi_d, u_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)

    dt_loc = dt

    cuda_result = cudaMemcpyAsync(dt_loc_d, dt, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_estdt<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, u, u_lo_d, u_hi_d, dx_d, dt_loc_d)

    cuda_result = cudaMemcpyAsync(dt_loc, dt_loc_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize(stream)

    dt = min(dt, dt_loc)

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
                                 flux1, flux1_lo, flux1_hi, &
                                 flux2, flux2_lo, flux2_hi, &
                                 flux3, flux3_lo, flux3_hi, &
                                 area1, area1_lo, area1_hi, &
                                 area2, area2_lo, area2_hi, &
                                 area3, area3_lo, area3_hi, &
                                 vol, vol_lo, vol_hi, &
                                 courno, verbose, idx) bind(C, name="ca_mol_single_stage")

    use mol_module, only: mol_single_stage
    use advection_util_module, only: ht, allocate_ht, deallocate_ht
#ifdef CUDA
    use cuda_interfaces_module, only: cuda_mol_single_stage
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose, idx
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: updt_lo(3), updt_hi(3)
    integer,  intent(in   ) :: flux1_lo(3), flux1_hi(3)
    integer,  intent(in   ) :: flux2_lo(3), flux2_hi(3)
    integer,  intent(in   ) :: flux3_lo(3), flux3_hi(3)
    integer,  intent(in   ) :: area1_lo(3), area1_hi(3)
    integer,  intent(in   ) :: area2_lo(3), area2_hi(3)
    integer,  intent(in   ) :: area3_lo(3), area3_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
    real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
    real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
    real(rt), intent(in   ) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
    real(rt), intent(in   ) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
    real(rt), intent(in   ) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in   ) :: dx(3), dt, time
    real(rt), intent(inout) :: courno

    integer :: ngf
    integer :: It_lo(3), It_hi(3)
    integer :: st_lo(3), st_hi(3)
    integer :: g_lo(3), g_hi(3)
    integer :: gd_lo(3), gd_hi(3)

    type(ht) :: h

    ngf = 1

    It_lo = lo - 1
    It_hi = hi + 1

    st_lo = lo - 2
    st_hi = hi + 2

    gd_lo = lo
    gd_hi = hi + 1

    g_lo = lo - ngf
    g_hi = hi + ngf

    ! Allocate all the temporaries we will need.
    call allocate_ht(h, lo, hi, flux1_lo, flux1_hi, flux2_lo, flux2_hi, &
                     flux3_lo, flux3_hi, st_lo, st_hi, It_lo, It_hi, &
                     g_lo, g_hi, gd_lo, gd_hi, q_lo, q_hi)

#ifdef CUDA

    attributes(device) :: uin, uout, q, qaux, update, flux1, flux2, flux3, area1, area2, area3, vol

    integer                   :: cuda_result
    integer(cuda_stream_kind) :: stream
    type(dim3)                :: numThreads, numBlocks

    real(rt), device :: time_d
    integer,  device :: lo_d(3), hi_d(3), domlo_d(3), domhi_d(3)
    integer,  device :: uin_lo_d(3), uin_hi_d(3)
    integer,  device :: uout_lo_d(3), uin_hi_d(3)
    integer,  device :: q_lo_d(3), q_hi_d(3)
    integer,  device :: qa_lo_d(3), qa_hi_d(3)
    integer,  device :: updt_lo_d(3), updt_hi_d(3)
    real(rt), device :: dx_d(3), dt_d
    integer,  device :: flux1_lo_d(3), flux1_hi_d(3)
    integer,  device :: flux2_lo_d(3), flux2_hi_d(3)
    integer,  device :: flux3_lo_d(3), flux3_hi_d(3)
    integer,  device :: area1_lo_d(3), area1_hi_d(3)
    integer,  device :: area2_lo_d(3), area2_hi_d(3)
    integer,  device :: area3_lo_d(3), area3_hi_d(3)
    integer,  device :: vol_lo_d(3), vol_hi_d(3)
    real(rt), device :: courno_d
    integer,  device :: verbose_d

    real(rt) :: courno_loc

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(domlo_d, domlo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(domhi_d, domhi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uin_lo_d, uin_lo, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(uout_lo_d, uout_lo, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(q_lo_d, q_lo, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(qa_lo_d, qa_lo, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(updt_lo_d, updt_lo, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(flux1_lo_d, flux1_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(flux1_hi_d, flux1_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(flux2_lo_d, flux2_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(flux2_hi_d, flux2_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(flux3_lo_d, flux3_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(flux3_hi_d, flux3_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(area1_lo_d, area1_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(area1_hi_d, area1_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(area2_lo_d, area2_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(area2_hi_d, area2_hi, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(area3_lo_d, area3_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(area3_hi_d, area3_hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(vol_lo_d, vol_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(vol_hi_d, vol_hi, 3, cudaMemcpyHostToDevice, stream)

    courno_loc = courno

    cuda_result = cudaMemcpyAsync(courno_d, courno_loc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(verbose_d, verbose, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_mol_single_stage<<<numBlocks, numThreads, 0, stream>>>(time_d, &
                                                                     lo_d, hi_d, domlo_d, domhi_d, &
                                                                     uin, uin_lo_d, uin_hi_d, &
                                                                     uout, uout_lo_d, uout_hi_d, &
                                                                     q, q_lo_d, q_hi_d, &
                                                                     qaux, qa_lo_d, qa_hi_d, &
                                                                     update, updt_lo_d, updt_hi_d, &
                                                                     dx_d, dt_d, h, &
                                                                     flux1, flux1_lo_d, flux1_hi_d, &
                                                                     flux2, flux2_lo_d, flux2_hi_d, &
                                                                     flux3, flux3_lo_d, flux3_hi_d, &
                                                                     area1, area1_lo_d, area1_hi_d, &
                                                                     area2, area2_lo_d, area2_hi_d, &
                                                                     area3, area3_lo_d, area3_hi_d, &
                                                                     vol, vol_lo_d, vol_hi_d, &
                                                                     courno_d, verbose_d)

    cuda_result = cudaMemcpyAsync(courno_loc, courno_d, 1, cudaMemcpyDeviceToHost, stream)

    cuda_result = cudaStreamSynchronize()

    courno = max(courno, courno_loc)

#else

    call mol_single_stage(time, &
                          lo, hi, domlo, domhi, &
                          uin, uin_lo, uin_hi, &
                          uout, uout_lo, uout_hi, &
                          q, q_lo, q_hi, &
                          qaux, qa_lo, qa_hi, &
                          update, updt_lo, updt_hi, &
                          dx, dt, h, &
                          flux1, flux1_lo, flux1_hi, &
                          flux2, flux2_lo, flux2_hi, &
                          flux3, flux3_lo, flux3_hi, &
                          area1, area1_lo, area1_hi, &
                          area2, area2_lo, area2_hi, &
                          area3, area3_lo, area3_hi, &
                          vol, vol_lo, vol_hi, &
                          courno, verbose)

#endif

    call deallocate_ht(h)

  end subroutine ca_mol_single_stage

end module c_interface_modules
