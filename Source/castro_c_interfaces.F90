module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
  use mempool_module, only: bl_allocate, bl_deallocate
#ifdef CUDA
  use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, &
                     cudaStreamSynchronize, cudaDeviceSynchronize, dim3, cuda_stream_kind
  use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams, stream_from_index, stream_index
#endif

contains

  subroutine ca_initdata(level, lo, hi, &
                         state, s_lo, s_hi, dx, &
                         xlo, xhi) bind(C, name="ca_initdata")

#ifdef CUDA
    use cuda_interfaces_module, only: cuda_initdata
#endif
    use initdata_module, only: initdata

    implicit none

    integer,  intent(in   ), value :: level
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: xlo(3), xhi(3), dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: state, lo, hi, s_lo, s_hi, dx, xlo, xhi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_initdata &
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
         (level, lo, hi, state, s_lo, s_hi, dx, xlo, xhi)

#else

    call initdata(level, lo, hi, state, s_lo(1), s_lo(2), s_lo(3), s_hi(1), s_hi(2), s_hi(3), dx, xlo, xhi)

#endif

  end subroutine ca_initdata


  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi) &
                                     bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: state, lo, hi, s_lo, s_hi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call enforce_consistent_e &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
#endif
         (lo, hi, state, s_lo, s_hi)

  end subroutine ca_enforce_consistent_e


  subroutine ca_compute_temp(lo, hi, state, s_lo, s_hi) &
                             bind(C, name="ca_compute_temp")

    use castro_util_module, only: compute_temp

    implicit none

    integer,  intent(in   ) :: lo(3),hi(3)
    integer,  intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: state, lo, hi, s_lo, s_hi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call compute_temp &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
#endif
         (lo, hi, state, s_lo, s_hi)

  end subroutine ca_compute_temp


  subroutine ca_reset_internal_e(lo, hi, u, u_lo, u_hi, verbose) &
                                 bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: verbose
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call reset_internal_e &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
#endif
         (lo, hi, u, u_lo, u_hi, verbose)

  end subroutine ca_reset_internal_e


  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi) &
                                  bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call normalize_species &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
#endif
         (u, u_lo, u_hi, lo, hi)

  end subroutine ca_normalize_species


  subroutine ca_enforce_minimum_density(uin, uin_lo, uin_hi, &
                                        uout, uout_lo, uout_hi, &
                                        vol, vol_lo, vol_hi, &
                                        lo, hi, frac_change, verbose) &
                                        bind(C, name="ca_enforce_minimum_density")

    use advection_util_module, only: enforce_minimum_density

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: verbose
    integer,  intent(in   ) ::  uin_lo(3),  uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) ::  vol_lo(3),  vol_hi(3)

    real(rt), intent(in   ) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in   ) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt), intent(inout) :: frac_change

#ifdef CUDA
    attributes(managed) :: uin, uout, vol, lo, hi, uin_lo, uin_hi, uout_lo, uout_hi, vol_lo, vol_hi

    real(rt), managed, pointer :: frac_change_d(:)

    integer :: cuda_result
    type(dim3) :: numThreads, numBlocks

    real(rt) :: local_frac_change

    call bl_allocate(frac_change_d, 1, 1)

    local_frac_change = frac_change

    cuda_result = cudaMemcpyAsync(frac_change_d, local_frac_change, 1, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call enforce_minimum_density &
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
         (uin, uin_lo, uin_hi, &
          uout, uout_lo, uout_hi, &
          vol, vol_lo, vol_hi, &
          lo, hi, frac_change_d(1), verbose)

    cuda_result = cudaMemcpyAsync(local_frac_change, frac_change_d, 1, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(stream_index)))

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    frac_change = min(frac_change, local_frac_change)

    call bl_deallocate(frac_change_d)

#else

    call enforce_minimum_density(uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 vol, vol_lo, vol_hi, &
                                 lo, hi, frac_change, verbose)

#endif

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

    call check_initial_species(lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species


  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi) bind(C, name = "ca_ctoprim")

    use advection_util_module, only: ctoprim

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: uin_lo(3), uin_hi(3)
    integer,  intent(in) :: q_lo(3), q_hi(3)
    integer,  intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

#ifdef CUDA
    attributes(managed) :: uin, q, qaux, lo, hi, uin_lo, uin_hi, q_lo, q_hi, qa_lo, qa_hi

    type(dim3) :: numThreads, numBlocks

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call ctoprim &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>> &
#endif
         (lo, hi, uin, uin_lo, uin_hi, q, q_lo, q_hi, qaux, qa_lo, qa_hi)

  end subroutine ca_ctoprim


  subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,dx,xlo,time,dt,bc,level,grid_no) &
                       bind(c, name="ca_dervel")

    use castro_util_module, only: dervel

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3), nv
    integer,  intent(in   ) :: d_lo(3), d_hi(3), nc
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(3,2,nc)
    real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ) :: level, grid_no

#ifdef CUDA
    attributes(managed) :: vel, dat, lo, hi, v_lo, v_hi, d_lo, d_hi, dx, xlo, domlo, domhi

    integer,  managed, pointer :: nv_d(:), nc_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)
    real(rt), managed, pointer :: time_d(:), dt_d(:)
    integer,  managed, pointer :: level_d(:), grid_no_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    ! Note that this stream calculation is not ideal because there are
    ! potentially multiple tiles per box.

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    call bl_allocate(nv_d, 1, 1)
    call bl_allocate(nc_d, 1, 1)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, nc)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(dt_d, 1, 1)
    call bl_allocate(level_d, 1, 1)
    call bl_allocate(grid_no_d, 1, 1)

    cuda_result = cudaMemcpyAsync(nv_d, nv, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call dervel<<<numBlocks, numThreads, 0, stream>>>(vel,v_lo,v_hi,nv_d(1), &
                                                      dat,d_lo,d_hi,nc_d(1), &
                                                      lo,hi,domlo,domhi, &
                                                      dx,xlo,time_d(1),dt_d(1), &
                                                      bc_d,level_d(1),grid_no_d(1))

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(nv_d)
    call bl_deallocate(nc_d)
    call bl_deallocate(bc_d)
    call bl_deallocate(time_d)
    call bl_deallocate(dt_d)
    call bl_deallocate(level_d)
    call bl_deallocate(grid_no_d)

#else

    call dervel(vel,v_lo,v_hi,nv, &
                dat,d_lo,d_hi,nc,lo,hi,domlo, &
                domhi,dx,xlo,time,dt,bc,level,grid_no)

#endif

  end subroutine ca_dervel



  subroutine ca_derpres(p,p_lo,p_hi,np, &
                        u,u_lo,u_hi,nc,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(c, name="ca_derpres")

    use castro_util_module, only: derpres

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
    attributes(managed) :: p, u, lo, hi, p_lo, p_hi, u_lo, u_hi, dx, xlo, domlo, domhi

    integer,  managed, pointer :: np_d(:), nc_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)
    real(rt), managed, pointer :: time_d(:), dt_d(:)
    integer,  managed, pointer :: level_d(:), grid_no_d(:)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(grid_no, max_cuda_streams) + 1)

    call bl_allocate(np_d, 1, 1)
    call bl_allocate(nc_d, 1, 1)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, nc)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(dt_d, 1, 1)
    call bl_allocate(level_d, 1, 1)
    call bl_allocate(grid_no_d, 1, 1)

    cuda_result = cudaMemcpyAsync(np_d, np, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(nc_d, nc, 1, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6 * nc, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(dt_d, dt, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(level_d, level, 1, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(grid_no_d, grid_no, 1, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call derpres<<<numBlocks, numThreads, 0, stream>>>(p,p_lo,p_hi,np_d(1), &
                                                       u,u_lo,u_hi,nc_d(1), &
                                                       lo,hi,domlo,domhi, &
                                                       dx,xlo,time_d(1),dt_d(1), &
                                                       bc_d,level_d(1),grid_no_d(1))

    cuda_result = cudaStreamSynchronize(stream)

    call bl_deallocate(np_d)
    call bl_deallocate(nc_d)
    call bl_deallocate(bc_d)
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

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="ca_estdt")

    use timestep_module, only: estdt

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi, dx

    real(rt), managed, pointer :: dt_loc_d(:)

    integer                   :: cuda_result
    type(dim3)                :: numThreads, numBlocks

    real(rt) :: dt_loc

    call bl_allocate(dt_loc_d, 1, 1)

    dt_loc = dt

    cuda_result = cudaMemcpyAsync(dt_loc_d, dt, 1, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call estdt<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(lo, hi, u, u_lo, u_hi, dx, dt_loc_d(1))

    cuda_result = cudaMemcpyAsync(dt_loc, dt_loc_d, 1, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(stream_index)))

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    dt = min(dt, dt_loc)

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
                                 flatn, flatn_lo, flatn_hi, &
                                 div, div_lo, div_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 update, updt_lo, updt_hi, &
                                 dx, dt, &
                                 q1, q1_lo, q1_hi, &
                                 q2, q2_lo, q2_hi, &
                                 q3, q3_lo, q3_hi, &
                                 qm, qm_lo, qp_hi, &
                                 qp, qp_lo, qm_hi, &
                                 sxm, sxm_lo, sxm_hi, &
                                 sxp, sxp_lo, sxp_hi, &
                                 sym, sym_lo, sym_hi, &
                                 syp, syp_lo, syp_hi, &
                                 szm, szm_lo, szm_hi, &
                                 szp, szp_lo, szp_hi, &
                                 flux1, f1_lo, f1_hi, &
                                 flux2, f2_lo, f2_hi, &
                                 flux3, f3_lo, f3_hi, &
                                 area1, a1_lo, a1_hi, &
                                 area2, a2_lo, a2_hi, &
                                 area3, a3_lo, a3_hi, &
                                 vol, vol_lo, vol_hi, &
                                 verbose) bind(C, name="ca_mol_single_stage")

    use mol_module, only: prepare_for_fluxes, prepare_profile, construct_flux
    use advection_util_module, only: construct_hydro_update

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: verbose
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: flatn_lo(3), flatn_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: updt_lo(3), updt_hi(3)
    integer,  intent(in   ) :: q1_lo(3), q1_hi(3)
    integer,  intent(in   ) :: q2_lo(3), q2_hi(3)
    integer,  intent(in   ) :: q3_lo(3), q3_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: sxm_lo(3), sxm_hi(3)
    integer,  intent(in   ) :: sxp_lo(3), sxp_hi(3)
    integer,  intent(in   ) :: sym_lo(3), sym_hi(3)
    integer,  intent(in   ) :: syp_lo(3), syp_hi(3)
    integer,  intent(in   ) :: szm_lo(3), szm_hi(3)
    integer,  intent(in   ) :: szp_lo(3), szp_hi(3)
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
    real(rt), intent(inout) :: flatn(flatn_lo(1):flatn_hi(1), flatn_lo(2):flatn_hi(2), flatn_lo(3):flatn_hi(3), NQ)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(inout) :: q1(q1_lo(1):q1_hi(1), q1_lo(2):q1_hi(2), q1_lo(3):q1_hi(3), NGDNV)
    real(rt), intent(inout) :: q2(q2_lo(1):q2_hi(1), q2_lo(2):q2_hi(2), q2_lo(3):q2_hi(3), NGDNV)
    real(rt), intent(inout) :: q3(q3_lo(1):q3_hi(1), q3_lo(2):q3_hi(2), q3_lo(3):q3_hi(3), NGDNV)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ)
    real(rt), intent(inout) :: sxm(sxm_lo(1):sxm_hi(1), sxm_lo(2):sxm_hi(2), sxm_lo(3):sxm_hi(3), NQ)
    real(rt), intent(inout) :: sxp(sxp_lo(1):sxp_hi(1), sxp_lo(2):sxp_hi(2), sxp_lo(3):sxp_hi(3), NQ)
    real(rt), intent(inout) :: sym(sym_lo(1):sym_hi(1), sym_lo(2):sym_hi(2), sym_lo(3):sym_hi(3), NQ)
    real(rt), intent(inout) :: syp(syp_lo(1):syp_hi(1), syp_lo(2):syp_hi(2), syp_lo(3):syp_hi(3), NQ)
    real(rt), intent(inout) :: szm(szm_lo(1):szm_hi(1), szm_lo(2):szm_hi(2), szm_lo(3):szm_hi(3), NQ)
    real(rt), intent(inout) :: szp(szp_lo(1):szp_hi(1), szp_lo(2):szp_hi(2), szp_lo(3):szp_hi(3), NQ)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1), f1_lo(2):f1_hi(2), f1_lo(3):f1_hi(3), NVAR)
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1), f2_lo(2):f2_hi(2), f2_lo(3):f2_hi(3), NVAR)
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1), f3_lo(2):f3_hi(2), f3_lo(3):f3_hi(3), NVAR)
    real(rt), intent(in   ) :: area1(a1_lo(1):a1_hi(1), a1_lo(2):a1_hi(2), a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: area2(a2_lo(1):a2_hi(1), a2_lo(2):a2_hi(2), a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: area3(a3_lo(1):a3_hi(1), a3_lo(2):a3_hi(2), a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time

    integer :: k_lo(3), k_hi(3)
    integer :: idir

#ifdef CUDA
    attributes(managed) :: uin, uout, q, qaux, update, flux1, flux2, flux3, area1, area2, area3, vol, &
                           lo, hi, uin_lo, uin_hi, uout_lo, uout_hi, q_lo, q_hi, qa_lo, qa_hi, &
                           updt_lo, updt_hi, dx, f1_lo, f1_hi, f2_lo, f2_hi, f3_lo, f3_hi, &
                           a1_lo, a1_hi, a2_lo, a2_hi, a3_lo, a3_hi, vol_lo, vol_hi, &
                           q1, q1_lo, q1_hi, q2, q2_lo, q2_hi, q3, q3_lo, q3_hi, &
                           qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, domlo, domhi, &
                           sxm, sxm_lo, sxm_hi, sxp, sxp_lo, sxp_hi, &
                           sym, sym_lo, sym_hi, syp, syp_lo, syp_hi, &
                           szm, szm_lo, szm_hi, szp, szp_lo, szp_hi, &
                           flatn, flatn_lo, flatn_hi, div, div_lo, div_hi

    integer                   :: cuda_result
    type(dim3)                :: numThreads, numBlocks

    integer,  managed, pointer :: k_lo_d(:), k_hi_d(:)
#endif

#ifdef CUDA

    call bl_allocate(k_lo_d, 1, 3)
    call bl_allocate(k_hi_d, 1, 3)

    ! Construct edge states as inputs to flux construction

    k_lo = div_lo
    k_hi = div_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call prepare_for_fluxes<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, dt, dx, &
                                                                  q, flatn, q_lo, q_hi, &
                                                                  div, div_lo, div_hi, &
                                                                  qaux, qa_lo, qa_hi, &
                                                                  sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi, &
                                                                  qm, qp, qm_lo, qm_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    call prepare_profile<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, &
                                                               q, flatn, q_lo, q_hi, &
                                                               sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi, &
                                                               qm, qp, qm_lo, qm_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    ! Compute F^x

    idir = 1

    k_lo = f1_lo
    k_hi = f1_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call construct_flux<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, domlo, domhi, dx, dt, idir, &
                                                              div, div_lo, div_hi, &
                                                              uin, uin_lo, uin_hi, &
                                                              qm, qp, qm_lo, qm_hi, &
                                                              flux1, q1, f1_lo, f1_hi, &
                                                              area1, a1_lo, a1_hi, &
                                                              qaux, qa_lo, qa_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    ! Compute F^y

    idir = 2

    k_lo = f2_lo
    k_hi = f2_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call construct_flux<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, domlo, domhi, dx, dt, idir, &
                                                              div, div_lo, div_hi, &
                                                              uin, uin_lo, uin_hi, &
                                                              qm, qp, qm_lo, qm_hi, &
                                                              flux2, q2, f2_lo, f2_hi, &
                                                              area2, a2_lo, a2_hi, &
                                                              qaux, qa_lo, qa_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    ! Compute F^z

    idir = 3

    k_lo = f3_lo
    k_hi = f3_hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call construct_flux<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, domlo, domhi, dx, dt, idir, &
                                                              div, div_lo, div_hi, &
                                                              uin, uin_lo, uin_hi, &
                                                              qm, qp, qm_lo, qm_hi, &
                                                              flux3, q3, f3_lo, f3_hi, &
                                                              area3, a3_lo, a3_hi, &
                                                              qaux, qa_lo, qa_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    ! Create an update source term based on the flux divergence.

    k_lo = lo
    k_hi = hi

    cuda_result = cudaMemcpyAsync(k_lo_d, k_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))
    cuda_result = cudaMemcpyAsync(k_hi_d, k_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(k_lo, k_hi, numBlocks, numThreads)

    call construct_hydro_update<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(k_lo_d, k_hi_d, dx, dt, &
                                                                      flux1, q1, f1_lo, f1_hi, &
                                                                      flux2, q2, f2_lo, f2_hi, &
                                                                      flux3, q3, f3_lo, f3_hi, &
                                                                      area1, a1_lo, a1_hi, &
                                                                      area2, a2_lo, a2_hi, &
                                                                      area3, a3_lo, a3_hi, &
                                                                      vol, vol_lo, vol_hi, &
                                                                      update, updt_lo, updt_hi)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

#else

    ! Construct edge states as inputs to flux construction

    k_lo = div_lo
    k_hi = div_hi
    call prepare_for_fluxes(k_lo, k_hi, dt, dx, &
                            q, flatn, q_lo, q_hi, &
                            div, div_lo, div_hi, &
                            qaux, qa_lo, qa_hi, &
                            sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi, &
                            qm, qp, qm_lo, qm_hi)

    call prepare_profile(k_lo, k_hi, &
                         q, flatn, q_lo, q_hi, &
                         sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi, &
                         qm, qp, qm_lo, qm_hi)

    ! Compute F^x

    idir = 1
    k_lo = f1_lo
    k_hi = f1_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, div_lo, div_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, qm_lo, qm_hi, &
                        flux1, q1, f1_lo, f1_hi, &
                        area1, a1_lo, a1_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^y

    idir = 2
    k_lo = f2_lo
    k_hi = f2_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, div_lo, div_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, qm_lo, qm_hi, &
                        flux2, q2, f2_lo, f2_hi, &
                        area2, a2_lo, a2_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^z

    idir = 3
    k_lo = f3_lo
    k_hi = f3_hi
    call construct_flux(k_lo, k_hi, domlo, domhi, dx, dt, idir, &
                        div, div_lo, div_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, qm_lo, qm_hi, &
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

#ifdef CUDA

    call bl_deallocate(k_lo_d)
    call bl_deallocate(k_hi_d)

#endif

  end subroutine ca_mol_single_stage



  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx, &
                        vol,v_lo,v_hi,mass) bind(c, name='ca_summass')

    use bl_constants_module, only: ZERO
    use amrex_fort_module, only: rt => amrex_real
    use castro_util_module, only: summass

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(inout) :: mass

#ifdef CUDA
    attributes(managed) :: rho, vol, lo, hi, r_lo, r_hi, v_lo, v_hi, dx

    integer                   :: cuda_result
    type(dim3)                :: numThreads, numBlocks

    real(rt)         :: mass_loc
    real(rt), managed, pointer :: mass_d(:)

    call bl_allocate(mass_d, 1, 1)

    mass_loc = ZERO

    cuda_result = cudaMemcpyAsync(mass_d, mass_loc, 1, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(stream_index)))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call summass<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(stream_index))>>>(lo, hi, rho, r_lo, r_hi, &
                                                       dx, vol, v_lo, v_hi, mass_d(1))

    cuda_result = cudaMemcpyAsync(mass_loc, mass_d, 1, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(stream_index)))

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(stream_index)))

    mass = mass + mass_loc

    call bl_deallocate(mass_d)

#else

    mass = ZERO

    call summass(lo,hi,rho,r_lo,r_hi,dx, &
                 vol,v_lo,v_hi,mass)

#endif

  end subroutine ca_summass



  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_hypfill")

    use bc_fill_module, only: hypfill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

#ifdef CUDA

    attributes(device) :: adv, domlo, domhi

    integer                   :: cuda_result
    type(dim3)                :: numThreads, numBlocks

    integer,  managed, pointer :: adv_lo_d(:), adv_hi_d(:)
    real(rt), managed, pointer :: dx_d(:), xlo_d(:), time_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)

#endif

    integer :: adv_lo(3), adv_hi(3)

    adv_lo(1) = adv_l1
    adv_lo(2) = adv_l2
    adv_lo(3) = adv_l3
    adv_hi(1) = adv_h1
    adv_hi(2) = adv_h2
    adv_hi(3) = adv_h3

#ifdef CUDA

    call bl_allocate(adv_lo_d, 1, 3)
    call bl_allocate(adv_hi_d, 1, 3)
    call bl_allocate(dx_d, 1, 3)
    call bl_allocate(xlo_d, 1, 3)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, NVAR)

    cuda_result = cudaMemcpyAsync(adv_lo_d, adv_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))
    cuda_result = cudaMemcpyAsync(adv_hi_d, adv_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6*NVAR, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    call threads_and_blocks(adv_lo, adv_hi, numBlocks, numThreads)

    call hypfill<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(0))>>>(adv_lo_d, adv_hi_d, adv, adv_lo_d, adv_hi_d, domlo, domhi, &
                                                       dx_d, xlo_d, time_d(1), bc_d)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(0)))

    call bl_deallocate(adv_lo_d)
    call bl_deallocate(adv_hi_d)
    call bl_deallocate(dx_d)
    call bl_deallocate(xlo_d)
    call bl_deallocate(time_d)
    call bl_deallocate(bc_d)

#else

    call hypfill(adv_lo, adv_hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc)

#endif

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real
    use bc_fill_module, only: denfill

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
#ifdef CUDA

    attributes(device) :: adv, domlo, domhi

    integer                   :: cuda_result
    type(dim3)                :: numThreads, numBlocks

    integer,  managed, pointer :: adv_lo_d(:), adv_hi_d(:)
    real(rt), managed, pointer :: dx_d(:), xlo_d(:), time_d(:)
    integer,  managed, pointer :: bc_d(:,:,:)

#endif

    integer :: adv_lo(3), adv_hi(3)

    adv_lo(1) = adv_l1
    adv_lo(2) = adv_l2
    adv_lo(3) = adv_l3
    adv_hi(1) = adv_h1
    adv_hi(2) = adv_h2
    adv_hi(3) = adv_h3

#ifdef CUDA

    call bl_allocate(adv_lo_d, 1, 3)
    call bl_allocate(adv_hi_d, 1, 3)
    call bl_allocate(dx_d, 1, 3)
    call bl_allocate(xlo_d, 1, 3)
    call bl_allocate(time_d, 1, 1)
    call bl_allocate(bc_d, 1, 3, 1, 2, 1, NVAR)

    cuda_result = cudaMemcpyAsync(adv_lo_d, adv_lo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))
    cuda_result = cudaMemcpyAsync(adv_hi_d, adv_hi, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(dx_d, dx, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))
    cuda_result = cudaMemcpyAsync(xlo_d, xlo, 3, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(time_d, time, 1, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    cuda_result = cudaMemcpyAsync(bc_d, bc, 6*NVAR, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(0)))

    call threads_and_blocks(adv_lo, adv_hi, numBlocks, numThreads)

    call denfill<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(0))>>>(adv_lo_d, adv_hi_d, adv, adv_lo_d, adv_hi_d, domlo, domhi, &
                                                       dx_d, xlo_d, time_d(1), bc_d)

    cuda_result = cudaStreamSynchronize(cuda_streams(stream_from_index(0)))

    call bl_deallocate(adv_lo_d)
    call bl_deallocate(adv_hi_d)
    call bl_deallocate(dx_d)
    call bl_deallocate(xlo_d)
    call bl_deallocate(time_d)
    call bl_deallocate(bc_d)

#else

    call denfill(adv_lo, adv_hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc)

#endif

  end subroutine ca_denfill

end module c_interface_modules
