module cuda_interfaces_module

  implicit none

contains

  attributes(global) &
  subroutine cuda_enforce_consistent_e(lo,hi,state,s_lo,s_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call enforce_consistent_e(idx, idx, state, s_lo, s_hi)

  end subroutine cuda_enforce_consistent_e

  

  attributes(global) &
  subroutine cuda_initdata(level, time, lo, hi, ns, &
                           state, s_lo, s_hi, xlo, xhi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: level, ns
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call initdata(level, time, idx, idx, ns, state, s_lo, s_hi, xlo, xhi)

  end subroutine cuda_initdata
  


  attributes(global) &
  subroutine cuda_enforce_minimum_density(uin, uin_lo, uin_hi, &
                                          uout, uout_lo, uout_hi, &
                                          vol, vol_lo, vol_hi, &
                                          lo, hi, frac_change, verbose)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use advection_util_module, only: enforce_minimum_density

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: frac_change
    integer,  intent(in   ) :: verbose

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call enforce_minimum_density(uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 vol, vol_lo, vol_hi, &
                                 idx, idx, frac_change, verbose)

  end subroutine cuda_enforce_minimum_density



  attributes(global) &
  subroutine cuda_normalize_species(u, u_lo, u_hi, lo, hi)

    use amrex_fort_module, only: rt => amrex_real
    use castro_util_module, only: normalize_species
    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call normalize_species(u, u_lo, u_hi, idx, idx)

  end subroutine cuda_normalize_species



  attributes(global) &
  subroutine cuda_ctoprim(lo, hi, &
                          uin, uin_lo, uin_hi, &
                          q,     q_lo,   q_hi, &
                          qaux, qa_lo,  qa_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NQ, NQAUX
    use advection_util_module, only: ctoprim

    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(in   ) :: uin_lo(3), uin_hi(3)
    integer, intent(in   ) :: q_lo(3), q_hi(3)
    integer, intent(in   ) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call ctoprim(idx, idx, &
                 uin, uin_lo, uin_hi, &
                 q,     q_lo,   q_hi, &
                 qaux, qa_lo,  qa_hi)

  end subroutine cuda_ctoprim



  attributes(global) &
  subroutine cuda_reset_internal_e(lo,hi,u,u_lo,u_hi,verbose)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use castro_util_module, only: reset_internal_e

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call reset_internal_e(idx, idx, u, u_lo, u_hi, verbose)

  end subroutine cuda_reset_internal_e



  attributes(global) &
  subroutine cuda_dervel(vel,v_lo,v_hi,nv, &
                         dat,d_lo,d_hi,nc,lo,hi,domlo, &
                         domhi,delta,xlo,time,dt,bc,level,grid_no)

    use amrex_fort_module, only: rt => amrex_real
    use castro_util_module, only: dervel

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

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call dervel(vel,v_lo,v_hi,nv, &
                dat,d_lo,d_hi,nc,idx,idx,domlo, &
                domhi,delta,xlo,time,dt,bc,level,grid_no)

  end subroutine cuda_dervel



  attributes(global) &
  subroutine cuda_derpres(p,p_lo,p_hi,ncomp_p, &
                          u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no)

    use amrex_fort_module, only: rt => amrex_real
    use castro_util_module, only: derpres

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3), ncomp_p
    integer,  intent(in   ) :: u_lo(3), u_hi(3), ncomp_u
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
    integer,  intent(in   ) :: bc(3,2,ncomp_u), level, grid_no

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call derpres(p,p_lo,p_hi,ncomp_p, &
                 u,u_lo,u_hi,ncomp_u,idx,idx,domlo, &
                 domhi,dx,xlo,time,dt,bc,level,grid_no)

  end subroutine cuda_derpres



  attributes(global) &
  subroutine cuda_compute_temp(lo,hi,state,s_lo,s_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use castro_util_module, only: compute_temp

    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call compute_temp(idx, idx, state, s_lo, s_hi)

  end subroutine cuda_compute_temp



  attributes(global) &
  subroutine cuda_estdt(lo,hi,u,u_lo,u_hi,dx,dt)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use timestep_module, only: estdt

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call estdt(idx,idx,u,u_lo,u_hi,dx,dt)

  end subroutine cuda_estdt



  attributes(global) &
  subroutine cuda_prepare_for_fluxes(lo, hi, dt, dx, courno, &
                                     q, flatn, q_lo, q_hi, &
                                     div, g_lo, g_hi, &
                                     qaux, qa_lo, qa_hi, &
                                     sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                                     qm, qp, It_lo, It_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use mol_module, only: prepare_for_fluxes

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: g_lo(3), g_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: st_lo(3), st_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: flatn(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3))
    real(rt), intent(inout) :: div(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: sxm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sxp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sym(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: syp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(in   ) :: dx(3), dt
    real(rt), intent(inout) :: courno

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call prepare_for_fluxes(idx, idx, dt, dx, courno, &
                            q, flatn, q_lo, q_hi, &
                            div, g_lo, g_hi, &
                            qaux, qa_lo, qa_hi, &
                            sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                            qm, qp, It_lo, It_hi)

  end subroutine cuda_prepare_for_fluxes



  attributes(global) &
  subroutine cuda_prepare_profile(lo, hi, &
                                  q, flatn, q_lo, q_hi, &
                                  sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                                  qm, qp, It_lo, It_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use mol_module, only: prepare_profile

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: st_lo(3), st_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: flatn(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3))
    real(rt), intent(inout) :: sxm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sxp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sym(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: syp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call prepare_profile(idx, idx, &
                         q, flatn, q_lo, q_hi, &
                         sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                         qm, qp, It_lo, It_hi)

  end subroutine cuda_prepare_profile



  attributes(global) &
  subroutine cuda_construct_flux(lo, hi, domlo, domhi, dx, dt, idir, &
                                 div, g_lo, g_hi, &
                                 uin, uin_lo, uin_hi, &
                                 qm, qp, It_lo, It_hi, &
                                 flux, qint, f_lo, f_hi, &
                                 area, a_lo, a_hi, &
                                 qaux, qa_lo, qa_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use mol_module, only: construct_flux

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: g_lo(3), g_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)

    real(rt), intent(in   ) :: div(g_lo(1):g_hi(1), g_lo(2):g_hi(2), g_lo(3):g_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: dx(3), dt

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call construct_flux(idx, idx, domlo, domhi, dx, dt, idir, &
                        div, g_lo, g_hi, &
                        uin, uin_lo, uin_hi, &
                        qm, qp, It_lo, It_hi, &
                        flux, qint, f_lo, f_hi, &
                        area, a_lo, a_hi, &
                        qaux, qa_lo, qa_hi)

  end subroutine cuda_construct_flux



  attributes(global) &
  subroutine cuda_construct_hydro_update(lo, hi, dx, dt, &
                                         f1, q1, f1_lo, f1_hi, &
                                         f2, q2, f2_lo, f2_hi, &
                                         f3, q3, f3_lo, f3_hi, &
                                         a1, a1_lo, a1_hi, &
                                         a2, a2_lo, a2_hi, &
                                         a3, a3_lo, a3_hi, &
                                         vol, vol_lo, vol_hi, &
                                         update, u_lo, u_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NQ
    use advection_util_module, only: construct_hydro_update

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)
    integer,  intent(in   ) :: a1_lo(3), a1_hi(3)
    integer,  intent(in   ) :: a2_lo(3), a2_hi(3)
    integer,  intent(in   ) :: a3_lo(3), a3_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: dx(3), dt

    real(rt), intent(in   ) :: q1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NQ)
    real(rt), intent(in   ) :: q2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NQ)
    real(rt), intent(in   ) :: q3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NQ)
    real(rt), intent(in   ) :: f1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in   ) :: f2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(in   ) :: f3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(in   ) :: a1(a1_lo(1):a1_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: a2(a2_lo(1):a2_hi(1),a2_lo(2):a2_hi(2),a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: a3(a3_lo(1):a3_hi(1),a3_lo(2):a3_hi(2),a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: update(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call construct_hydro_update(idx, idx, dx, dt, &
                                f1, q1, f1_lo, f1_hi, &
                                f2, q2, f2_lo, f2_hi, &
                                f3, q3, f3_lo, f3_hi, &
                                a1, a1_lo, a1_hi, &
                                a2, a2_lo, a2_hi, &
                                a3, a3_lo, a3_hi, &
                                vol, vol_lo, vol_hi, &
                                update, u_lo, u_hi)

  end subroutine cuda_construct_hydro_update



#ifdef CUDA
  attributes(global) &
#endif
  subroutine cuda_summass(lo,hi,rho,r_lo,r_hi,dx, &
                          vol,v_lo,v_hi,mass)

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

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call summass(idx,idx,rho,r_lo,r_hi,dx,vol,v_lo,v_hi,mass)

  end subroutine cuda_summass

end module cuda_interfaces_module
