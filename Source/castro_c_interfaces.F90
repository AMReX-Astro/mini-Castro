module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

contains

  subroutine ca_initdata(level, lo, hi, &
                         state, s_lo, s_hi, dx, &
                         xlo, xhi) bind(C, name="ca_initdata")

    use initdata_module, only: initdata

    implicit none

    integer,  intent(in   ), value :: level
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: xlo(3), xhi(3), dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(device) :: lo, hi, s_lo, s_hi, xlo, xhi, dx, state
#endif

    call initdata &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (level, lo, hi, state, s_lo, s_hi, dx, xlo, xhi)

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
#endif

    call enforce_consistent_e &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
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
#endif

    call compute_temp &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, state, s_lo, s_hi)

  end subroutine ca_compute_temp



  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi) &
                                  bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi
#endif

    call normalize_species &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
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
    attributes(managed) :: uin, uout, vol, lo, hi, uin_lo, uin_hi, uout_lo, uout_hi, vol_lo, vol_hi, frac_change
#endif

    call enforce_minimum_density &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (uin, uin_lo, uin_hi, &
          uout, uout_lo, uout_hi, &
          vol, vol_lo, vol_hi, &
          lo, hi, frac_change, verbose)

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: lo, hi, state, state_lo, state_hi
#endif

    call check_initial_species &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species



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
    attributes(device) :: vel, dat, lo, hi, v_lo, v_hi, d_lo, d_hi, dx, xlo, domlo, domhi, nv, nc, bc, time, dt, level, grid_no
#endif

    call dervel &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (vel,v_lo,v_hi,nv, &
          dat,d_lo,d_hi,nc, &
          lo,hi,domlo,domhi, &
          dx,xlo,time,dt, &
          bc,level,grid_no)

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
    attributes(managed) :: p, u, lo, hi, p_lo, p_hi, u_lo, u_hi, dx, xlo, domlo, domhi, np, nc, bc, time, dt, level, grid_no
#endif

    call derpres &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (p,p_lo,p_hi,np, &
          u,u_lo,u_hi,nc, &
          lo,hi,domlo,domhi, &
          dx,xlo,time,dt, &
          bc,level,grid_no)

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
    attributes(device) :: u, lo, hi, u_lo, u_hi, dx, dt
#endif

    call estdt &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, u, u_lo, u_hi, dx, dt)

  end subroutine ca_estdt



  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx, &
                        vol,v_lo,v_hi,mass) bind(c, name='ca_summass')

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
    attributes(device) :: rho, vol, lo, hi, r_lo, r_hi, v_lo, v_hi, dx, mass
#endif

    call summass &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, rho, r_lo, r_hi, dx, vol, v_lo, v_hi, mass)

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
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

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
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

  end subroutine ca_denfill

end module c_interface_modules
