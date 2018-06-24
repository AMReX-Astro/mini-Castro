module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
#ifdef AMREX_USE_CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

contains

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

#ifdef AMREX_USE_CUDA
    attributes(device) :: vel, dat, lo, hi, v_lo, v_hi, d_lo, d_hi, dx, xlo, domlo, domhi, nv, nc, bc, time, dt, level, grid_no
#endif

    call dervel &
#ifdef AMREX_USE_CUDA
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

#ifdef AMREX_USE_CUDA
    attributes(managed) :: p, u, lo, hi, p_lo, p_hi, u_lo, u_hi, dx, xlo, domlo, domhi, np, nc, bc, time, dt, level, grid_no
#endif

    call derpres &
#ifdef AMREX_USE_CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (p,p_lo,p_hi,np, &
          u,u_lo,u_hi,nc, &
          lo,hi,domlo,domhi, &
          dx,xlo,time,dt, &
          bc,level,grid_no)

  end subroutine ca_derpres



  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_hypfill")

    use bc_fill_module, only: hypfill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

#ifdef AMREX_USE_CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef AMREX_USE_CUDA
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

#ifdef AMREX_USE_CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef AMREX_USE_CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

  end subroutine ca_denfill

end module c_interface_modules
