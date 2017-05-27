module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.
  
  subroutine ca_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR
    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_denfill")

    use prob_params_module, only: dim  

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill


  
#ifdef GRAVITY
  subroutine ca_phigravfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_phigravfill")

    use prob_params_module, only: dim  

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_lo(3),phi_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call filcc_nd(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill
  

  
  subroutine ca_gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravxfill")

    use prob_params_module, only: dim  

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravyfill")

    use prob_params_module, only: dim  

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravzfill")

    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill
#endif

  

end module bc_fill_module
