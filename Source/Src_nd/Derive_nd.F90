module derive_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="ca_dervel")
    !
    ! This routine will derive the velocity from the momentum.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine ca_dervel



  subroutine ca_derpres(p,p_lo,p_hi,ncomp_p, &
                        u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_derpres")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: p_lo(3), p_hi(3), ncomp_p
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho  = u(i,j,k,URHO)
             eos_state % T    = u(i,j,k,UTEMP)
             eos_state % e    = u(i,j,k,UEINT) * rhoInv
             eos_state % xn   = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux  = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             p(i,j,k,1) = eos_state % p
          enddo
       enddo
    enddo

  end subroutine ca_derpres

end module derive_module
