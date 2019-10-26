module timestep_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  ! Courant-condition limited timestep

  AMREX_CUDA_FORT_DEVICE subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name='ca_estdt')

    use network, only: nspec
    use castro_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS
    use amrex_constants_module, only: ONE
    use eos_module, only: eos_t, eos_input_re, eos
    use reduction_module, only: reduce_min

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt

    real(rt) :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3, dt_tmp
    integer  :: i, j, k

    type (eos_t) :: eos_state

    ! Call EOS for the purpose of computing sound speed

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho = u(i,j,k,URHO )
             eos_state % T   = u(i,j,k,UTEMP)
             eos_state % e   = u(i,j,k,UEINT) * rhoInv
             eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv

             call eos(eos_input_re, eos_state)

             ! Compute velocity and then calculate CFL timestep.

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

             c = eos_state % cs

             dt1 = dx(1)/(c + abs(ux))
             dt2 = dx(2)/(c + abs(uy))
             dt3 = dx(3)/(c + abs(uz))

             dt_tmp = ONE / dt1 + ONE / dt2 + ONE / dt3

             call reduce_min(dt, ONE / dt_tmp)

          enddo
       enddo
    enddo

  end subroutine ca_estdt

end module timestep_module
