module timestep_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  ! Courant-condition limited timestep

#ifdef CUDA
  attributes(global) &
#endif
  subroutine estdt(lo,hi,u,u_lo,u_hi,dx,dt)

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use bl_constants_module, only: ONE
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds, amrex_min
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX
    use prob_params_module, only: dim

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt

    real(rt) :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3, dt_tmp
    integer  :: i, j, k
    integer  :: blo(3), bhi(3)

    type (eos_t) :: eos_state

    call get_loop_bounds(blo, bhi, lo, hi)

    ! Call EOS for the purpose of computing sound speed

    do k = blo(3), bhi(3)
       do j = blo(2), bhi(2)
          do i = blo(1), bhi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho = u(i,j,k,URHO )
             eos_state % T   = u(i,j,k,UTEMP)
             eos_state % e   = u(i,j,k,UEINT) * rhoInv
             eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             ! Compute velocity and then calculate CFL timestep.

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

             c = eos_state % cs

             dt1 = dx(1)/(c + abs(ux))
             if (dim >= 2) then
                dt2 = dx(2)/(c + abs(uy))
             else
                dt2 = dt1
             endif
             if (dim == 3) then
                dt3 = dx(3)/(c + abs(uz))
             else
                dt3 = dt1
             endif

             dt_tmp = ONE/dt1
             if (dim >= 2) then
                dt_tmp = dt_tmp + ONE/dt2
             endif
             if (dim == 3) then
                dt_tmp = dt_tmp + ONE/dt3
             endif

             call amrex_min(dt, ONE / dt_tmp)

          enddo
       enddo
    enddo

  end subroutine estdt

end module timestep_module
