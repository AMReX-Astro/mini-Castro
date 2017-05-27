module timestep_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! Courant-condition limited timestep

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="ca_estdt")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX, do_ctu
    use prob_params_module, only: dim
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt)         :: dx(3), dt, dt_tmp

    real(rt)         :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3
    integer          :: i, j, k

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

             if (do_ctu == 1) then
                dt  = min(dt,dt1,dt2,dt3)
             else
                ! method of lines constraint is tougher
                dt_tmp = ONE/dt1
                if (dim >= 2) then
                   dt_tmp = dt_tmp + ONE/dt2
                endif
                if (dim == 3) then
                   dt_tmp = dt_tmp + ONE/dt3
                endif
                dt = min(dt, ONE/dt_tmp)
             endif

          enddo
       enddo
    enddo

  end subroutine ca_estdt

  ! Reactions-limited timestep

  ! Diffusion-limited timestep

  ! Check whether the last timestep violated any of our stability criteria.
  ! If so, suggest a new timestep which would not.

  subroutine ca_check_timestep(s_old, so_lo, so_hi, &
                               s_new, sn_lo, sn_hi, &
                               lo, hi, &
                               dx, dt_old, dt_new) &
                               bind(C, name="ca_check_timestep")

    use bl_constants_module, only: HALF, ONE
    use meth_params_module, only: NVAR, URHO, UTEMP, UEINT, UFS, UFX, UMX, UMZ, &
                                  cfl, do_hydro
    use prob_params_module, only: dim
    use network, only: nspec, naux
    use extern_probin_module, only: small_x
    use eos_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: so_lo(3), so_hi(3)
    integer          :: sn_lo(3), sn_hi(3)
    real(rt)         :: s_old(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),NVAR)
    real(rt)         :: s_new(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
    real(rt)         :: dx(3), dt_old, dt_new

    integer          :: i, j, k
    real(rt)         :: rhooinv, rhoninv
    real(rt)         :: tau_CFL


    real(rt)         :: v(3), c
    type (eos_t)     :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhooinv = ONE / s_old(i,j,k,URHO)
             rhoninv = ONE / s_new(i,j,k,URHO)

             ! CFL hydrodynamic stability criterion

             ! If the timestep violated (v+c) * dt / dx > 1,
             ! suggest a new timestep such that (v+c) * dt / dx <= CFL,
             ! where CFL is the user's chosen timestep constraint.
             ! We don't use the CFL choice in the test for violation
             ! because if we did, then even a small increase in velocity
             ! over the timestep would be enough to trigger a retry
             ! (and empirically we have found that this can happen
             ! basically every timestep, which can greatly increase
             ! the cost of a simulation for not much benefit);
             ! but these types of issues are exactly why a user picks a
             ! safe buffer value like CFL = 0.8 or CFL = 0.5. We only
             ! want to trigger a retry if the timestep strongly violated
             ! the stability criterion.

             if (do_hydro .eq. 1) then

                eos_state % rho = s_new(i,j,k,URHO )
                eos_state % T   = s_new(i,j,k,UTEMP)
                eos_state % e   = s_new(i,j,k,UEINT) * rhoninv
                eos_state % xn  = s_new(i,j,k,UFS:UFS+nspec-1) * rhoninv
                eos_state % aux = s_new(i,j,k,UFX:UFX+naux-1) * rhoninv

                call eos(eos_input_re, eos_state)

                v = HALF * (s_old(i,j,k,UMX:UMZ) * rhooinv + s_new(i,j,k,UMX:UMZ) * rhoninv)

                c = eos_state % cs

                tau_CFL = minval(dx(1:dim) / (c + abs(v(1:dim))))

                if (dt_old > tau_CFL) then
                   dt_new = min(dt_new, cfl * tau_CFL)
                endif

             endif

          enddo
       enddo
    enddo

  end subroutine ca_check_timestep

end module timestep_module
