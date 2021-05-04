
module initdata_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream

  implicit none

  ! Problem setup parameters

  real(rt), parameter :: p_ambient = 1.e21_rt        ! ambient pressure (in erg/cc)
  real(rt), parameter :: dens_ambient = 1.e4_rt      ! ambient density (in g/cc)
  real(rt), parameter :: exp_energy = 1.e52_rt       ! absolute energy of the explosion (in erg)
  real(rt), parameter :: r_init = 1.25e8_rt          ! initial radius of the explosion (in cm)
  integer, parameter  :: nsub = 10                   ! subgrid zones in the initial model

contains

  CASTRO_FORT_DEVICE subroutine initdata(lo, hi, &
                                         u, u_lo, u_hi, &
                                         dx, problo, probhi) &
                                         bind(C, name='initdata')

    use amrex_constants_module, only: M_PI, FOUR3RD
    use castro_module , only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
    use eos_module, only: eos_t, eos_input_rp, eos
    use network, only: aion, zion

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3), problo(3), probhi(3)

    real(rt) :: xmin, ymin, zmin
    real(rt) :: xx, yy, zz
    real(rt) :: dist
    real(rt) :: vctr, e_exp, eint

    integer :: i,j,k, ii, jj, kk
    integer :: npert, nambient

    real(rt) :: center(3)
    real(rt) :: e_ambient

    type(eos_t) :: eos_state

    ! Convert the ambient pressure into an ambient energy

    eos_state % rho = dens_ambient
    eos_state % p = p_ambient
    eos_state % T = 1.e4   ! an initial guess -- needed for iteration
    eos_state % abar = aion(1)
    eos_state % zbar = zion(1)

    call eos(eos_input_rp, eos_state)

    e_ambient = eos_state % e

    ! Set explosion energy -- we will convert the point-explosion energy into
    ! a corresponding energy distributed throughout the perturbed volume.
    ! Note that this is done to avoid EOS calls in the initialization.

    vctr  = FOUR3RD*M_PI*r_init**3

    e_exp = exp_energy / vctr / dens_ambient

    center(:) = (problo(:)+probhi(:)) / 2.e0_rt

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(u)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(u)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             xmin = problo(1) + dx(1) * dble(i)
             ymin = problo(2) + dx(2) * dble(j)
             zmin = problo(3) + dx(3) * dble(k)

             npert = 0
             nambient = 0

             do kk = 0, nsub-1
                zz = zmin + (dx(3)/dble(nsub))*(kk + 0.5e0_rt)

                do jj = 0, nsub-1
                   yy = ymin + (dx(2)/dble(nsub))*(jj + 0.5e0_rt)

                   do ii = 0, nsub-1
                      xx = xmin + (dx(1)/dble(nsub))*(ii + 0.5e0_rt)

                      dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2

                      if(dist <= r_init**2) then
                         npert = npert + 1
                      else
                         nambient = nambient + 1
                      endif

                   enddo
                enddo
             enddo

             eint = (npert * e_exp + nambient * e_ambient) / nsub**3

             u(i,j,k,URHO) = dens_ambient
             u(i,j,k,UMX) = 0.e0_rt
             u(i,j,k,UMY) = 0.e0_rt
             u(i,j,k,UMZ) = 0.e0_rt

             u(i,j,k,UTEMP) = 1.d9 ! Arbitrary temperature that will be overwritten on a computeTemp call.

             u(i,j,k,UEDEN) = u(i,j,k,URHO) * eint

             u(i,j,k,UEINT) = u(i,j,k,URHO) * eint

             ! initialize species
             u(i,j,k,UFS:) = 0.e0_rt
             u(i,j,k,UFS) = u(i,j,k,URHO)

          enddo
       enddo
    enddo

  end subroutine initdata

end module initdata_module



subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name='amrex_probinit')

  use initdata_module, only: p_ambient, dens_ambient, exp_energy, r_init, nsub
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

end subroutine amrex_probinit
