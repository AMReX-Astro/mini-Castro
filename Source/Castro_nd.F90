module castro_module

  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec
  use amrex_acc_module, only: acc_stream

  implicit none

  !---------------------------------------------------------------------
  ! conserved state components
  !---------------------------------------------------------------------

  ! We use these to index into the state "U"
  integer, parameter :: URHO = 1
  integer, parameter :: UMX = 2
  integer, parameter :: UMY = 3
  integer, parameter :: UMZ = 4
  integer, parameter :: UEDEN = 5
  integer, parameter :: UEINT = 6
  integer, parameter :: UTEMP = 7
  integer, parameter :: NTHERM = UTEMP
  integer, parameter :: UFS = NTHERM + 1
  integer, parameter :: NVAR = NTHERM + nspec

  !---------------------------------------------------------------------
  ! primitive state components
  !---------------------------------------------------------------------

  ! We use these to index into the state "Q"
  integer, parameter :: QRHO = 1
  integer, parameter :: QU = 2
  integer, parameter :: QV = 3
  integer, parameter :: QW = 4
  integer, parameter :: QGAME = 5
  integer, parameter :: QPRES = 6
  integer, parameter :: QREINT = 7
  integer, parameter :: QTEMP = 8
  integer, parameter :: QTHERM = QTEMP
  integer, parameter :: QFS = QTHERM + 1
  integer, parameter :: QVAR = QTHERM + nspec

  ! The NQAUX here are auxiliary quantities for the primitive variable
  ! call but that do not need to participate in tracing.
  integer, parameter :: QGAMC = 1
  integer, parameter :: QC    = 2
  integer, parameter :: QDPDR = 3
  integer, parameter :: QDPDE = 4
  integer, parameter :: NQAUX = 4

  ! These are used for the Godunov state
  integer, parameter :: GDRHO = 1
  integer, parameter :: GDU = 2
  integer, parameter :: GDV = 3
  integer, parameter :: GDW = 4
  integer, parameter :: GDPRES = 5
  integer, parameter :: GDGAME = 6
  integer, parameter :: NGDNV = 6

  real(rt), parameter :: small_dens = 1.0d-12
  real(rt), parameter :: small_temp = 1.0d3
  real(rt), parameter :: small_pres = 1.e-200_rt
  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt
  real(rt), parameter :: dual_energy_eta1 = 1.e0_rt
  real(rt), parameter :: cfl = 0.5d0

contains

  CASTRO_FORT_DEVICE subroutine enforce_minimum_density(lo, hi, u, u_lo, u_hi) bind(C, name='enforce_minimum_density')

    use amrex_constants_module, only: ZERO, ONE
    use network, only: nspec, aion_inv, zion
    use eos_module, only: eos_t, eos_input_rt, eos

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer      :: i, j, k
    integer      :: n, ispec
    type (eos_t) :: eos_state

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) private(eos_state) deviceptr(u) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) private(eos_state) is_device_ptr(u)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (u(i,j,k,URHO) < small_dens) then

                do ispec = 1, nspec
                   n = UFS + ispec - 1
                   u(i,j,k,n) = u(i,j,k,n) * (small_dens / u(i,j,k,URHO))
                end do

                eos_state % rho = small_dens
                eos_state % T   = small_temp
                eos_state % abar = ONE / (sum(u(i,j,k,UFS:UFS+nspec-1) * aion_inv(:)) / small_dens)
                eos_state % zbar = eos_state % abar * (sum(u(i,j,k,UFS:UFS+nspec-1) * zion(:) * aion_inv(:)) / small_dens)

                call eos(eos_input_rt, eos_state)

                u(i,j,k,URHO ) = eos_state % rho
                u(i,j,k,UTEMP) = eos_state % T

                u(i,j,k,UMX  ) = ZERO
                u(i,j,k,UMY  ) = ZERO
                u(i,j,k,UMZ  ) = ZERO

                u(i,j,k,UEINT) = eos_state % rho * eos_state % e
                u(i,j,k,UEDEN) = u(i,j,k,UEINT)

             endif

          enddo
       enddo
    enddo

  end subroutine enforce_minimum_density



  CASTRO_FORT_DEVICE subroutine normalize_species(lo, hi, u, u_lo, u_hi) bind(C, name='normalize_species')

    use network, only: nspec
    use amrex_constants_module, only: ONE

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer  :: i, j, k

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(u) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(u)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             u(i,j,k,UFS:UFS+nspec-1) = max(1.0d-30 * u(i,j,k,URHO), min(u(i,j,k,URHO), u(i,j,k,UFS:UFS+nspec-1)))

             u(i,j,k,UFS:UFS+nspec-1) = u(i,j,k,UFS:UFS+nspec-1) / sum(u(i,j,k,UFS:UFS+nspec-1))

          enddo
       enddo
    enddo

  end subroutine normalize_species



  CASTRO_FORT_DEVICE subroutine reset_internal_e(lo, hi, u, u_lo, u_hi) bind(C, name='reset_internal_e')

    use eos_module, only: eos_t, eos_input_re, eos_input_rt, eos
    use network, only: nspec, aion_inv, zion
    use amrex_constants_module, only: ZERO, HALF, ONE

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer  :: i,j,k
    real(rt) :: Up, Vp, Wp, ke, rho_eint, eden, small_e, eint_new, rhoInv

    real(rt), parameter :: dual_energy_eta2 = 1.e-4_rt

    type (eos_t) :: eos_state

    ! Reset internal energy

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) private(eos_state) deviceptr(u) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) private(eos_state) is_device_ptr(u)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE/u(i,j,k,URHO)
             Up = u(i,j,k,UMX) * rhoInv
             Vp = u(i,j,k,UMY) * rhoInv
             Wp = u(i,j,k,UMZ) * rhoInv
             ke = HALF * (Up**2 + Vp**2 + Wp**2)

             if (u(i,j,k,UEDEN) < ZERO) then

                if (u(i,j,k,UEINT) < ZERO) then

                   eos_state % rho   = u(i,j,k,URHO)
                   eos_state % T     = small_temp
                   eos_state % abar = ONE / (sum(u(i,j,k,UFS:UFS+nspec-1) * aion_inv(:)) * rhoInv)
                   eos_state % zbar = eos_state % abar * (sum(u(i,j,k,UFS:UFS+nspec-1) * zion(:) * aion_inv(:)) * rhoInv)

                   call eos(eos_input_rt, eos_state)

                   u(i,j,k,UEINT) = u(i,j,k,URHO) * eos_state % e

                endif

                u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

             else

                rho_eint = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

                ! Reset (e from e) if it's greater than eta * E.
                if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,k,UEDEN) .gt. dual_energy_eta2) then

                   u(i,j,k,UEINT) = rho_eint

                ! If not resetting and little e is negative ...
                else if (u(i,j,k,UEINT) .le. ZERO) then

                   eos_state % rho   = u(i,j,k,URHO)
                   eos_state % T     = small_temp
                   eos_state % abar = ONE / (sum(u(i,j,k,UFS:UFS+nspec-1) * aion_inv(:)) * rhoInv)
                   eos_state % zbar = eos_state % abar * (sum(u(i,j,k,UFS:UFS+nspec-1) * zion(:) * aion_inv(:)) * rhoInv)

                   call eos(eos_input_rt, eos_state)

                   eint_new = eos_state % e

                   u(i,j,k,UEINT) = u(i,j,k,URHO) * eint_new

                endif

             end if
          enddo
       enddo
    enddo

  end subroutine reset_internal_e



  CASTRO_FORT_DEVICE subroutine compute_temp(lo, hi, u, u_lo, u_hi) bind(C, name='compute_temp')

    use network, only: nspec, aion_inv, zion
    use eos_module, only: eos_input_re, eos_t, eos
    use amrex_constants_module, only: ZERO, ONE

    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    integer , intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer  :: i,j,k
    real(rt) :: rhoInv

    type (eos_t) :: eos_state

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) private(eos_state) deviceptr(u) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) private(eos_state) is_device_ptr(u)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho = u(i,j,k,URHO)
             eos_state % T   = u(i,j,k,UTEMP) ! Initial guess for the EOS
             eos_state % e   = u(i,j,k,UEINT) * rhoInv
             eos_state % abar = ONE / (sum(u(i,j,k,UFS:UFS+nspec-1) * aion_inv(:)) * rhoInv)
             eos_state % zbar = eos_state % abar * (sum(u(i,j,k,UFS:UFS+nspec-1) * zion(:) * aion_inv(:)) * rhoInv)

             call eos(eos_input_re, eos_state)

             u(i,j,k,UTEMP) = eos_state % T
             u(i,j,k,UEINT) = u(i,j,k,URHO) * eos_state % e

          enddo
       enddo
    enddo

  end subroutine compute_temp



  CASTRO_FORT_DEVICE subroutine denerror(lo, hi, &
                                         tag, taglo, taghi, &
                                         den, denlo, denhi, &
                                         set, clear) &
                                         bind(C, name="denerror")

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: denlo(3), denhi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3))
    integer(1), intent(in   ), value :: set, clear

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    real(rt), parameter :: dengrad_rel = 0.25d0

    ! Tag on regions of high density gradient

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(tag, den) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(tag, den)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ax = ABS(den(i+1,j,k) - den(i,j,k))
             ay = ABS(den(i,j+1,k) - den(i,j,k))
             az = ABS(den(i,j,k+1) - den(i,j,k))
             ax = MAX(ax,ABS(den(i,j,k) - den(i-1,j,k)))
             ay = MAX(ay,ABS(den(i,j,k) - den(i,j-1,k)))
             az = MAX(az,ABS(den(i,j,k) - den(i,j,k-1)))
             if (MAX(ax,ay,az) .ge. ABS(dengrad_rel * den(i,j,k))) then
                tag(i,j,k) = set
             end if
          end do
       end do
    end do

  end subroutine denerror



  CASTRO_FORT_DEVICE subroutine calculate_blast_radius(lo, hi, &
                                                       u, u_lo, u_hi, &
                                                       dx, problo, probhi, &
                                                       blast_mass, blast_radius, &
                                                       max_density) &
                                                       bind(C, name="calculate_blast_radius")

    use amrex_constants_module, only: HALF, TWO
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3), problo(3), probhi(3)
    real(rt), intent(inout) :: blast_mass, blast_radius
    real(rt), intent(in   ), value :: max_density

    integer  :: i, j, k
    real(rt) :: x, y, z
    real(rt) :: center(3)

    real(rt), parameter :: density_tolerance = 0.1d0

    ! Add to the (mass-weighted) blast radius if the density of this zone
    ! is within density_tolerance of the maximum.

    center = (probhi - problo) / TWO

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(u) reduction(+:blast_mass, blast_radius) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(u) reduction(+:blast_mass, blast_radius)
#endif
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = problo(1) + (i + HALF) * dx(1) - center(1)
             y = problo(2) + (j + HALF) * dx(2) - center(2)
             z = problo(3) + (k + HALF) * dx(3) - center(3)

             if (abs(u(i,j,k,URHO) - max_density) / max_density <= density_tolerance) then
                call reduce_add(blast_mass, u(i,j,k,URHO))
                call reduce_add(blast_radius, u(i,j,k,URHO) * sqrt(x**2 + y**2 + z**2))
             end if
          end do
       end do
    end do

  end subroutine calculate_blast_radius
  
end module castro_module
