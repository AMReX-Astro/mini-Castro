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

  CASTRO_FORT_DEVICE subroutine enforce_minimum_density(lo, hi, state, s_lo, s_hi) bind(C, name='enforce_minimum_density')

    use amrex_constants_module, only: ZERO
    use network, only: nspec
    use eos_module, only: eos_t, eos_input_rt, eos

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer      :: i, ii, j, jj, k, kk
    integer      :: i_set, j_set, k_set
    real(rt)     :: max_dens
    integer      :: n, ispec
    type (eos_t) :: eos_state

    max_dens = ZERO

    !$acc parallel loop gang vector collapse(3) deviceptr(state) private(eos_state) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) < small_dens) then

                ! Reset to the characteristics of the adjacent state with the highest density.

                max_dens = state(i,j,k,URHO)

                i_set = i
                j_set = j
                k_set = k

                do kk = -1, 1
                   do jj = -1, 1
                      do ii = -1, 1

                         if (i+ii >= s_lo(1) .and. j+jj <= s_lo(2) .and. k+kk .ge. lo(3) .and. &
                             i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk .le. hi(3)) then

                            if (state(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then
                               i_set = i + ii
                               j_set = j + jj
                               k_set = k + kk
                               max_dens = state(i_set,j_set,k_set,URHO)
                            end if

                         end if

                      end do
                   end do
                end do

                if (max_dens < small_dens) then

                   ! We could not find any nearby zones with sufficient density.
                   ! Our only recourse is to set the density equal to small_dens,
                   ! and the temperature equal to small_temp. We set the velocities
                   ! to zero, though any choice here would be arbitrary.

                   do ispec = 1, nspec
                      n = UFS + ispec - 1
                      state(i,j,k,n) = state(i,j,k,n) * (small_dens / state(i,j,k,URHO))
                   end do

                   eos_state % rho = small_dens
                   eos_state % T   = small_temp
                   eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / small_dens

                   call eos(eos_input_rt, eos_state)

                   state(i,j,k,URHO ) = eos_state % rho
                   state(i,j,k,UTEMP) = eos_state % T

                   state(i,j,k,UMX  ) = ZERO
                   state(i,j,k,UMY  ) = ZERO
                   state(i,j,k,UMZ  ) = ZERO

                   state(i,j,k,UEINT) = eos_state % rho * eos_state % e
                   state(i,j,k,UEDEN) = state(i,j,k,UEINT)

                else

                   ! Reset to the characteristics of the target zone.

                   state(i,j,k,:) = state(i_set,j_set,k_set,:)

                endif

             endif

          enddo
       enddo
    enddo

  end subroutine enforce_minimum_density



  CASTRO_FORT_DEVICE subroutine normalize_species(lo, hi, state, s_lo, s_hi) bind(C, name='normalize_species')

    use network, only: nspec
    use amrex_constants_module, only: ONE

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i, j, k

    !$acc parallel loop gang vector collapse(3) deviceptr(state) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             state(i,j,k,UFS:UFS+nspec-1) = max(1.0d-30 * state(i,j,k,URHO), min(state(i,j,k,URHO), state(i,j,k,UFS:UFS+nspec-1)))

             state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,UFS:UFS+nspec-1) / sum(state(i,j,k,UFS:UFS+nspec-1))

          enddo
       enddo
    enddo

  end subroutine normalize_species



  CASTRO_FORT_DEVICE subroutine reset_internal_e(lo, hi, u, u_lo, u_hi) bind(C, name='reset_internal_e')

    use eos_module, only: eos_t, eos_input_re, eos_input_rt, eos
    use network, only: nspec
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

    !$acc parallel loop gang vector collapse(3) deviceptr(u) private(eos_state) async(acc_stream)
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
                   eos_state % xn(:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv

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
                   eos_state % xn(:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv

                   call eos(eos_input_rt, eos_state)

                   eint_new = eos_state % e

                   u(i,j,k,UEINT) = u(i,j,k,URHO) * eint_new

                endif

             end if
          enddo
       enddo
    enddo

  end subroutine reset_internal_e



  CASTRO_FORT_DEVICE subroutine compute_temp(lo, hi, state, s_lo, s_hi) bind(C, name='compute_temp')

    use network, only: nspec
    use eos_module, only: eos_input_re, eos_t, eos
    use amrex_constants_module, only: ZERO, ONE

    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i,j,k
    real(rt) :: rhoInv

    type (eos_t) :: eos_state

    !$acc parallel loop gang vector collapse(3) deviceptr(state) private(eos_state) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)

             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
             eos_state % e   = state(i,j,k,UEINT) * rhoInv
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv

             call eos(eos_input_re, eos_state)

             state(i,j,k,UTEMP) = eos_state % T
             state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

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

    !$acc parallel loop gang vector collapse(3) deviceptr(tag, den) async(acc_stream)
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
                                                       state, s_lo, s_hi, &
                                                       dx, problo, probhi, &
                                                       blast_mass, blast_radius, &
                                                       max_density) &
                                                       bind(C, name="calculate_blast_radius")

    use amrex_constants_module, only: HALF, TWO
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
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

    !$acc parallel loop gang vector collapse(3) deviceptr(state) reduction(+:blast_mass, blast_radius) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = problo(1) + (i + HALF) * dx(1) - center(1)
             y = problo(2) + (j + HALF) * dx(2) - center(2)
             z = problo(3) + (k + HALF) * dx(3) - center(3)

             if (abs(state(i,j,k,URHO) - max_density) / max_density <= density_tolerance) then
                call reduce_add(blast_mass, state(i,j,k,URHO))
                call reduce_add(blast_radius, state(i,j,k,URHO) * sqrt(x**2 + y**2 + z**2))
             end if
          end do
       end do
    end do

  end subroutine calculate_blast_radius
  
end module castro_module
