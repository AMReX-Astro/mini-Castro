module castro_module

  use amrex_fort_module, only: rt => amrex_real
  use actual_network, only: nspec

  implicit none

  !---------------------------------------------------------------------
  ! conserved state components
  !---------------------------------------------------------------------

  ! NTHERM: number of thermodynamic variables (rho, 3 momenta, rho*e, rho*E, T)
  integer, parameter :: NTHERM = 7

  ! NVAR  : number of total variables in initial system  
  integer, parameter :: NVAR = NTHERM + nspec

  ! We use these to index into the state "U"
  integer, parameter :: URHO = 1
  integer, parameter :: UMX = 2
  integer, parameter :: UMY = 3
  integer, parameter :: UMZ = 4
  integer, parameter :: UEDEN = 5
  integer, parameter :: UEINT = 6
  integer, parameter :: UTEMP = 7 ! == NTHERM
  integer, parameter :: UFS = NTHERM + 1

  !---------------------------------------------------------------------
  ! primitive state components
  !---------------------------------------------------------------------

  ! QTHERM: number of primitive variables: rho, p, (rho e), T + 3 velocity components 
  integer, parameter :: QTHERM = NTHERM + 1 ! the + 1 is for QGAME which is always defined in primitive mode

  ! QVAR  : number of total variables in primitive form
  integer, parameter :: QVAR = QTHERM + nspec

  ! We use these to index into the state "Q"
  integer, parameter :: QRHO = 1
  integer, parameter :: QU = 2
  integer, parameter :: QV = 3
  integer, parameter :: QW = 4
  integer, parameter :: QGAME = 5
  integer, parameter :: QPRES = 6
  integer, parameter :: QREINT = 7
  integer, parameter :: QTEMP = 8 ! == QTHERM
  integer, parameter :: QFS = QTHERM + 1

  ! The NQAUX here are auxiliary quantities (game, gamc, c, csml, dpdr, dpde)
  ! that we create in the primitive variable call but that do not need to
  ! participate in tracing.
  integer, parameter :: NQAUX = 5
  integer, parameter :: QGAMC = 1
  integer, parameter :: QC    = 2
  integer, parameter :: QDPDR = 3
  integer, parameter :: QDPDE = 4

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, parameter :: NGDNV = 6
  integer, parameter :: GDRHO = 1
  integer, parameter :: GDU = 2
  integer, parameter :: GDV = 3
  integer, parameter :: GDW = 4
  integer, parameter :: GDPRES = 5
  integer, parameter :: GDGAME = 6

  real(rt), parameter :: small_dens = 1.0d-12
  real(rt), parameter :: small_temp = 1.0d3
  real(rt), parameter :: cfl = 0.5d0

contains

  subroutine ca_reset_internal_e(lo,hi,u,u_lo,u_hi) bind(c,name='ca_reset_internal_e')

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re, eos_input_rt
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

    !$gpu

    ! Reset internal energy

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

  end subroutine ca_reset_internal_e



  subroutine ca_compute_temp(lo,hi,state,s_lo,s_hi) bind(c,name='ca_compute_temp')

    use network, only: nspec
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use amrex_constants_module, only: ZERO, ONE

    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i,j,k
    real(rt) :: rhoInv

    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! First check the inputs for validity.

#ifndef AMREX_USE_CUDA
             if (state(i,j,k,URHO) <= ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_util.F90::ca_compute_temp ',i,j,k
                print *,'>>> ... negative density ',state(i,j,k,URHO)
                print *,'    '
                call bl_error("Error:: compute_temp_nd.F90")
             end if

             if (state(i,j,k,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_util.F90::ca_compute_temp ',i,j,k
                print *,'>>> ... negative (rho e) ',state(i,j,k,UEINT)
                print *,'   '
                call bl_error("Error:: compute_temp_nd.F90")
             end if
#endif

             rhoInv = ONE / state(i,j,k,URHO)

             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
             eos_state % e   = state(i,j,k,UEINT) * rhoInv
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv

             call eos(eos_input_re, eos_state)

             state(i,j,k,UTEMP) = eos_state % T

             ! In case we've floored, or otherwise allowed the energy to change, update the energy accordingly.

             state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

          enddo
       enddo
    enddo

  end subroutine ca_compute_temp
  


  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi) bind(c,name='ca_normalize_species')

    use network, only: nspec
    use amrex_constants_module, only: ONE

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k
    real(rt) :: xn(nspec)

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             xn = u(i,j,k,UFS:UFS+nspec-1)

             xn = max(1.0d-30 * u(i,j,k,URHO), min(u(i,j,k,URHO), xn))

             xn = u(i,j,k,URHO) * (xn / sum(xn))

             u(i,j,k,UFS:UFS+nspec-1) = xn

          enddo
       enddo
    enddo

  end subroutine ca_normalize_species



  subroutine ca_denerror(lo, hi, &
                         tag, taglo, taghi, &
                         den, denlo, denhi, &
                         set, clear) &
                         bind(C, name="ca_denerror")

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

    !$gpu

    ! Tag on regions of high density gradient
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

  end subroutine ca_denerror
  
end module castro_module



subroutine ca_network_init() bind(C, name="ca_network_init")

  use network, only: network_init
  call network_init()

end subroutine ca_network_init



subroutine ca_get_num_spec(nspec_out) bind(C, name="ca_get_num_spec")

  use network, only: nspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(out) :: nspec_out

  nspec_out = nspec

end subroutine ca_get_num_spec



subroutine ca_get_spec_names(spec_names,ispec,len) &
     bind(C, name="ca_get_spec_names")

  use network, only: nspec, short_spec_names
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in   ) :: ispec
  integer, intent(inout) :: len
  integer, intent(inout) :: spec_names(len)

  ! Local variables
  integer   :: i

  len = len_trim(short_spec_names(ispec+1))

  do i = 1,len
     spec_names(i) = ichar(short_spec_names(ispec+1)(i:i))
  end do

end subroutine ca_get_spec_names



subroutine ca_get_qvar(qvar_in) bind(C, name="ca_get_qvar")

  use castro_module, only: QVAR
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(inout) :: qvar_in

  qvar_in = QVAR

end subroutine ca_get_qvar



subroutine ca_get_nqaux(nqaux_in) bind(C, name="ca_get_nqaux")

  use castro_module, only: NQAUX
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(inout) :: nqaux_in

  nqaux_in = NQAUX

end subroutine ca_get_nqaux



subroutine ca_get_ngdnv(ngdnv_in) bind(C, name="ca_get_ngdnv")

  use castro_module, only: NGDNV

  implicit none

  integer, intent(inout) :: ngdnv_in

  ngdnv_in = NGDNV

end subroutine ca_get_ngdnv



subroutine ca_set_method_params() bind(C, name="ca_set_method_params")

  use eos_module, only: eos_init

  implicit none

  integer :: ispec

  integer :: i
  integer :: ioproc

  call eos_init()

end subroutine ca_set_method_params
