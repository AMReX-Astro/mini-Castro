module eos_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: M_PI
  use network, only: nspec

  implicit none

  public eos_t, eos_init, eos_finalize, eos

  integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter :: eos_input_re = 2  ! rho, e are inputs
  integer, parameter :: eos_input_rp = 3  ! rho, p are inputs

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  type :: eos_t

    real(rt) :: rho       ! mass density (g/cm**3)
    real(rt) :: T         ! temperature (K)
    real(rt) :: p         ! pressure (dyn/cm**2)
    real(rt) :: e         ! internal energy (erg/g)
    real(rt) :: h         ! enthalpy (erg/g)
    real(rt) :: s         ! entropy (erg/g/K)
    real(rt) :: abar      ! mean nucleon number
    real(rt) :: zbar      ! mean proton number

    real(rt) :: dpdT      ! d pressure / d temperature
    real(rt) :: dpdr      ! d pressure / d density
    real(rt) :: dedT      ! d energy / d temperature
    real(rt) :: dedr      ! d energy / d density
    real(rt) :: dhdT      ! d enthalpy / d temperature
    real(rt) :: dhdr      ! d enthalpy / d density
    real(rt) :: dsdT      ! d entropy / d temperature
    real(rt) :: dsdr      ! d entropy / d density
    real(rt) :: dpde      ! d pressure / d energy
    real(rt) :: dpdr_e    ! d pressure / d density at constant energy

    real(rt) :: cv        ! specific heat at constant volume
    real(rt) :: cp        ! specific heat at constant pressure
    real(rt) :: xne       ! electron number density
    real(rt) :: xnp       ! positron number density
    real(rt) :: eta       ! electron degeneracy parameter
    real(rt) :: pele      ! electron pressure
    real(rt) :: ppos      ! positron pressure
    real(rt) :: mu        ! mean molecular weight
    real(rt) :: mu_e      ! mean number of nucleons per electron
    real(rt) :: y_e       ! electron fraction == 1 / mu_e
    real(rt) :: gam1      ! first adiabatic index (d log P/ d log rho) |_s
    real(rt) :: cs        ! sound speed

    real(rt) :: dpdA      ! d pressure / d abar
    real(rt) :: dpdZ      ! d pressure / d zbar
    real(rt) :: dedA      ! d energy / d abar
    real(rt) :: dedZ      ! d energy / d zbar

  end type eos_t

  integer, parameter, private :: imax = 541, jmax = 201

  real(rt), parameter :: tlo = 3.0d0, thi = 13.0d0
  real(rt), parameter :: dlo = -12.0d0, dhi = 15.0d0

  real(rt), parameter :: tstp = (thi - tlo) / float(jmax-1)
  real(rt), parameter :: tstpi = 1.0d0 / tstp

  real(rt), parameter :: dstp = (dhi - dlo) / float(imax-1)
  real(rt), parameter :: dstpi = 1.0d0/dstp

  real(rt), parameter :: mintemp = 10.0d0**tlo
  real(rt), parameter :: mindens = 10.0d0**dlo

  real(rt), parameter :: ttol = 1.d-8

  ! Density and temperature
  real(rt), allocatable :: d(:), t(:)
  real(rt), allocatable :: dt(:), dt2(:), dti(:), dt2i(:)
  real(rt), allocatable :: dd(:), dd2(:), ddi(:), dd2i(:)

  ! Helmholtz free energy and derivatives
  real(rt), allocatable :: f(:,:)
  real(rt), allocatable :: fd(:,:), fdd(:,:)
  real(rt), allocatable :: ft(:,:), ftt(:,:)
  real(rt), allocatable :: fdt(:,:), fddt(:,:), fdtt(:,:), fddtt(:,:)

  ! Pressure derivatives
  real(rt), allocatable :: dpdf(:,:), dpdfd(:,:), dpdft(:,:), dpdfdt(:,:)

  ! Chemical potential and derivatives
  real(rt), allocatable :: ef(:,:), efd(:,:), eft(:,:), efdt(:,:)

  ! Number density and derivatives
  real(rt), allocatable :: xf(:,:), xfd(:,:), xft(:,:), xfdt(:,:)

#if (defined(AMREX_USE_CUDA) && !(defined(AMREX_USE_ACC) || defined(AMREX_USE_OMP_OFFLOAD)))
  attributes(managed) :: d, t
  attributes(managed) :: dt, dt2, dti, dt2i
  attributes(managed) :: dd, dd2, ddi, dd2i
  attributes(managed) :: f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt
  attributes(managed) :: dpdf, dpdfd, dpdft, dpdfdt
  attributes(managed) :: ef, efd, eft, efdt
  attributes(managed) :: xf, xfd, xft, xfdt
#endif

#ifdef AMREX_USE_ACC
  !$acc declare create(d, t)
  !$acc declare create(dt, dt2, dti, dt2i)
  !$acc declare create(dd, dd2, ddi, dd2i)
  !$acc declare create(f, fd, fdd, ft, ftt, fdt, fddt, fdtt, fddtt)
  !$acc declare create(dpdf, dpdfd, dpdft, dpdfdt)
  !$acc declare create(ef, efd, eft, efdt)
  !$acc declare create(xf, xfd, xft, xfdt)
#endif

#ifdef AMREX_USE_OMP_OFFLOAD
  !$omp declare target(d, t)
  !$omp declare target(dt, dt2, dti, dt2i)
  !$omp declare target(dd, dd2, ddi, dd2i)
  !$omp declare target(f, fd, fdd, ft, ftt, fdt, fddt, fdtt, fddtt)
  !$omp declare target(dpdf, dpdfd, dpdft, dpdfdt)
  !$omp declare target(ef, efd, eft, efdt)
  !$omp declare target(xf, xfd, xft, xfdt)
#endif

  integer, parameter :: max_newton = 100

  ! Physical constants
  real(rt), parameter :: h       = 6.6260689633d-27
  real(rt), parameter :: avo_eos = 6.0221417930d23
  real(rt), parameter :: clight  = 2.99792458d10
  real(rt), parameter :: kerg    = 1.380650424d-16
  real(rt), parameter :: amu     = 1.66053878283d-24
  real(rt), parameter :: asol    = 4.0d0 * 5.67051d-5 / clight
  real(rt), parameter :: sioncon = (2.0d0 * M_PI * amu * kerg) / (h * h)
  real(rt), parameter :: kergavo = kerg * avo_eos
  real(rt), parameter :: asoli3  = asol / 3.0d0

contains

  CASTRO_FORT_DEVICE subroutine eos(input, state)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    use amrex_constants_module, only: ZERO, HALF, ONE, TWO

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    logical :: converged
    integer :: iter
    real(rt) :: temp_old, v_want, v, dvdx, error

    real(rt) :: temp, den, din, deni, tempi, abar, zbar, ytot1, ye
    real(rt) :: pres, ener, entr, dpresdd, dpresdt, denerdd, denerdt, dentrdd, dentrdt
    real(rt) :: pele, dpepdt, dpepdd, eele, deepdt, deepdd, sele, dsepdd, dsepdt
    real(rt) :: prad, dpraddd, dpraddt, erad, deraddd, deraddt, srad, dsraddd, dsraddt
    real(rt) :: pion, dpiondd, dpiondt, eion, deiondd, deiondt, sion, dsiondd, dsiondt
    real(rt) :: s, x, y, z, zz, zzi, chit, chid

    integer  :: iat, jat
    real(rt) :: free, df_d, df_t, df_tt, df_dt
    real(rt) :: xt, xd, mxt, mxd
    real(rt) :: fi(36)
    real(rt) :: si0t, si1t, si2t, si0mt, si1mt, si2mt
    real(rt) :: si0d, si1d, si2d, si0md, si1md, si2md
    real(rt) :: dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt
    real(rt) :: dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md
    real(rt) :: ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    ! Initial setup for iterations

    if (input .eq. eos_input_rt) then
       converged = .true.
    else
       converged = .false.
    end if

    if (input .eq. eos_input_re) then
       v_want = state % e
    else if (input .eq. eos_input_rp) then
       v_want = state % p
    end if

    do iter = 1, max_newton

       temp  = state % T
       den   = state % rho

       ytot1 = 1.0d0 / state % abar
       ye    = state % zbar / state % abar
       din   = ye * den

       !..initialize
       deni    = 1.0d0 / den
       tempi   = 1.0d0 / temp

       !..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad * tempi

       erad    = 3.0d0 * prad * deni
       deraddd = -erad * deni
       deraddt = 3.0d0 * dpraddt * deni

       srad    = (prad * deni + erad) * tempi
       dsraddd = (dpraddd * deni - prad * deni * deni + deraddd) * tempi
       dsraddt = (dpraddt * deni + deraddt - srad) * tempi

       !..ion section:
       pion    = kergavo * ytot1 * den * temp
       dpiondd = kergavo * ytot1 * temp
       dpiondt = kergavo * ytot1 * den

       eion    = 1.5d0 * pion * deni
       deiondd = (1.5d0 * dpiondd - eion) * deni
       deiondt = 1.5d0 * dpiondt * deni

       x       = state % abar * state % abar * sqrt(state % abar) * deni / avo_eos
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion * deni + eion) * tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd * deni - pion * deni * deni + deiondd) * tempi &
                 - kergavo * deni * ytot1
       dsiondt = (dpiondt * deni + deiondt) * tempi - &
                 (pion*deni + eion) * tempi * tempi &
                 + 1.5d0 * kergavo * tempi * ytot1

       !..electron-positron section:

       !..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))

       !..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)

       !..various differences
       xt  = max( (temp - t(jat))*dti(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

       !..the six density and six temperature basis functions
       si0t  =  psi0(xt)
       si1t  =  psi1(xt)*dt(jat)
       si2t  =  psi2(xt)*dt2(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt(jat)
       si2mt =  psi2(mxt)*dt2(jat)

       si0d  =  psi0(xd)
       si1d  =  psi1(xd)*dd(iat)
       si2d  =  psi2(xd)*dd2(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd(iat)
       si2md =  psi2(mxd)*dd2(iat)

       !..derivatives of the weight functions
       dsi0t  =  dpsi0(xt)*dti(jat)
       dsi1t  =  dpsi1(xt)
       dsi2t  =  dpsi2(xt)*dt(jat)

       dsi0mt = -dpsi0(mxt)*dti(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt(jat)

       dsi0d  =  dpsi0(xd)*ddi(iat)
       dsi1d  =  dpsi1(xd)
       dsi2d  =  dpsi2(xd)*dd(iat)

       dsi0md = -dpsi0(mxd)*ddi(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd(iat)

       !..second derivatives of the weight functions
       ddsi0t  =  ddpsi0(xt)*dt2i(jat)
       ddsi1t  =  ddpsi1(xt)*dti(jat)
       ddsi2t  =  ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i(jat)
       ddsi1mt = -ddpsi1(mxt)*dti(jat)
       ddsi2mt =  ddpsi2(mxt)

       !..the free energy
       free = h5(fi, &
                 si0t, si1t, si2t, si0mt, si1mt, si2mt, &
                 si0d, si1d, si2d, si0md, si1md, si2md)

       !..derivative with respect to density
       df_d = h5(fi, &
                 si0t,  si1t,  si2t,  si0mt,  si1mt,  si2mt, &
                 dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md)

       !..derivative with respect to temperature
       df_t = h5(fi, &
                 dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
                 si0d,  si1d,  si2d,  si0md,  si1md,  si2md)

       !..derivative with respect to temperature**2
       df_tt = h5(fi, &
                  ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                  si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

       !..derivative with respect to temperature and density
       df_dt = h5(fi, &
                  dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
                  dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md)

       !..now get the pressure derivative with density, chemical potential, and
       !..electron positron number densities
       !..get the interpolation weight functions
       si0t  = xpsi0(xt)
       si1t  = xpsi1(xt) * dt(jat)

       si0mt = xpsi0(mxt)
       si1mt = -xpsi1(mxt) * dt(jat)

       si0d  = xpsi0(xd)
       si1d  = xpsi1(xd) * dd(iat)

       si0md = xpsi0(mxd)
       si1md = -xpsi1(mxd) * dd(iat)

       !..derivatives of weight functions
       dsi0t  = xdpsi0(xt) * dti(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt) * dti(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd) * ddi(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd) * ddi(iat)
       dsi1md = xdpsi1(mxd)

       !..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

       !..pressure derivative with density
       dpepdd = h3(fi, &
                   si0t, si1t, si0mt, si1mt, &
                   si0d, si1d, si0md, si1md)
       dpepdd  = max(ye * dpepdd, 0.0d0)

       !..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)

       !..electron chemical potential eta
       state % eta = h3(fi, &
                        si0t, si1t, si0mt, si1mt, &
                        si0d, si1d, si0md, si1md)

       !..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

       !..electron + positron number densities
       state % xne = h3(fi, &
                        si0t, si1t, si0mt, si1mt, &
                        si0d, si1d, si0md, si1md)
       state % xnp = 0.0d0

       !..the desired electron-positron thermodynamic quantities

       !..dpepdd at high temperatures and low densities is below the
       !..floating point limit of the subtraction of two large terms.
       !..since dpresdd doesn't enter the maxwell relations at all, use the
       !..bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt

       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x

       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd

       !..sum all the components
       pres    = prad + pion + pele
       ener    = erad + eion + eele
       entr    = srad + sion + sele

       dpresdd = dpraddd + dpiondd + dpepdd
       dpresdt = dpraddt + dpiondt + dpepdt

       denerdd = deraddd + deiondd + deepdd
       denerdt = deraddt + deiondt + deepdt

       dentrdd = dsraddd + dsiondd + dsepdd
       dentrdt = dsraddt + dsiondt + dsepdt

       zz    = pres * deni
       zzi   = den / pres
       chit  = temp / pres * dpresdt
       chid  = dpresdd * zzi
       state % cv = denerdt
       state % gam1 = chit * zz * chit / (temp * state % cv) + chid
       state % cp = state % cv * state % gam1 / chid

       state % p = pres
       state % dpdT = dpresdt
       state % dpdr = dpresdd
       state % dpde = dpresdt / denerdt
       state % dpdr_e = dpresdd - dpresdt * denerdd / denerdt

       state % e = ener
       state % dedT = denerdt
       state % dedr = denerdd

       state % s = entr
       state % dsdT = dentrdt
       state % dsdr = dentrdd

       state % h = ener + pres / den
       state % dhdr = denerdd + dpresdd / den - pres / den**2
       state % dhdT = denerdt + dpresdt / den

       state % pele = pele
       state % ppos = 0.0d0

       if (converged) then

          exit

       else

          temp_old = temp

          if (input .eq. eos_input_re) then
             v    = state % e
             dvdx = state % dedT
          else ! input .eq. eos_input_rp
             v    = state % p
             dvdx = state % dpdT
          end if

          ! Now do the calculation for the next guess for T
          temp = temp - (v - v_want) / dvdx

          ! Don't let the temperature change by more than a factor of two
          temp = max(0.5 * temp_old, min(temp, 2.0 * temp_old))

          ! Don't let us freeze
          temp = max(mintemp, temp)

          state % T = temp

          ! Compute the error from the last iteration
          error = abs((temp - temp_old) / temp_old)

          if (error .lt. ttol) converged = .true.

       end if

    enddo

    state % cs = sqrt(state % gam1 * state % p / state % rho)

    state % dpdA = 0.0d0
    state % dpdZ = 0.0d0

    state % dedA = 0.0d0
    state % dedZ = 0.0d0

    ! Ensure the inputs don't change as a result of the EOS call.

    if (input .eq. eos_input_re) then

       state % e = v_want

    else if (input .eq. eos_input_rp) then

       state % p = v_want

    endif

  end subroutine eos



  subroutine eos_init() bind(C, name='eos_init')

    use amrex_error_module, only: amrex_error
    use amrex_paralleldescriptor_module, only: parallel_bcast => amrex_pd_bcast, amrex_pd_ioprocessor

    implicit none

    integer :: i, j
    integer :: status

    ! Allocate managed module variables

    allocate(d(imax))
    allocate(t(jmax))
    allocate(f(imax,jmax))
    allocate(fd(imax,jmax))
    allocate(ft(imax,jmax))
    allocate(fdd(imax,jmax))
    allocate(ftt(imax,jmax))
    allocate(fdt(imax,jmax))
    allocate(fddt(imax,jmax))
    allocate(fdtt(imax,jmax))
    allocate(fddtt(imax,jmax))
    allocate(dpdf(imax,jmax))
    allocate(dpdfd(imax,jmax))
    allocate(dpdft(imax,jmax))
    allocate(dpdfdt(imax,jmax))
    allocate(ef(imax,jmax))
    allocate(efd(imax,jmax))
    allocate(eft(imax,jmax))
    allocate(efdt(imax,jmax))
    allocate(xf(imax,jmax))
    allocate(xfd(imax,jmax))
    allocate(xft(imax,jmax))
    allocate(xfdt(imax,jmax))
    allocate(dt(jmax))
    allocate(dt2(jmax))
    allocate(dti(jmax))
    allocate(dt2i(jmax))
    allocate(dd(imax))
    allocate(dd2(imax))
    allocate(ddi(imax))
    allocate(dd2i(imax))

    ! Read the table

    do j = 1, jmax
       t(j) = 10.0d0**(tlo + (j-1)*tstp)
    end do

    do i = 1, imax
       d(i) = 10.0d0**(dlo + (i-1)*dstp)
    end do

    if (amrex_pd_ioprocessor()) then

       ! Open the table
       open(unit=2, file='helm_table.dat', status='old', iostat=status, action='read')

       if (status > 0) then
          call amrex_error('eos_init: Failed to open helm_table.dat')
       endif

       ! Read the free energy and derivatives
       do j=1,jmax
          do i=1,imax
             read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                  fddt(i,j),fdtt(i,j),fddtt(i,j)
          end do
       end do

       ! Read the pressure derivatives
       do j = 1, jmax
          do i = 1, imax
             read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
          end do
       end do

       ! Read the chemical potential and derivatives
       do j = 1, jmax
          do i = 1, imax
             read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
          end do
       end do

       ! Read the number density and derivatives
       do j = 1, jmax
          do i = 1, imax
             read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
          end do
       end do
    end if

    if (amrex_pd_ioprocessor()) then
       close(unit=2)
    endif

    call parallel_bcast(f)
    call parallel_bcast(fd)
    call parallel_bcast(ft)
    call parallel_bcast(fdd)
    call parallel_bcast(ftt)
    call parallel_bcast(fdt)
    call parallel_bcast(fddt)
    call parallel_bcast(fdtt)
    call parallel_bcast(fddtt)
    call parallel_bcast(dpdf)
    call parallel_bcast(dpdfd)
    call parallel_bcast(dpdft)
    call parallel_bcast(dpdfdt)
    call parallel_bcast(ef)
    call parallel_bcast(efd)
    call parallel_bcast(eft)
    call parallel_bcast(efdt)
    call parallel_bcast(xf)
    call parallel_bcast(xfd)
    call parallel_bcast(xft)
    call parallel_bcast(xfdt)

    ! Construct the temperature and density deltas and their inverses
    do j = 1, jmax-1
       dt(j)   = t(j+1) - t(j)
       dt2(j)  = dt(j) * dt(j)
       dti(j)  = 1.0d0 / dt(j)
       dt2i(j) = 1.0d0 / dt2(j)
    end do

    do i = 1, imax-1
       dd(i)   = d(i+1) - d(i)
       dd2(i)  = dd(i) * dd(i)
       ddi(i)  = 1.0d0 / dd(i)
       dd2i(i) = 1.0d0 / dd2(i)
    end do

#ifdef AMREX_USE_ACC
    !$acc update device(d, t)
    !$acc update device(dt, dt2, dti, dt2i)
    !$acc update device(dd, dd2, ddi, dd2i)
    !$acc update device(f, fd, fdd, ft, ftt, fdt, fddt, fdtt, fddtt)
    !$acc update device(dpdf, dpdfd, dpdft, dpdfdt)
    !$acc update device(ef, efd, eft, efdt)
    !$acc update device(xf, xfd, xft, xfdt)
#endif

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target update to(d, t)
    !$omp target update to(dt, dt2, dti, dt2i)
    !$omp target update to(dd, dd2, ddi, dd2i)
    !$omp target update to(f, fd, fdd, ft, ftt, fdt, fddt, fdtt, fddtt)
    !$omp target update to(dpdf, dpdfd, dpdft, dpdfdt)
    !$omp target update to(ef, efd, eft, efdt)
    !$omp target update to(xf, xfd, xft, xfdt)
#endif

  end subroutine eos_init



  ! quintic hermite polynomial functions
  ! psi0 and its derivatives
  CASTRO_FORT_DEVICE pure function psi0(z) result(psi0r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: psi0r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    psi0r = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0

  end function psi0

  CASTRO_FORT_DEVICE pure function dpsi0(z) result(dpsi0r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: dpsi0r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    dpsi0r = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)

  end function dpsi0

  CASTRO_FORT_DEVICE pure function ddpsi0(z) result(ddpsi0r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: ddpsi0r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    ddpsi0r = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

  end function ddpsi0

  ! psi1 and its derivatives
  CASTRO_FORT_DEVICE pure function psi1(z) result(psi1r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: psi1r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    psi1r = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)

  end function psi1

  CASTRO_FORT_DEVICE pure function dpsi1(z) result(dpsi1r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: dpsi1r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    dpsi1r = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0

  end function dpsi1

  CASTRO_FORT_DEVICE pure function ddpsi1(z) result(ddpsi1r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: ddpsi1r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    ddpsi1r = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

  end function ddpsi1

  ! psi2  and its derivatives
  CASTRO_FORT_DEVICE pure function psi2(z) result(psi2r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: psi2r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    psi2r = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)

  end function psi2

  CASTRO_FORT_DEVICE pure function dpsi2(z) result(dpsi2r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: dpsi2r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    dpsi2r = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)

  end function dpsi2

  CASTRO_FORT_DEVICE pure function ddpsi2(z) result(ddpsi2r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: ddpsi2r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    ddpsi2r = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

  end function ddpsi2


  ! biquintic hermite polynomial function
  CASTRO_FORT_DEVICE pure function h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md) result(h5r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: fi(36)
    real(rt), intent(in) :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md
    real(rt) :: h5r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    h5r =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
         + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
         + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
         + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
         + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
         + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
         + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
         + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
         + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
         + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
         + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
         + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
         + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
         + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
         + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt

  end function h5


  ! cubic hermite polynomial functions
  ! psi0 & derivatives
  CASTRO_FORT_DEVICE pure function xpsi0(z) result(xpsi0r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: xpsi0r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    xpsi0r = z * z * (2.0d0*z - 3.0d0) + 1.0

  end function xpsi0

  CASTRO_FORT_DEVICE pure function xdpsi0(z) result(xdpsi0r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: xdpsi0r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    xdpsi0r = z * (6.0d0*z - 6.0d0)

  end function xdpsi0


  ! psi1 & derivatives
  CASTRO_FORT_DEVICE pure function xpsi1(z) result(xpsi1r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: xpsi1r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    xpsi1r = z * ( z * (z - 2.0d0) + 1.0d0)

  end function xpsi1

  CASTRO_FORT_DEVICE pure function xdpsi1(z) result(xdpsi1r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: z
    real(rt) :: xdpsi1r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    xdpsi1r = z * (3.0d0*z - 4.0d0) + 1.0d0

  end function xdpsi1

  ! bicubic hermite polynomial function
  CASTRO_FORT_DEVICE pure function h3(fi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) result(h3r)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    real(rt), intent(in) :: fi(36)
    real(rt), intent(in) :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md
    real(rt) :: h3r

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

    h3r =  fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
         + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
         + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
         + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
         + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

  end function h3



  subroutine eos_finalize() bind(C, name='eos_finalize')

    implicit none

    ! Deallocate managed module variables

    deallocate(d)
    deallocate(t)
    deallocate(f)
    deallocate(fd)
    deallocate(ft)
    deallocate(fdd)
    deallocate(ftt)
    deallocate(fdt)
    deallocate(fddt)
    deallocate(fdtt)
    deallocate(fddtt)
    deallocate(dpdf)
    deallocate(dpdfd)
    deallocate(dpdft)
    deallocate(dpdfdt)
    deallocate(ef)
    deallocate(efd)
    deallocate(eft)
    deallocate(efdt)
    deallocate(xf)
    deallocate(xfd)
    deallocate(xft)
    deallocate(xfdt)
    deallocate(dt)
    deallocate(dt2)
    deallocate(dti)
    deallocate(dt2i)
    deallocate(dd)
    deallocate(dd2)
    deallocate(ddi)
    deallocate(dd2i)

  end subroutine eos_finalize

end module eos_module
