module riemann_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                QPRES, QGAME, QREINT, QFS, &
                                QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                UFS, UFX, &
                                NGDNV, GDRHO, GDPRES, GDGAME, &
                                QC, QCSML, QGAMC, &
                                small_dens, small_temp, &
                                npassive, upass_map, qpass_map
  use advection_util_module, only: ht

  implicit none

  real(rt), parameter :: smallu = 1.e-12_rt

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine cmpflx(flx, flx_lo, flx_hi, &
                    qaux, qa_lo, qa_hi, &
                    h, idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use network, only: nspec, naux
    use amrex_fort_module, only: rt => amrex_real
    use bl_constants_module, only: ZERO, HALF, ONE

    integer,  intent(in   ) :: flx_lo(3), flx_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: idir,ilo,ihi,jlo,jhi,kc,kflux,k3d
    integer,  intent(in   ) :: domlo(3),domhi(3)

    type(ht), intent(inout) :: h

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)

    ! qaux come in dimensioned as the full box, so we use k3d here to
    ! index it in z

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! local variables

    integer      :: i, j
    integer      :: is_shock
    real(rt)     :: cl, cr
    type (eos_t) :: eos_state
    real(rt)     :: rhoInv

    if (idir == 1) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             h%smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i-1,j,k3d,QCSML) )
             h%cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i-1,j,k3d,QC) )
             h%gamcm(i,j) = qaux(i-1,j,k3d,QGAMC)
             h%gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    elseif (idir == 2) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             h%smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j-1,k3d,QCSML) )
             h%cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j-1,k3d,QC) )
             h%gamcm(i,j) = qaux(i,j-1,k3d,QGAMC)
             h%gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    else
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             h%smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j,k3d-1,QCSML) )
             h%cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j,k3d-1,QC) )
             h%gamcm(i,j) = qaux(i,j,k3d-1,QGAMC)
             h%gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    endif

    ! The following section is temporarily disabled when CUDA
    ! is used. It will be turned on later, once the hydro update
    ! is device-ready.

#ifndef CUDA
    ! recompute the thermodynamics on the interface to make it
    ! all consistent

    ! we want to take the edge states of rho, p, and X, and get
    ! new values for gamc and (rho e) on the edges that are
    ! thermodynamically consistent.
    
    do j = jlo, jhi
       do i = ilo, ihi

          ! this is an initial guess for iterations, since we
          ! can't be certain that temp is on interfaces
          eos_state % T = 10000.0e0_rt

          ! minus state
          eos_state % rho = h%qm(i,j,kc,QRHO,idir)
          eos_state % p   = h%qm(i,j,kc,QPRES,idir)
          eos_state % e   = h%qm(i,j,kc,QREINT,idir)/h%qm(i,j,kc,QRHO,idir)
          eos_state % xn  = h%qm(i,j,kc,QFS:QFS+nspec-1,idir)
          eos_state % aux = h%qm(i,j,kc,QFX:QFX+naux-1,idir)

          call eos(eos_input_re, eos_state)

          h%qm(i,j,kc,QREINT,idir) = eos_state % e * eos_state % rho
          h%qm(i,j,kc,QPRES,idir)  = eos_state % p
          h%gamcm(i,j)             = eos_state % gam1

       enddo
    enddo

    ! plus state
    do j = jlo, jhi
       do i = ilo, ihi
          rhoInv = ONE / h%qp(i,j,kc,QRHO,idir)

          eos_state % rho = h%qp(i,j,kc,QRHO,idir)
          eos_state % p   = h%qp(i,j,kc,QPRES,idir)
          eos_state % e   = h%qp(i,j,kc,QREINT,idir)/h%qp(i,j,kc,QRHO,idir)
          eos_state % xn  = h%qp(i,j,kc,QFS:QFS+nspec-1,idir) * rhoInv
          eos_state % aux = h%qp(i,j,kc,QFX:QFX+naux-1,idir) * rhoInv

          call eos(eos_input_re, eos_state)

          h%qp(i,j,kc,QREINT,idir) = eos_state % e * eos_state % rho
          h%qp(i,j,kc,QPRES,idir)  = eos_state % p
          h%gamcp(i,j)             = eos_state % gam1

       enddo
    enddo
#endif

    call riemannus(flx, flx_lo, flx_hi, &
                   h, idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

  end subroutine cmpflx



#ifdef CUDA
  attributes(device) &
#endif
  subroutine riemannus(uflx,uflx_lo,uflx_hi, &
                       h,idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use amrex_fort_module, only: rt => amrex_real
    use bl_constants_module, only: ZERO, HALF, ONE
    use prob_params_module, only: physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall

    real(rt), parameter :: small = 1.e-8_rt
    real(rt), parameter :: small_pres = 1.e-200_rt

    integer,  intent(in   ) :: uflx_lo(3),uflx_hi(3)
    integer,  intent(in   ) :: idir,ilo,ihi,jlo,jhi
    integer,  intent(in   ) :: domlo(3),domhi(3)

    type(ht), intent(inout) :: h
    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nqp, ipassive

    real(rt) :: regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, rhoetot, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall,qavg

    real(rt) :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt) :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    real(rt) :: wwinv, roinv, co2inv

    if (idir .eq. 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ

    else if (idir .eq. 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ

    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
    end if

    special_bnd_lo = (physbc_lo(idir) .eq. Symmetry &
         .or.         physbc_lo(idir) .eq. SlipWall &
         .or.         physbc_lo(idir) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) .eq. Symmetry &
         .or.         physbc_hi(idir) .eq. SlipWall &
         .or.         physbc_hi(idir) .eq. NoSlipWall)

    if (idir .eq. 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = ONE
    if (idir.eq.3) then
       if ( k3d .eq. domlo(3)   .and. special_bnd_lo .or. &
            k3d .eq. domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = ZERO
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = ONE
       if (idir .eq. 2) then
          if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
               j .eq. domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = ZERO
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(h%qm(i,j,kc,QRHO,idir), small_dens)

          ! pick left velocities based on direction
          ul  = h%qm(i,j,kc,iu,idir)
          v1l = h%qm(i,j,kc,iv1,idir)
          v2l = h%qm(i,j,kc,iv2,idir)

          pl  = max(h%qm(i,j,kc,QPRES ,idir), small_pres)
          rel =     h%qm(i,j,kc,QREINT,idir)
          rr = max(h%qp(i,j,kc,QRHO,idir), small_dens)

          ! pick right velocities based on direction
          ur  = h%qp(i,j,kc,iu,idir)
          v1r = h%qp(i,j,kc,iv1,idir)
          v2r = h%qp(i,j,kc,iv2,idir)

          pr  = max(h%qp(i,j,kc,QPRES,idir), small_pres)
          rer =     h%qp(i,j,kc,QREINT,idir)
          csmall = h%smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(h%gamcm(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(h%gamcp(i,j)*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar,small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = h%gamcm(i,j)
          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = h%gamcp(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(h%gamcm(i,j)+h%gamcp(i,j))
          endif

          ro = max(small_dens,ro)

          roinv = ONE/ro

          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall,co)
          co2inv = ONE/(co*co)

          drho = (pstar - po)*co2inv
          rstar = ro + drho
          rstar = max(small_dens,rstar)

          entho = (reo + po)*roinv*co2inv
          estar = reo + (pstar - po)*entho
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif

          if (spout-spin == ZERO) then
             scr = small*h%cavg(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar > ZERO) then
             h%qint(i,j,kc,iv1) = v1l
             h%qint(i,j,kc,iv2) = v2l
          else if (ustar < ZERO) then
             h%qint(i,j,kc,iv1) = v1r
             h%qint(i,j,kc,iv2) = v2r
          else
             h%qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             h%qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif
          h%qint(i,j,kc,GDRHO) = frac*rstar + (ONE - frac)*ro
          h%qint(i,j,kc,iu  ) = frac*ustar + (ONE - frac)*uo

          h%qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          regdnv = frac*estar + (ONE - frac)*reo
          if (spout < ZERO) then
             h%qint(i,j,kc,GDRHO) = ro
             h%qint(i,j,kc,iu  ) = uo
             h%qint(i,j,kc,GDPRES) = po
             regdnv = reo
          endif

          if (spin >= ZERO) then
             h%qint(i,j,kc,GDRHO) = rstar
             h%qint(i,j,kc,iu  ) = ustar
             h%qint(i,j,kc,GDPRES) = pstar
             regdnv = estar
          endif


          h%qint(i,j,kc,GDGAME) = h%qint(i,j,kc,GDPRES)/regdnv + ONE
          h%qint(i,j,kc,GDPRES) = max(h%qint(i,j,kc,GDPRES),small_pres)
          u_adv = h%qint(i,j,kc,iu)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z


          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = h%qint(i,j,kc,GDRHO)*u_adv

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*h%qint(i,j,kc,iu ) + h%qint(i,j,kc,GDPRES)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*h%qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*h%qint(i,j,kc,iv2)

          rhoetot = regdnv + HALF*h%qint(i,j,kc,GDRHO)*(h%qint(i,j,kc,iu)**2 + h%qint(i,j,kc,iv1)**2 + h%qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + h%qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = u_adv*regdnv
          ! store this for vectorization
          h%us1d(i) = ustar

       end do

       ! passively advected quantities
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (h%us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*h%qm(i,j,kc,nqp,idir)

             else if (h%us1d(i) < ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*h%qp(i,j,kc,nqp,idir)

             else
                qavg = HALF * (h%qm(i,j,kc,nqp,idir) + h%qp(i,j,kc,nqp,idir))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

  end subroutine riemannus

end module riemann_module
