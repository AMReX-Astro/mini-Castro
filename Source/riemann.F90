module riemann_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine compute_flux(lo, hi, &
                                             ql, ql_lo, ql_hi, &
                                             qr, qr_lo, qr_hi, &
                                             flx, flx_lo, flx_hi, &
                                             qint, q_lo, q_hi, &
                                             qgdnv, qg_lo, qg_hi, &
                                             qaux, qa_lo, qa_hi, &
                                             idir) bind(C, name="compute_flux")

    use amrex_constants_module, only: ZERO, HALF, ONE
    use castro_module, only: QVAR, QRHO, QU, QV, QW, QPRES, QC, QGAMC, QGAME, QFS, QREINT, &
                             NQAUX, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS, &
                             NGDNV, GDRHO, GDPRES, GDGAME, GDRHO, GDU, GDV, GDW, &
                             small, small_dens, smallu, small_pres
    use network, only: nspec

    implicit none

    ! Solve Riemann problem with the Colella, Glaz, and Ferguson solver.
    ! This is a 2-shock solver that uses a very simple approximation for the
    ! star state, and carries an auxiliary jump condition for (rho e) to
    ! deal with a real gas.

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)

    integer, intent(in), value :: idir

    real(rt), intent(in   ) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),QVAR)
    real(rt), intent(in   ) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),QVAR)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    integer :: i, j, k
    integer :: n, nqp

    real(rt) :: regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall, qavg
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, im1, im2, im3
    real(rt) :: wwinv, roinv, co2inv
    real(rt) :: u_adv, rhoeint, rhoetot

    integer :: ispec

    ! set integer pointers for the normal and transverse velocity and
    ! momentum

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ
    else if (idir == 2) then
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

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! ------------------------------------------------------------------
             ! set the left and right states for this interface
             ! ------------------------------------------------------------------

             rl = max(ql(i,j,k,QRHO), small_dens)

             ! pick left velocities based on direction
             ul  = ql(i,j,k,iu)
             v1l = ql(i,j,k,iv1)
             v2l = ql(i,j,k,iv2)
             pl  = max(ql(i,j,k,QPRES), small_pres)
             rel = ql(i,j,k,QREINT)

             rr = max(qr(i,j,k,QRHO), small_dens)

             ! pick right velocities based on direction
             ur  = qr(i,j,k,iu)
             v1r = qr(i,j,k,iv1)
             v2r = qr(i,j,k,iv2)
             pr  = max(qr(i,j,k,QPRES), small_pres)
             rer = qr(i,j,k,QREINT)

             ! ------------------------------------------------------------------
             ! estimate the star state: pstar, ustar
             ! ------------------------------------------------------------------

             if (idir == 1) then
                csmall = max( small, max( small * qaux(i,j,k,QC) , small * qaux(i-1,j,k,QC))  )
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i-1,j,k,QC))
                gamcl = qaux(i-1,j,k,QGAMC)
                gamcr = qaux(i,j,k,QGAMC)
             else if (idir == 2) then
                csmall = max( small, max( small * qaux(i,j,k,QC) , small * qaux(i,j-1,k,QC))  )
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i,j-1,k,QC))
                gamcl = qaux(i,j-1,k,QGAMC)
                gamcr = qaux(i,j,k,QGAMC)
             else
                csmall = max( small, max( small * qaux(i,j,k,QC) , small * qaux(i,j,k-1,QC))  )
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i,j,k-1,QC))
                gamcl = qaux(i,j,k-1,QGAMC)
                gamcr = qaux(i,j,k,QGAMC)
             end if

             wsmall = small_dens*csmall

             ! this is Castro I: Eq. 33
             wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
             wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

             wwinv = ONE/(wl + wr)
             pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
             ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

             pstar = max(pstar, small_pres)

             ! for symmetry preservation, if ustar is really small, then we
             ! set it to zero
             if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
                ustar = ZERO
             endif

             ! ------------------------------------------------------------------
             ! look at the contact to determine which region we are in
             ! ------------------------------------------------------------------

             ! this just determines which of the left or right states is still
             ! in play.  We still need to look at the other wave to determine
             ! if the star state or this state is on the interface.

             if (ustar > ZERO) then
                ro = rl
                uo = ul
                po = pl
                reo = rel
                gamco = gamcl
             else if (ustar < ZERO) then
                ro = rr
                uo = ur
                po = pr
                reo = rer
                gamco = gamcr
             else
                ro = HALF*(rl + rr)
                uo = HALF*(ul + ur)
                po = HALF*(pl + pr)
                reo = HALF*(rel + rer)
                gamco = HALF*(gamcl + gamcr)
             endif

             ro = max(small_dens, ro)

             roinv = ONE/ro

             co = sqrt(abs(gamco*po*roinv))
             co = max(csmall,co)
             co2inv = ONE/(co*co)

             ! we can already deal with the transverse velocities -- they
             ! only jump across the contact
             if (ustar > ZERO) then
                qint(i,j,k,iv1) = v1l
                qint(i,j,k,iv2) = v2l

             else if (ustar < ZERO) then
                qint(i,j,k,iv1) = v1r
                qint(i,j,k,iv2) = v2r

             else
                qint(i,j,k,iv1) = HALF*(v1l+v1r)
                qint(i,j,k,iv2) = HALF*(v2l+v2r)
             endif


             ! ------------------------------------------------------------------
             ! compute the rest of the star state
             ! ------------------------------------------------------------------

             drho = (pstar - po)*co2inv
             rstar = ro + drho
             rstar = max(small_dens, rstar)

             entho = (reo + po)*roinv*co2inv
             estar = reo + (pstar - po)*entho

             cstar = sqrt(abs(gamco*pstar/rstar))
             cstar = max(cstar, csmall)

             ! ------------------------------------------------------------------
             ! finish sampling the solution
             ! ------------------------------------------------------------------

             ! look at the remaining wave to determine if the star state or the
             ! 'o' state above is on the interface

             sgnm = sign(ONE, ustar)

             ! the values of u +/- c on either side of the non-contact
             ! wave
             spout = co - sgnm*uo
             spin = cstar - sgnm*ustar

             ! a simple estimate of the shock speed
             ushock = HALF*(spin + spout)

             if (pstar-po > ZERO) then
                spin = ushock
                spout = ushock
             endif

             if (spout-spin == ZERO) then
                scr = small*cavg
             else
                scr = spout-spin
             endif

             ! interpolate for the case that we are in a rarefaction
             frac = (ONE + (spout + spin)/scr)*HALF
             frac = max(ZERO, min(ONE, frac))

             qint(i,j,k,QRHO) = frac*rstar + (ONE - frac)*ro
             qint(i,j,k,iu  ) = frac*ustar + (ONE - frac)*uo

             qint(i,j,k,QPRES) = frac*pstar + (ONE - frac)*po
             regdnv = frac*estar + (ONE - frac)*reo

             ! as it stands now, we set things assuming that the rarefaction
             ! spans the interface.  We overwrite that here depending on the
             ! wave speeds

             ! look at the speeds on either side of the remaining wave
             ! to determine which region we are in
             if (spout < ZERO) then
                ! the l or r state is on the interface
                qint(i,j,k,QRHO) = ro
                qint(i,j,k,iu  ) = uo
                qint(i,j,k,QPRES) = po
                regdnv = reo
             endif

             if (spin >= ZERO) then
                ! the star state is on the interface
                qint(i,j,k,QRHO) = rstar
                qint(i,j,k,iu  ) = ustar
                qint(i,j,k,QPRES) = pstar
                regdnv = estar
             endif

             qint(i,j,k,QGAME) = qint(i,j,k,QPRES)/regdnv + ONE
             qint(i,j,k,QPRES) = max(qint(i,j,k,QPRES),small_pres)
             qint(i,j,k,QREINT) = regdnv

             ! passively advected quantities
             do ispec = 1, nspec

                n  = UFS + ispec - 1
                nqp = QFS + ispec - 1

                if (ustar > ZERO) then
                   qint(i,j,k,nqp) = ql(i,j,k,nqp)
                else if (ustar < ZERO) then
                   qint(i,j,k,nqp) = qr(i,j,k,nqp)
                else
                   qavg = HALF * (ql(i,j,k,nqp) + qr(i,j,k,nqp))
                   qint(i,j,k,nqp) = qavg
                end if

             end do

             ! Store results in the Godunov state

             qgdnv(i,j,k,GDRHO) = qint(i,j,k,QRHO)
             qgdnv(i,j,k,GDU) = qint(i,j,k,QU)
             qgdnv(i,j,k,GDV) = qint(i,j,k,QV)
             qgdnv(i,j,k,GDW) = qint(i,j,k,QW)
             qgdnv(i,j,k,GDPRES) = qint(i,j,k,QPRES)
             qgdnv(i,j,k,GDGAME) = qint(i,j,k,QGAME)

             ! Compute fluxes, order as conserved state (not q)

             u_adv = qint(i,j,k,iu)
             rhoeint = qint(i,j,k,QREINT)

             flx(i,j,k,URHO) = qint(i,j,k,QRHO) * u_adv

             flx(i,j,k,im1) = flx(i,j,k,URHO) * qint(i,j,k,iu)
             flx(i,j,k,im2) = flx(i,j,k,URHO) * qint(i,j,k,iv1)
             flx(i,j,k,im3) = flx(i,j,k,URHO) * qint(i,j,k,iv2)

             rhoetot = rhoeint + HALF * qint(i,j,k,QRHO) * &
                                 (qint(i,j,k,iu)**2 + &
                                  qint(i,j,k,iv1)**2 + &
                                  qint(i,j,k,iv2)**2)

             flx(i,j,k,UEDEN) = u_adv * (rhoetot + qint(i,j,k,QPRES))
             flx(i,j,k,UEINT) = u_adv * rhoeint

             flx(i,j,k,UTEMP) = ZERO

             ! passively advected quantities
             do ispec = 1, nspec
                n  = UFS + ispec - 1
                nqp = QFS + ispec - 1
                flx(i,j,k,n) = flx(i,j,k,URHO) * qint(i,j,k,nqp)
             end do

          end do
       end do
    end do

  end subroutine compute_flux

end module riemann_module
