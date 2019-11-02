module ppm_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, SIXTH, HALF, TWO3RD, ONE, TWO, THREE, SIX

  implicit none

contains

  subroutine ppm_reconstruct(s, flatn, sm, sp) bind(C, name='ppm_reconstruct')
    ! This routine does the reconstruction of the zone data into a parabola.

    implicit none

    real(rt), intent(in   ) :: s(-2:2), flatn
    real(rt), intent(inout) :: sm, sp

    ! local
    real(rt) :: dsl, dsr, dsc, dsvl_l, dsvl_r

    ! Compute van Leer slopes

    dsl = TWO * (s(-1) - s(-2))
    dsr = TWO * (s(0) - s(-1))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(0) - s(-2))
       dsvl_l = sign(ONE, dsc) * min(abs(dsc), abs(dsl), abs(dsr))
    else
       dsvl_l = ZERO
    end if

    dsl = TWO * (s(0) - s(-1))
    dsr = TWO * (s(1) - s(0))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(1) - s(-1))
       dsvl_r = sign(ONE, dsc) * min(abs(dsc), abs(dsl), abs(dsr))
    else
       dsvl_r = ZERO
    end if

    ! Interpolate s to edges

    sm = HALF*(s(0) + s(-1)) - SIXTH*(dsvl_r - dsvl_l)

    ! Make sure sedge lies in between adjacent cell-centered values

    sm = max(sm, min(s(0), s(-1)))
    sm = min(sm, max(s(0), s(-1)))



    ! Compute van Leer slopes

    dsl = TWO  * (s(0) - s(-1))
    dsr = TWO  * (s(1) - s(0))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(1) - s(-1))
       dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_l = ZERO
    end if

    dsl = TWO  * (s(1) - s(0))
    dsr = TWO  * (s(2) - s(1))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(2) - s(0))
       dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_r = ZERO
    end if

    ! Interpolate s to edges

    sp = HALF*(s(1) + s(0)) - SIXTH*(dsvl_r - dsvl_l)

    ! Make sure sedge lies in between adjacent cell-centered values

    sp = max(sp, min(s(1), s(0)))
    sp = min(sp, max(s(1), s(0)))



    ! Flatten the parabola

    sm = flatn * sm + (ONE - flatn) * s(0)
    sp = flatn * sp + (ONE - flatn) * s(0)

    ! Modify using quadratic limiters -- note this version of the limiting comes
    ! from Colella and Sekora (2008), not the original PPM paper.

    if ((sp - s(0)) * (s(0) - sm) .le. ZERO) then

       sp = s(0)
       sm = s(0)

    else if (abs(sp - s(0)) .ge. TWO * abs(sm - s(0))) then

       sp = THREE * s(0) - TWO * sm

    else if (abs(sm - s(0)) .ge. TWO * abs(sp - s(0))) then

       sm = THREE * s(0) - TWO * sp

    end if

  end subroutine ppm_reconstruct



  subroutine ppm_int_profile(sm, sp, sc, u, c, dtdx, Ip, Im) bind(C, name='ppm_int_profile')
    ! Integrate the parabolic profile to the edge of the cell.

    implicit none

    real(rt), intent(in   ) :: sm, sp, sc, u, c, dtdx
    real(rt), intent(inout) :: Ip(1:3), Im(1:3)

    ! local
    real(rt) :: speed, sigma, s6

    ! compute x-component of Ip and Im
    s6 = SIX * sc - THREE * (sm + sp)

    ! Ip/m is the integral under the parabola for the extent
    ! that a wave can travel over a timestep
    !
    ! Ip integrates to the right edge of a cell
    ! Im integrates to the left edge of a cell

    ! u-c wave
    speed = u - c
    sigma = abs(speed) * dtdx

    ! if speed == ZERO, then either branch is the same
    if (speed <= ZERO) then
       Ip(1) = sp
       Im(1) = sm + HALF*sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(1) = sp - HALF*sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(1) = sm
    endif

    ! u wave
    speed = u
    sigma = abs(speed) * dtdx

    if (speed <= ZERO) then
       Ip(2) = sp
       Im(2) = sm + HALF * sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(2) = sp - HALF * sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(2) = sm
    endif

    ! u+c wave
    speed = u + c
    sigma = abs(speed) * dtdx

    if (speed <= ZERO) then
       Ip(3) = sp
       Im(3) = sm + HALF * sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(3) = sp - HALF * sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(3) = sm
    endif

  end subroutine ppm_int_profile



  subroutine trace_ppm(lo, hi, &
                       idir, &
                       q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       flatn, f_lo, f_hi, &
                       qm, qm_lo, qm_hi, &
                       qp, qp_lo, qp_hi, &
                       vlo, vhi, domlo, domhi, &
                       dx, dt) bind(C, name='trace_ppm')

    use network, only: nspec
    use castro_module, only: QVAR, NQAUX, QRHO, QU, QV, QW, QC, QGAMC, QGAME, &
                             QREINT, QTEMP, QFS, QPRES, small_dens, small_pres

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) ::  flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)

    real(rt), intent(in) :: dt, dx(3)

    ! Local variables

    integer :: n, i, j, k, ispec

    logical :: reconstruct_state(QVAR)

    real(rt) :: hdt, dtdx

    real(rt) :: sm, sp

    real(rt) :: s(-2:2)
    real(rt) :: Ip(1:3,QVAR), Im(1:3,QVAR)
    real(rt) :: Ip_gc(1:3,1), Im_gc(1:3,1)

    integer :: QUN, QUT, QUTT

    real(rt) :: cc, csq
    real(rt) :: rho, un

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm

    real(rt) :: rho_ref, rho_ref_inv, un_ref, p_ref, rhoe_g_ref, h_g_ref
    real(rt) :: cc_ref, cc_ref_inv, csq_ref, gam_g_ref

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp

    hdt = HALF * dt
    dtdx = dt / dx(idir)

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! We don't need to reconstruct all of the QVAR state variables.
    reconstruct_state(:) = .true.
    reconstruct_state(QGAME) = .false.
    reconstruct_state(QTEMP) = .false.

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2

             un = q(i,j,k,QUN)

             ! do the parabolic reconstruction and compute the
             ! integrals under the characteristic waves
             do n = 1, QVAR
                if (.not. reconstruct_state(n)) cycle

                if (idir == 1) then
                   s(:) = q(i-2:i+2,j,k,n)
                else if (idir == 2) then
                   s(:) = q(i,j-2:j+2,k,n)
                else
                   s(:) = q(i,j,k-2:k+2,n)
                end if

                call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip(:,n), Im(:,n))

             end do


             if (idir == 1) then
                s(:) = qaux(i-2:i+2,j,k,QGAMC)
             else if (idir == 2) then
                s(:) = qaux(i,j-2:j+2,k,QGAMC)
             else
                s(:) = qaux(i,j,k-2:k+2,QGAMC)
             end if

             call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

             call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip_gc, Im_gc)

             ! do the passives separately
             do ispec = 1, nspec
                n = QFS + ispec - 1

                ! Plus state on face i
                if ((idir == 1 .and. i >= vlo(1)) .or. &
                    (idir == 2 .and. j >= vlo(2)) .or. &
                    (idir == 3 .and. k >= vlo(3))) then

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1 wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   qp(i,j,k,n) = Im(2,n)

                end if

                ! Minus state on face i+1
                if (idir == 1 .and. i <= vhi(1)) then
                   qm(i+1,j,k,n) = Ip(2,n)
                else if (idir == 2 .and. j <= vhi(2)) then
                   qm(i,j+1,k,n) = Ip(2,n)
                else if (idir == 3 .and. k <= vhi(3)) then
                   qm(i,j,k+1,n) = Ip(2,n)
                end if

             end do



             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(1,QRHO)
                un_ref    = Im(1,QUN)

                p_ref    = Im(1,QPRES)
                rhoe_g_ref = Im(1,QREINT)

                gam_g_ref  = Im_gc(1,1)

                rho_ref = max(rho_ref, small_dens)
                rho_ref_inv = ONE/rho_ref
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref*rho_ref_inv
                cc_ref = sqrt(csq_ref)
                cc_ref_inv = ONE/cc_ref
                h_g_ref = (p_ref + rhoe_g_ref)*rho_ref_inv

                ! *m are the jumps carried by un-c
                ! *p are the jumps carried by un+c

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)

                dum = un_ref - Im(1,QUN)
                dptotm = p_ref - Im(1,QPRES)

                drho = rho_ref - Im(2,QRHO)
                dptot = p_ref - Im(2,QPRES)
                drhoe_g = rhoe_g_ref - Im(2,QREINT)

                dup = un_ref - Im(3,QUN)
                dptotp = p_ref - Im(3,QPRES)

                ! (rho, u, p, (rho e) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv
                alphap = HALF*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv
                alpha0r = drho - dptot/csq_ref
                alpha0e_g = drhoe_g - dptot*h_g_ref/csq_ref

                if (un-cc > ZERO) then
                   alpham = ZERO
                else
                   alpham = -alpham
                end if

                if (un+cc > ZERO) then
                   alphap = ZERO
                else
                   alphap = -alphap
                end if

                if (un > ZERO) then
                   alpha0r = ZERO
                else
                   alpha0r = -alpha0r
                end if

                if (un > ZERO) then
                   alpha0e_g = ZERO
                else
                   alpha0e_g = -alpha0e_g
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                qp(i,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv
                qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g
                qp(i,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(2,QUT)
                qp(i,j,k,QUTT) = Im(2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(3,QRHO)
                un_ref    = Ip(3,QUN)

                p_ref    = Ip(3,QPRES)
                rhoe_g_ref = Ip(3,QREINT)

                gam_g_ref  = Ip_gc(3,1)

                rho_ref = max(rho_ref, small_dens)
                rho_ref_inv = ONE/rho_ref
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref*rho_ref_inv
                cc_ref = sqrt(csq_ref)
                cc_ref_inv = ONE/cc_ref
                h_g_ref = (p_ref + rhoe_g_ref)*rho_ref_inv

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(1,QUN)
                dptotm  = p_ref - Ip(1,QPRES)

                drho = rho_ref - Ip(2,QRHO)
                dptot = p_ref - Ip(2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(2,QREINT)

                dup = un_ref - Ip(3,QUN)
                dptotp = p_ref - Ip(3,QPRES)

                ! (rho, u, p, (rho e)) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv
                alphap = HALF*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv
                alpha0r = drho - dptot/csq_ref
                alpha0e_g = drhoe_g - dptot*h_g_ref/csq_ref

                if (un-cc > ZERO) then
                   alpham = -alpham
                else
                   alpham = ZERO
                end if

                if (un+cc > ZERO) then
                   alphap = -alphap
                else
                   alphap = ZERO
                end if

                if (un > ZERO) then
                   alpha0r = -alpha0r
                else
                   alpha0r = ZERO
                end if

                if (un > ZERO) then
                   alpha0e_g = -alpha0e_g
                else
                   alpha0e_g = ZERO
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   qm(i+1,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv
                   qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref)

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT)

                else if (idir == 2) then
                   qm(i,j+1,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv
                   qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref)

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT)

                else if (idir == 3) then
                   qm(i,j,k+1,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv
                   qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref)

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT)
                endif

             end if

          end do
       end do
    end do

  end subroutine trace_ppm



  subroutine ctu_ppm_states(lo, hi, &
                            vlo, vhi, &
                            q, qd_lo, qd_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            qxm, qxm_lo, qxm_hi, &
                            qxp, qxp_lo, qxp_hi, &
                            qym, qym_lo, qym_hi, &
                            qyp, qyp_lo, qyp_hi, &
                            qzm, qzm_lo, qzm_hi, &
                            qzp, qzp_lo, qzp_hi, &
                            dx, dt, &
                            domlo, domhi) bind(C, name="ctu_ppm_states")
    ! Compute the normal interface states by reconstructing
    ! the primitive variables using the piecewise parabolic method
    ! and doing characteristic tracing.  We do not apply the
    ! transverse terms here.

    use castro_module, only: QVAR, NVAR, &
                             QFS, QTEMP, QREINT, &
                             QC, QGAMC, NQAUX, QGAME, QREINT, &
                             NGDNV, GDU, GDV, GDW, GDPRES

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
    real(rt), intent(in) :: dx(3)   ! grid spacing in X, Y, Z direction
    real(rt), intent(in), value :: dt    ! time stepsize
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)   ! input state, primitives
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)   ! auxiliary hydro data
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))   ! flattening parameter

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), QVAR)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), QVAR)
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), QVAR)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), QVAR)
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), QVAR)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), QVAR)

    real(rt) :: hdt
    integer :: i, j, k, n, idir

    hdt = HALF * dt

    ! compute the interface states

    call trace_ppm(lo, hi, &
                   1, &
                   q, qd_lo, qd_hi, &
                   qaux, qa_lo, qa_hi, &
                   flatn, f_lo, f_hi, &
                   qxm, qxm_lo, qxm_hi, &
                   qxp, qxp_lo, qxp_hi, &
                   vlo, vhi, domlo, domhi, &
                   dx, dt)

    call trace_ppm(lo, hi, &
                   2, &
                   q, qd_lo, qd_hi, &
                   qaux, qa_lo, qa_hi, &
                   flatn, f_lo, f_hi, &
                   qym, qym_lo, qym_hi, &
                   qyp, qyp_lo, qyp_hi, &
                   vlo, vhi, domlo, domhi, &
                   dx, dt)

    call trace_ppm(lo, hi, &
                   3, &
                   q, qd_lo, qd_hi, &
                   qaux, qa_lo, qa_hi, &
                   flatn, f_lo, f_hi, &
                   qzm, qzm_lo, qzm_hi, &
                   qzp, qzp_lo, qzp_hi, &
                   vlo, vhi, domlo, domhi, &
                   dx, dt)

  end subroutine ctu_ppm_states

end module ppm_module
