module hydro_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, SIXTH, FOURTH, HALF, TWO3RD, ONE, TWO, THREE, SIX

  implicit none

  real(rt), parameter :: small = 1.e-8_rt
  real(rt), parameter :: small_pres = 1.e-200_rt
  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: dual_energy_eta1 = 1.e0_rt

contains

  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi) bind(c,name='ca_ctoprim')

    use network, only: nspec
    use eos_module, only: eos_t, eos_input_re, eos
    use castro_module, only: NVAR, URHO, UMX, UMZ, &
                             UEDEN, UEINT, UTEMP, &
                             QRHO, QU, QV, QW, UFS, &
                             QREINT, QPRES, QTEMP, QGAME, QFS, &
                             QVAR, QC, QGAMC, QDPDR, QDPDE, NQAUX, &
                             small_dens

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), parameter :: small = 1.e-8_rt

    integer  :: i, j, k, g
    integer  :: n, iq, ispec
    real(rt) :: kineng, rhoinv
    real(rt) :: vel(3)

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = uin(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Get the internal energy, which we'll use for
             ! determining the pressure.  We use a dual energy
             ! formalism. If (E - K) < eta1 and eta1 is suitably
             ! small, then we risk serious numerical truncation error
             ! in the internal energy.  Therefore we'll use the result
             ! of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv
             endif

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

             ! Load passively advected quatities into q
             do ispec = 1, nspec
                n  = UFS + ispec - 1
                iq = QFS + ispec - 1
                q(i,j,k,iq) = uin(i,j,k,n) * rhoinv
             enddo

             ! get gamc, p, T, c, csml using q state
             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
             q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde
             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs

          enddo
       enddo
    enddo

  end subroutine ca_ctoprim



  subroutine ppm_reconstruct(s, flatn, sm, sp)
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



  subroutine ppm_int_profile(sm, sp, sc, u, c, dtdx, Ip, Im)
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



  subroutine trace_ppm_species(i, j, k, &
                               idir, &
                               q, qd_lo, qd_hi, &
                               Ip, Im, &
                               qm, qm_lo, qm_hi, &
                               qp, qp_lo, qp_hi, &
                               vlo, vhi, domlo, domhi, &
                               dx, dt)
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only: nspec
    use castro_module, only: QVAR, QU, QV, QW, UFS

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)

    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)

    real(rt), intent(in) :: Ip(1:3,QVAR)
    real(rt), intent(in) :: Im(1:3,QVAR)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
    real(rt), intent(in) :: dt, dx(3)

    integer, intent(in) :: i, j, k

    integer :: ispec, n

    ! the passive stuff is the same regardless of the tracing
    do ispec = 1, nspec
       n = UFS + ispec - 1

       ! Plus state on face i
       if ((idir == 1 .and. i >= vlo(1)) .or. &
           (idir == 2 .and. j >= vlo(2)) .or. &
           (idir == 3 .and. k >= vlo(3))) then

          ! We have
          !
          ! q_l = q_ref - Proj{(q_ref - I)}
          !
          ! and Proj{} represents the characteristic projection.
          ! But for these, there is only 1-wave that matters, the u
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

  end subroutine trace_ppm_species



  subroutine trace_ppm_rhoe(lo, hi, &
                            idir, &
                            q, qd_lo, qd_hi, &
                            qaux, qa_lo, qa_hi, &
                            flatn, f_lo, f_hi, &
                            qm, qm_lo, qm_hi, &
                            qp, qp_lo, qp_hi, &
                            vlo, vhi, domlo, domhi, &
                            reconstruct_state, &
                            dx, dt)

    use network, only: nspec
    use castro_module, only: QVAR, NQAUX, QRHO, QU, QV, QW, &
                             QREINT, QPRES, QGAME, QC, QGAMC, &
                             small_dens

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

    logical, intent(in) :: reconstruct_state(QVAR)

    ! Local variables
    integer :: i, j, k, n

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
             call trace_ppm_species(i, j, k, &
                                    idir, &
                                    q, qd_lo, qd_hi, &
                                    Ip, Im, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    vlo, vhi, domlo, domhi, &
                                    dx, dt)

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


  end subroutine trace_ppm_rhoe



  subroutine trace_ppm(lo, hi, &
                       idir, &
                       q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       flatn, f_lo, f_hi, &
                       qm, qm_lo, qm_hi, &
                       qp, qp_lo, qp_hi, &
                       vlo, vhi, domlo, domhi, &
                       dx, dt)
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only: nspec
    use castro_module, only: QVAR, NQAUX, QU, QV, QW, QGAME, QREINT, QTEMP

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

    integer :: n, i, j, k

    logical :: reconstruct_state(QVAR)

    ! we don't need to reconstruct all of the QVAR state variables
    reconstruct_state(:) = .true.
    reconstruct_state(QGAME) = .false.
    reconstruct_state(QTEMP) = .false.

    call trace_ppm_rhoe(lo, hi, &
                        idir, &
                        q, qd_lo, qd_hi, &
                        qaux, qa_lo, qa_hi, &
                        flatn, f_lo, f_hi, &
                        qm, qm_lo, qm_hi, &
                        qp, qp_lo, qp_hi, &
                        vlo, vhi, domlo, domhi, &
                        reconstruct_state, &
                        dx, dt)

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

    hdt = HALF*dt

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



  subroutine compute_flux_q(lo, hi, &
                            qint, q_lo, q_hi, &
                            F, F_lo, F_hi, &
                            idir)
    ! given a primitive state, compute the flux in direction idir

    use castro_module, only : QVAR, NVAR, NQAUX, &
                              URHO, UMX, UMY, UMZ, &
                              UEDEN, UEINT, UTEMP, UFS, &
                              QRHO, QU, QV, QW, &
                              QPRES, QGAME, QREINT, &
                              QC, QGAMC, QFS
    use network, only: nspec

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: F_lo(3), F_hi(3)

    real(rt), intent(in) :: qint(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), QVAR)
    real(rt), intent(out) :: F(F_lo(1):F_hi(1), F_lo(2):F_hi(2), F_lo(3):F_hi(3), NVAR)

    integer :: iu, iv1, iv2, im1, im2, im3
    integer :: g, n, ispec, nqp
    real(rt) :: u_adv, rhoetot, rhoeint
    integer :: i, j, k

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

             u_adv = qint(i,j,k,iu)

             rhoeint = qint(i,j,k,QREINT)

             ! Compute fluxes, order as conserved state (not q)
             F(i,j,k,URHO) = qint(i,j,k,QRHO)*u_adv

             F(i,j,k,im1) = F(i,j,k,URHO)*qint(i,j,k,iu)
             F(i,j,k,im2) = F(i,j,k,URHO)*qint(i,j,k,iv1)
             F(i,j,k,im3) = F(i,j,k,URHO)*qint(i,j,k,iv2)

             rhoetot = rhoeint + &
                  HALF*qint(i,j,k,QRHO)*(qint(i,j,k,iu)**2 + &
                  qint(i,j,k,iv1)**2 + &
                  qint(i,j,k,iv2)**2)

             F(i,j,k,UEDEN) = u_adv*(rhoetot + qint(i,j,k,QPRES))
             F(i,j,k,UEINT) = u_adv*rhoeint

             F(i,j,k,UTEMP) = ZERO

             ! passively advected quantities
             do ispec = 1, nspec
                n  = UFS + ispec - 1
                nqp = QFS + ispec - 1
                F(i,j,k,n) = F(i,j,k,URHO)*qint(i,j,k,nqp)
             end do

          end do
       end do
    end do

  end subroutine compute_flux_q



  subroutine ca_store_godunov_state(lo, hi, &
                                    qint, qi_lo, qi_hi, &
                                    qgdnv, qg_lo, qg_hi) bind(C, name="ca_store_godunov_state")
    ! this copies the full interface state (NQ -- one for each primitive
    ! variable) over to a smaller subset of size NGDNV for use later in the
    ! hydro advancement.

    use castro_module, only: QVAR, NVAR, NQAUX, &
                             URHO, &
                             QRHO, QU, QV, QW, &
                             QPRES, QGAME, &
                             NGDNV, GDRHO, GDPRES, GDGAME, &
                             GDRHO, GDU, GDV, GDW

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: qi_lo(3), qi_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)
    real(rt), intent(in) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), QVAR)
    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! the hybrid routine uses the Godunov indices, not the full NQ state
             qgdnv(i,j,k,GDRHO) = qint(i,j,k,QRHO)
             qgdnv(i,j,k,GDU) = qint(i,j,k,QU)
             qgdnv(i,j,k,GDV) = qint(i,j,k,QV)
             qgdnv(i,j,k,GDW) = qint(i,j,k,QW)
             qgdnv(i,j,k,GDPRES) = qint(i,j,k,QPRES)
             qgdnv(i,j,k,GDGAME) = qint(i,j,k,QGAME)

          end do
       end do
    end do

  end subroutine ca_store_godunov_state


  
  subroutine cmpflx_plus_godunov(lo, hi, &
                                 qm, qm_lo, qm_hi, &
                                 qp, qp_lo, qp_hi, nc, comp, &
                                 flx, flx_lo, flx_hi, &
                                 qint, q_lo, q_hi, &
                                 qgdnv, qg_lo, qg_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 idir, domlo, domhi) bind(C, name="cmpflx_plus_godunov")

    use castro_module, only: NVAR, QVAR, NQAUX, NGDNV
    
    implicit none

    ! note: lo, hi necessarily the limits of the valid (no ghost
    ! cells) domain, but could be hi+1 in some dimensions.  We rely on
    ! the caller to specify the interfaces over which to solve the
    ! Riemann problems

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)

    integer, intent(in), value :: idir

    integer, intent(in) :: domlo(3),domhi(3)
    integer, intent(in), value :: nc, comp

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR,nc)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR,nc)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    call cmpflx(lo, hi, &
                qm, qm_lo, qm_hi, &
                qp, qp_lo, qp_hi, nc, comp, &
                flx, flx_lo, flx_hi, &
                qint, q_lo, q_hi, &
                qaux, qa_lo, qa_hi, &
                idir, domlo, domhi)

    call ca_store_godunov_state(lo, hi, &
                                qint, q_lo, q_hi, &
                                qgdnv, qg_lo, qg_hi)

  end subroutine cmpflx_plus_godunov

  subroutine cmpflx(lo, hi, &
                    qm, qm_lo, qm_hi, &
                    qp, qp_lo, qp_hi, nc, comp, &
                    flx, flx_lo, flx_hi, &
                    qint, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    idir, domlo, domhi)

    use castro_module, only: NVAR, QVAR, NQAUX

    implicit none

    ! note: lo, hi necessarily the limits of the valid (no ghost
    ! cells) domain, but could be hi+1 in some dimensions.  We rely on
    ! the caller to specific the interfaces over which to solve the
    ! Riemann problems

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: idir

    integer, intent(in) :: domlo(3),domhi(3)
    integer, intent(in) :: nc, comp

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR,nc)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR,nc)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! Solve Riemann problem to get the fluxes

    call riemann_state(qm, qm_lo, qm_hi, &
                       qp, qp_lo, qp_hi, nc, comp, &
                       qint, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       idir, lo, hi, &
                       domlo, domhi)

    call compute_flux_q(lo, hi, &
                        qint, q_lo, q_hi, &
                        flx, flx_lo, flx_hi, &
                        idir)

  end subroutine cmpflx




  subroutine riemann_state(qm, qm_lo, qm_hi, &
                           qp, qp_lo, qp_hi, nc, comp, &
                           qint, q_lo, q_hi, &
                           qaux, qa_lo, qa_hi, &
                           idir, lo, hi, domlo, domhi)
    ! just compute the hydrodynamic state on the interfaces
    ! don't compute the fluxes

    use castro_module, only: QVAR, NQAUX

    implicit none

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: idir
    ! note: lo, hi are not necessarily the limits of the valid (no
    ! ghost cells) domain, but could be hi+1 in some dimensions.  We
    ! rely on the caller to specific the interfaces over which to
    ! solve the Riemann problems
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: nc, comp

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR,nc)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR,nc)

    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! Solve Riemann problem

    call riemannus(qm, qm_lo, qm_hi, &
                   qp, qp_lo, qp_hi, nc, comp, &
                   qaux, qa_lo, qa_hi, &
                   qint, q_lo, q_hi, &
                   idir, lo, hi, &
                   domlo, domhi)

  end subroutine riemann_state



  subroutine riemannus(ql, ql_lo, ql_hi, &
                       qr, qr_lo, qr_hi, nc, comp, &
                       qaux, qa_lo, qa_hi, &
                       qint, q_lo, q_hi, &
                       idir, lo, hi, &
                       domlo, domhi)
    ! Colella, Glaz, and Ferguson solver
    !
    ! this is a 2-shock solver that uses a very simple approximation for the
    ! star state, and carries an auxiliary jump condition for (rho e) to
    ! deal with a real gas

    use castro_module, only: QVAR, QRHO, QU, QV, QW, QPRES, QC, QGAMC, QGAME, QFS, QREINT, &
                             NQAUX, URHO, UMX, UMY, UMZ, UFS, small_dens
    use network, only: nspec

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3),domhi(3)
    integer, intent(in) :: nc, comp

    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),QVAR,nc)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),QVAR,nc)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

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

             rl = max(ql(i,j,k,QRHO,comp), small_dens)

             ! pick left velocities based on direction
             ul  = ql(i,j,k,iu,comp)
             v1l = ql(i,j,k,iv1,comp)
             v2l = ql(i,j,k,iv2,comp)
             pl  = max(ql(i,j,k,QPRES,comp), small_pres)
             rel = ql(i,j,k,QREINT,comp)

             rr = max(qr(i,j,k,QRHO,comp), small_dens)

             ! pick right velocities based on direction
             ur  = qr(i,j,k,iu,comp)
             v1r = qr(i,j,k,iv1,comp)
             v2r = qr(i,j,k,iv2,comp)
             pr  = max(qr(i,j,k,QPRES,comp), small_pres)
             rer = qr(i,j,k,QREINT,comp)

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
                   qint(i,j,k,nqp) = ql(i,j,k,nqp,comp)
                else if (ustar < ZERO) then
                   qint(i,j,k,nqp) = qr(i,j,k,nqp,comp)
                else
                   qavg = HALF * (ql(i,j,k,nqp,comp) + qr(i,j,k,nqp,comp))
                   qint(i,j,k,nqp) = qavg
                end if

             end do

          end do
       end do
    end do

  end subroutine riemannus



  subroutine apply_av(lo, hi, idir, dx, &
                      div, div_lo, div_hi, &
                      uin, uin_lo, uin_hi, &
                      flux, f_lo, f_hi) bind(c, name="apply_av")

    use castro_module, only: NVAR, UTEMP

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1

    real(rt), parameter :: difmag = 0.1d0

    do n = 1, NVAR

       if ( n == UTEMP ) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i,j+1,k) + &
                                    div(i,j,k+1) + div(i,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i+1,j,k) + &
                                    div(i,j,k+1) + div(i+1,j,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j-1,k,n))

                else

                   div1 = FOURTH * (div(i,j  ,k) + div(i+1,j  ,k) + &
                                    div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-1,n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1

             end do
          end do
       end do

    end do

  end subroutine apply_av

  
  subroutine normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")
    ! Normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179.

    use network, only: nspec
    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: NVAR, URHO, UFS

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n
    real(rt) :: sum, fac

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             sum = ZERO

             do n = UFS, UFS+nspec-1
                sum = sum + flux(i,j,k,n)
             end do

             if (sum .ne. ZERO) then
                fac = flux(i,j,k,URHO) / sum
             else
                fac = ONE
             end if

             do n = UFS, UFS+nspec-1
                flux(i,j,k,n) = flux(i,j,k,n) * fac
             end do

          end do
       end do
    end do

  end subroutine normalize_species_fluxes



  subroutine calc_pdivu(lo, hi, &
                        q1, q1_lo, q1_hi, &
                        area1, a1_lo, a1_hi, &
                        q2, q2_lo, q2_hi, &
                        area2, a2_lo, a2_hi, &
                        q3, q3_lo, q3_hi, &
                        area3, a3_lo, a3_hi, &
                        vol, v_lo, v_hi, &
                        dx, pdivu, div_lo, div_hi)
    ! this computes the cell-centered p div(U) term from the
    ! edge-centered Godunov state.  This is used in the internal energy
    ! update

    use castro_module, only: NGDNV, GDPRES, GDU, GDV, GDW
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: pdivu(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: a1_lo(3), a1_hi(3)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt), intent(in) :: area1(a1_lo(1):a1_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: a2_lo(3), a2_hi(3)
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(in) :: area2(a2_lo(1):a2_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: a3_lo(3), a3_hi(3)
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt), intent(in) :: area3(a3_lo(1):a3_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: v_lo(3), v_hi(3)
    real(rt), intent(in) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             pdivu(i,j,k) = &
                  HALF*(q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                  (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                  HALF*(q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                  (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                  HALF*(q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)

          enddo
       enddo
    enddo

  end subroutine calc_pdivu

  

  subroutine ctu_consup(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q, q_lo, q_hi, &
                        update, updt_lo, updt_hi, &
                        flux1, flux1_lo, flux1_hi, &
                        flux2, flux2_lo, flux2_hi, &
                        flux3, flux3_lo, flux3_hi, &
                        qx, qx_lo, qx_hi, &
                        qy, qy_lo, qy_hi, &
                        qz, qz_lo, qz_hi, &
                        area1, area1_lo, area1_hi, &
                        area2, area2_lo, area2_hi, &
                        area3, area3_lo, area3_hi, &
                        vol, vol_lo, vol_hi, &
                        pdivu, pdivu_lo, pdivu_hi, &
                        dx, dt) bind(C, name="ctu_consup")

    use castro_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                             UEINT, UTEMP, NGDNV, QVAR, GDPRES

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)
    integer, intent(in) ::   pdivu_lo(3),   pdivu_hi(3)

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)
    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: pdivu(pdivu_lo(1):pdivu_hi(1),pdivu_lo(2):pdivu_hi(2),pdivu_lo(3):pdivu_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    integer :: i, j, g, k, n
    integer :: domlo(3), domhi(3)
    real(rt) :: volInv

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
                    qz, qz_lo, qz_hi, &
                    area3, area3_lo, area3_hi, &
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, pdivu_lo, pdivu_hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine ctu_consup


  subroutine scale_flux(lo, hi, &
                        flux, f_lo, f_hi, &
                        area, a_lo, a_hi, dt) bind(c, name="scale_flux")

    use castro_module, only: NVAR, GDPRES, UMX, NGDNV

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: dt

    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k, n

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k,n) = dt * flux(i,j,k,n) * area(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine scale_flux

end module hydro_module
