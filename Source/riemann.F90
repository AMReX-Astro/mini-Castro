module riemann_module

  use bl_types
  use bl_constants_module
  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                 UFS, UFX, &
                                 NGDNV, GDRHO, GDPRES, GDGAME, &
                                 QC, QCSML, QGAMC, &
                                 small_dens, small_temp, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, hybrid_riemann, &
                                 allow_negative_energy
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public cmpflx, shock

  real(rt)        , parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    shk, s_lo, s_hi, &
                    idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use eos_module
    use network, only: nspec, naux

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)

    integer, intent(in) :: idir,ilo,ihi,jlo,jhi,kc,kflux,k3d
    integer, intent(in) :: domlo(3),domhi(3)

    ! note: qm, qp, q come in as planes (all of x,y
    ! zones but only 2 elements in the z dir) instead of being
    ! dimensioned the same as the full box.  We index these with kc
    ! flux either comes in as planes (like qm, qp, ... above), or
    ! comes in dimensioned as the full box.  We index the flux with
    ! kflux -- this will be set correctly for the different cases.

    real(rt)        , intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)        , intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt)        , intent(inout) ::    flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt)        , intent(inout) ::   qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! qaux come in dimensioned as the full box, so we use k3d here to
    ! index it in z

    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)        , intent(in) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    ! local variables

    integer i, j
    integer :: gd_lo(2), gd_hi(2)
    real(rt)        , pointer :: smallc(:,:), cavg(:,:)
    real(rt)        , pointer :: gamcm(:,:), gamcp(:,:)
    integer :: is_shock
    real(rt)         :: cl, cr
    type (eos_t) :: eos_state

    real(rt)         :: rhoInv

    gd_lo = (/ ilo, jlo /)
    gd_hi = (/ ihi, jhi /)

    call bl_allocate ( smallc, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (   cavg, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcm, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcp, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    if (idir == 1) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i-1,j,k3d,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i-1,j,k3d,QC) )
             gamcm(i,j) = qaux(i-1,j,k3d,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    elseif (idir == 2) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j-1,k3d,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j-1,k3d,QC) )
             gamcm(i,j) = qaux(i,j-1,k3d,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    else
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j,k3d-1,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j,k3d-1,QC) )
             gamcm(i,j) = qaux(i,j,k3d-1,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
          enddo
       enddo
    endif

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
          eos_state % rho = qm(i,j,kc,QRHO)
          eos_state % p   = qm(i,j,kc,QPRES)
          eos_state % e   = qm(i,j,kc,QREINT)/qm(i,j,kc,QRHO)
          eos_state % xn  = qm(i,j,kc,QFS:QFS+nspec-1)
          eos_state % aux = qm(i,j,kc,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qm(i,j,kc,QREINT) = eos_state % e * eos_state % rho
          qm(i,j,kc,QPRES)  = eos_state % p
          gamcm(i,j)        = eos_state % gam1

       enddo
    enddo

    ! plus state
    do j = jlo, jhi
       do i = ilo, ihi
          rhoInv = ONE / qp(i,j,kc,QRHO)

          eos_state % rho = qp(i,j,kc,QRHO)
          eos_state % p   = qp(i,j,kc,QPRES)
          eos_state % e   = qp(i,j,kc,QREINT)/qp(i,j,kc,QRHO)
          eos_state % xn  = qp(i,j,kc,QFS:QFS+nspec-1) * rhoInv
          eos_state % aux = qp(i,j,kc,QFX:QFX+naux-1) * rhoInv

          call eos(eos_input_re, eos_state)

          qp(i,j,kc,QREINT) = eos_state % e * eos_state % rho
          qp(i,j,kc,QPRES)  = eos_state % p
          gamcp(i,j)        = eos_state % gam1

       enddo
    enddo

    call riemannus(qm, qp, qpd_lo, qpd_hi, &
                   gamcm, gamcp, cavg, smallc, gd_lo, gd_hi, &
                   flx, flx_lo, flx_hi, &
                   qint, q_lo, q_hi, &
                   idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)


    call bl_deallocate(smallc)
    call bl_deallocate(  cavg)
    call bl_deallocate( gamcm)
    call bl_deallocate( gamcp)
  end subroutine cmpflx


  subroutine shock(q,qd_lo,qd_hi,shk,s_lo,s_hi,lo,hi,dx)

    use prob_params_module, only : coord_type
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt)        , intent(in) :: dx(3)
    real(rt)        , intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)        , intent(inout) :: shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: i, j, k

    real(rt)         :: dxinv, dyinv, dzinv
    real(rt)         :: divU
    real(rt)         :: px_pre, px_post, py_pre, py_post, pz_pre, pz_post
    real(rt)         :: e_x, e_y, e_z, d
    real(rt)         :: p_pre, p_post, pjump

    real(rt)        , parameter :: small = 1.e-10_rt
    real(rt)        , parameter :: eps = 0.33e0_rt

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! construct div{U}
             divU = HALF*(q(i+1,j,k,QU) - q(i-1,j,k,QU))*dxinv + &
                    HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))*dyinv + &
                    HALF*(q(i,j,k+1,QW) - q(i,j,k-1,QW))*dzinv

             ! find the pre- and post-shock pressures in each direction
             if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < ZERO) then
                px_pre  = q(i+1,j,k,QPRES)
                px_post = q(i-1,j,k,QPRES)
             else
                px_pre  = q(i-1,j,k,QPRES)
                px_post = q(i+1,j,k,QPRES)
             endif

             if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < ZERO) then
                py_pre  = q(i,j+1,k,QPRES)
                py_post = q(i,j-1,k,QPRES)
             else
                py_pre  = q(i,j-1,k,QPRES)
                py_post = q(i,j+1,k,QPRES)
             endif

             if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < ZERO) then
                pz_pre  = q(i,j,k+1,QPRES)
                pz_post = q(i,j,k-1,QPRES)
             else
                pz_pre  = q(i,j,k-1,QPRES)
                pz_post = q(i,j,k+1,QPRES)
             endif

             ! use compression to create unit vectors for the shock direction
             e_x = (q(i+1,j,k,QU) - q(i-1,j,k,QU))**2
             e_y = (q(i,j+1,k,QV) - q(i,j-1,k,QV))**2
             e_z = (q(i,j,k+1,QW) - q(i,j,k-1,QW))**2
             d = ONE/(e_x + e_y + e_z + small)

             e_x = e_x*d
             e_y = e_y*d
             e_z = e_z*d

             ! project the pressures onto the shock direction
             p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre
             p_post = e_x*px_post + e_y*py_post + e_z*pz_post

             ! test for compression + pressure jump to flag a shock
             if (p_pre == ZERO) then
                ! this can happen if U = 0, so e_x, ... = 0
                pjump = ZERO
             else
                pjump = eps - (p_post - p_pre)/p_pre
             endif

             if (pjump < ZERO .and. divU < ZERO) then
                shk(i,j,k) = ONE
             else
                shk(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine shock



  subroutine riemannus(ql,qr,qpd_lo,qpd_hi, &
                       gamcl,gamcr,cav,smallc,gd_lo,gd_hi, &
                       uflx,uflx_lo,uflx_hi, &
                       qint,q_lo,q_hi, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use amrex_fort_module, only : rt => amrex_real

    real(rt), parameter :: small = 1.e-8_rt
    real(rt), parameter :: small_pres = 1.e-200_rt

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

    real(rt)         :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)         :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt)         ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt)         ::    qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nqp, ipassive

    real(rt)         :: regdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, rhoetot, scr
    real(rt)         :: rstar, cstar, estar, pstar, ustar
    real(rt)         :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt)         :: sgnm, spin, spout, ushock, frac
    real(rt)         :: wsmall, csmall,qavg

    real(rt)        , pointer :: us1d(:)

    real(rt)         :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt)         :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    real(rt)         :: wwinv, roinv, co2inv

    call bl_allocate(us1d,ilo,ihi)

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

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel =     ql(i,j,kc,QREINT)
          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer =     qr(i,j,kc,QREINT)
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

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
             gamco = gamcl(i,j)
          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
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
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar > ZERO) then
             qint(i,j,kc,iv1) = v1l
             qint(i,j,kc,iv2) = v2l
          else if (ustar < ZERO) then
             qint(i,j,kc,iv1) = v1r
             qint(i,j,kc,iv2) = v2r
          else
             qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif
          qint(i,j,kc,GDRHO) = frac*rstar + (ONE - frac)*ro
          qint(i,j,kc,iu  ) = frac*ustar + (ONE - frac)*uo

          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          regdnv = frac*estar + (ONE - frac)*reo
          if (spout < ZERO) then
             qint(i,j,kc,GDRHO) = ro
             qint(i,j,kc,iu  ) = uo
             qint(i,j,kc,GDPRES) = po
             regdnv = reo
          endif

          if (spin >= ZERO) then
             qint(i,j,kc,GDRHO) = rstar
             qint(i,j,kc,iu  ) = ustar
             qint(i,j,kc,GDPRES) = pstar
             regdnv = estar
          endif


          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE
          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES),small_pres)
          u_adv = qint(i,j,kc,iu)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z


          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*u_adv

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu ) + qint(i,j,kc,GDPRES)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

          rhoetot = regdnv + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = u_adv*regdnv
          ! store this for vectorization
          us1d(i) = ustar

       end do

       ! passively advected quantities
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nqp)

             else if (us1d(i) < ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nqp)

             else
                qavg = HALF * (ql(i,j,kc,nqp) + qr(i,j,kc,nqp))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call bl_deallocate(us1d)

  end subroutine riemannus



  pure function bc_test(idir, i, j, domlo, domhi) result (f)

    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall

    integer, intent(in) :: idir, i, j, domlo(*), domhi(*)
    integer :: f

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    f = 1

    if (idir == 1) then
       if (i == domlo(1) .and. &
            (physbc_lo(1) == Symmetry .or. &
             physbc_lo(1) == SlipWall .or. &
             physbc_lo(1) == NoSlipWall) ) then
          f = 0
       endif

       if (i == domhi(1)+1 .and. &
            (physbc_hi(1) == Symmetry .or. &
             physbc_hi(1) == SlipWall .or. &
             physbc_hi(1) == NoSlipWall) ) then
          f = 0
       endif
    end if

    if (idir == 2) then
       if (j == domlo(2) .and. &
            (physbc_lo(2) == Symmetry .or. &
             physbc_lo(2) == SlipWall .or. &
             physbc_lo(2) == NoSlipWall) ) then
          f = 0
       endif

       if (j == domhi(2)+1 .and. &
            (physbc_hi(2) == Symmetry .or. &
             physbc_hi(2) == SlipWall .or. &
             physbc_hi(2) == NoSlipWall) ) then
          f = 0
       end if
    endif

  end function bc_test


  pure subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    ! compute the lagrangian wave speeds.

    real(rt)        , intent(in) :: p,v,gam,gdot,pstar,csq,gmin,gmax
    real(rt)        , intent(out) :: wsq, gstar

    real(rt)        , parameter :: smlp1 = 1.e-10_rt
    real(rt)        , parameter :: small = 1.e-7_rt

    real(rt)         :: alpha, beta

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar = (pstar-p)*gdot/(pstar+p) + gam
    gstar = max(gmin, min(gmax, gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! this is CG Eq. 34
    alpha = pstar - (gstar-ONE)*p/(gam-ONE)
    if (alpha == ZERO) alpha = smlp1*(pstar + p)

    beta = pstar + HALF*(gstar-ONE)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq = max(wsq, (HALF*(gam-ONE)/gam)*csq)

    return
  end subroutine wsqge


  pure subroutine pstar_bisection(pstar_lo, pstar_hi, &
                                  ul, pl, taul, gamel, clsql, &
                                  ur, pr, taur, gamer, clsqr, &
                                  gdot, gmin, gmax, &
                                  pstar, gamstar, converged, pstar_hist_extra)

    ! we want to zero                                                                     
    ! f(p*) = u*_l(p*) - u*_r(p*)                                                         
    ! we'll do bisection                                                                  
                                                  
    use meth_params_module, only : cg_maxiter, cg_tol

    real(rt)        , intent(inout) :: pstar_lo, pstar_hi
    real(rt)        , intent(in) :: ul, pl, taul, gamel, clsql
    real(rt)        , intent(in) :: ur, pr, taur, gamer, clsqr
    real(rt)        , intent(in) :: gdot, gmin, gmax
    real(rt)        , intent(out) :: pstar, gamstar
    logical, intent(out) :: converged
    real(rt)        , intent(out) :: pstar_hist_extra(:)

    real(rt)         :: pstar_c, ustar_l, ustar_r, f_lo, f_hi, f_c
    real(rt)         :: wl, wr, wlsq, wrsq

    integer :: iter


    ! lo bounds
    call wsqge(pl, taul, gamel, gdot,  &
               gamstar, pstar_lo, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
               gamstar, pstar_lo, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_lo - pstar)*wl
    ustar_r = ur + (pstar_lo - pstar)*wr

    f_lo = ustar_l - ustar_r


    ! hi bounds
    call wsqge(pl, taul, gamel, gdot,  &
               gamstar, pstar_hi, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
               gamstar, pstar_hi, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_hi - pstar)*wl
    ustar_r = ur + (pstar_hi - pstar)*wr

    f_hi = ustar_l - ustar_r

    ! bisection
    iter = 1
    do while (iter <= cg_maxiter)

       pstar_c = HALF * (pstar_lo + pstar_hi)

       pstar_hist_extra(iter) = pstar_c

       call wsqge(pl, taul, gamel, gdot,  &
                  gamstar, pstar_c, wlsq, clsql, gmin, gmax)

       call wsqge(pr, taur, gamer, gdot,  &
                  gamstar, pstar_c, wrsq, clsqr, gmin, gmax)

       wl = ONE / sqrt(wlsq)
       wr = ONE / sqrt(wrsq)

       ustar_l = ul - (pstar_c - pl)*wl
       ustar_r = ur - (pstar_c - pr)*wr

       f_c = ustar_l - ustar_r

       if ( HALF * abs(pstar_lo - pstar_hi) < cg_tol * pstar_c ) then
          converged = .true.
          exit
       endif

       if (f_lo * f_c < ZERO) then
          ! root is in the left half
          pstar_hi = pstar_c
          f_hi = f_c
       else
          pstar_lo = pstar_c
          f_lo = f_c
       endif
    enddo

    pstar = pstar_c

  end subroutine pstar_bisection


  subroutine HLL(ql, qr, cl, cr, idir, ndim, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map
    use prob_params_module, only : mom_flux_has_p

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    real(rt)        , intent(inout) :: f(NVAR)
    integer, intent(in) :: idir, ndim

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    real(rt)         :: a1, a4, bd, bl, bm, bp, br
    real(rt)         :: cavg, uavg
    real(rt)         :: fl_tmp, fr_tmp
    real(rt)         :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    real(rt)        , parameter :: small = 1.e-10_rt

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)



    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )


    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg


    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd


    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux.  Note for 1-d and 2-d non cartesian
    ! r-coordinate, we leave off the pressure term and handle that
    ! separately in the update, to accommodate different geometries
    fl_tmp = ql(QRHO)*ql(ivel)**2
    fr_tmp = qr(QRHO)*qr(ivel)**2
    if (mom_flux_has_p(idir)%comp(UMX-1+idir)) then
       fl_tmp = fl_tmp + ql(QPRES)
       fr_tmp = fr_tmp + qr(QPRES)
    endif

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL


  pure subroutine cons_state(q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    real(rt)        , intent(in)  :: q(QVAR)
    real(rt)        , intent(out) :: U(NVAR)

    integer :: ipassive, n, nq

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = q(QRHO)*q(QU)
    U(UMY)  = q(QRHO)*q(QV)
    U(UMZ)  = q(QRHO)*q(QW)

    U(UEDEN) = q(QREINT) + HALF*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = ZERO

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = q(QRHO)*q(nq)
    enddo

  end subroutine cons_state


  pure subroutine HLLC_state(idir, S_k, S_c, q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, QPRES, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    integer, intent(in) :: idir
    real(rt)        , intent(in)  :: S_k, S_c
    real(rt)        , intent(in)  :: q(QVAR)
    real(rt)        , intent(out) :: U(NVAR)

    real(rt)         :: hllc_factor, u_k
    integer :: ipassive, n, nq

    if (idir == 1) then
       u_k = q(QU)
    elseif (idir == 2) then
       u_k = q(QV)
    elseif (idir == 3) then
       u_k = q(QW)
    endif

    hllc_factor = q(QRHO)*(S_k - u_k)/(S_k - S_c)
    U(URHO) = hllc_factor
    if (idir == 1) then
       U(UMX)  = hllc_factor*S_c
       U(UMY)  = hllc_factor*q(QV)
       U(UMZ)  = hllc_factor*q(QW)
    elseif (idir == 2) then
       U(UMX)  = hllc_factor*q(QU)
       U(UMY)  = hllc_factor*S_c
       U(UMZ)  = hllc_factor*q(QW)
    elseif (idir == 3) then
       U(UMX)  = hllc_factor*q(QU)
       U(UMY)  = hllc_factor*q(QV)
       U(UMZ)  = hllc_factor*S_c
    endif

    U(UEDEN) = hllc_factor*(q(QREINT)/q(QRHO) + &
         HALF*(q(QU)**2 + q(QV)**2 + q(QW)**2) + &
         (S_c - u_k)*(S_c + q(QPRES)/(q(QRHO)*(S_k - u_k))))
    U(UEINT) = hllc_factor*q(QREINT)/q(QRHO)

    U(UTEMP) = ZERO  ! we don't evolve T

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = hllc_factor*q(nq)
    enddo

  end subroutine HLLC_state


  pure subroutine compute_flux(idir, ndim, bnd_fac, U, p, F)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map
    use prob_params_module, only : mom_flux_has_p

    integer, intent(in) :: idir, ndim, bnd_fac
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(in) :: p
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = U(UMX)/U(URHO)
    elseif (idir == 2) then
       u_flx = U(UMY)/U(URHO)
    elseif (idir == 3) then
       u_flx = U(UMZ)/U(URHO)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(URHO) = U(URHO)*u_flx

    F(UMX) = U(UMX)*u_flx
    F(UMY) = U(UMY)*u_flx
    F(UMZ) = U(UMZ)*u_flx

    if (mom_flux_has_p(idir)%comp(UMX-1+idir)) then
       ! we do not include the pressure term in any non-Cartesian
       ! coordinate directions
       F(UMX-1+idir) = F(UMX-1+idir) + p
    endif

    F(UEINT) = U(UEINT)*u_flx
    F(UEDEN) = (U(UEDEN) + p)*u_flx

    F(UTEMP) = ZERO

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       F(n) = U(n)*u_flx
    enddo

  end subroutine compute_flux

end module riemann_module
