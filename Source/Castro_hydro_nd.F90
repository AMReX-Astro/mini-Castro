module hydro_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, HALF, ONE

  implicit none

contains

  subroutine divu(lo, hi, &
                  q, q_lo, q_hi, &
                  dx, div, div_lo, div_hi) bind(C, name='divu')
    ! this computes the *node-centered* divergence

    use amrex_constants_module, only: FOURTH
    use castro_module, only: QU, QV, QW, QVAR

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE / dx(1)
    dyinv = ONE / dx(2)
    dzinv = ONE / dx(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ux = FOURTH*( &
                  + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                  + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                  + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                  + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                  + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                  + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                  + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                  + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                  + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                  + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                  + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                  + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu



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
                             small_dens, dual_energy_eta1

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

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
                                 qp, qp_lo, qp_hi, &
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

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    call cmpflx(lo, hi, &
                qm, qm_lo, qm_hi, &
                qp, qp_lo, qp_hi, &
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
                    qp, qp_lo, qp_hi, &
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

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! Solve Riemann problem to get the fluxes

    call riemann_state(qm, qm_lo, qm_hi, &
                       qp, qp_lo, qp_hi, &
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
                           qp, qp_lo, qp_hi, &
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

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)

    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! Solve Riemann problem

    call riemannus(qm, qm_lo, qm_hi, &
                   qp, qp_lo, qp_hi, &
                   qaux, qa_lo, qa_hi, &
                   qint, q_lo, q_hi, &
                   idir, lo, hi, &
                   domlo, domhi)

  end subroutine riemann_state



  subroutine riemannus(ql, ql_lo, ql_hi, &
                       qr, qr_lo, qr_hi, &
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
                             NQAUX, URHO, UMX, UMY, UMZ, UFS, small, small_dens, smallu, small_pres
    use network, only: nspec

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3),domhi(3)

    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),QVAR)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),QVAR)

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

          end do
       end do
    end do

  end subroutine riemannus



  subroutine apply_av(lo, hi, idir, dx, &
                      div, div_lo, div_hi, &
                      uin, uin_lo, uin_hi, &
                      flux, f_lo, f_hi) bind(c, name="apply_av")

    use amrex_constants_module, only: FOURTH
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
                        dx, dt) bind(C, name="ctu_consup")

    use castro_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                             UEINT, UTEMP, NGDNV, QVAR, &
                             GDU, GDV, GDW, GDPRES

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
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    integer :: i, j, g, k, n
    integer :: domlo(3), domhi(3)
    real(rt) :: volInv
    real(rt) :: pdivu

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
                   pdivu = HALF * (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) * &
                                  (qx(i+1,j,k,GDU) - qx(i,j,k,GDU)) / dx(1) + &
                           HALF * (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) * &
                                  (qy(i,j+1,k,GDV) - qy(i,j,k,GDV)) / dx(2) + &
                           HALF * (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) * &
                                  (qz(i,j,k+1,GDW) - qz(i,j,k,GDW)) / dx(3)

                   update(i,j,k,n) = update(i,j,k,n) - pdivu
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine ctu_consup



  subroutine store_flux(lo, hi, &
                        flux_out, fo_lo, fo_hi, &
                        flux_in, fi_lo, fi_hi, &
                        area, a_lo, a_hi, &
                        dt) bind(C, name="store_flux")

    use castro_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fo_lo(3), fo_hi(3)
    integer,  intent(in   ) :: fi_lo(3), fi_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(inout) :: flux_out(fo_lo(1):fo_hi(1),fo_lo(2):fo_hi(2),fo_lo(3):fo_hi(3),NVAR)
    real(rt), intent(in   ) :: flux_in(fi_lo(1):fi_hi(1),fi_lo(2):fi_hi(2),fi_lo(3):fi_hi(3),NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    real(rt), intent(in), value :: dt

    integer :: i, j, k, n

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux_out(i,j,k,n) = dt * flux_in(i,j,k,n) * area(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine store_flux

end module hydro_module
