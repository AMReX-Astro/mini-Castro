module hydro_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, HALF, ONE
  use castro_module, only: small_pres, small

  implicit none

contains

  subroutine trans1_on_2states(lo, hi, &
                               idir1, idir2, &
                               q2m, q2m_lo, q2m_hi, &
                               q2mo, q2mo_lo, q2mo_hi, &
                               q2p, q2p_lo, q2p_hi, &
                               q2po, q2po_lo, q2po_hi, &
                               qaux, qa_lo, qa_hi, &
                               f1, f1_lo, f1_hi, &
                               q1, q1_lo, q1_hi, &
                               cdtdx) bind(C, name="trans1_on_2states")

    use network, only: nspec
    use castro_module, only: QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                             QPRES, QREINT, QGAME, QFS, &
                             QC, QGAMC, &
                             URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME

    integer, intent(in) :: q2m_lo(3), q2m_hi(3)
    integer, intent(in) :: q2p_lo(3), q2p_hi(3)
    integer, intent(in) :: q2mo_lo(3), q2mo_hi(3)
    integer, intent(in) :: q2po_lo(3), q2po_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: cdtdx
    integer,  intent(in), value :: idir1, idir2

    real(rt), intent(in) :: q2m(q2m_lo(1):q2m_hi(1),q2m_lo(2):q2m_hi(2),q2m_lo(3):q2m_hi(3),QVAR)
    real(rt), intent(in) :: q2p(q2p_lo(1):q2p_hi(1),q2p_lo(2):q2p_hi(2),q2p_lo(3):q2p_hi(3),QVAR)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: f1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)

    real(rt), intent(out) :: q2mo(q2mo_lo(1):q2mo_hi(1),q2mo_lo(2):q2mo_hi(2),q2mo_lo(3):q2mo_hi(3),QVAR)
    real(rt), intent(out) :: q2po(q2po_lo(1):q2po_hi(1),q2po_lo(2):q2po_hi(2),q2po_lo(3):q2po_hi(3),QVAR)

    integer :: d, il, jl, kl, ir, jr, kr

    integer  :: i, j, k, n, nqp, ispec

    real(rt) :: lq2(QVAR), lq2o(QVAR)

    real(rt) :: rhoinv
    real(rt) :: rrnew
    real(rt) :: rrr2, rrl2
    real(rt) :: rur2, rul2
    real(rt) :: rvr2, rvl2
    real(rt) :: rwr2, rwl2
    real(rt) :: ekenr2, ekenl2
    real(rt) :: rer2, rel2
    real(rt) :: rrnewr2, rrnewl2
    real(rt) :: runewr2, runewl2
    real(rt) :: rvnewr2, rvnewl2
    real(rt) :: rwnewr2, rwnewl2
    real(rt) :: renewr2, renewl2
    real(rt) :: pnewr2, pnewl2
    real(rt) :: rhoekenr2, rhoekenl2
    real(rt) :: pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt) :: compu
    real(rt) :: gamc

    logical :: reset_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !       qm|qp
             !         |
             ! --------+--------
             !   i-1       i
             !        i-1/2
             !
             ! the qm state will see the transverse flux in zone i-1

             ! Loop over plus and minus states

             do d = -1, 0

                ! We are handling the states at the interface of
                ! (i, i+1) in the x-direction, and similarly for
                ! the y- and z- directions.

                il = i
                jl = j
                kl = k

                if (idir1 == 1) then
                   ir = i+1
                   jr = j
                   kr = k
                else if (idir1 == 2) then
                   ir = i
                   jr = j+1
                   kr = k
                else
                   ir = i
                   jr = j
                   kr = k+1
                end if

                ! We're handling both the plus and minus states;
                ! for the minus state we're shifting one zone to
                ! the left in our chosen direction.

                if (idir2 == 1) then
                   il = il+d
                   ir = ir+d
                else if (idir2 == 2) then
                   jl = jl+d
                   jr = jr+d
                else
                   kl = kl+d
                   kr = kr+d
                end if

                if (d == -1) then
                   lq2(:) = q2m(i,j,k,:)
                else
                   lq2(:) = q2p(i,j,k,:)
                end if

                !-------------------------------------------------------------------------
                ! update all of the passively-advected quantities with the
                ! transverse term and convert back to the primitive quantity
                !-------------------------------------------------------------------------

                do ispec = 1, nspec
                   n  = UFS + ispec - 1
                   nqp = QFS + ispec - 1

                   rrnew = lq2(QRHO) - cdtdx*(f1(ir,jr,kr,URHO) - f1(il,jl,kl,URHO))
                   compu = lq2(QRHO)*lq2(nqp) - cdtdx*(f1(ir,jr,kr,n) - f1(il,jl,kl,n))
                   lq2o(nqp) = compu/rrnew
                end do

                !-------------------------------------------------------------------
                ! add the transverse flux difference in the 1-direction to 2-states
                ! for the fluid variables
                !-------------------------------------------------------------------

                pgp  = q1(ir,jr,kr,GDPRES)
                pgm  = q1(il,jl,kl,GDPRES)
                ugp  = q1(ir,jr,kr,GDU+idir1-1)
                ugm  = q1(il,jl,kl,GDU+idir1-1)
                gegp = q1(ir,jr,kr,GDGAME)
                gegm = q1(il,jl,kl,GDGAME)

                ! we need to augment our conserved system with either a p
                ! equation or gammae (if we have ppm_predict_gammae = 1) to
                ! be able to deal with the general EOS

                dup = pgp*ugp - pgm*ugm
                du = ugp-ugm
                pav = HALF*(pgp+pgm)
                uav = HALF*(ugp+ugm)
                geav = HALF*(gegp+gegm)
                dge = gegp-gegm

                ! this is the gas gamma_1
                gamc = qaux(il,jl,kl,QGAMC)

                ! Convert to conservation form
                rrl2 = lq2(QRHO)
                rul2 = rrl2*lq2(QU)
                rvl2 = rrl2*lq2(QV)
                rwl2 = rrl2*lq2(QW)
                ekenl2 = HALF*rrl2*sum(lq2(QU:QW)**2)
                rel2 = lq2(QREINT) + ekenl2

                ! Add transverse predictor
                rrnewl2 = rrl2 - cdtdx*(f1(ir,jr,kr,URHO) - f1(il,jl,kl,URHO))
                runewl2 = rul2 - cdtdx*(f1(ir,jr,kr,UMX) - f1(il,jl,kl,UMX))
                rvnewl2 = rvl2 - cdtdx*(f1(ir,jr,kr,UMY) - f1(il,jl,kl,UMY))
                rwnewl2 = rwl2 - cdtdx*(f1(ir,jr,kr,UMZ) - f1(il,jl,kl,UMZ))
                renewl2 = rel2 - cdtdx*(f1(ir,jr,kr,UEDEN) - f1(il,jl,kl,UEDEN))

                ! Reset to original value if adding transverse terms made density negative
                reset_state = .false.
                if (rrnewl2 < ZERO) then
                   rrnewl2 = rrl2
                   runewl2 = rul2
                   rvnewl2 = rvl2
                   rwnewl2 = rwl2
                   renewl2 = rel2
                   reset_state = .true.
                endif

                ! Convert back to primitive form
                lq2o(QRHO) = rrnewl2
                rhoinv = ONE/rrnewl2
                lq2o(QU) = runewl2*rhoinv
                lq2o(QV) = rvnewl2*rhoinv
                lq2o(QW) = rwnewl2*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenl2 = HALF*(runewl2**2 + rvnewl2**2 + rwnewl2**2)*rhoinv
                lq2o(QREINT) = renewl2 - rhoekenl2

                if (.not. reset_state) then
                   pnewl2 = lq2(QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   lq2o(QPRES) = max(pnewl2,small_pres)
                else
                   lq2o(QPRES) = lq2(QPRES)
                   lq2o(QGAME) = lq2(QGAME)
                endif

                if (d == -1) then
                   q2mo(i,j,k,:) = lq2o(:)
                else
                   q2po(i,j,k,:) = lq2o(:)
                end if

             end do

          end do
       end do
    end do

  end subroutine trans1_on_2states



  ! Add the transverse corrections from directions 2 and 3
  ! to the states in direction 1.

  subroutine trans2(lo, hi, &
                    idir1, idir2, idir3, &
                    qm1, qm1_lo, qm1_hi, &
                    qm1o, qm1o_lo, qm1o_hi, &
                    qp1, qp1_lo, qp1_hi, &
                    qp1o, qp1o_lo, qp1o_hi, &
                    qaux, qa_lo, qa_hi, &
                    f2, f2_lo, f2_hi, &
                    f3, f3_lo, f3_hi, &
                    q2, q2_lo, q2_hi, &
                    q3, q3_lo, q3_hi, &
                    cdtdx1, cdtdx2, cdtdx3) bind(C, name="trans2")

    use network, only: nspec
    use castro_module, only: QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                             QPRES, QREINT, QGAME, QFS, &
                             QC, QGAMC, &
                             URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME

    integer, intent(in) :: qm1_lo(3), qm1_hi(3)
    integer, intent(in) :: qm1o_lo(3), qm1o_hi(3)
    integer, intent(in) :: qp1_lo(3), qp1_hi(3)
    integer, intent(in) :: qp1o_lo(3), qp1o_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: f2_lo(3), f2_hi(3)
    integer, intent(in) :: f3_lo(3), f3_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: idir1, idir2, idir3

    real(rt), intent(in), value :: cdtdx1, cdtdx2, cdtdx3

    real(rt), intent(in) :: qm1(qm1_lo(1):qm1_hi(1),qm1_lo(2):qm1_hi(2),qm1_lo(3):qm1_hi(3),QVAR)
    real(rt), intent(in) :: qp1(qp1_lo(1):qp1_hi(1),qp1_lo(2):qp1_hi(2),qp1_lo(3):qp1_hi(3),QVAR)
    real(rt), intent(out) :: qm1o(qm1o_lo(1):qm1o_hi(1),qm1o_lo(2):qm1o_hi(2),qm1o_lo(3):qm1o_hi(3),QVAR)
    real(rt), intent(out) :: qp1o(qp1o_lo(1):qp1o_hi(1),qp1o_lo(2):qp1o_hi(2),qp1o_lo(3):qp1o_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: f2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(in) :: f3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)

    integer :: i, j, k, n, nqp, ispec
    integer :: d
    integer :: il1, jl1, kl1, ir1, jr1, kr1, il2, jl2, kl2, ir2, jr2, kr2, il3, jl3, kl3, ir3, jr3, kr3

    real(rt) :: lqo(QVAR), lq(QVAR)

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pg2p, pg2m, ug2p, ug2m, geg2p, geg2m, du2p, p2av, du2, p2new, ge2new
    real(rt) :: pg3p, pg3m, ug3p, ug3m, geg3p, geg3m, du3p, p3av, du3, p3new, ge3new
    real(rt) :: u2av, ge2av, dge2, u3av, ge3av, dge3
    real(rt) :: compr, compl, compnr, compnl

    logical :: reset_state

    !-------------------------------------------------------------------
    ! add the transverse differences to the states for the fluid variables
    ! the states we're updating are determined by the 1-index, while the
    ! transverse differences come from the 2 and 3 indices
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             do d = -1, 0

                il1 = i
                jl1 = j
                kl1 = k

                il2 = i
                jl2 = j
                kl2 = k

                il3 = i
                jl3 = j
                kl3 = k

                if (idir1 == 1) then
                   ir2 = i+d
                   jr2 = j+1
                   kr2 = k

                   ir3 = i+d
                   jr3 = j
                   kr3 = k+1

                   il1 = i+d
                   il2 = i+d
                   il3 = i+d
                else if (idir1 == 2) then
                   ir2 = i+1
                   jr2 = j+d
                   kr2 = k

                   ir3 = i
                   jr3 = j+d
                   kr3 = k+1

                   jl1 = j+d
                   jl2 = j+d
                   jl3 = j+d
                else
                   ir2 = i+1
                   jr2 = j
                   kr2 = k+d

                   ir3 = i
                   jr3 = j+1
                   kr3 = k+d

                   kl1 = k+d
                   kl2 = k+d
                   kl3 = k+d
                end if

                if (d == -1) then
                   lq(:) = qm1(i,j,k,:)
                else
                   lq(:) = qp1(i,j,k,:)
                end if

                !-------------------------------------------------------------------------
                ! update all of the passively-advected quantities with the
                ! transerse term and convert back to the primitive quantity
                !-------------------------------------------------------------------------

                do ispec = 1, nspec
                   n  = UFS + ispec - 1
                   nqp = QFS + ispec - 1

                   rrr = lq(QRHO)
                   compr = rrr*lq(nqp)
                   rrnewr = rrr - cdtdx2*(f2(ir2,jr2,kr2,URHO) - f2(il2,jl2,kl2,URHO)) &
                                - cdtdx3*(f3(ir3,jr3,kr3,URHO) - f3(il3,jl3,kl3,URHO))
                   compnr = compr - cdtdx2*(f2(ir2,jr2,kr2,n) - f2(il2,jl2,kl2,n)) &
                                  - cdtdx3*(f3(ir3,jr3,kr3,n) - f3(il3,jl3,kl3,n))

                   lqo(nqp) = compnr/rrnewr
                end do

                pg2p  = q2(ir2,jr2,kr2,GDPRES)
                pg2m  = q2(il2,jl2,kl2,GDPRES)
                ug2p  = q2(ir2,jr2,kr2,GDU+idir2-1)
                ug2m  = q2(il2,jl2,kl2,GDU+idir2-1)
                geg2p = q2(ir2,jr2,kr2,GDGAME)
                geg2m = q2(il2,jl2,kl2,GDGAME)

                du2p = pg2p*ug2p - pg2m*ug2m
                p2av = HALF*(pg2p+pg2m)
                u2av = HALF*(ug2p+ug2m)
                ge2av = HALF*(geg2p+geg2m)
                du2 = ug2p-ug2m
                dge2 = geg2p-geg2m

                p2new = cdtdx2*(du2p + p2av*du2*(qaux(il1,jl1,kl1,QGAMC) - ONE))
                ge2new = cdtdx2*( (ge2av-ONE)*(ge2av - qaux(il1,jl1,kl1,QGAMC))*du2 - u2av*dge2 )

                pg3p  = q3(ir3,jr3,kr3,GDPRES)
                pg3m  = q3(il3,jl3,kl3,GDPRES)
                ug3p  = q3(ir3,jr3,kr3,GDU+idir3-1)
                ug3m  = q3(il3,jl3,kl3,GDU+idir3-1)
                geg3p = q3(ir3,jr3,kr3,GDGAME)
                geg3m = q3(il3,jl3,kl3,GDGAME)

                du3p = pg3p*ug3p - pg3m*ug3m
                p3av = HALF*(pg3p+pg3m)
                u3av = HALF*(ug3p+ug3m)
                ge3av = HALF*(geg3p+geg3m)
                du3 = ug3p-ug3m
                dge3 = geg3p-geg3m

                p3new = cdtdx3*(du3p + p3av*du3*(qaux(il1,jl1,kl1,QGAMC) - ONE))
                ge3new = cdtdx3*( (ge3av-ONE)*(ge3av - qaux(il1,jl1,kl1,QGAMC))*du3 - u3av*dge3 )

                ! Convert to conservation form
                rrr = lq(QRHO)
                rur = rrr*lq(QU)
                rvr = rrr*lq(QV)
                rwr = rrr*lq(QW)
                ekenr = HALF*rrr*sum(lq(QU:QW)**2)
                rer = lq(QREINT) + ekenr

                ! Add transverse predictor
                rrnewr = rrr - cdtdx2*(f2(ir2,jr2,kr2,URHO) - f2(il2,jl2,kl2,URHO)) &
                             - cdtdx3*(f3(ir3,jr3,kr3,URHO) - f3(il3,jl3,kl3,URHO))
                runewr = rur - cdtdx2*(f2(ir2,jr2,kr2,UMX) - f2(il2,jl2,kl2,UMX)) &
                             - cdtdx3*(f3(ir3,jr3,kr3,UMX) - f3(il3,jl3,kl3,UMX))
                rvnewr = rvr - cdtdx2*(f2(ir2,jr2,kr2,UMY) - f2(il2,jl2,kl2,UMY)) &
                             - cdtdx3*(f3(ir3,jr3,kr3,UMY) - f3(il3,jl3,kl3,UMY))
                rwnewr = rwr - cdtdx2*(f2(ir2,jr2,kr2,UMZ) - f2(il2,jl2,kl2,UMZ)) &
                             - cdtdx3*(f3(ir3,jr3,kr3,UMZ) - f3(il3,jl3,kl3,UMZ))
                renewr = rer - cdtdx2*(f2(ir2,jr2,kr2,UEDEN) - f2(il2,jl2,kl2,UEDEN)) &
                             - cdtdx3*(f3(ir3,jr3,kr3,UEDEN) - f3(il3,jl3,kl3,UEDEN))

                ! Reset to original value if adding transverse terms
                ! made density negative
                reset_state = .false.
                if (rrnewr < ZERO) then
                   rrnewr = rrr
                   runewr = rur
                   rvnewr = rvr
                   rwnewr = rwr
                   renewr = rer
                   reset_state = .true.
                end if

                lqo(QRHO  ) = rrnewr
                lqo(QU    ) = runewr/rrnewr
                lqo(QV    ) = rvnewr/rrnewr
                lqo(QW    ) = rwnewr/rrnewr

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                lqo(QREINT) = renewr - rhoekenr

                if (.not. reset_state) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = lq(QPRES) - p2new - p3new
                   lqo(QPRES) = pnewr
                else
                   lqo(QPRES) = lq(QPRES)
                   lqo(QGAME) = lq(QGAME)
                endif

                lqo(QPRES) = max(lqo(QPRES), small_pres)

                if (d == -1) then
                   qm1o(i,j,k,:) = lqo(:)
                else
                   qp1o(i,j,k,:) = lqo(:)
                end if

             end do

          end do
       end do
    end do

  end subroutine trans2



  subroutine transxz(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxz, fxz_lo, fxz_hi, &
                     fzx, fzx_lo, fzx_hi, &
                     qx, qx_lo, qx_hi, &
                     qz, qz_lo, qz_hi, &
                     hdt, cdtdx, cdtdz) bind(C, name="transxz")

    use network, only: nspec
    use castro_module, only: QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                             QPRES, QREINT, QGAME, QFS, &
                             QC, QGAMC, &
                             URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fxz_lo(3), fxz_hi(3)
    integer, intent(in) :: fzx_lo(3), fzx_hi(3)
    integer, intent(in) :: qx_lo(3),qx_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdz

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),QVAR)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxz(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3),NVAR)
    real(rt), intent(in) :: fzx(fzx_lo(1):fzx_hi(1),fzx_lo(2):fzx_hi(2),fzx_lo(3):fzx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    integer  :: i, j, k, n, nqp, ispec

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt) :: pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt) :: uxav, gexav, dgex, uzav, gezav, dgez
    real(rt) :: compr, compl, compnr, compnl

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nqp = QFS + ispec - 1

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp(i,j,k,QRHO)
                compr = rrr*qp(i,j,k,nqp)
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                               - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                qpo(i,j  ,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                             - cdtdz*(fzx(i  ,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
                compnl = compl - cdtdx*(fxz(i+1,j-1,k,n) - fxz(i,j-1,k,n)) &
                               - cdtdz*(fzx(i  ,j-1,k+1,n) - fzx(i,j-1,k,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse xz and zx differences to the y-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qypo state
             !-------------------------------------------------------------------

             pgxp  = qx(i+1,j,k,GDPRES)
             pgxm  = qx(i,j,k,GDPRES)
             ugxp  = qx(i+1,j,k,GDU)
             ugxm  = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)

             pgzp  = qz(i,j,k+1,GDPRES)
             pgzm  = qz(i,j,k,GDPRES)
             ugzp  = qz(i,j,k+1,GDW)
             ugzm  = qz(i,j,k,GDW)
             gegzp = qz(i,j,k+1,GDGAME)
             gegzm = qz(i,j,k,GDGAME)

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm

             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm

             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                          - cdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
             runewr = rur - cdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                          - cdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
             rvnewr = rvr - cdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                          - cdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
             rwnewr = rwr - cdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                          - cdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
             renewr = rer - cdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                          - cdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qpo(i,j,k,QRHO  ) = rrnewr
             qpo(i,j,k,QU    ) = runewr/rrnewr
             qpo(i,j,k,QV    ) = rvnewr/rrnewr
             qpo(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,k,QPRES) - pxnew - pznew
                qpo(i,j,k,QPRES) = pnewr
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)



             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             pgxp  = qx(i+1,j-1,k,GDPRES)
             pgxm  = qx(i,j-1,k,GDPRES)
             ugxp  = qx(i+1,j-1,k,GDU)
             ugxm  = qx(i,j-1,k,GDU)
             gegxp = qx(i+1,j-1,k,GDGAME)
             gegxm = qx(i,j-1,k,GDGAME)

             pgzp  = qz(i,j-1,k+1,GDPRES)
             pgzm  = qz(i,j-1,k,GDPRES)
             ugzp  = qz(i,j-1,k+1,GDW)
             ugzm  = qz(i,j-1,k,GDW)
             gegzp = qz(i,j-1,k+1,GDGAME)
             gegzm = qz(i,j-1,k,GDGAME)

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm

             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j-1,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j-1,k,QGAMC))*dux - uxav*dgex )

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm

             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j-1,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j-1,k,QGAMC))*duz - uzav*dgez )

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm(i,j,k,QU:QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                          - cdtdz*(fzx(i,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j-1,k,UMX) - fxz(i,j-1,k,UMX)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMX) - fzx(i,j-1,k,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j-1,k,UMY) - fxz(i,j-1,k,UMY)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMY) - fzx(i,j-1,k,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j-1,k,UMZ) - fxz(i,j-1,k,UMZ)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMZ) - fzx(i,j-1,k,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j-1,k,UEDEN) - fxz(i,j-1,k,UEDEN)) &
                          - cdtdz*(fzx(i,j-1,k+1,UEDEN) - fzx(i,j-1,k,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,k,QPRES) - pxnew - pznew
                qmo(i,j,k,QPRES) = pnewl
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

          end do
       end do
    end do

  end subroutine transxz



  subroutine transxy(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxy, fxy_lo, fxy_hi, &
                     fyx, fyx_lo, fyx_hi, &
                     qx, qx_lo, qx_hi, &
                     qy, qy_lo, qy_hi, &
                     hdt, cdtdx, cdtdy) bind(C, name="transxy")

    use network, only: nspec
    use castro_module, only: QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                             QPRES, QREINT, QGAME, QFS, &
                             QC, QGAMC, &
                             URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fxy_lo(3), fxy_hi(3)
    integer, intent(in) :: fyx_lo(3), fyx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdy

    real(rt), intent(in) ::   qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(in) ::   qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),QVAR)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),QVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxy(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3),NVAR)
    real(rt), intent(in) :: fyx(fyx_lo(1):fyx_hi(1),fyx_lo(2):fyx_hi(2),fyx_lo(3):fyx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    integer  :: i, j, k, n, nqp, ispec

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt) :: pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt) :: uxav, gexav, dgex, uyav, geyav, dgey
    real(rt) :: compr, compl, compnr, compnl

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nqp = QFS + ispec - 1

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp(i,j,k,QRHO)
                compr = rrr*qp(i,j,k,nqp)
                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                               - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                qpo(i,j,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                             - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
                compnl = compl - cdtdx*(fxy(i+1,j,k-1,n) - fxy(i,j,k-1,n)) &
                               - cdtdy*(fyx(i,j+1,k-1,n) - fyx(i,j,k-1,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse xy and yx differences to the z-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

             pgxp = qx(i+1,j,k,GDPRES)
             pgxm = qx(i,j,k,GDPRES)
             ugxp = qx(i+1,j,k,GDU)
             ugxm = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)

             pgyp = qy(i,j+1,k,GDPRES)
             pgym = qy(i,j,k,GDPRES)
             ugyp = qy(i,j+1,k,GDV)
             ugym = qy(i,j,k,GDV)
             gegyp = qy(i,j+1,k,GDGAME)
             gegym = qy(i,j,k,GDGAME)

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm

             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym

             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                          - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
             runewr = rur - cdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                          - cdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
             rvnewr = rvr - cdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                          - cdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
             rwnewr = rwr - cdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
             renewr = rer - cdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qpo(i,j,k,QRHO  ) = rrnewr
             qpo(i,j,k,QU    ) = runewr/rrnewr
             qpo(i,j,k,QV    ) = rvnewr/rrnewr
             qpo(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                qpo(i,j,k,QPRES) = pnewr
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)



             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pgxp = qx(i+1,j,k-1,GDPRES)
             pgxm = qx(i,j,k-1,GDPRES)
             ugxp = qx(i+1,j,k-1,GDU)
             ugxm = qx(i,j,k-1,GDU)
             gegxp = qx(i+1,j,k-1,GDGAME)
             gegxm = qx(i,j,k-1,GDGAME)

             pgyp = qy(i,j+1,k-1,GDPRES)
             pgym = qy(i,j,k-1,GDPRES)
             ugyp = qy(i,j+1,k-1,GDV)
             ugym = qy(i,j,k-1,GDV)
             gegyp = qy(i,j+1,k-1,GDGAME)
             gegym = qy(i,j,k-1,GDGAME)

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm

             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k-1,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k-1,QGAMC))*dux - uxav*dgex )

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym

             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k-1,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k-1,QGAMC))*duy - uyav*dgey )

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm(i,j,k,QU:QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                          - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
             runewl = rul - cdtdx*(fxy(i+1,j,k-1,UMX) - fxy(i,j,k-1,UMX)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMX) - fyx(i,j,k-1,UMX))
             rvnewl = rvl - cdtdx*(fxy(i+1,j,k-1,UMY) - fxy(i,j,k-1,UMY)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMY) - fyx(i,j,k-1,UMY))
             rwnewl = rwl - cdtdx*(fxy(i+1,j,k-1,UMZ) - fxy(i,j,k-1,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMZ) - fyx(i,j,k-1,UMZ))
             renewl = rel - cdtdx*(fxy(i+1,j,k-1,UEDEN) - fxy(i,j,k-1,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k-1,UEDEN) - fyx(i,j,k-1,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,k,QPRES) - pxnew - pynew
                qmo(i,j,k,QPRES) = pnewl
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

          end do
       end do
    end do

  end subroutine transxy



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
                             NQAUX, URHO, UMX, UMY, UMZ, UFS, small_dens, smallu
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

    use castro_module, only: NVAR

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



  subroutine store_flux(lo, hi, &
                        flux_out, fo_lo, fo_hi, &
                        flux_in, fi_lo, fi_hi) bind(C, name="store_flux")

    use castro_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fo_lo(3), fo_hi(3)
    integer,  intent(in   ) :: fi_lo(3), fi_hi(3)

    real(rt), intent(inout) :: flux_out(fo_lo(1):fo_hi(1),fo_lo(2):fo_hi(2),fo_lo(3):fo_hi(3),NVAR)
    real(rt), intent(inout) :: flux_in(fi_lo(1):fi_hi(1),fi_lo(2):fi_hi(2),fi_lo(3):fi_hi(3),NVAR)

    integer :: i, j, k, n

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux_out(i,j,k,n) = flux_in(i,j,k,n)
             enddo
          enddo
       enddo
    enddo

  end subroutine store_flux

end module hydro_module
