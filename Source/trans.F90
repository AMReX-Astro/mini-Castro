module transverse_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use amrex_constants_module, only: ZERO, HALF, ONE

  implicit none

contains

  ! Add the transverse corrections from directions 2 and 3
  ! to the states in direction 1.

  CASTRO_FORT_DEVICE subroutine trans1(lo, hi, &
                                       idir1, idir2, &
                                       q2m, q2m_lo, q2m_hi, &
                                       q2mo, q2mo_lo, q2mo_hi, &
                                       q2p, q2p_lo, q2p_hi, &
                                       q2po, q2po_lo, q2po_hi, &
                                       qaux, qa_lo, qa_hi, &
                                       f1, f1_lo, f1_hi, &
                                       q1, q1_lo, q1_hi, &
                                       cdtdx) bind(C, name="trans1")

    use network, only: nspec
    use castro_module, only: QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                             QPRES, QREINT, QGAME, QFS, &
                             QC, QGAMC, &
                             URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                             small_pres

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

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(q2m, q2p, q2mo, q2po, qaux, f1, q1) private(lq2, lq2o) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(q2m, q2p, q2mo, q2po, qaux, f1, q1) private(lq2, lq2o)
#endif
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
                   lq2o(QPRES) = max(pnewl2, small_pres)
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

  end subroutine trans1



  ! Add the transverse corrections from directions 2 and 3
  ! to the states in direction 1.

  CASTRO_FORT_DEVICE subroutine trans2(lo, hi, &
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
                             NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, small_pres

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

#ifdef AMREX_USE_ACC
    !$acc parallel loop gang vector collapse(3) deviceptr(qm1, qp1, qm1o, qp1o, qaux, f2, f3, q2, q3) private(lqo, lq) async(acc_stream)
#endif
#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) is_device_ptr(qm1, qp1, qm1o, qp1o, qaux, f2, f3, q2, q3) private(lqo, lq)
#endif
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

end module transverse_module
