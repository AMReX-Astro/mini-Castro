module flatten_module

  use amrex_acc_module, only: acc_stream

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine ca_uflatten(lo, hi, q, q_lo, q_hi, flatn, f_lo, f_hi) bind(C, name='ca_uflatten')

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: QVAR, QU, QV, QW, QPRES

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    integer :: i, j, k, ishft, idir

    real(rt), parameter :: small_pres = 1.e-200_rt

    real(rt) :: denom, zeta, tst, tmp
    real(rt) :: pl, pr, dp, dp2, du, z, z2, chi, chi2

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    ! x-direction flattening coef

    idir = 1
    
    !$acc parallel loop gang vector collapse(3) deviceptr(q, flatn) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             flatn(i,j,k) = ONE

             do idir = 1, 3

                if (idir == 1) then
                   pr = q(i+1,j,k,QPRES)
                   pl = q(i-1,j,k,QPRES)
                   dp2 = q(i+2,j,k,QPRES) - q(i-2,j,k,QPRES)
                   du = q(i+1,j,k,QU) - q(i-1,j,k,QU)
                else if (idir == 2) then
                   pr = q(i,j+1,k,QPRES)
                   pl = q(i,j-1,k,QPRES)
                   dp2 = q(i,j+2,k,QPRES) - q(i,j-2,k,QPRES)
                   du = q(i,j+1,k,QV) - q(i,j-1,k,QV)
                else
                   pr = q(i,j,k+1,QPRES)
                   pl = q(i,j,k-1,QPRES)
                   dp2 = q(i,j,k+2,QPRES) - q(i,j,k-2,QPRES)
                   du = q(i,j,k+1,QW) - q(i,j,k-1,QW)
                end if
                dp = pr - pl

                if (pr - pl .gt. ZERO) then
                   ishft = 1
                else
                   ishft = -1
                endif

                denom = max(small_pres, abs(dp2))
                zeta = abs(dp) / denom
                z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

                if (du .le. ZERO) then
                   tst = ONE
                else
                   tst = ZERO
                endif

                tmp = min(pr, pl)

                if ((abs(dp)/tmp) .gt. shktst) then
                   chi = tst
                else
                   chi = ZERO
                endif

                if (idir == 1) then
                   pr = q(i+1-ishft,j,k,QPRES)
                   pl = q(i-1-ishft,j,k,QPRES)
                   dp2 = q(i+2-ishft,j,k,QPRES) - q(i-2-ishft,j,k,QPRES)
                   du = q(i+1-ishft,j,k,QU) - q(i-1-ishft,j,k,QU)
                else if (idir == 2) then
                   pr = q(i,j+1-ishft,k,QPRES)
                   pl = q(i,j-1-ishft,k,QPRES)
                   dp2 = q(i,j+2-ishft,k,QPRES) - q(i,j-2-ishft,k,QPRES)
                   du = q(i,j+1-ishft,k,QV) - q(i,j-1-ishft,k,QV)
                else
                   pr = q(i,j,k+1-ishft,QPRES)
                   pl = q(i,j,k-1-ishft,QPRES)
                   dp2 = q(i,j,k+2-ishft,QPRES) - q(i,j,k-2-ishft,QPRES)
                   du = q(i,j,k+1-ishft,QW) - q(i,j,k-1-ishft,QW)
                end if
                dp = pr - pl

                denom = max(small_pres, abs(dp2))
                zeta = abs(dp) / denom
                z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

                if (du .le. ZERO) then
                   tst = ONE
                else
                   tst = ZERO
                endif

                tmp = min(pr, pl)

                if ((abs(dp)/tmp) .gt. shktst) then
                   chi2 = tst
                else
                   chi2 = ZERO
                endif

                flatn(i,j,k) = min(flatn(i,j,k), ONE - max(chi2 * z2, chi * z))

             end do

          end do
       end do
    end do

  end subroutine ca_uflatten

end module flatten_module
