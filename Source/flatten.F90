module flatten_module

  use amrex_acc_module, only: acc_stream

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine uflatten(i, j, k, q, q_lo, q_hi, flatn) bind(C, name='uflatten')

    !$acc routine seq

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: QVAR, QU, QV, QW, QPRES, small_pres

    implicit none

    integer,  intent(in   ) :: i, j, k
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: flatn

    integer :: ishft, idir

    real(rt) :: denom, denominv, zeta, tst, tmp
    real(rt) :: pl, pr, dp, dp2, du, z, z2, chi, chi2

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE / (zcut2 - zcut1)

    flatn = ONE

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

       if (dp .gt. ZERO) then
          ishft = 1
       else
          ishft = -1
       endif

       denom = max(small_pres, abs(dp2))
       denominv = ONE / denom
       zeta = abs(dp) * denominv
       z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

       if (du .le. ZERO) then
          tst = ONE
       else
          tst = ZERO
       endif

       tmp = min(pr, pl)

       if (abs(dp) .gt. shktst * tmp) then
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
       denominv = ONE / denom
       zeta = abs(dp) * denominv
       z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

       if (du .le. ZERO) then
          tst = ONE
       else
          tst = ZERO
       endif

       tmp = min(pr, pl)

       if (abs(dp) .gt. shktst * tmp) then
          chi2 = tst
       else
          chi2 = ZERO
       endif

       flatn = min(flatn, ONE - max(chi2 * z2, chi * z))

    end do

  end subroutine uflatten

end module flatten_module
