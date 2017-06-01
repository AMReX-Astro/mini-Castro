module flatten_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine uflaten(lo, hi, q, q_lo, q_hi, h)

    use bl_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dg
    use advection_util_module, only: ht
    use meth_params_module, only: NQ, QU, QV, QW, QPRES

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    type(ht), intent(inout) :: h

    integer :: i, j, k, ishft

    real(rt), parameter :: small_pres = 1.e-200_rt

    real(rt) :: denom, zeta, tst, tmp, ftmp

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    ! x-direction flattening coef
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          !dir$ ivdep
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             h%dp(i,j,k) = q(i+1*dg(1),j,k,QPRES) - q(i-1*dg(1),j,k,QPRES)
             denom = max(small_pres,abs(q(i+2*dg(1),j,k,QPRES)-q(i-2*dg(1),j,k,QPRES)))
             zeta = abs(h%dp(i,j,k))/denom
             h%z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i-1*dg(1),j,k,QU)-q(i+1*dg(1),j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i+1*dg(1),j,k,QPRES),q(i-1*dg(1),j,k,QPRES))
             if ((abs(h%dp(i,j,k))/tmp).gt.shktst) then
                h%chi(i,j,k) = tst
             else
                h%chi(i,j,k) = ZERO
             endif
          enddo
          do i = lo(1),hi(1)
             if(h%dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             h%flatn(i,j,k) = ONE - &
                  max(h%chi(i-ishft*dg(1),j,k)*h%z(i-ishft*dg(1),j,k),h%chi(i,j,k)*h%z(i,j,k))
          enddo
       enddo
    enddo

    ! y-direction flattening coef
    do k = lo(3),hi(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          !dir$ ivdep
          do i = lo(1),hi(1)
             h%dp(i,j,k) = q(i,j+1*dg(2),k,QPRES) - q(i,j-1*dg(2),k,QPRES)
             denom = max(small_pres,abs(q(i,j+2*dg(2),k,QPRES)-q(i,j-2*dg(2),k,QPRES)))
             zeta = abs(h%dp(i,j,k))/denom
             h%z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i,j-1*dg(2),k,QV)-q(i,j+1*dg(2),k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i,j+1*dg(2),k,QPRES),q(i,j-1*dg(2),k,QPRES))
             if ((abs(h%dp(i,j,k))/tmp).gt.shktst) then
                h%chi(i,j,k) = tst
             else
                h%chi(i,j,k) = ZERO
             endif
          enddo
       end do
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if(h%dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(h%chi(i,j-ishft*dg(2),k)*h%z(i,j-ishft*dg(2),k),h%chi(i,j,k)*h%z(i,j,k))
             h%flatn(i,j,k) = min( h%flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

    ! z-direction flattening coef
    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          !dir$ ivdep
          do i = lo(1),hi(1)
             h%dp(i,j,k) = q(i,j,k+1*dg(3),QPRES) - q(i,j,k-1*dg(3),QPRES)
             denom = max(small_pres,abs(q(i,j,k+2*dg(3),QPRES)-q(i,j,k-2*dg(3),QPRES)))
             zeta = abs(h%dp(i,j,k))/denom
             h%z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i,j,k-1*dg(3),QW)-q(i,j,k+1*dg(3),QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i,j,k+1*dg(3),QPRES),q(i,j,k-1*dg(3),QPRES))
             if ((abs(h%dp(i,j,k))/tmp).gt.shktst) then
                h%chi(i,j,k) = tst
             else
                h%chi(i,j,k) = ZERO
             endif
          enddo
       enddo
    enddo
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if(h%dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(h%chi(i,j,k-ishft*dg(3))*h%z(i,j,k-ishft*dg(3)),h%chi(i,j,k)*h%z(i,j,k))
             h%flatn(i,j,k) = min( h%flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

  end subroutine uflaten

end module flatten_module
