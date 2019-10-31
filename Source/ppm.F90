module ppm_module

  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use amrex_constants_module, only: ZERO, SIXTH, HALF, ONE, TWO, THREE
  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream
  use castro_module, only: QVAR

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine ca_ppm_reconstruct(lo, hi, &
                                                   s, s_lo, s_hi, &
                                                   flatn, f_lo, f_hi, &
                                                   qm, qm_lo, qm_hi, &
                                                   qp, qp_lo, qp_hi) bind(C, name='ca_ppm_reconstruct')

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)

    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3),QVAR)
    real(rt), intent(in   ) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR,3)

    ! local
    integer :: i, j, k, n, idir

    real(rt) :: q(-2:2)
    real(rt) :: dsl, dsr, dsc
    real(rt) :: dsvl_l, dsvl_r
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    !$acc parallel loop gang vector collapse(4) deviceptr(s, flatn, qm, qp) private(q) async(acc_stream)
    do n = 1, QVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                do idir = 1, 3

                   ! Compute van Leer slopes
                   
                   if (idir == 1) then
                      q(-2:2) = s(i-2:i+2,j,k,n)
                   else if (idir == 2) then
                      q(-2:2) = s(i,j-2:j+2,k,n)
                   else
                      q(-2:2) = s(i,j,k-2:k+2,n)
                   end if

                   dsc = HALF * (q( 0) - q(-2))
                   dsl = TWO  * (q(-1) - q(-2))
                   dsr = TWO  * (q( 0) - q(-1))
                   if (dsl*dsr .gt. ZERO) then
                      dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                   else
                      dsvl_l = ZERO
                   end if

                   dsc = HALF * (q(1) - q(-1))
                   dsl = TWO  * (q(0) - q(-1))
                   dsr = TWO  * (q(1) - q( 0))
                   if (dsl*dsr .gt. ZERO) then
                      dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                   else
                      dsvl_r = ZERO
                   end if

                   ! Interpolate s to x-edges

                   sm = HALF*(q(0) + q(-1)) - SIXTH*(dsvl_r - dsvl_l)

                   ! Make sure sedge lies in between adjacent cell-centered values

                   sm = max(sm, min(q(0), q(-1)))
                   sm = min(sm, max(q(0), q(-1)))

                   ! Compute van Leer slopes

                   dsc = HALF * (q(1) - q(-1))
                   dsl = TWO  * (q(0) - q(-1))
                   dsr = TWO  * (q(1) - q( 0))
                   if (dsl*dsr .gt. ZERO) then
                      dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                   else
                      dsvl_l = ZERO
                   end if

                   dsc = HALF * (q(2) - q(0))
                   dsl = TWO  * (q(1) - q(0))
                   dsr = TWO  * (q(2) - q(1))
                   if (dsl*dsr .gt. ZERO) then
                      dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                   else
                      dsvl_r = ZERO
                   end if

                   ! Interpolate s to x-edges

                   sp = HALF*(q(1) + q(0)) - SIXTH*(dsvl_r - dsvl_l)

                   ! Make sure sedge lies in between adjacent cell-centered values

                   sp = max(sp, min(q(1), q(0)))
                   sp = min(sp, max(q(1), q(0)))

                   ! Flatten the parabola
                   sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*q(0)
                   sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*q(0)

                   ! Modify using quadratic limiters -- note this version of the limiting comes
                   ! from Colella and Sekora (2008), not the original PPM paper.
                   if ((sp-q(0))*(q(0)-sm) .le. ZERO) then

                      sp = q(0)
                      sm = q(0)

                   else if (abs(sp-q(0)) .ge. TWO*abs(sm-q(0))) then

                      sp = THREE*q(0) - TWO*sm

                   else if (abs(sm-q(0)) .ge. TWO*abs(sp-q(0))) then

                      sm = THREE*q(0) - TWO*sp

                   end if

                   if (idir == 1) then
                      qp(i  ,j,k,n,1) = sm
                      qm(i+1,j,k,n,1) = sp
                   else if (idir == 2) then
                      qp(i,j  ,k,n,2) = sm
                      qm(i,j+1,k,n,2) = sp
                   else
                      qp(i,j,k  ,n,3) = sm
                      qm(i,j,k+1,n,3) = sp
                   end if

                end do

             end do
          end do
       end do
    end do

  end subroutine ca_ppm_reconstruct

end module ppm_module
