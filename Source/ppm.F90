module ppm_module

  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use amrex_constants_module, only: ZERO, SIXTH, HALF, ONE, TWO, THREE
  use amrex_fort_module, only: rt => amrex_real
  use castro_module, only: QVAR

  implicit none

contains

  AMREX_CUDA_FORT_DEVICE subroutine ca_ppm_reconstruct(lo, hi, &
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
    integer :: i, j, k, n

    real(rt) :: dsl, dsr, dsc
    real(rt) :: dsvl_l, dsvl_r
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, QVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsc = HALF * (s(i  ,j,k,n) - s(i-2,j,k,n))
                dsl = TWO  * (s(i-1,j,k,n) - s(i-2,j,k,n))
                dsr = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to x-edges

                sm = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i-1,j,k,n)))
                sm = min(sm, max(s(i,j,k,n),s(i-1,j,k,n)))

                ! Compute van Leer slopes

                dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i+2,j,k,n) - s(i  ,j,k,n))
                dsl = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                dsr = TWO  * (s(i+2,j,k,n) - s(i+1,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to x-edges

                sp = HALF*(s(i+1,j,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i+1,j,k,n),s(i,j,k,n)))
                sp = min(sp, max(s(i+1,j,k,n),s(i,j,k,n)))

                ! Flatten the parabola
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters -- note this version of the limiting comes
                ! from Colella and Sekora (2008), not the original PPM paper.
                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i  ,j,k,n,1) = sm
                qm(i+1,j,k,n,1) = sp

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, QVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsc = HALF * (s(i,j  ,k,n) - s(i,j-2,k,n))
                dsl = TWO  * (s(i,j-1,k,n) - s(i,j-2,k,n))
                dsr = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to y-edges

                sm = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i,j-1,k,n)))
                sm = min(sm, max(s(i,j,k,n),s(i,j-1,k,n)))

                ! Compute van Leer slopes

                dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i,j+2,k,n) - s(i,j  ,k,n))
                dsl = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                dsr = TWO  * (s(i,j+2,k,n) - s(i,j+1,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to y-edges

                sp = HALF*(s(i,j+1,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i,j+1,k,n),s(i,j,k,n)))
                sp = min(sp, max(s(i,j+1,k,n),s(i,j,k,n)))

                ! Flatten the parabola

                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i,j  ,k,n,2) = sm
                qm(i,j+1,k,n,2) = sp

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, QVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsc = HALF * (s(i,j,k  ,n) - s(i,j,k-2,n))
                dsl = TWO  * (s(i,j,k-1,n) - s(i,j,k-2,n))
                dsr = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to z-edges

                sm = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i,j,k-1,n)))
                sm = min(sm, max(s(i,j,k,n),s(i,j,k-1,n)))

                ! Compute van Leer slopes

                dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsc = HALF * (s(i,j,k+2,n) - s(i,j,k  ,n))
                dsl = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                dsr = TWO  * (s(i,j,k+2,n) - s(i,j,k+1,n))
                if (dsl*dsr .gt. ZERO) then
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to z-edges

                sp = HALF*(s(i,j,k+1,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i,j,k+1,n),s(i,j,k,n)))
                sp = min(sp, max(s(i,j,k+1,n),s(i,j,k,n)))

                ! Flatten the parabola

                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i,j,k  ,n,3) = sm
                qm(i,j,k+1,n,3) = sp

             end do
          end do
       end do
    end do

  end subroutine ca_ppm_reconstruct

end module ppm_module
