module ppm_module

  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use bl_constants_module, only: ZERO, SIXTH, HALF, ONE, TWO, THREE
  use amrex_fort_module, only: rt => amrex_real
  use advection_util_module, only: ht
  use meth_params_module, only: NQ

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_reconstruct(s, s_lo, s_hi, h, lo, hi)

    implicit none

    integer,  intent(in) :: s_lo(3), s_hi(3)
    integer,  intent(in) :: lo(3), hi(3)

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)

    call ppm_type1(s, s_lo, s_hi, h, lo, hi)

  end subroutine ppm_reconstruct

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_type1(s, s_lo, s_hi, h, lo, hi)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)

    ! local
    integer :: i, j, k, n

    real(rt) :: dsl, dsr, dsc
    real(rt) :: dsvl_l, dsvl_r
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

#ifndef CUDA
    if (s_lo(1) .gt. lo(1)-3 .or. s_lo(2) .gt. lo(2)-3 .or. s_lo(3) .gt. lo(3)-3) then
         print *,'Low bounds of array: ',s_lo(1), s_lo(2),s_lo(3)
         print *,'Low bounds of  loop: ',lo(1),lo(2),lo(3)
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_hi(1) .lt. hi(1)+3 .or. s_hi(2) .lt. hi(2)+3 .or. s_hi(3) .lt. hi(3)+3) then
         print *,'Hi  bounds of array: ',s_hi(1), s_hi(2), s_hi(3)
         print *,'Hi  bounds of  loop: ',hi(1),hi(2),hi(3)
         call bl_error("Need more ghost cells on array in ppm_type1")
      end if
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

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
                sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
                sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

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

                h%sxp(i,j,k,n) = sp
                h%sxm(i,j,k,n) = sm

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

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

                sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
                sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                h%syp(i,j,k,n) = sp
                h%sym(i,j,k,n) = sm

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

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

                sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
                sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                h%szp(i,j,k,n) = sp
                h%szm(i,j,k,n) = sm

             end do
          end do
       end do
    end do

  end subroutine ppm_type1

end module ppm_module
