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
  subroutine ppm_reconstruct(s, s_lo, s_hi, h, lo, hi, dx, n)

    implicit none

    integer,  intent(in) :: s_lo(3), s_hi(3)
    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: n

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)
    real(rt), intent(in   ) :: dx(3)

    call ppm_type1(s, s_lo, s_hi, h, lo, hi, dx, n)

  end subroutine ppm_reconstruct

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_type1(s, s_lo, s_hi, h, lo, hi, dx, n)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: n

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)
    real(rt), intent(in   ) :: dx(3)

    ! local
    integer :: i, j, k

    real(rt) :: dsl, dsr, dsc
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

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-2, hi(1)+2
             dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
             dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
             dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
             if (dsl*dsr .gt. ZERO) then
                h%dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                h%dsvl(i,j,k) = ZERO
             end if
          end do
       end do
    end do

    ! interpolate s to x-edges
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          !dir$ ivdep
          do i = lo(1)-1, hi(1)+2
             h%sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) &
                            - SIXTH*(h%dsvl(i,j,k)-h%dsvl(i-1,j,k))
             ! make sure sedge lies in between adjacent cell-centered values
             h%sedge(i,j,k) = max(h%sedge(i,j,k),min(s(i,j,k,n),s(i-1,j,k,n)))
             h%sedge(i,j,k) = min(h%sedge(i,j,k),max(s(i,j,k,n),s(i-1,j,k,n)))
          end do
       end do
    end do

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sm = h%sedge(i  ,j,k)
             sp = h%sedge(i+1,j,k)

             ! flatten the parabola
             sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
             sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters -- note this version of the limiting comes
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-2, hi(2)+2
          do i = lo(1)-1, hi(1)+1
             dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
             dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
             dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
             if (dsl*dsr .gt. ZERO) then
                h%dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                h%dsvl(i,j,k) = ZERO
             end if
          end do
       end do
    end do

    ! interpolate s to y-edges
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+2
          !dir$ ivdep
          do i = lo(1)-1, hi(1)+1
             h%sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) &
                            - SIXTH*(h%dsvl(i,j,k)-h%dsvl(i,j-1,k))
             ! make sure sedge lies in between adjacent cell-centered values
             h%sedge(i,j,k) = max(h%sedge(i,j,k),min(s(i,j,k,n),s(i,j-1,k,n)))
             h%sedge(i,j,k) = min(h%sedge(i,j,k),max(s(i,j,k,n),s(i,j-1,k,n)))
          end do
       end do
    end do

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sm = h%sedge(i,j  ,k)
             sp = h%sedge(i,j+1,k)

             ! flatten the parabola
             sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
             sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    do k = lo(3)-2, hi(3)+2
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
             dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
             dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
             if (dsl*dsr .gt. ZERO) then
                h%dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                h%dsvl(i,j,k) = ZERO
             end if
          end do
       end do
    end do

    ! interpolate s to z-edges
    do k = lo(3)-1, hi(3)+2
       do j = lo(2)-1, hi(2)+1
          !dir$ ivdep
          do i = lo(1)-1, hi(1)+1

             h%sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) &
                            - SIXTH*(h%dsvl(i,j,k)-h%dsvl(i,j,k-1))
             ! make sure sedge lies in between adjacent cell-centered values
             h%sedge(i,j,k) = max(h%sedge(i,j,k),min(s(i,j,k,n),s(i,j,k-1,n)))
             h%sedge(i,j,k) = min(h%sedge(i,j,k),max(s(i,j,k,n),s(i,j,k-1,n)))

          end do
       end do
    end do
       

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sm = h%sedge(i,j,k  )
             sp = h%sedge(i,j,k+1)

             ! flatten the parabola
             sm = h%flatn(i,j,k)*sm + (ONE-h%flatn(i,j,k))*s(i,j,k,n)
             sp = h%flatn(i,j,k)*sp + (ONE-h%flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters
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

  end subroutine ppm_type1

end module ppm_module
