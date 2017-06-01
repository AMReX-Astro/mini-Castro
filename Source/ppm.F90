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
  subroutine ppm_reconstruct(s, s_lo, s_hi, h, ilo1, ilo2, ihi1, ihi2, dx, k3d, kc, n)

    implicit none

    integer,  intent(in) ::  s_lo(3),  s_hi(3)
    integer,  intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer,  intent(in) :: k3d, kc, n

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)
    real(rt), intent(in   ) :: dx(3)

    call ppm_type1(s, s_lo, s_hi, h, ilo1, ilo2, ihi1, ihi2, dx, k3d, kc, n)

  end subroutine ppm_reconstruct

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_type1(s, s_lo, s_hi, h, ilo1, ilo2, ihi1, ihi2, dx, k3d, kc, n)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: k3d, kc, n

    type(ht), intent(inout) :: h
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)
    real(rt), intent(in   ) :: dx(3)

    ! local
    integer :: i,j,k

    real(rt) :: dsl, dsr, dsc
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    real(rt) :: dsvlm, dsvl0, dsvlp

#ifndef CUDA
    if (s_lo(1) .gt. ilo1-3 .or. s_lo(2) .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_lo(1), s_lo(2)
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_hi(1) .lt. ihi1+3 .or. s_hi(2) .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_hi(1), s_hi(2)
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call bl_error("Need more ghost cells on array in ppm_type1")
      end if
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc = HALF * (s(i+1,j,k3d,n) - s(i-1,j,k3d,n))
          dsl = TWO  * (s(i  ,j,k3d,n) - s(i-1,j,k3d,n))
          dsr = TWO  * (s(i+1,j,k3d,n) - s(i  ,j,k3d,n))
          if (dsl*dsr .gt. ZERO) then
             h%dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             h%dsvl(i,j) = ZERO
          end if
       end do
    end do

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       !dir$ ivdep
       do i=ilo1-1,ihi1+2
          h%sedge(i,j) = HALF*(s(i,j,k3d,n)+s(i-1,j,k3d,n)) &
               - SIXTH*(h%dsvl(i,j)-h%dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          h%sedge(i,j) = max(h%sedge(i,j),min(s(i,j,k3d,n),s(i-1,j,k3d,n)))
          h%sedge(i,j) = min(h%sedge(i,j),max(s(i,j,k3d,n),s(i-1,j,k3d,n)))
       end do
    end do

    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sm = h%sedge(i  ,j)
          sp = h%sedge(i+1,j)

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = h%flatn(i,j,k3d)*sm + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)
          sp = h%flatn(i,j,k3d)*sp + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)

          ! modify using quadratic limiters -- note this version of the limiting comes
          ! from Colella and Sekora (2008), not the original PPM paper.
          if ((sp-s(i,j,k3d,n))*(s(i,j,k3d,n)-sm) .le. ZERO) then
             sp = s(i,j,k3d,n)
             sm = s(i,j,k3d,n)

          else if (abs(sp-s(i,j,k3d,n)) .ge. TWO*abs(sm-s(i,j,k3d,n))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d,n) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d,n) - TWO*sm

          else if (abs(sm-s(i,j,k3d,n)) .ge. TWO*abs(sp-s(i,j,k3d,n))) then
          !else if ((sp-sm)*(s(i,j,k3d,n) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d,n) - TWO*sp
          end if

          h%sxp(i,j,kc,n) = sp
          h%sxm(i,j,kc,n) = sm

       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc = HALF * (s(i,j+1,k3d,n) - s(i,j-1,k3d,n))
          dsl = TWO  * (s(i,j  ,k3d,n) - s(i,j-1,k3d,n))
          dsr = TWO  * (s(i,j+1,k3d,n) - s(i,j  ,k3d,n))
          if (dsl*dsr .gt. ZERO) then
             h%dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             h%dsvl(i,j) = ZERO
          end if
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       !dir$ ivdep
       do i=ilo1-1,ihi1+1
          h%sedge(i,j) = HALF*(s(i,j,k3d,n)+s(i,j-1,k3d,n)) &
               - SIXTH*(h%dsvl(i,j)-h%dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          h%sedge(i,j) = max(h%sedge(i,j),min(s(i,j,k3d,n),s(i,j-1,k3d,n)))
          h%sedge(i,j) = min(h%sedge(i,j),max(s(i,j,k3d,n),s(i,j-1,k3d,n)))
       end do
    end do

    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sm = h%sedge(i,j  )
          sp = h%sedge(i,j+1)

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = h%flatn(i,j,k3d)*sm + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)
          sp = h%flatn(i,j,k3d)*sp + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)

          ! modify using quadratic limiters
          if ((sp-s(i,j,k3d,n))*(s(i,j,k3d,n)-sm) .le. ZERO) then
             sp = s(i,j,k3d,n)
             sm = s(i,j,k3d,n)

          else if (abs(sp-s(i,j,k3d,n)) .ge. TWO*abs(sm-s(i,j,k3d,n))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d,n) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d,n) - TWO*sm

          else if (abs(sm-s(i,j,k3d,n)) .ge. TWO*abs(sp-s(i,j,k3d,n))) then
          !else if ((sp-sm)*(s(i,j,k3d,n) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d,n) - TWO*sp
          end if

          h%syp(i,j,kc,n) = sp
          h%sym(i,j,kc,n) = sm

       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction

    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! compute on slab below
          k = k3d-1
          dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
          dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          if (dsl*dsr .gt. ZERO) then
             dsvlm = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvlm = ZERO
          end if

          ! compute on slab above
          k = k3d+1
          dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
          dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          if (dsl*dsr .gt. ZERO) then
             dsvlp = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvlp = ZERO
          end if

          ! compute on current slab
          k = k3d
          dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
          dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          if (dsl*dsr .gt. ZERO) then
             dsvl0 = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl0 = ZERO
          end if

          ! interpolate to lo face
          k = k3d
          sm = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvl0-dsvlm)
          ! make sure sedge lies in between adjacent cell-centered values
          sm = max(sm,min(s(i,j,k,n),s(i,j,k-1,n)))
          sm = min(sm,max(s(i,j,k,n),s(i,j,k-1,n)))

          ! interpolate to hi face
          k = k3d+1
          sp = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvlp-dsvl0)

          ! make sure sedge lies in between adjacent cell-centered values
          sp = max(sp,min(s(i,j,k,n),s(i,j,k-1,n)))
          sp = min(sp,max(s(i,j,k,n),s(i,j,k-1,n)))

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = h%flatn(i,j,k3d)*sm + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)
          sp = h%flatn(i,j,k3d)*sp + (ONE-h%flatn(i,j,k3d))*s(i,j,k3d,n)

          ! modify using quadratic limiters
          if ((sp-s(i,j,k3d,n))*(s(i,j,k3d,n)-sm) .le. ZERO) then
             sp = s(i,j,k3d,n)
             sm = s(i,j,k3d,n)

          else if (abs(sp-s(i,j,k3d,n)) .ge. TWO*abs(sm-s(i,j,k3d,n))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d,n) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d,n) - TWO*sm

          else if (abs(sm-s(i,j,k3d,n)) .ge. TWO*abs(sp-s(i,j,k3d,n))) then
          !else if ((sp-sm)*(s(i,j,k3d,n) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d,n) - TWO*sp
          end if

          h%szp(i,j,kc,n) = sp
          h%szm(i,j,kc,n) = sm

       end do
    end do

  end subroutine ppm_type1

end module ppm_module
