module ppm_module

  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use bl_constants_module, only: ZERO, SIXTH, HALF, ONE, TWO, THREE
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_reconstruct(s, s_lo, s_hi, &
                             flatn, f_lo, f_hi, &
                             sxm, sxp, sym, syp, szm, szp, sd_lo, sd_hi, &
                             ilo1, ilo2, ihi1, ihi2, dx, k3d, kc)

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  sd_lo(3),  sd_hi(3)
    integer, intent(in) ::  f_lo(3),  f_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: k3d, kc

    real(rt)        , intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3))
    real(rt)        , intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))
    real(rt)        , intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(in) :: dx(3)

    call ppm_type1(s, s_lo, s_hi, &
                   flatn, f_lo, f_hi, &
                   sxm, sxp, sym, syp, szm, szp, sd_lo, sd_hi, &
                   ilo1, ilo2, ihi1, ihi2, dx, k3d, kc)


  end subroutine ppm_reconstruct

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

#ifdef CUDA
  attributes(device) &
#endif
  subroutine ppm_type1(s, s_lo, s_hi, &
                       flatn, f_lo, f_hi, &
                       sxm, sxp, sym, syp, szm, szp, sd_lo, sd_hi, &
                       ilo1, ilo2, ihi1, ihi2, dx, k3d, kc)

    use mempool_module, only: bl_allocate, bl_deallocate
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: sd_lo(3), sd_hi(3)
    integer, intent(in) ::  f_lo(3),  f_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: k3d, kc

    real(rt)        , intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3))
    real(rt)        , intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))
    real(rt)        , intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(in) :: dx(3)

    ! local
    integer :: i,j,k

    real(rt) :: dsl, dsr, dsc
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    ! \delta s_{\ib}^{vL}
    real(rt), pointer :: dsvl(:,:)

    real(rt) :: dsvlm, dsvl0, dsvlp

    ! s_{i+\half}^{H.O.}
    real(rt), pointer :: sedge(:,:)

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

    ! cell-centered indexing w/extra ghost cell
    call bl_allocate(dsvl,ilo1-2,ihi1+2,ilo2-2,ihi2+2)

    ! edge-centered indexing
    call bl_allocate(sedge,ilo1-1,ihi1+2,ilo2-1,ihi2+2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc = HALF * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl = TWO  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr = TWO  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl*dsr .gt. ZERO) then
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl(i,j) = ZERO
          end if
       end do
    end do

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       !dir$ ivdep
       do i=ilo1-1,ihi1+2
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do
    end do

    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sm = sedge(i  ,j)
          sp = sedge(i+1,j)

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(i,j,k3d)*sm + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          sp = flatn(i,j,k3d)*sp + (ONE-flatn(i,j,k3d))*s(i,j,k3d)

          ! modify using quadratic limiters -- note this version of the limiting comes
          ! from Colella and Sekora (2008), not the original PPM paper.
          if ((sp-s(i,j,k3d))*(s(i,j,k3d)-sm) .le. ZERO) then
             sp = s(i,j,k3d)
             sm = s(i,j,k3d)

          else if (abs(sp-s(i,j,k3d)) .ge. TWO*abs(sm-s(i,j,k3d))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d) - TWO*sm

          else if (abs(sm-s(i,j,k3d)) .ge. TWO*abs(sp-s(i,j,k3d))) then
          !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d) - TWO*sp
          end if

          sxp(i,j,kc) = sp
          sxm(i,j,kc) = sm

       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc = HALF * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl = TWO  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr = TWO  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl*dsr .gt. ZERO) then
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl(i,j) = ZERO
          end if
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       !dir$ ivdep
       do i=ilo1-1,ihi1+1
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sm = sedge(i,j  )
          sp = sedge(i,j+1)

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(i,j,k3d)*sm + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          sp = flatn(i,j,k3d)*sp + (ONE-flatn(i,j,k3d))*s(i,j,k3d)

          ! modify using quadratic limiters
          if ((sp-s(i,j,k3d))*(s(i,j,k3d)-sm) .le. ZERO) then
             sp = s(i,j,k3d)
             sm = s(i,j,k3d)

          else if (abs(sp-s(i,j,k3d)) .ge. TWO*abs(sm-s(i,j,k3d))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d) - TWO*sm

          else if (abs(sm-s(i,j,k3d)) .ge. TWO*abs(sp-s(i,j,k3d))) then
          !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d) - TWO*sp
          end if

          syp(i,j,kc) = sp
          sym(i,j,kc) = sm

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
          dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
          dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
          dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. ZERO) then
             dsvlm = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvlm = ZERO
          end if

          ! compute on slab above
          k = k3d+1
          dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
          dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
          dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. ZERO) then
             dsvlp = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvlp = ZERO
          end if

          ! compute on current slab
          k = k3d
          dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
          dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
          dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. ZERO) then
             dsvl0 = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl0 = ZERO
          end if

          ! interpolate to lo face
          k = k3d
          sm = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl0-dsvlm)
          ! make sure sedge lies in between adjacent cell-centered values
          sm = max(sm,min(s(i,j,k),s(i,j,k-1)))
          sm = min(sm,max(s(i,j,k),s(i,j,k-1)))

          ! interpolate to hi face
          k = k3d+1
          sp = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvlp-dsvl0)

          ! make sure sedge lies in between adjacent cell-centered values
          sp = max(sp,min(s(i,j,k),s(i,j,k-1)))
          sp = min(sp,max(s(i,j,k),s(i,j,k-1)))

          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(i,j,k3d)*sm + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          sp = flatn(i,j,k3d)*sp + (ONE-flatn(i,j,k3d))*s(i,j,k3d)

          ! modify using quadratic limiters
          if ((sp-s(i,j,k3d))*(s(i,j,k3d)-sm) .le. ZERO) then
             sp = s(i,j,k3d)
             sm = s(i,j,k3d)

          else if (abs(sp-s(i,j,k3d)) .ge. TWO*abs(sm-s(i,j,k3d))) then
          !else if (-(sp-sm)**2/SIX > &
          !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
             sp = THREE*s(i,j,k3d) - TWO*sm

          else if (abs(sm-s(i,j,k3d)) .ge. TWO*abs(sp-s(i,j,k3d))) then
          !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
          !     (sp - sm)**2/SIX) then
             sm = THREE*s(i,j,k3d) - TWO*sp
          end if

          szp(i,j,kc) = sp
          szm(i,j,kc) = sm

       end do
    end do

    call bl_deallocate(dsvl)
    call bl_deallocate(sedge)

  end subroutine ppm_type1

end module ppm_module
