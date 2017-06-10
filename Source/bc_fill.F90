module bc_fill_module

  implicit none

  public :: hypfill, denfill

  private

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine hypfill(blo, bhi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc)

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real
    use filcc_module, only: filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: blo(3), bhi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    real(rt) :: state(NVAR)
    real(rt) :: staten(NVAR)

    integer  :: i, j, k, n, lo(3), hi(3)
    real(rt) :: x, y, z
    logical  :: rho_only

    do n = 1,NVAR
       call filccn(blo, bhi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, dx, xlo, bc, n)
    enddo

    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the proper
    ! upwinding.  However, this decision makes this routine look
    ! somewhat non-orthodox, in that we need to set external values in
    ! either case....how do we know it's Outflow?  We have to assume
    ! that the setup routines converted Outflow to FOEXTRAP.

    ! Set flag for bc function
    rho_only = .FALSE.

    !     XLO
    if ( (bc(1,1,1) == EXT_DIR .or. bc(1,1,1) == FOEXTRAP) .and. (blo(1) < domlo(1)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (i < domlo(1)) then
                   do n = 1, NVAR
                      state(n) = adv(domlo(1),j,k,n)
                   end do
                   call bcnormal(state,staten,1,+1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   end do
                end if
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1) == EXT_DIR .or. bc(1,2,1) == FOEXTRAP) .and. (bhi(1) > domhi(1)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (i > domhi(1)) then
                   do n = 1, NVAR
                      state(n) = adv(domhi(1),j,k,n)
                   end do
                   call bcnormal(state,staten,1,-1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   end do
                end if
             end do
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1,1) == EXT_DIR .or. bc(2,1,1) == FOEXTRAP) .and. (blo(2) < domlo(2)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (j < domlo(2)) then
                   do n = 1, NVAR
                      state(n) = adv(i,domlo(2),k,n)
                   end do
                   call bcnormal(state,staten,2,+1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   end do
                end if
             end do
          end do
       end do
    end if

    !     YHI
    if ( (bc(2,2,1) == EXT_DIR .or. bc(2,2,1) == FOEXTRAP) .and. (bhi(2) > domhi(2)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (j > domhi(2)) then
                   do n = 1, NVAR
                      state(n) = adv(i,domhi(2),k,n)
                   end do
                   call bcnormal(state,staten,2,-1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   enddo
                end if
             end do
          end do
       end do
    end if

    !     ZLO
    if ( (bc(3,1,1) == EXT_DIR .or. bc(3,1,1) == FOEXTRAP) .and. (blo(3) < domlo(3)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (k < domlo(3)) then
                   do n = 1, NVAR
                      state(n) = adv(i,j,domlo(3),n)
                   end do
                   call bcnormal(state,staten,3,+1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   end do
                end if
             end do
          end do
       end do
    end if

    !     ZHI
    if ( (bc(3,2,1) == EXT_DIR .or. bc(3,2,1) == FOEXTRAP) .and. (bhi(3) > domhi(3)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (k > domhi(3)) then
                   do n = 1, NVAR
                      state(n) = adv(i,j,domhi(3),n)
                   end do
                   call bcnormal(state,staten,3,-1,rho_only)
                   do n = 1, NVAR
                      adv(i,j,k,n) = staten(n)
                   end do
                end if
             end do
          end do
       end do
    end if

  end subroutine hypfill



#ifdef CUDA
  attributes(device) &
#endif
  subroutine denfill(blo, bhi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real
    use filcc_module, only: filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: blo(3), bhi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    logical :: rho_only
    integer :: i, j, k

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filccn(blo, bhi, adv, adv_lo, adv_hi, 1, domlo, domhi, dx, xlo, bc, 1)

    rho_only = .TRUE.

    !     XLO
    if ( (bc(1,1,1) == EXT_DIR .or. bc(1,1,1) == FOEXTRAP) .and. (blo(1) < domlo(1)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (i < domlo(1)) then
                   call bcnormal(adv(domlo(1),j,k),adv(i,j,k),1,+1,rho_only)
                end if
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1) == EXT_DIR .or. bc(1,2,1) == FOEXTRAP) .and. (bhi(1) > domhi(1)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (i > domhi(1)) then
                   call bcnormal(adv(domhi(1),j,k),adv(i,j,k),1,-1,rho_only)
                end if
             end do
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1,1) == EXT_DIR .or. bc(2,1,1) == FOEXTRAP) .and. (blo(2) < domlo(2)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (j < domlo(2)) then
                   call bcnormal(adv(i,domlo(2),k),adv(i,j,k),2,+1,rho_only)
                end if
             end do
          end do
       end do
    end if

    !     YHI
    if ( (bc(2,2,1) == EXT_DIR .or. bc(2,2,1) == FOEXTRAP) .and. (bhi(2) > domhi(2)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (j > domhi(2)) then
                   call bcnormal(adv(i,domhi(2),k),adv(i,j,k),2,-1,rho_only)
                end if
             end do
          end do
       end do
    end if

    !     ZLO
    if ( (bc(3,1,1) == EXT_DIR .or. bc(3,1,1) == FOEXTRAP) .and. (blo(3) < domlo(3)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (k < domlo(3)) then
                   call bcnormal(adv(i,j,domlo(3)),adv(i,j,k),3,+1,rho_only)
                end if
             end do
          end do
       end do
    end if

    !     ZHI
    if ( (bc(3,2,1) == EXT_DIR .or. bc(3,2,1) == FOEXTRAP) .and. (bhi(3) > domhi(3)) ) then
       do i = blo(1), bhi(1)
          do j = blo(2), bhi(2)
             do k = blo(3), bhi(3)
                if (k > domhi(3)) then
                   call bcnormal(adv(i,j,domhi(3)),adv(i,j,k),3,-1,rho_only)
                end if
             end do
          end do
       end do
    end if
  end subroutine denfill



#ifdef CUDA
  attributes(device) &
#endif
  subroutine bcnormal(u_int,u_ext,dir,sgn,rho_only)

    use probdata_module, only: dens_ambient, p_ambient
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in   ) :: u_int(*)
    real(rt), intent(inout) :: u_ext(*)
    logical,  intent(in   ) :: rho_only
    integer,  intent(in   ) :: dir, sgn

    real(rt) :: rho, rhou(3), eden, T, Y
    integer  :: n, t1, t2, i

    ! for the Sedov problem, we will always set the state to the ambient conditions

    if (rho_only) then

       u_ext(1) = dens_ambient

    else

       ! First set everything from internal data (this is probably a bad
       ! thing to do...)  That is, we should have explicit boundary data
       ! for advected fields and species

       do i = 1, NVAR
          u_ext(i) = u_int(i)
       enddo

       u_ext(URHO)   = dens_ambient
       u_ext(UMX)    = 0.e0_rt
       u_ext(UMY)    = 0.e0_rt
       u_ext(UMZ)    = 0.e0_rt
       !u_ext(UEDEN)  = p_ambient/(gamma_const-1.e0_rt)
       !u_ext(UEINT)  = u_ext(UEDEN)

    endif

  end subroutine bcnormal

end module bc_fill_module
