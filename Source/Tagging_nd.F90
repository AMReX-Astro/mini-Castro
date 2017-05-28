module tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), save ::    denerr,   dengrad
  real(rt), save ::    enterr,   entgrad
  real(rt), save ::    velerr,   velgrad
  real(rt), save ::   temperr,  tempgrad
  real(rt), save ::  presserr, pressgrad
  real(rt), save ::    raderr,   radgrad
  integer,  save ::  max_denerr_lev,   max_dengrad_lev
  integer,  save ::  max_enterr_lev,   max_entgrad_lev
  integer,  save ::  max_velerr_lev,   max_velgrad_lev
  integer,  save ::  max_temperr_lev,  max_tempgrad_lev
  integer,  save ::  max_presserr_lev, max_pressgrad_lev
  integer,  save ::  max_raderr_lev,   max_radgrad_lev

contains

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the density
  ! ::: -----------------------------------------------------------

  subroutine ca_denerror(tag,taglo,taghi, &
                         set,clear, &
                         den,denlo,denhi, &
                         lo,hi,nd,domlo,domhi, &
                         delta,xlo,problo,time,level) &
                         bind(C, name="ca_denerror")

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: set, clear, nd, level
    integer,  intent(in   ) :: taglo(3), taghi(3)
    integer,  intent(in   ) :: denlo(3), denhi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    integer,  intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt), intent(in   ) :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3),nd)
    real(rt), intent(in   ) :: delta(3), xlo(3), problo(3), time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    ! Tag on regions of high density
    if (level .lt. max_denerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (den(i,j,k,1) .ge. denerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    ! Tag on regions of high density gradient
    if (level .lt. max_dengrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(den(i+1*dg(1),j,k,1) - den(i,j,k,1))
                ay = ABS(den(i,j+1*dg(2),k,1) - den(i,j,k,1))
                az = ABS(den(i,j,k+1*dg(3),1) - den(i,j,k,1))
                ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. dengrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_denerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the temperature
  ! ::: -----------------------------------------------------------

  subroutine ca_temperror(tag,taglo,taghi, &
                          set,clear, &
                          temp,templo,temphi, &
                          lo,hi,np,domlo,domhi, &
                          delta,xlo,problo,time,level) &
                          bind(C, name="ca_temperror")

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: set, clear, np, level
    integer,  intent(in   ) :: taglo(3), taghi(3)
    integer,  intent(in   ) :: templo(3), temphi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    integer,  intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt), intent(in   ) :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
    real(rt), intent(in   ) :: delta(3), xlo(3), problo(3), time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    ! Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (temp(i,j,k,1) .ge. temperr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    ! Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(temp(i+1*dg(1),j,k,1) - temp(i,j,k,1))
                ay = ABS(temp(i,j+1*dg(2),k,1) - temp(i,j,k,1))
                az = ABS(temp(i,j,k+1*dg(3),1) - temp(i,j,k,1))
                ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. tempgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_temperror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the pressure
  ! ::: -----------------------------------------------------------

  subroutine ca_presserror(tag,taglo,taghi, &
                           set,clear, &
                           press,presslo,presshi, &
                           lo,hi,np,domlo,domhi, &
                           delta,xlo,problo,time,level) &
                           bind(C, name="ca_presserror")

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: set, clear, np, level
    integer,  intent(in   ) :: taglo(3), taghi(3)
    integer,  intent(in   ) :: presslo(3), presshi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    integer,  intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt), intent(in   ) :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
    real(rt), intent(in   ) :: delta(3), xlo(3), problo(3), time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    ! Tag on regions of high pressure
    if (level .lt. max_presserr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (press(i,j,k,1) .ge. presserr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    ! Tag on regions of high pressure gradient
    if (level .lt. max_pressgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
                ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
                az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
                ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. pressgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_presserror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the velocity
  ! ::: -----------------------------------------------------------

  subroutine ca_velerror(tag,taglo,taghi, &
                         set,clear, &
                         vel,vello,velhi, &
                         lo,hi,nv,domlo,domhi, &
                         delta,xlo,problo,time,level) &
                         bind(C, name="ca_velerror")

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: set, clear, nv, level
    integer,  intent(in   ) :: taglo(3), taghi(3)
    integer,  intent(in   ) :: vello(3), velhi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    integer,  intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt), intent(in   ) :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),nv)
    real(rt), intent(in   ) :: delta(3), xlo(3), problo(3), time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    !     Tag on regions of high velocity
    if (level .lt. max_velerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (vel(i,j,k,1) .ge. velerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high velocity gradient
    if (level .lt. max_velgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(vel(i+1*dg(1),j,k,1) - vel(i,j,k,1))
                ay = ABS(vel(i,j+1*dg(2),k,1) - vel(i,j,k,1))
                az = ABS(vel(i,j,k+1*dg(3),1) - vel(i,j,k,1))
                ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. velgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_velerror

end module tagging_module
