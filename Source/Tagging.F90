module tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_denerror(tag,taglo,taghi, &
                         set,clear, &
                         den,denlo,denhi, &
                         lo,hi,nd,domlo,domhi, &
                         delta,xlo,problo,time,level) &
                         bind(C, name="ca_denerror")

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

    real(rt), parameter :: dengrad_rel = 0.25d0

    ! Tag on regions of high density gradient
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
             ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
             az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
             ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
             ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
             az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
             if (MAX(ax,ay,az) .ge. ABS(dengrad_rel * den(i,j,k,1))) then
                tag(i,j,k) = set
             end if
          end do
       end do
    end do

  end subroutine ca_denerror

end module tagging_module
