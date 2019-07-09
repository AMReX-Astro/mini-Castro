module tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_denerror(lo, hi, &
                         tag, taglo, taghi, &
                         den, denlo, denhi, &
                         set, clear) &
                         bind(C, name="ca_denerror")

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: denlo(3), denhi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3))
    integer(1), intent(in   ), value :: set, clear

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    real(rt), parameter :: dengrad_rel = 0.25d0

    !$gpu

    ! Tag on regions of high density gradient
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ax = ABS(den(i+1,j,k) - den(i,j,k))
             ay = ABS(den(i,j+1,k) - den(i,j,k))
             az = ABS(den(i,j,k+1) - den(i,j,k))
             ax = MAX(ax,ABS(den(i,j,k) - den(i-1,j,k)))
             ay = MAX(ay,ABS(den(i,j,k) - den(i,j-1,k)))
             az = MAX(az,ABS(den(i,j,k) - den(i,j,k-1)))
             if (MAX(ax,ay,az) .ge. ABS(dengrad_rel * den(i,j,k))) then
                tag(i,j,k) = set
             end if
          end do
       end do
    end do

  end subroutine ca_denerror

end module tagging_module
