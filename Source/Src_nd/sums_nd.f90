module castro_sums_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains
   
  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx,&
                        vol,v_lo,v_hi,mass) bind(C, name="ca_summass")

    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k

    mass = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mass = mass + rho(i,j,k) * vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_summass

end module castro_sums_module
