module bc_fill_module

  implicit none

  public :: hypfill, denfill

  private

contains

  AMREX_LAUNCH subroutine hypfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, &
                                  domlo, domhi, dx, xlo, time, bc)

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    real(rt) :: state(NVAR)
    real(rt) :: staten(NVAR)

    integer  :: i, j, k, n
    integer  :: blo(3), bhi(3), lo(3), hi(3)
    real(rt) :: x, y, z
    logical  :: rho_only

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = adv_l3
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = adv_h3

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, NVAR, domlo, domhi, dx, xlo, bc)

  end subroutine hypfill



  AMREX_LAUNCH subroutine denfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, &
                                  domlo, domhi, dx, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    logical :: rho_only
    integer :: i, j, k
    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = adv_l3
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = adv_h3

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, 1, domlo, domhi, dx, xlo, bc)

  end subroutine denfill

end module bc_fill_module
