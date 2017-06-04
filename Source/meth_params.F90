
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_castro_method_params().

module meth_params_module

  use bl_error_module, only: bl_error
  use amrex_fort_module, only: rt => amrex_real

implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter :: NHYP = 4

  ! NTHERM: number of thermodynamic variables
  integer, allocatable, save :: NTHERM, NVAR
  integer, allocatable, save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS, UFX

  ! QTHERM: number of primitive variables
  integer, allocatable, save :: QTHERM, QVAR
  integer, allocatable, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  integer, allocatable, save :: NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE
  integer, allocatable, save :: QFA, QFS, QFX

  integer, allocatable, save :: nadv

  ! NQ will be the total number of primitive variables, hydro + radiation
  integer, allocatable, save :: NQ

  integer, allocatable, save :: npassive
  integer, allocatable, save :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, allocatable, save :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME

  integer, save :: numpts_1d

  real(rt), save, allocatable :: outflow_data_old(:,:)
  real(rt), save, allocatable :: outflow_data_new(:,:)
  real(rt), save :: outflow_data_old_time
  real(rt), save :: outflow_data_new_time
  logical,  save :: outflow_data_allocated
  real(rt), save :: max_dist

  ! Create versions of these variables on the GPU
  ! the device update is then done in Castro_nd.F90

  !$acc declare &
  !$acc create(NTHERM, NVAR) &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX) &
  !$acc create(QTHERM, QVAR) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QGAMC, QGAME) &
  !$acc create(NQ) &
  !$acc create(QFA, QFS, QFX)

#ifdef CUDA
  attributes(managed) :: NTHERM, NVAR, nadv
  attributes(managed) :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  attributes(managed) :: QTHERM, QVAR
  attributes(managed) :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  attributes(managed) :: NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE
  attributes(managed) :: QFA, QFS, QFX
  attributes(managed) :: NQ
  attributes(managed) :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
  attributes(managed) :: upass_map, qpass_map, npassive
#endif

  ! Begin the declarations of the ParmParse parameters

  real(rt), allocatable, save :: small_dens
  real(rt), allocatable, save :: small_temp
  real(rt), allocatable, save :: cfl

#ifdef CUDA
  attributes(managed) :: small_dens
  attributes(managed) :: small_temp
  attributes(managed) :: cfl
#endif

  !$acc declare &
  !$acc create(small_dens, small_temp, cfl)

  ! End the declarations of the ParmParse parameters

contains

  subroutine ca_set_castro_method_params() bind(C, name="ca_set_castro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse
    use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync
#endif

    implicit none

#ifdef CUDA
    integer :: istat
#endif

    type (amrex_parmparse) :: pp

    call amrex_parmparse_build(pp, "castro")

    allocate(small_dens)
    small_dens = -1.d200
    allocate(small_temp)
    small_temp = -1.d200
    allocate(cfl)
    cfl = 0.8d0

    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("cfl", cfl)

    !$acc update &
    !$acc device(small_dens, small_temp, cfl)

    call amrex_parmparse_destroy(pp)

  end subroutine ca_set_castro_method_params

end module meth_params_module
