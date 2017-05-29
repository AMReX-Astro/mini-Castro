
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
  integer, save :: NTHERM, NVAR
  integer, save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS, UFX

  ! QTHERM: number of primitive variables
  integer, save :: QTHERM, QVAR
  integer, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  integer, save :: NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE
  integer, save :: QFA, QFS, QFX

  integer, save :: nadv

  ! NQ will be the total number of primitive variables, hydro + radiation
  integer, save :: NQ         

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME

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
  integer, device :: NTHERM_d, NVAR_d
  integer, device :: URHO_d, UMX_d, UMY_d, UMZ_d, UMR_d, UML_d, UMP_d, UEDEN_d, UEINT_d, UTEMP_d, UFA_d, UFS_d, UFX_d
  integer, device :: QTHERM_d, QVAR_d
  integer, device :: QRHO_d, QU_d, QV_d, QW_d, QPRES_d, QREINT_d, QTEMP_d, QGAME_d
  integer, device :: NQAUX_d, QGAMC_d, QC_d, QCSML_d, QDPDR_d, QDPDE_d
  integer, device :: QFA_d, QFS_d, QFX_d
  integer, device :: NQ_d
  integer, device :: npassive_d
  integer, allocatable, device :: upass_map_d(:), qpass_map_d(:)
#endif

  ! Begin the declarations of the ParmParse parameters

  real(rt), save :: small_dens
  real(rt), save :: small_temp
  real(rt), save :: cfl

#ifdef CUDA
  real(rt), device :: small_dens_d
  real(rt), device :: small_temp_d
  real(rt), device :: cfl_d
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

    small_dens = -1.d200;
    small_temp = -1.d200;
    cfl = 0.8d0;

    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("cfl", cfl)

#ifdef CUDA
  istat = cudaMemcpyAsync(small_dens_d, small_dens, 1)
  istat = cudaMemcpyAsync(small_temp_d, small_temp, 1)
  istat = cudaMemcpyAsync(cfl_d, cfl, 1)
#endif

    !$acc update &
    !$acc device(small_dens, small_temp, cfl)

    call amrex_parmparse_destroy(pp)

  end subroutine ca_set_castro_method_params

end module meth_params_module
