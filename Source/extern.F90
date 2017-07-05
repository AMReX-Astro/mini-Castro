module extern_probin_module

  use bl_types
  use bl_space

  implicit none

  private

  logical, allocatable, public :: use_eos_coulomb
  logical, allocatable, public :: eos_input_is_constant
  real (kind=dp_t), allocatable, public :: small_x
#ifdef CUDA
  attributes(managed) :: small_x
  attributes(managed) :: eos_input_is_constant
  attributes(managed) :: use_eos_coulomb
#endif

end module extern_probin_module

subroutine runtime_init(name,namlen)

  use extern_probin_module

#ifdef CUDA
  use cudafor, only: cudaMemAdvise, cudaMemAdviseSetPreferredLocation
  use cuda_module, only: cuda_device_id
#endif

  implicit none

#ifdef CUDA
  integer :: cuda_result
#endif

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /extern/ use_eos_coulomb
  namelist /extern/ eos_input_is_constant
  namelist /extern/ small_x

  allocate(use_eos_coulomb)
  allocate(eos_input_is_constant)
  allocate(small_x)

  use_eos_coulomb = .true.
  eos_input_is_constant = .false.
  small_x = 1.d-30


  ! create the filename
  if (namlen > maxlen) then
     print *, 'probin file name too long'
     stop
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=extern, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     print *, 'ERROR: problem in the extern namelist'
     stop
  endif

  close (unit=un)

  !$acc update &
  !$acc device(use_eos_coulomb, eos_input_is_constant, small_x)

#ifdef CUDA
  cuda_result = cudaMemAdvise(use_eos_coulomb, 1, cudaMemAdviseSetPreferredLocation, cuda_device_id)
  cuda_result = cudaMemAdvise(eos_input_is_constant, 1, cudaMemAdviseSetPreferredLocation, cuda_device_id)
  cuda_result = cudaMemAdvise(small_x, 1, cudaMemAdviseSetPreferredLocation, cuda_device_id)
#endif

end subroutine runtime_init

