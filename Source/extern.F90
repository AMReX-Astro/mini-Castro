module extern_probin_module

  use bl_types
  use bl_space

  implicit none

  private

  real (kind=dp_t), allocatable, save, public :: eos_gamma
  !$acc declare create(eos_gamma)
#ifdef CUDA
  attributes(managed) :: eos_gamma
#endif
  logical, allocatable, save, public :: eos_assume_neutral
  !$acc declare create(eos_assume_neutral)
#ifdef CUDA
  attributes(managed) :: eos_assume_neutral
#endif
  real (kind=dp_t), allocatable, save, public :: small_x
  !$acc declare create(small_x)
#ifdef CUDA
  attributes(managed) :: small_x
#endif

end module extern_probin_module

subroutine runtime_init(name,namlen)

  use extern_probin_module

  implicit none

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /extern/ eos_gamma
  namelist /extern/ eos_assume_neutral
  namelist /extern/ small_x

  allocate(eos_gamma)
  allocate(eos_assume_neutral)
  allocate(small_x)

  eos_gamma = 5.d0/3.d0
  eos_assume_neutral = .true.
  small_x = 1.d-3


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
  !$acc device(eos_gamma, eos_assume_neutral, small_x)

end subroutine runtime_init

