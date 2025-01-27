!================================================================================!
! This file is part of libpvol.
!
! Copyright (C) 2023 Felix Zeller, Tim Neudecker, Philipp Pracht
!
! libpvol is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! libpvol is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with libpvol. If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

program libpvol_main_tester
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use omp_lib
  use libpvol_interface
  use pvol_type_timer
  use xyzreader
  implicit none

  integer :: nat
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  integer :: i,j,k,l
  character(len=40) :: atmp
!========================================================================================!
  real(wp) :: p,probe
  real(wp) :: energy
  real(wp),allocatable :: gradient(:,:)
  real(wp) :: gnorm
  logical :: fail,pr
  integer :: io,grdpts,radii
  type(libpvol_calculator) :: xhcff
  type(libpvol_calculator) :: pv
  type(pvol_timer) :: timer

  character(len=1028) :: inputfile
  integer :: threads,calctype
  real(wp),parameter :: autoaa = 0.529177249_wp
!========================================================================================!

  call print_header()

  !> DEFAULTS
  pr = .true.

  calctype = 0 !> 0 = XHCFF, 1 = XHCFF(+PV)
  p = 1.0_wp
  probe = 1.2_wp !> 1.2 is a typical value for water
  grdpts = 230
  radii = 1 !> 0 = D3, 1 = Bondi
  energy = 0.0_wp
  gnorm = 0.0_wp

!=======================================================================================!
!> parse optional command line args and their defaults
  threads = 1
  call ParseCommandLineArgs(threads,inputfile,calctype,grdpts,radii,probe,p)

!> update setting based on input
#ifdef WITH_OpenMP
  call OMP_Set_Num_Threads(threads)
#ifdef WITH_MKL
  call MKL_Set_Num_Threads(threads)
#endif
  call ompprint_intern(atmp)
  call sleep(1)
#endif

!> read xyz file
  if (.not.exists(inputfile)) then
    error stop 'Input file not found!'
  end if
  call readxyz(inputfile,nat,at,xyz)
  allocate (gradient(3,nat),source=0.0_wp)

  write (*,*) 'Input coords:'
  call printxyz(stdout,nat,at,xyz)
  write (*,*)

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  call timer%new(1,.true.)
  call timer%measure(1,atmp)
  select case (calctype)
  case (0)
    !> XHCFF
    call xhcff%init(nat,at,xyz,p,'XHCFF',proberad=probe,gridpts=grdpts, &
    &    vdwSet=radii,verbose=pr,printlevel=2)
    call xhcff%singlepoint(nat,at,xyz,energy,gradient)

  case (1)
    !> PV
    call pv%init(nat,at,xyz,p,'PV',proberad=probe,gridpts=grdpts, &
    &    vdwSet=radii,verbose=pr,printlevel=2)
    call pv%singlepoint(nat,at,xyz,energy,gradient)

  end select
  write (*,*)
  call timer%measure(1)
  call timer%write_timing(stdout,1)

!=======================================================================================!
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program libpvol_main_tester

!=======================================================================================!
subroutine ompprint_intern(str)
  use omp_lib
  implicit none
  integer :: nproc,TID
  character(len=*) :: str
!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  IF (TID .EQ. 0) THEN
    nproc = OMP_GET_NUM_THREADS()
    write (str,'(a,i0,a)') 'PVol runtime (',nproc,' threads)'
  END IF
!$OMP END PARALLEL
end subroutine ompprint_intern

!=======================================================================================!
subroutine ParseCommandLineArgs(threads,inputfile,calctype,grdpts,radii,probe,pressure)
  use iso_fortran_env,only:wp => real64
  use xyzreader
  use xhcff_surface_lebedev,only:gridSize
  implicit none
  character(len=256) :: arg,arg2 ! Buffer to hold each argument
  integer :: numArgs,i,j     ! Variables to store argument count and loop index
  integer :: io,dumi
  real(wp) :: dum
  !> IN/OUTPUTS
  integer,intent(inout) :: threads
  integer,intent(inout) :: calctype
  integer,intent(inout) :: grdpts
  integer,intent(inout) :: radii
  character(len=*),intent(inout) :: inputfile
  real(wp),intent(inout) :: probe,pressure

  ! Get the number of command-line arguments
  numArgs = COMMAND_ARGUMENT_COUNT()

  ! Check if there are any arguments passed
  if (numArgs == 0) then
    return
  end if

  ! Loop through the command-line arguments
  do i = 1,numArgs
    ! Fetch each argument and store it in 'arg'
    call GET_COMMAND_ARGUMENT(i,arg)

    !> First argument can be input file
    if (i == 1) then
      if (exists(arg)) then
        inputfile = trim(arg)
      end if
    end if

    select case (trim(arg))
    case ('-T','--threads')
      !> parallelization
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dumi
      if (io == 0) threads = dumi

    case ('-i','--input')
      !> input file
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      if (exists(arg2)) inputfile = trim(arg2)

    case ('-xhcff','--xhcff')
      !> switch from XHCFF to XHCFF+(PV)
      calctype = 0

    case ('-pr','--proberad')
      !> Probe radius
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dum
      if (io == 0) probe = dum
      probe = max(probe,0.0_wp)

    case ('-p','--pressure')
      !> pressure in GPa
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dum
      if (io == 0) pressure = dum

    case ('-g','--gridpoints')
      !> Lebedev grid points per atom (the closest match to
      !> predifined Lebedev-Laikov grids is chosen automatically)
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dumi
      if (io == 0) grdpts = dumi
      do j = 1,32
        if (grdpts >= gridSize(j)) dumi = gridSize(j)
      end do
      grdpts = dumi

    case ('-r','--radii')
      !> input file
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      select case (trim(arg2))
      case ('d3')
        radii = 0
      case ('bondi')
        radii = 1
      end select

    case ('-h','--help')
      !> display info about these flags and program usage
      call printhelp()

    end select
  end do

end subroutine ParseCommandLineArgs

subroutine printhelp()
  !> Print a brief description of the program
  write (*,"(1x,a)") "This program performs a single PVol or XHCFF evaluation for a molecule."
  write (*,"(1x,a)") "Usage: pvol <inputfile>.xyz [options]"
  write (*,*)
  write (*,"(1x,a)") "Options (all optional):"
  !> Describing each command-line argument
  write (*,"(1x,a)") " -T, --threads <int>      Number of threads for parallelization."
  write (*,"(1x,a)") " -i, --input <str>        Specify input file (xyz format)."
  write (*,"(1x,a)") "                          Alternatively, the first argument can be the input file."
  write (*,"(1x,a)") " -xhcff, --xchff               Switch from PV to XHCFF(+PV)."
  write (*,"(1x,a)") " -pr, --proberad <real>   Set the probe radius In Angström."
  write (*,"(1x,a)") " -p, --pressure <real>    Set the pressure in GPa."
  write (*,"(1x,a)") " -g, --gridpoints <int>   Number of Lebedev grid points per atom."
  write (*,"(1x,a)") ' -r, --radii <str>        Choose reference radii. Choice of "bondi" (default) or "d3"'
  write (*,"(1x,a)") " -h, --help               Display this help message and exit."

  !> Examples of common usage
  write (*,*)
  write (*,"(1x,a)") "Examples:"
  write (*,"(1x,a)") " xhcff inputfile.xyz -p 2.0"
  write (*,"(1x,a)") " xhcff -i inputfile.xyz -T 4"
  write (*,"(1x,a)") " xhcff --pressure 1.0 --proberad 1.2"

  stop
end subroutine printhelp

subroutine print_header()
  implicit none
  write (*,*)
  write (*,"(10x,a)") repeat('+',42)
  write (*,"(10x,a)") "+"//repeat(' ',12)//"PVol App v0.0.1"//repeat(' ',12)//'+'
  write (*,"(10x,a)") repeat('+',42)
  write (*,"(10x,a)") "Authors:  F.Zeller, T.Neudecker, P.Pracht"
  write (*,*)

end subroutine print_header
