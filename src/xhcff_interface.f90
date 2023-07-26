!================================================================================!
! This file is part of xhcfflib.
!
! Copyright (C) 2023 Felix Zeller, Tim Neudecker, Philipp Pracht
!
! xhcfflib is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xhcfflib is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xhcfflib. If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> The typical usage is:
!> 1. declare "use xhcff_interface" in your code
!> 2. Create a xhcff_calculator
!> 3. Initialize it with the %init procedure (returns the first XHCFF grad)
!> 4. Obtain the XHCFF grad with teh %singlepoint procedure

module xhcff_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use xhcff_engrad
  use xhcff_surface_module
  use xhcff_surface_lebedev, only: gridSize
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: xhcff_calculator

!> Main class for interface
  type :: xhcff_calculator

    !> jobdata
    integer :: nat        !> number of atoms
    integer, allocatable :: at(:)   !> atom types
    real(wp), allocatable :: xyz(:,:) !> Cartesian coordinates in Bohr
    real(wp) :: pressure_au !> pressure in a.u.
    real(wp) :: pressure_gpa !> pressure in GPa


    ! IO stuff
    logical :: verbose
    integer :: myunit !> filehandling unit

    !> Output
    real(wp) :: energy
    real(wp), allocatable :: gradient(:,:)

    !> controle variables
    logical :: is_initialized = .false.
    logical :: Bondi !> use Bond radii instead of D3

    !> Errorcode
    integer :: io = 1

!> TODO add vdw rad selection Bondi

    type(surface_calculator),allocatable :: surf

  contains
    procedure :: init => xhcff_initialize
    procedure :: deallocate => xhcff_data_deallocate
    procedure :: reset => xhcff_data_deallocate
    procedure :: info => print_xhcff_results
    procedure :: singlepoint => xhcff_singlepoint
  end type xhcff_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine xhcff_singlepoint(self, nat,at,xyz, energy, gradient, iostat)
    implicit none

    !> DATA CONTAINER
    class(xhcff_calculator),intent(inout) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(:,:)
    integer,intent(out),optional  :: iostat

    !> reset output data elements
    self%energy = 0.0_wp
    self%gradient(:,:) = 0.0_wp

    !> Error handling if not initialized
    if (.not.self%is_initialized) then
      if(self%io == 0) then
        self%io = 1
      end if
      if(present(iostat)) then
        iostat = 1
      end if
    end if

    !> update coordinates
    if(nat /= self%nat .or. any(at.ne.self%at) )then
      if(self%io == 0) then
        self%io = 1
      end if
      if(present(iostat)) then
        iostat = 1
      end if
      return
    endif

    if(self%io /= 0) then
      if(self%verbose) then
        call print_error(self%myunit, self%io)
      end if
      return
    end if

    self%xyz(:,:) = xyz(:,:)

    !> update surface calculator
    call self%surf%update( at, xyz )

    !> singlpoint + gradient calculation
    call xhcff_eg(self%nat, self%at, self%xyz, self%pressure_au, self%surf, self%energy, self%gradient)

    if (self%verbose) then
      call print_xhcff_results(self)
    end if

    !> return singlepoint results
    energy = self%energy
    gradient = self%gradient

  end subroutine xhcff_singlepoint

!========================================================================================!

  subroutine print_xhcff_results(self,iunit)
    class(xhcff_calculator),intent(in) :: self
    integer,intent(in),optional :: iunit ! file handle (usually output_unit=6)
    integer :: myunit, i
    character(len=*),parameter :: outfmt = '(2x,a,f23.12,1x,a)'

    if (present(iunit)) then
      myunit = iunit
    else
      myunit = self%myunit
    end if

    write(myunit,*) '================================================================'
    write(myunit,*) '====================== XHCFF Results ==========================='
    write(myunit,*) '================================================================'

    !> print Errormessage if Error was detected and quit
    if(self%io /= 0) then
      call print_error(self%myunit, self%io)
      return
    end if

    write (myunit,'(2x, a, t40, f14.4, 1x, a)') "Pressure   ",self%pressure_gpa,"/ GPa   "

    !> surface printout
    if (allocated(self%surf)) then
      call self%surf%info(myunit)
    end if

    write (*,*)
    print *, 'XHCFF Gradient:'
    do i=1,self%nat
      write (*,'(2x,i3,3x,3f16.6)') , i, self%gradient(1:3,i)
    end do

  end subroutine print_xhcff_results

!========================================================================================!

  subroutine xhcff_initialize(self,nat,at,xyz,pressure, &
  &                 gridpts,proberad,verbose,iunit,vdwSet, iostat)
    character(len=*),parameter :: source = 'xhcff_initialize'
    class(xhcff_calculator),intent(inout) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: pressure !> pressure in GPa
    integer,intent(in),optional :: gridpts !> gridpoints per atom to construct lebedev grid
    real(wp),intent(in),optional :: proberad !> proberadius for sas calculation
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer, intent(in), optional :: vdwSet !> Set of vdwRad to use: 0 -> D3, 1 -> Bondi
    integer,intent(out),optional :: iostat
    !> LOCAL
    integer :: ich,io,myunit, surferr
    logical :: ex,okbas,pr,pr2
    logical :: exitRun
    real(wp),parameter :: gpatoau = 3.3989309735473356e-05_wp



    !> Reset datatypes
    call self%reset()


    io = 0

    !> check input variables
    if(present(gridpts)) then
      if(ALL(gridSize/= gridpts)) then
        io = 2
      end if
    end if

    if(present(proberad) .and. (proberad .lt. 0.0_wp)) then
      io = 3
    end if

    if(present(vdwSet) .and. (vdwSet == 1)) then
      if(ANY(at > 18)) then
        io = 5
      end if
      else if( present(vdwSet) .and. (vdwSet /= 0)) then
            io = 6
    end if
    write(*,*) 'calling surface calculator'

  !> mapping of optional instuctions
    if (present(verbose)) then
      self%verbose = verbose
    else
      self%verbose = .false.
    end if

    if (present(iunit)) then
      self%myunit = iunit
    else
      self%myunit = stdout
    end if

    if(present(vdwSet) .and. (vdwSet == 1)) then
        self%bondi = .true.
    else
      self%bondi = .false.
    end if


    !> surface calculator
    allocate (self%surf)
    call self%surf%setup(nat,at,xyz,.false., surferr, ngrid=gridpts, probe=proberad, bondi=self%bondi)


    if(surferr /= 0) then
      io = 4
    end if

    !> save input data
    self%pressure_gpa = pressure
    self%pressure_au = pressure * 3.3989309735473356e-05
    self%nat = nat
    allocate (self%at(nat))
    self%at = at
    allocate(self%xyz(3,nat))
    self%xyz = xyz

    !> init calc storage
    self%energy = 0.0_wp
    allocate(self%gradient(3, nat))
    self%gradient = 0.0_wp

    if ((io /= 0).and.pr) then
      call print_error(myunit, io)
      write (myunit,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if

    if (io == 0) then
      self%is_initialized = .True.
    end if

    self%io = io


  end subroutine xhcff_initialize

!========================================================================================!
  subroutine xhcff_data_deallocate(self)
    implicit none
    class(xhcff_calculator) :: self

    self%energy = 0.0_Wp
    self%nat = 0
    self%pressure_au = 0
    self%pressure_gpa = 0
    self%io = 1
    self%myunit = 6
    self%is_initialized = .false.
    self%bondi = .false.

    if (allocated(self%surf)) deallocate (self%surf)
    if (allocated(self%gradient)) deallocate(self%gradient)
    if (allocated(self%at)) deallocate(self%at)
    if (allocated(self%xyz)) deallocate(self%xyz)
  end subroutine xhcff_data_deallocate

  !======================================================================================!

  subroutine print_error(myunit, errcode)
    integer, intent(in) :: myunit, errcode
    integer :: i

    write(myunit,*) 'Error in XHCFF module:'
    select case(errcode)
      case(1)
      write(myunit,*) 'was not initialized before calculation call!'

      case(2)
      write(myunit,*) 'Demanded illegal gridsize!'
      write(myunit,*) 'allowed gridpts per atom:'
      do i=0,2
      write(myunit,*) gridSize(i*8+1:(i+1)*8)
      end do
      case(3)
      write(myunit,*) 'Proberadius cannot be negative!'

      case(4)
      write(myunit,*) 'could not create surface grid!'

      case(5)
      write(myunit,*) 'Bondi VDW radii only impemented for H-Ar, use D3 radii instead!'

      case(6)
      write(myunit,*) 'Demanded set of radii not implemented!'

      case(7)
      write(myunit,*) 'Passed geometry does not match previous one!'
    end select

  end subroutine

!========================================================================================!
end module xhcff_interface

