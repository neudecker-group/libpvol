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
module xhcff_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use xhcff_engrad
  use xhcff_surface_module
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: xhcff_calculator
  !public :: xhcff_initialize
  !public :: xhcff_singlepoint
  !public :: print_xhcff_results

!> Main class for interface
  type :: xhcff_calculator

    !> INPUT
    integer :: nat        !> number of atoms
    integer, allocatable :: at(:)   !> atom types
    real(wp), allocatable :: xyz(:,:) !> Cartesian coordinates in Bohr
    real(wp) :: pressure_au !> pressure in a.u.
    real(wp) :: pressure_gpa !> pressure in GPa
    logical :: verbose
    integer :: myunit !> filehandling unit

    real(wp) :: energy
    real(wp), allocatable :: gradient(:,:)

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

  subroutine xhcff_singlepoint(self, energy, gradient, iostat)
    implicit none

    !> DATA CONTAINER
    class(xhcff_calculator),intent(inout) :: self
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(:,:)
    integer,intent(out),optional  :: iostat

    ! data elements
    self%energy = 0.0_wp
    self%gradient(:,:) = 0.0_wp

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

    write (myunit,'(2x, a, t40, f14.4, 1x, a)') "Pressure   ",self%pressure_gpa,"/ GPa   "

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
  &                 gridsize,proberad,verbose,iunit,iostat)
    character(len=*),parameter :: source = 'xhcff_initialize'
    class(xhcff_calculator),intent(inout) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: pressure !> pressure in GPa
    integer,intent(in),optional :: gridsize
    real(wp),intent(in),optional :: proberad !> proberadius for sas calculation
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(out),optional :: iostat
    !> LOCAL
    ! TODO add inputchecks
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun

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

  !> Reset datatypes
    call self%reset()


    !> surface calculator
    allocate (self%surf)
    io = 0
    ! call surface calculator setup with otional parameters
    if (present(gridsize) .and. (present(proberad))) then
      call self%surf%setup(nat,at,xyz,.false., io, ngrid=gridsize, probe=proberad)

    elseif (present(gridsize)) then
      call self%surf%setup(nat,at,xyz,.false., io, ngrid=gridsize)

    elseif (present(proberad)) then
      call self%surf%setup(nat,at,xyz,.false., io, probe=proberad)

    else
      call self%surf%setup(nat,at,xyz,.false., io)
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
      write (myunit,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine xhcff_initialize

!========================================================================================!
  subroutine xhcff_data_deallocate(self)
    implicit none
    class(xhcff_calculator) :: self

    self%energy = 0.0_Wp
    self%nat = 0
    self%pressure_au = 0
    self%pressure_gpa = 0
    self%myunit = 6

    if (allocated(self%surf)) deallocate (self%surf)
    if (allocated(self%gradient)) deallocate(self%gradient)
    if (allocated(self%at)) deallocate(self%at)
    if (allocated(self%xyz)) deallocate(self%xyz)
  end subroutine xhcff_data_deallocate

  ! TODO write update routine

!========================================================================================!
end module xhcff_interface

