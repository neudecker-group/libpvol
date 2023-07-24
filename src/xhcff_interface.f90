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
  public :: xhcff_data
  public :: xhcff_initialize
  public :: xhcff_singlepoint
  public :: print_xhcff_results

!> this type bundles together most of the
!> data required for a XHCFF calculation
  !TODO rewrite as class
  type :: xhcff_data

    !> INPUT
    integer :: nat        !> number of atoms
    integer, allocatable :: at(:)   !> atom types
    real(wp), allocatable :: xyz(:,:) !> Cartesian coordinates in Bohr
    real(wp) :: pressure_au !> pressure in a.u.
    real(wp) :: pressure_gpa !> pressure in GPa
    logical :: verbose

    real(wp) :: energy
    real(wp), allocatable :: gradient(:,:)

!> TODO add vdw rad selection Bondi

    type(surface_calculator),allocatable :: surf

  contains
    procedure :: deallocate => xhcff_data_deallocate
    procedure :: reset => xhcff_data_deallocate
    procedure :: info => print_xhcff_results
  end type xhcff_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine xhcff_singlepoint(dat, energy, gradient, iostat)
    implicit none

    !> DATA CONTAINER
    type(xhcff_data),intent(inout) :: dat  !> collection of xhcff datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(:,:)
    integer,intent(out),optional  :: iostat

    ! data elements
    dat%energy = 0.0_wp
    dat%gradient(:,:) = 0.0_wp

    !> singlpoint + gradient call
    call xhcff_eg(dat%nat, dat%at, dat%xyz, dat%pressure_au, dat%surf, dat%energy, dat%gradient)

    !TODO add printout of gradient if verbose
    !> return singlepoint results from data container
    energy = dat%energy
    gradient = dat%gradient


  end subroutine xhcff_singlepoint
!========================================================================================!

  subroutine print_xhcff_results(dat,iunit)
    class(xhcff_data),intent(in) :: dat
    integer,intent(in),optional :: iunit ! file handle (usually output_unit=6)
    integer :: myunit
    character(len=*),parameter :: outfmt = '(2x,a,f23.12,1x,a)'
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if
!TODO also add printout for selected pressure, e.g.
!    write (myunit,outfmt) "Pressure   ",dat%pressure,"GPa   "
    if (allocated(dat%surf)) then
      call dat%surf%info(myunit)
    end if
  end subroutine print_xhcff_results
!========================================================================================!

  !TODO set defaults and do error handling
  subroutine xhcff_initialize(nat,at,xyz,pressure,dat, &
  &                 gridsize,proberad,print,verbose,iunit,iostat)
    character(len=*),parameter :: source = 'xhcff_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: pressure !> pressure in GPa
    integer,intent(in),optional :: gridsize
    real(wp),intent(in),optional :: proberad !> proberadius for sas calculation
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(out),optional :: iostat
    !> OUTPUT
    type(xhcff_data),intent(inout) :: dat
    !> LOCAL
    ! TODO this seems like Error handling to me, shall I include that as well
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun

  !> mapping of optional instuctions
    if (present(print)) then
      dat%verbose = print
    else
      dat%verbose = .false.
    end if

    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if

  !> Reset datatypes
    call dat%reset()


    !> surface calculator
    allocate (dat%surf)
    io = 0
    ! call surface calculator setup with otional parameters
    if (present(gridsize) .and. (present(proberad))) then
      call dat%surf%setup(nat,at,xyz,pr, io, ngrid=gridsize, probe=proberad)

    elseif (present(gridsize)) then
      call dat%surf%setup(nat,at,xyz,pr, io, ngrid=gridsize)

    elseif (present(proberad)) then
      call dat%surf%setup(nat,at,xyz,pr, io, probe=proberad)

    else
      call dat%surf%setup(nat,at,xyz,pr, io)
    end if

    !> save input data
    dat%pressure_gpa = pressure
    dat%pressure_au = pressure * 3.3989309735473356e-05
    dat%nat = nat
    allocate (dat%at(nat))
    dat%at = at
    allocate(dat%xyz(3,nat))
    dat%xyz = xyz

    !> init calc storage
    dat%energy = 0.0_wp
    allocate(dat%gradient(3, nat))
    dat%gradient = 0.0_wp

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
    class(xhcff_data) :: self

    self%energy = 0.0_Wp
    self%nat = 0
    self%pressure_au = 0
    self%pressure_gpa = 0

    if (allocated(self%surf)) deallocate (self%surf)
    if (allocated(self%gradient)) deallocate(self%gradient)
    if (allocated(self%at)) deallocate(self%at)
    if (allocated(self%xyz)) deallocate(self%xyz)
  end subroutine xhcff_data_deallocate

  ! TODO write update routine

!========================================================================================!
end module xhcff_interface

