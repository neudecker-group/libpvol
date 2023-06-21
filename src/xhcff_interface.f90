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
  type :: xhcff_data

    real(wp) :: xhcff_energy

!> TODO variables go here:
    ! gridsize
    ! probe rad
    ! vdW rad selection (?)
    ! input pressure
    ! allocatable surface calculator
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

  subroutine xhcff_singlepoint(nat,at,xyz,dat,proberad,pressure,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    real(wp),intent(in) :: proberad !> proberadius for sas calculation
    real(wp),intent(in) :: pressure !> pressure in au
    logical,intent(in),optional    :: verbose  !> printout activation
    type(xhcff_data),intent(inout) :: dat  !> collection of xhcff datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io
    logical :: pr

    !> printout activation via verbosity
    if (present(verbose)) then
      pr = verbose
    else
      pr = .false. !> (there is close to no printout anyways)
    end if

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0

    !> singlpoint + gradient call goes here (best would be another module)
    !TODO rewrite to pass only the xhcff_data object
    call xhcff_eg(nat,at,xyz,proberad,pressure,energy,gradient,verbose,iostat)

    if (present(iostat)) then
      iostat = io
    end if

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

!TODO the xhcff_initialize: still needs pressure, gridsize, proberad
  subroutine xhcff_initialize(nat,at,xyz,dat, &
  &                 print,verbose,iunit,iostat)
    character(len=*),parameter :: source = 'xhcff_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(out),optional :: iostat
    !> OUTPUT
    type(xhcff_data),intent(inout) :: dat
    !> LOCAL
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun

!> mapping of optional instuctions
    if (present(print)) then
      pr = print
    else
      pr = .false.
    end if
    if (present(verbose)) then
      pr2 = verbose
    else
      pr2 = .false.
    end if
    if (pr2) pr = pr2
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if

!> Reset datatypes
    call dat%reset()

!> TODO XHCFF calculator setup goes here
    ! allocate surface calculator
    allocate (dat%surf)
    ! call surface calculator setup
    !call surf%setup( ... )

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
    self%xhcff_energy = 0.0_Wp
    ! TODO reset pressure, probe rad, grid size

    if (allocated(self%surf)) deallocate (self%surf)
  end subroutine xhcff_data_deallocate

!========================================================================================!
end module xhcff_interface

