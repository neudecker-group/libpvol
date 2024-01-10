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
!> 1. declare "use pv_interface" in your code
!> 2. Create a pv_calculator
!> 3. Initialize it with the %init procedure (returns the first energy and grad)
!> 4. Obtain the energy and grad with the %singlepoint procedure

module pv_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use pv_engrad
  use xhcff_surface_module
  use xhcff_surface_lebedev,only:gridSize
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: pv_calculator

!> Main class for interface
  type :: pv_calculator

    !> System data
    integer :: nat                   !> number of atoms
    integer,allocatable :: at(:)     !> atom types
    real(wp),allocatable :: xyz(:,:) !> Cartesian coordinates in Bohr
    real(wp) :: pressure_au          !> pressure in a.u.
    real(wp) :: pressure_gpa         !> pressure in GPa

    ! IO stuff
    logical :: verbose
    integer :: printlevel !> amount of printout
    integer :: myunit !> filehandling unit

    !> Output
    real(wp) :: energy
    real(wp) :: volume !> volume in bohr ** 3
    real(wp),allocatable :: gradient(:,:)

    !> controle variables
    logical :: is_initialized = .false.
    logical :: Bondi !> use Bond radii instead of D3

    !> Errorcode
    integer :: io = 1

    type(surface_calculator),allocatable :: surf

  contains
    procedure :: init => pv_initialize
    procedure :: deallocate => pv_data_deallocate
    procedure :: reset => pv_data_deallocate
    procedure :: info => print_pv_results
    procedure :: singlepoint => pv_singlepoint
  end type pv_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine pv_singlepoint(self,nat,at,xyz,energy,gradient,iostat)
    implicit none

    !> DATA CONTAINER
    class(pv_calculator),intent(inout) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(:,:)
    integer,intent(out),optional  :: iostat
    
    real(wp),parameter :: geotol = 1.0e-7_wp


    !> reset output data elements
    self%energy = 0.0_wp
    self%gradient(:,:) = 0.0_wp

    !> Error handling if not initialized
    if (.not.self%is_initialized) then
      if (self%io == 0) then
        self%io = 1
      end if
      if (present(iostat)) then
        iostat = 1
      end if
    end if

    !> update coordinates
    if (nat /= self%nat.or.any(at .ne. self%at)) then
      if (self%io == 0) then
        self%io = 1
      end if
      if (present(iostat)) then
        iostat = 1
      end if
      return
    end if

    if (self%io /= 0) then
      if (self%verbose) then
        call print_error(self%myunit,self%io)
      end if
      return
    end if

    !> check if we need to update the geometry?
    if(any(abs(self%xyz(:,:) - xyz(:,:)) > geotol))then
      self%xyz(:,:) = xyz(:,:)

      !> update surface calculator
      call self%surf%update(at,xyz)
    endif

    !> singlpoint + gradient calculation
    call pv_eg(self%nat,self%at,self%xyz,self%pressure_au,self%surf,self%energy,self%gradient, self%volume)

    if (self%verbose) then
      call print_pv_results(self)
    end if

    !> return singlepoint results
    energy = self%energy
    gradient = self%gradient

  end subroutine pv_singlepoint

!========================================================================================!
  subroutine print_pv_results(self,iunit)
    class(pv_calculator),intent(in) :: self
    integer,intent(in),optional :: iunit !> file handle (usually output_unit=6)
    integer :: myunit,i
    character(len=*),parameter :: outfmt = '(2x,a,f23.12,1x,a)'

    if (present(iunit)) then
      myunit = iunit
    else
      myunit = self%myunit
    end if

    !> print Errormessage if Error was detected and quit
    if (self%io /= 0) then
      call print_error(self%myunit,self%io)
      return
    end if

    if(self%printlevel >= 1) then
      write (myunit,*) '================================================================'
      write (myunit,*) '==================== PV MODULE Results ========================='
      write (myunit,*) '================================================================'

      write (myunit,'(2x, a, t40, f14.4, 1x, a)') "Pressure   ",self%pressure_gpa,"/ GPa   "

      !> surface printout
      if (allocated(self%surf)) then
        call self%surf%info(myunit)
      end if
    end if

    !> always print Volume and Energy
    write (myunit,'(2x, a, t40, f14.4, 1x, a)') "Volume   ",self%volume ,"/ Bohr ** 3  "
    write (myunit,'(2x, a, t40, f14.4, 1x, a)') "Energy   ",self%energy ,"/ Eh "
    if(self%printlevel >= 2) then
      write (myunit,*)
      write (myunit,'(a)') '> PV Gradient ( Eh/a0 ):'
      do i = 1,self%nat
        write (myunit,'(2x,i3,3x,3f16.6)'),i,self%gradient(1:3,i)
      end do
    end if
    end subroutine print_pv_results

!========================================================================================!
  subroutine pv_initialize(self,nat,at,xyz,pressure, &
  &                 gridpts,proberad,scaling,verbose,iunit,vdwSet,printlevel,iostat)
    character(len=*),parameter :: source = 'xhcff_initialize'
    class(pv_calculator),intent(inout) :: self
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)        !> coordinates in Bohr
    real(wp),intent(in) :: pressure          !> pressure in GPa
    integer,intent(in),optional :: gridpts   !> gridpoints per atom to construct lebedev grid
    real(wp),intent(in),optional :: proberad !> proberadius for sas calculation in angstrom
    real(wp),intent(in),optional :: scaling  !> scaling of vdw radii to simulate sas
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(in),optional  :: printlevel !> printlevel, > 1 full printout
    integer,intent(in),optional :: vdwSet !> Set of vdwRad to use: 0 -> D3, 1 -> Bondi
    integer,intent(out),optional :: iostat
    !> LOCAL
    integer :: ich,io,myunit,surferr
    logical :: ex,okbas,pr,pr2
    logical :: exitRun
    real(wp),parameter :: gpatoau = 3.3989309735473356e-05_wp

    !> Reset datatypes
    call self%reset()
    io = 0

    !>
    !> mapping of optional instuctions
    !>

    if (present(verbose)) then
      self%verbose = verbose
    else
      self%verbose = .false.
    end if

    if(present(printlevel)) then
      self%printlevel = printlevel
      else
        self%printlevel = 0
    end if

    if (present(iunit)) then
      self%myunit = iunit
    else
      self%myunit = stdout
    end if

    !>
    !> check input variables
    !>

    if (present(gridpts).and.(io == 0)) then
      if (ALL(gridSize /= gridpts)) then
        io = 2
      end if
    end if

    if (present(proberad).and.(io == 0)) then
      if (proberad .lt. 0.0_wp) then
        io = 3
      end if
    end if

    !> (optional) selection of vdW radii
    if (present(vdwSet).and.(io == 0)) then
      if (vdwSet == 1) then
        if (ANY(at > 18)) then
          io = 5
        end if
        self%bondi = .true.
      elseif (vdwSet /= 0) then
        io = 6
        self%bondi = .false.
      end if
    end if

    if (present(scaling) .and. (io == 0)) then
      if (scaling .lt. 0.0_wp) then
        io = 8
      end if
    end if

    if (self%verbose) write (self%myunit,'(a)') '> XHCFF: calling surface calculator'

    !> surface calculator
    if (io == 0) then
      allocate (self%surf)
      call self%surf%setup(nat,at,xyz,.false.,surferr,ngrid=gridpts, &
      &    probe=proberad,scaling=scaling,bondi=self%bondi)

      if (surferr /= 0) then
        io = 4
      end if
    end if

    !> save input data
    self%pressure_gpa = pressure
    self%pressure_au = pressure*gpatoau
    self%nat = nat
    allocate (self%at(nat))
    self%at = at
    allocate (self%xyz(3,nat))
    self%xyz = xyz

    !> init calc storage
    self%energy = 0.0_wp
    self%volume = 0.0_wp
    allocate (self%gradient(3,nat))
    self%gradient = 0.0_wp

    if ((io /= 0).and.pr) then
      call print_error(myunit,io)
      write (myunit,'("Could not create PV calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if

    if (io == 0) then
      self%is_initialized = .True.
    end if

    self%io = io

  end subroutine pv_initialize

!========================================================================================!
  subroutine pv_data_deallocate(self)
    implicit none
    class(pv_calculator) :: self

    self%energy = 0.0_Wp
    self%volume = 0.0_wp
    self%nat = 0
    self%pressure_au = 0
    self%pressure_gpa = 0
    self%io = 1
    self%printlevel=0
    self%myunit = 6
    self%is_initialized = .false.
    self%bondi = .false.

    if (allocated(self%surf)) deallocate (self%surf)
    if (allocated(self%gradient)) deallocate (self%gradient)
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
  end subroutine pv_data_deallocate

!======================================================================================!
  !> TODO Move to own file
  subroutine print_error(myunit,errcode)
    integer,intent(in) :: myunit,errcode
    integer :: i

    write (myunit,*) 'Error in PV module:'
    select case (errcode)
    case (1)
      write (myunit,*) 'was not initialized before calculation call!'

    case (2)
      write (myunit,*) 'Demanded illegal gridsize!'
      write (myunit,*) 'allowed gridpts per atom:'
      do i = 0,2
        write (myunit,*) gridSize(i*8+1:(i+1)*8)
      end do
    case (3)
      write (myunit,*) 'Proberadius cannot be negative!'

    case (4)
      write (myunit,*) 'could not create surface grid!'

    case (5)
      write (myunit,*) 'Bondi VDW radii only impemented for H-Ar, use D3 radii instead!'

    case (6)
      write (myunit,*) 'Demanded set of radii not implemented!'

    case (7)
      write (myunit,*) 'Passed geometry does not match previous one!'

    case(8)
      write (myunit,*) 'Scaling cannot be negative!'
    end select

  end subroutine print_error

!========================================================================================!
!========================================================================================!
end module pv_interface

