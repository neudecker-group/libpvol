!================================================================================!
! This file is part of xhcfflib.
!
! Copyright (C) 2024 Felix Zeller, Tim Neudecker, Philipp Pracht
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

!> c bindings for the xhcff_lib

module xhcfflib_interface_c
  use iso_c_binding
  use iso_fortran_env,only:wp => real64,stdout=>output_unit,stderr=>error_unit
  use xhcfflib_interface
  implicit none
  private

  !> Public C-compatible interface
  public :: c_xhcfflib_calculator
  public :: c_xhcfflib_calculator_init
  ! public :: c_xhcfflib_calculator_singlepoint, &
  !           c_xhcfflib_calculator_deallocate, c_xhcfflib_calculator_print

  !> C-compatible type containing a pointer to the original Fortran type
  type :: c_xhcfflib_calculator
    !> C will understand fortran types as pointers
    type(xhcfflib_calculator),pointer :: ptr
  end type c_xhcfflib_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!>--- C-compatible initialization function
  function c_xhcfflib_calculator_init(c_nat,c_at,c_xyz,c_pressure,c_model, &
    &                                  c_gridpts,c_proberad,c_verbose, &
    &                                  c_vdwSet) &
    &                                  result(libptr) &
    &                                  bind(C,name="c_xhcfflib_calculator_init")
    implicit none
    type(c_ptr) :: libptr
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !> WARNING: row-first vs column-first difference  in Fortran and C!
    real(c_double),target,intent(in) :: c_xyz(3,*)
    !> We assume here that a 3-by-x matrix is passed, which in C corresponds
    !> to a vector of length nat for x, y and z coordinates respectively
    !> The Fortran interface expects nat vectors of length 3, however.
    !> This means, in the actual call later we have to pass the transpose!
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    real(c_double),value,intent(in) :: c_pressure
    !character(kind=c_char,len=1),intent(in) :: c_model(*)
    integer(c_int),value,intent(in) :: c_model !> NOTE: IT IS A SWITCHASE INTEGER HERE
    integer(c_int),value,intent(in) :: c_gridpts
    real(c_double),value,intent(in) :: c_proberad
    !real(c_double),value,intent(in) :: c_scaling
    logical(c_bool),value,intent(in) :: c_verbose
    !integer(c_int),value,intent(in) :: c_iunit
    integer(c_int),value,intent(in) :: c_vdwSet
    !integer(c_int),value,intent(in) :: c_printlevel

    type(c_xhcfflib_calculator),pointer :: calc_ptr
    type(xhcfflib_calculator),pointer :: calc

    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp) :: pressure
    character(len=:),allocatable :: model
    integer :: gridpts
    real(wp) :: proberad
    real(wp) :: scaling
    logical :: verbose
    integer :: iunit,vdwSet,printlevel,iostatus
    integer :: i,modelint
    
    !> Convert C arguments to Fortran types
    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [nat,3])
    pressure = c_pressure
    modelint = c_model
    select case (modelint)
    case (0)
      model = "XHCFF"
    case (1)
      model = "PV"
    case default
      model = "unkown"
    end select
    gridpts = c_gridpts
    proberad = c_proberad
    !scaling = c_scaling
    verbose = logical(c_verbose,kind=c_int)
    iunit = stdout !> use a default here. 
    vdwSet = c_vdwSet
    !printlevel = c_printlevel
    printlevel = 0

    !> Allocate and initialize the Fortran calculator
    allocate (calc)
    !> AGAIN, WARNING: There is a row-first vs column-first change in Fortran vs C -> Transpose
    call calc%init(nat,at,transpose(xyz),pressure,trim(model), &
    &              gridpts=gridpts,proberad=proberad,verbose=verbose,&
    &              iunit=iunit,vdwSet=vdwSet,iostat=iostatus)
    if (iostatus == 0) then
      !> Store the pointer in the C-compatible structure
      allocate (calc_ptr)
      calc_ptr%ptr => calc
      libptr = c_loc(calc_ptr)
    else
      write(stderr,'(a,i0)') 'Error initializing '//model//' calculator. code ',iostatus
      libptr = c_null_ptr
      deallocate (calc)
    end if
  end function c_xhcfflib_calculator_init

!========================================================================================!

!========================================================================================!
!========================================================================================!
end module xhcfflib_interface_c
