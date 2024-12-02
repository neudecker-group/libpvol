!================================================================================!
! This file is part of libpvol.
!
! Copyright (C) 2024 Felix Zeller, Tim Neudecker, Philipp Pracht
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

!> c bindings for libpvol

module libpvol_interface_c
  use iso_c_binding
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  use libpvol_interface
  implicit none
  private

  !> Public C-compatible interface
  public :: c_libpvol_calculator
  public :: c_libpvol_calculator_init
  public :: c_libpvol_calculator_deallocate
  public :: c_libpvol_calculator_singlepoint
  public :: c_libpvol_calculator_info

  !> C-compatible type containing a pointer to the original Fortran type
  type,bind(C) :: c_libpvol_calculator
    !> C will understand fortran types as pointers
    type(c_ptr) :: ptr
  end type c_libpvol_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!>--- C-compatible initialization function
  function c_libpvol_calculator_init(c_nat,c_at,c_xyz,c_pressure,c_model, &
    &                                  c_gridpts,c_proberad,c_verbose,c_printlevel,&
    &                                  c_vdwSet) &
    &                                  result(calculator) &
    &                                  bind(C,name="c_libpvol_calculator_init")
    implicit none
    type(c_libpvol_calculator) :: calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !> WARNING: row-first vs column-first difference  in Fortran and C!
    real(c_double),target,intent(in) :: c_xyz(3,*)
    !> We assume here that a 3-by-x elements are passed, which in C corresponds
    !> to a vector of length nat for x, y and z coordinates respectively
    !> when xyz[3][nat] was defined. 
    !> Hence, it should be defined as xyz[nat][3] in C in order for Fortran
    !> to handle everything correctly in the following!
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
    integer(c_int),value,intent(in) :: c_printlevel

    type(libpvol_calculator),pointer :: calc

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
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat]) !> assumes xyz[nat][3] in C
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
    printlevel = c_printlevel

    !> Allocate and initialize the Fortran calculator
    allocate (calc)
    !> AGAIN, WARNING: There is a row-first vs column-first change in Fortran vs C -> Transpose
    call calc%init(nat,at,transpose(xyz),pressure,trim(model), &
    &              gridpts=gridpts,proberad=proberad,verbose=verbose, printlevel=printlevel,&
    &              iunit=iunit,vdwSet=vdwSet,iostat=iostatus)
    if (iostatus == 0) then
      !> Store the pointer in the C-compatible structure
      calculator%ptr = c_loc(calc)
    else
      write (stderr,'(a,i0)') 'Error initializing '//model//' calculator. code ',iostatus
      calculator%ptr = c_null_ptr
      deallocate (calc)
    end if
  end function c_libpvol_calculator_init

!========================================================================================!

  subroutine c_libpvol_calculator_deallocate(calculator) &
    &     bind(C,name="c_libpvol_calculator_deallocate")
    type(c_libpvol_calculator),intent(inout) :: calculator
    type(libpvol_calculator),pointer :: calc_ptr

    !> Convert the C pointer to a Fortran pointer
    call c_f_pointer(calculator%ptr,calc_ptr)

    !> Deallocate the Fortran object
    if (associated(calc_ptr)) then
      call calc_ptr%deallocate()
      deallocate (calc_ptr)
    end if

    !> Nullify the C pointer
    calculator%ptr = c_null_ptr
  end subroutine c_libpvol_calculator_deallocate

!========================================================================================!

  subroutine c_libpvol_calculator_singlepoint(c_calculator,c_nat,c_at,c_xyz, &
    &                                           c_energy,c_gradient,c_iostat) &
    &                        bind(C,name="c_libpvol_calculator_singlepoint")
    implicit none
    !> Input arguments from C
    type(c_libpvol_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),target,intent(in) :: c_at(*)
    real(c_double),target,intent(in) :: c_xyz(3,*) !> NOTE Fortran/C matrix orders

    !> Output arguments to C
    real(c_double),intent(out) :: c_energy
    real(c_double),target,intent(out) :: c_gradient(3,*) !> NOTE Fortran/C matrix orders
    integer(c_int),intent(out) :: c_iostat

    !> Local Fortran variables
    type(libpvol_calculator),pointer :: calc_ptr
    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp),pointer :: grad(:,:)
    real(wp) :: energy
    integer :: iostat,i,j

    !> Convert C pointers to Fortran pointers
    call c_f_pointer(c_calculator%ptr,calc_ptr)
    call c_f_pointer(c_loc(c_at),at, [c_nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,c_nat]) !> Assumes xyz[nat][3] in C     
    call c_f_pointer(c_loc(c_gradient),grad, [3,c_nat])  !> Assumes grad[nat][3] in C

    !> Set the integer variable
    nat = c_nat

    !> Call the Fortran subroutine
    call calc_ptr%singlepoint(nat,at,xyz,energy,grad,iostat)

    !> Pass back the results to C variables
    c_energy = energy
    c_gradient(1:3,1:nat) = grad(1:3,1:nat)
    c_iostat = iostat

  end subroutine c_libpvol_calculator_singlepoint

!========================================================================================!

  subroutine c_libpvol_calculator_info(c_calculator,c_iunit) &
      & bind(C,name="c_libpvol_calculator_info")
    implicit none
    !> Input arguments from C
    type(c_libpvol_calculator),intent(in) :: c_calculator
    integer(c_int),value,intent(in) :: c_iunit
    !> Local Fortran variables
    type(libpvol_calculator),pointer :: calc_ptr
    integer :: myunit
    !> Convert C pointer to Fortran pointer
    call c_f_pointer(c_calculator%ptr,calc_ptr)
    myunit = c_iunit
    !> Call the Fortran subroutine
    call calc_ptr%info(myunit)
  end subroutine c_libpvol_calculator_info

!========================================================================================!
!========================================================================================!
end module libpvol_interface_c
