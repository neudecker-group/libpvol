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

module interface_c
  use iso_c_binding,only:c_int,c_double,c_ptr
  use xhcfflib_interface,only:xhcfflib_calculator,xhcfflib_calculator_init,xhcfflib_calculator_singlepoint
  implicit none
  private

  !> Public C-compatible interface
  public :: c_xhcfflib_calculator,c_xhcfflib_calculator_init,c_xhcfflib_calculator_singlepoint

  !> C-compatible type containing the original Fortran type
  type,bind(C) :: c_xhcfflib_calculator
    type(xhcfflib_calculator) :: calc
  end type c_xhcfflib_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  !> C-compatible initialization subroutine
  subroutine c_xhcfflib_calculator_init(c_calculator, &
  &          c_nat,c_at,c_xyz,c_pressure,c_model, &
  &          c_gridpts,c_proberad,c_scaling,c_verbose, &
  &          c_iunit,c_vdwSet,c_printlevel,c_iostat) bind(C,name="c_xhcfflib_calculator_init")
    use iso_c_binding,only:c_int,c_double,c_char,c_ptr,c_f_pointer
    type(c_xhcfflib_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),intent(in) :: c_at(*)
    real(c_double),intent(in) :: c_xyz(3,*)
    real(c_double),value,intent(in) :: c_pressure
    character(kind=c_char,len=*),intent(in) :: c_model
    integer(c_int),value,intent(in),optional :: c_gridpts
    real(c_double),value,intent(in),optional :: c_proberad
    real(c_double),value,intent(in),optional :: c_scaling
    logical(c_int),value,intent(in),optional :: c_verbose
    integer(c_int),value,intent(in),optional :: c_iunit
    integer(c_int),value,intent(in),optional :: c_vdwSet
    integer(c_int),value,intent(in),optional :: c_printlevel
    integer(c_int),intent(out),optional :: c_iostat

    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp) :: pressure
    character(len=*) :: model
    integer,optional :: gridpts
    real(wp),optional :: proberad
    real(wp),optional :: scaling
    logical,optional :: verbose
    integer,optional :: iunit
    integer,optional :: vdwSet
    integer,optional :: printlevel
    integer,optional :: iostat

    ! Convert C arguments to Fortran types
    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat])
    pressure = c_pressure
    model = trim(c_model)  ! Handle the conversion from C string to Fortran string

    if (present(c_gridpts)) gridpts = c_gridpts
    if (present(c_proberad)) proberad = c_proberad
    if (present(c_scaling)) scaling = c_scaling
    if (present(c_verbose)) verbose = logical(c_verbose,kind=c_int)
    if (present(c_iunit)) iunit = c_iunit
    if (present(c_vdwSet)) vdwSet = c_vdwSet
    if (present(c_printlevel)) printlevel = c_printlevel
    if (present(c_iostat)) iostat = c_iostat

    ! Call the original Fortran subroutine
    call c_calculator%calc%init(nat,at,xyz,pressure,model,gridpts, &
    &    proberad,scaling,verbose,iunit,vdwSet,printlevel,iostat)
  end subroutine c_xhcfflib_calculator_init

!=========================================================================================!

  subroutine c_xhcfflib_calculator_print(c_calculator, &
      &      c_iunit) bind(C,name="c_xhcfflib_calculator_print")
    use iso_c_binding,only:c_int
    type(c_xhcfflib_calculator),intent(in) :: c_calculator
    integer(c_int),value,intent(in),optional :: c_iunit

    integer,optional :: iunit

    ! Convert C arguments to Fortran types
    if (present(c_iunit)) iunit = c_iunit

    ! Call the original Fortran subroutine
    call c_calculator%calc%print_xhcff_results(iunit)
  end subroutine c_xhcfflib_calculator_print

!=========================================================================================!

  subroutine c_xhcfflib_calculator_deallocate(c_calculator) &
      & bind(C,name="c_xhcfflib_calculator_deallocate")
    use iso_c_binding,only:c_int
    type(c_xhcfflib_calculator),intent(inout) :: c_calculator

    ! Call the original Fortran deallocate subroutine
    call c_calculator%calc%deallocate()
  end subroutine c_xhcfflib_calculator_deallocate

!=========================================================================================!

  subroutine c_xhcfflib_calculator_singlepoint(c_calculator, &
      &  c_nat,c_at,c_xyz,c_energy,c_gradient,c_iostat) &
      &  bind(C,name="c_xhcfflib_calculator_singlepoint")
    use iso_c_binding,only:c_int,c_double,c_ptr,c_f_pointer
    type(c_xhcfflib_calculator),intent(inout) :: c_calculator
    integer(c_int),value,intent(in) :: c_nat
    integer(c_int),intent(in) :: c_at(*)
    real(c_double),intent(in) :: c_xyz(3,*)
    real(c_double),intent(out) :: c_energy
    real(c_double),intent(out) :: c_gradient(3,*)
    integer(c_int),intent(out),optional :: c_iostat

    integer :: nat
    integer,pointer :: at(:)
    real(wp),pointer :: xyz(:,:)
    real(wp) :: energy
    real(wp),pointer :: gradient(:,:)
    integer,optional :: iostat

    ! Convert C arguments to Fortran types
    nat = c_nat
    call c_f_pointer(c_loc(c_at),at, [nat])
    call c_f_pointer(c_loc(c_xyz),xyz, [3,nat])
    call c_f_pointer(c_loc(c_gradient),gradient, [3,nat])

    ! Call the original Fortran subroutine
    call c_calculator%calc%singlepoint(nat,at,xyz,energy,gradient,iostat)

    ! Return the results to C variables
    c_energy = energy
    if (present(c_iostat)) c_iostat = iostat

  end subroutine c_xhcfflib_calculator_singlepoint

!========================================================================================!
!========================================================================================!
end module interface_c

