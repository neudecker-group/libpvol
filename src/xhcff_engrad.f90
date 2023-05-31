!================================================================================!
! This file is part of xhcfflib.
!
! Copyright (C) 2023 
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
module xhcff_engrad
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: xhcff_eg


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine xhcff_eg(nat,at,xyz,params,surface,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    real(wp),intent(in) :: params(:)  
    real(wp),intent(in) :: surface(:,:)
    logical,intent(in),optional    :: verbose  !> printout activation 
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io
    logical :: pr

    !> printout activation via verbosity
    if(present(verbose))then
      pr = verbose
    else
      pr =.false. !> (there is close to no printout anyways)
    endif

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0

    !> singlpoint + gradient call goes here (best would be another module)



    if (present(iostat)) then
      iostat = io
    end if

  end subroutine xhcff_eg

!========================================================================================!
!========================================================================================!
end module xhcff_engrad

