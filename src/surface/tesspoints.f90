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
! along with xhcfflib.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

! Created by felix on 6/12/23.
module tesspoints
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none
  private
  public :: tesspts !, fpt, printout

  !> tesselation points of one atom
  type tesspts
    !> number of tesselation points of this atom
    integer :: n = 0

    !> surface normals
    real(wp),allocatable :: xyz(:,:)

    !> areas of corresponding tesselation points
    real(wp),allocatable :: ap(:)

    !> derivative of tesspoint areas wrt atoms
    real(wp),allocatable :: dadr(:,:,:)

  contains
    procedure :: allocate => allocate_tsspt

  end type tesspts

!=========================================================================================!
!=========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!

  subroutine allocate_tsspt(self,n,nat)
    implicit none
    class(tesspts) :: self
    integer,intent(in) :: n
    integer,intent(in) :: nat
    if (allocated(self%xyz).and.allocated(self%ap).and.n == self%n) then
      !> Don't re-allocate if tesspts is already set up.
      !> There can be quite a lot of tessalation points so this allocation
      !> could create some overhead when updating the surface frequently
      return

    else !> fallback if n changed somehow
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%ap)) deallocate (self%ap)
      if (allocated(self%dadr)) deallocate (self%dadr)
    end if

    self%n = n
    allocate (self%xyz(3,n),source=0.0_wp)
    allocate (self%ap(n),source=0.0_wp)
    allocate (self%dadr(3,nat,n),source=0.0_wp)
  end subroutine allocate_tsspt

!=========================================================================================!
!=========================================================================================!
end module tesspoints
