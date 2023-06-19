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

!> calculate xhcff gradient as desribed in 10.1063/5.0024671
module xhcff_engrad
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use tesspoints, only: tesspts
  use xhcff_surface_module
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: xhcff_eg


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine xhcff_eg(nat,at,xyz,proberad,pressure,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    real(wp),intent(in) :: proberad !> radius of solvent
    real(wp),intent(in) :: pressure !> external pressure in au
    !real(wp),intent(in) :: params(:)
    !real(wp),intent(in) :: surface(:,:)
    logical,intent(in),optional    :: verbose  !> printout activation 
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat

    !> LOCAL
    type(surface_calculator) :: surf
    real(wp), allocatable :: fvecs(:,:,:) !> matrix of force vectirs
    real(wp) :: xyzt(3), nvec(3), trcorr(3) !> container for tesselation coords, normalvector and translationalcorrection
    integer :: io, iat, it !> counters
    logical :: pr
    integer  :: ntess        !> number of tesspoints per atom

    !> printout activation via verbosity
    if(present(verbose))then
      pr = verbose
    else
      pr =.false. !> (there is close to no printout anyways)
    endif

    !init surface calculation
    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%normal, probe=proberad)
    ntess = surf%tess(1)%n
    allocate(fvecs(3, ntess, nat))

    !evaluate eq. 3
    do iat = 1, nat
            do it = 1, ntess
            ! calc force vectors, eq 3
              xyzt = surf%tess(iat)%xyz(:,it)
              nvec = xyz(:,iat) - xyzt
              nvec = nvec / sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)
              fvecs(:, it, iat) = nvec * surf%tess(iat)%ap(it) * pressure
            end do
    end do


    !calc gradient, eq. 4
    gradient = sum(fvecs, dim=2)

    ! correct for translations and rotations
    trcorr(1:3) = sum(gradient, dim=2) / nat
    do iat=1, nat
      gradient(1:3, iat) = gradient(1:3, iat) - trcorr(1:3)
    end do

    energy = 0.0_wp

    io = 0
    call surf%deallocate()
    !> singlpoint + gradient call goes here (best would be another module)



    if (present(iostat)) then
      iostat = io
    end if

  end subroutine xhcff_eg

!========================================================================================!
!========================================================================================!
end module xhcff_engrad

