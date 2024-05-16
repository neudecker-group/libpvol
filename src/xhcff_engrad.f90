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

!>
!> calculate xhcff gradient as desribed in 10.1063/5.0024671
!>

module xhcff_engrad
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use tesspoints,only:tesspts
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

  !> pure gradient calculation
  subroutine xhcff_eg(nat,at,xyz,pressure,surf,energy,gradient,volume)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    real(wp),intent(in) :: pressure   !> external pressure in au
    type(surface_calculator) :: surf  !> surface information

    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    real(wp),intent(out) :: volume !> volume in bohr ** 3

    !> LOCAL
    real(wp),allocatable :: fvecs(:,:,:)  !> matrix of force vectors
    real(wp) :: xyzt(3),nvec(3),trcorr(3) !> container for tesselation coords, normalvector and translationalcorrection
    integer :: iat,it !> counters
    integer  :: ntess  !> number of tesspoints per atom

    ntess = surf%tess(1)%n
    !allocate (fvecs(3,ntess,nat))
    gradient(:,:) = 0.0_wp
    volume = 0.0_wp
    !> Evaluate gradient via eq. 3 and eq. 4

    do iat = 1,nat
      do it = 1,ntess
        !> calc force vectors, eq 3
        xyzt = surf%tess(iat)%xyz(:,it)
        nvec = xyz(:,iat)-xyzt
        nvec = nvec/sqrt(nvec(1)**2+nvec(2)**2+nvec(3)**2)

        !> add tess point contribution to gradient
        gradient(:,iat) = gradient(:,iat) - nvec*surf%tess(iat)%ap(it)*pressure

        !> evaluate volume using the Gauß integral sentence
        nvec = nvec * -1 !> give normal vectors correct orientation
        volume = volume + (nvec(1)*xyzt(1) + nvec(2)*xyzt(2) + nvec(3) * xyzt(3)) * surf%tess(iat)%ap(it) / 3
      end do
    end do

    !> correct for translations and rotations
    trcorr(1:3) = sum(gradient,dim=2)/nat
    do iat = 1,nat
      gradient(1:3,iat) = gradient(1:3,iat)-trcorr(1:3)
    end do
    energy = 0.0_wp

  end subroutine xhcff_eg

!========================================================================================!
!========================================================================================!
end module xhcff_engrad

