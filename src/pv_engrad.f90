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
!> calculate PV term based on lebedev fused balls and analytical gradient
!>

module pv_engrad
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use tesspoints,only:tesspts
  use xhcff_surface_module
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: pv_eg

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  !> pure gradient calculation
  subroutine pv_eg(nat,at,xyz,pressure,surf,energy,gradient,volume)
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
    real(wp) :: xyzt(3),nvec(3) !> container for tesselation coords, normalvector
    integer :: iat, it !> counters
    integer  :: ntess  !> number of tesspoints per atom
    real(wp) :: rdotn  !> dot product of normal an coordinate vector of tessera

    ntess = surf%tess(1)%n
    gradient(:,:) = 0.0_wp
    volume = 0.0_wp

    !> Evaluate Volume via sentence of gauss

    do iat = 1,nat
      do it = 1,ntess
        if (surf%tess(iat)%ap(it) /= 0) then
        !> calc surface normals
        xyzt = surf%tess(iat)%xyz(:,it)
        nvec = xyzt - xyz(:,iat)
        nvec = nvec / sqrt(nvec(1)**2+nvec(2)**2+nvec(3)**2)
        !> dot product of position and normal vector
        rdotn = nvec(1)*xyzt(1) + nvec(2)*xyzt(2) + nvec(3) * xyzt(3)
        !> evaluate volume using the Gauß integral sentence
        volume = volume + rdotn * surf%tess(iat)%ap(it) / 3.0_wp
        
        !> evaluate gradient
        gradient(:,iat) = gradient(:,iat) + (nvec * surf%tess(iat)%ap(it)) / 3.0_wp
        gradient = gradient + (rdotn * surf%tess(iat)%dadr(:,:,it)) / 3.0_wp
        end if
      end do
    end do

    !> Energy is simply PV term
    energy = volume * pressure
    gradient = gradient * pressure

!> Gradient
    ! TODO implement

  end subroutine pv_eg

!========================================================================================!
!========================================================================================!
end module pv_engrad

