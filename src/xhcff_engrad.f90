!================================================================================!
! This file is part of libpvol.
!
! Copyright (C) 2023
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

!>
!> calculate xhcff gradient as desribed in 10.1063/5.0024671
!> calculate SAS implemented as eq 6-10 in 10.1021/acs.jctc.1c00471
!>

module xhcff_engrad
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none
  private

  !> Smoothing dielectric function parameters
  real(wp),parameter :: autoaa = 0.52917726_wp
  real(wp),parameter :: aatoau = 1.0_wp/autoaa
  real(wp),parameter :: w = 0.3_wp*aatoau
  real(wp),parameter :: w3 = w*(w*w)
  ! const part of eq 9
  real(wp),parameter :: ah0 = 0.5_wp
  ! const 2nd part of eq 9
  real(wp),parameter :: ah1 = 3.0_wp/(4.0_wp*w)
  ! const 3rd part of eq 6
  real(wp),parameter :: ah3 = -1.0_wp/(4.0_wp*w3)

  !> real space cut-offs
  real(wp),parameter :: tolsesp = 1.e-6_wp

  !> pi
  real(wp),parameter :: pi = 3.14159265359

!> routines/datatypes that can be seen outside the module
  public :: xhcff_eg

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  !> pure gradient calculation
  subroutine xhcff_eg(nat,nnsas,nnlists,xyz,vdwsa, &
    & wrp,trj2,angWeight,angGrid,pressure,area,volume,energy,grad)
    implicit none

    !> Number of atoms
    integer,intent(in) :: nat

    !> Number of neighbours to consider
    integer,intent(in) :: nnsas(:)

    !> Neighbourlist
    integer,intent(in) :: nnlists(:,:)

    !> Cartesian coordinates
    real(wp),intent(in) :: xyz(:,:)

    !> Van-der-Waals radii including probe radius of solvent
    real(wp),intent(in) :: vdwsa(:)

    !> Radial weights for the numerical integration
    real(wp),intent(in) :: wrp(:)

    !> Radial smoothing function
    real(wp),intent(in) :: trj2(:,:)

    !> Angular weights for the numerical integration
    real(wp),intent(in) :: angWeight(:)

    !> Angular grid for each atom
    real(wp),intent(in) :: angGrid(:,:)

    real(wp),intent(in) :: pressure

    !> Surface area
    real(wp),intent(out) :: area

    real(wp),intent(out) :: volume

    real(wp),intent(out) :: energy

    !> Derivative of surface area w.r.t. cartesian coordinates
    real(wp),intent(out) :: grad(:,:)

    integer :: iat,ip,nno
    real(wp) :: rsas,sasai,xyza(3),xyzp(3),sasap,wr,wsa,voli,rdotn,trcorr(3)
    real(wp),allocatable :: gradi(:,:)

    area = 0.0_wp
    volume = 0.0_wp
    allocate (gradi(3,nat),source=0.0_wp)

    !$omp parallel do default(none) reduction(+:area, volume, grad) &
    !$omp shared(nat, vdwsa, nnsas, xyz, wrp, angGrid, angWeight, nnlists, trj2) &
    !$omp private(iat, rsas, nno, sasai, xyza, wr, ip, xyzp, wsa, voli, gradi, sasap, rdotn)
    do iat = 1,nat

      rsas = vdwsa(iat)
      nno = nnsas(iat)

      !> initialize storage
      gradi = 0.0_wp
      sasai = 0.0_wp
      voli = 0.0_wp

      !> atomic position
      xyza(:) = xyz(:,iat)

      !> radial atomic weight
      wr = wrp(iat)

      !> loop over grid points
      do ip = 1,size(angGrid,2)
        !> grid point position
        xyzp(:) = xyza(:)+rsas*angGrid(1:3,ip)

        !> atomic surface function at the grid point
        call compute_w_sp(nat,nnlists(:nno,iat),trj2,vdwsa,xyz,nno,xyzp,sasap)
        if (sasap .gt. tolsesp) then

          !> numerical quadrature weight
          wsa = angWeight(ip)*wr*sasap

          !> accumulate the surface area, sum in eq 10
          sasai = sasai+wsa

          !> calculate and accumulate the volume (as for PV term) fraction, since it is basically free
          rdotn = xyzp(1)*angGrid(1,ip)+xyzp(2)*angGrid(2,ip)+xyzp(3)*angGrid(3,ip)
          voli = voli+rdotn*wsa

          !> eq 3, calculation of xhcff grad
          gradi(:,iat) = gradi(:,iat)+(angGrid(:,ip)*wsa)
        end if
      end do

      !!$omp critical
      !> finalize calculation here to save multipilications
      area = area+sasai
      volume = volume+voli/3.0_wp
      grad = grad+gradi
      !!$omp end critical
    end do
    !$omp end parallel do

    !> xhcff force field is not conservative and thus has no energy contribution
    energy = volume*pressure
    grad = grad*pressure

    !> correct for translations and rotations
    trcorr(1:3) = sum(grad,dim=2)/nat
    do iat = 1,nat
      grad(1:3,iat) = grad(1:3,iat)-trcorr(1:3)
    end do

  end subroutine xhcff_eg

!========================================================================================!

  !> evaluate switching function without gradient
  pure subroutine compute_w_sp(nat,nnlists,trj2,vdwsa,xyza,nno,xyzp,sasap)
    implicit none

    integer,intent(in)  :: nat
    integer,intent(in)  :: nnlists(nno)
    integer,intent(in)  :: nno
    real(wp),intent(in)  :: xyza(3,nat)
    real(wp),intent(in)  :: xyzp(3)
    real(wp),intent(out) :: sasap
    real(wp),intent(in)  :: trj2(2,nat)
    real(wp),intent(in)  :: vdwsa(nat)

    integer  :: i,ia
    real(wp) :: tj(3),tj2,sqtj
    real(wp) :: uj,ah3uj2
    real(wp) :: sasaij

    !> initialize storage
    sasap = 1.0_wp
    do i = 1,nno
      ia = nnlists(i)
      !> compute the distance to the atom
      tj(:) = xyzp(:)-xyza(:,ia)
      tj2 = tj(1)*tj(1)+tj(2)*tj(2)+tj(3)*tj(3)
      !> if within the outer cut-off compute
      if (tj2 .lt. trj2(2,ia)) then
        if (tj2 .le. trj2(1,ia)) then
          sasap = 0.0_wp
          return
        else
          sqtj = sqrt(tj2)
          !> r in eq. 9
          uj = sqtj-vdwsa(ia)
          ah3uj2 = ah3*uj*uj
          !> eq 9, evaluation of atomic volume exclusion function
          sasaij = ah0+(ah1+ah3uj2)*uj

          !> accumulate the molecular surface, product in eq. 10
          sasap = sasap*sasaij
        end if
      end if
    end do

  end subroutine compute_w_sp

!========================================================================================!
!========================================================================================!
end module xhcff_engrad

