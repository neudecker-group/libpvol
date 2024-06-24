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
 ! use tesspoints,only:tesspts
 ! use xhcff_surface_module
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
  public :: pv_egtest

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  !> pure gradient calculation

  subroutine pv_egtest(nat, nnsas, nnlists, xyz, vdwsa, &
        & wrp,trj2,angWeight,angGrid,pressure,area,volume,energy,grad)
  !$  use omp_lib
    
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

    real(wp), intent(out) :: energy

    !> Derivative of surface area w.r.t. cartesian coordinates
    real(wp),intent(out) :: grad(:,:)

    ! make output
    !type(tesspts),intent(out) :: tess(:)

    !real(wp) :: tj(3),tj2

    integer :: iat,ip,jj,nnj,nni,nno
    real(wp) :: rsas,sasai,xyza(3),xyzp(3),sasap,wr,wsa,drjj(3), voli, rdotn
    real(wp),allocatable :: grds(:,:),grads(:,:),xyzt(:,:),areas(:), gradi(:,:)
    integer,allocatable :: grdi(:)

    area = 0.0_wp
    volume =0.0_wp
    !> allocate space for the gradient storage
    allocate (grads(3,nat),source=0.0_wp)
    allocate (gradi(3,nat),source=0.0_wp)
    allocate (grds(3,maxval(nnsas)))
    allocate (grdi(maxval(nnsas)))
    allocate (xyzt(3,size(angGrid,2)))
    allocate (areas(size(angGrid,2)))
    !do iat = 1,nat
    !  call tess(iat)%allocate(size(angGrid,2),nat)
    !end do
    !$omp parallel do default(none) shared(area, volume, grad) &
    !$omp shared(nat, vdwsa, nnsas, xyz, wrp, angGrid, angWeight, nnlists, trj2) &
    !$omp private(iat, rsas, nno, grads, sasai, xyza, wr, ip, xyzp, wsa, xyzt, areas, &
    !$omp& voli, gradi, sasap, jj, nni, nnj, grdi, grds, drjj, rdotn)
    do iat = 1,nat

      rsas = vdwsa(iat)
      nno = nnsas(iat)

      !> initialize storage
      gradi = 0.0_wp
      sasai = 0.0_wp
      voli  = 0.0_wp

      !> reset areas and xyzt
      !areas = 0.0_wp
      !xyzt(:,:) = 0.0_wp

      !> atomic position
      xyza(:) = xyz(:,iat)
      !> radial atomic weight
      wr = wrp(iat)

      !> loop over grid points
      do ip = 1,size(angGrid,2)
        !> grid point position
        xyzp(:) = xyza(:)+rsas*angGrid(1:3,ip)
        
        !> reset area gradient storage for point
        grads = 0.0_wp
        !> save gridpoint position
        !xyzt(1:3,ip) = xyzp
        !> atomic surface function at the grid point
        !> compute the distance to the atom

        call compute_w_sp(nat,nnlists(:nno,iat),trj2,vdwsa,xyz,nno,xyzp, &
           & sasap,grds,nni,grdi)
        if (sasap .gt. tolsesp) then
          !> numerical quadrature weight
          wsa = angWeight(ip)*wr*sasap

          !> accumulate the surface area, sum in eq 10
          sasai = sasai+wsa

          !> calculate and accumulate the volume fraction
          rdotn = xyzp(1) * angGrid(1,ip) + xyzp(2) * angGrid(2,ip) + xyzp(3) * angGrid(3,ip)
          voli = voli + rdotn * wsa

          !> save area tesspoint
          !areas(ip) = wsa*4.0_wp*pi

          !> accumulate the surface gradient
          do jj = 1,nni
            nnj = grdi(jj)
            drjj(:) = wsa*grds(:,jj)
            !tess(iat)%dadr(:,iat,ip) = tess(iat)%dadr(:,iat,ip) + drjj
            !tess(iat)%dadr(:,nnj,ip) = tess(iat)%dadr(:,nnj,ip) - drjj
            grads(:,iat) = grads(:,iat)+drjj(:)
            grads(:,nnj) = grads(:,nnj)-drjj(:)
          end do
          gradi(:,iat) = gradi(:,iat) + (angGrid(:,ip) * wsa)
          gradi = gradi + (rdotn * grads)
        end if
      end do
      !> finalize eq 10
      area = area + sasai * 4.0_wp * pi
      volume = volume + voli * 4.0_wp * pi / 3.0_wp
      grad = grad + gradi * 4.0_wp * pi / 3.0_wp
      !dsdrt(:,:,iat) = grads
      !ess(iat)%n = size(angGrid,2)
      !tess(iat)%ap = areas
      !tess(iat)%xyz = xyzt
    end do
    !$omp end parallel do
    energy = volume * pressure
    grad = grad * pressure

  end subroutine pv_egtest

!========================================================================================!
  pure subroutine compute_w_sp(nat,nnlists,trj2,vdwsa,xyza,nno,xyzp,sasap,grds, &
        &                      nni,grdi)
    implicit none

    integer,intent(in)  :: nat
    integer,intent(in)  :: nnlists(nno)
    integer,intent(in)  :: nno
    integer,intent(out) :: nni
    real(wp),intent(in)  :: xyza(3,nat)
    real(wp),intent(in)  :: xyzp(3)
    real(wp),intent(out) :: sasap
    real(wp),intent(out) :: grds(3,nno)
    integer,intent(out) :: grdi(nno)
    real(wp),intent(in)  :: trj2(2,nat)
    real(wp),intent(in)  :: vdwsa(nat)

    integer  :: i,ia
    real(wp) :: tj(3),tj2,sqtj
    real(wp) :: uj,uj3,ah3uj2
    real(wp) :: sasaij,dsasaij

    !> initialize storage
    nni = 0
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
          dsasaij = ah1+3.0_wp*ah3uj2
          !> eq 9, evaluation of atomic volume exclusion function
          sasaij = ah0+(ah1+ah3uj2)*uj

          !> accumulate the molecular surface, product in eq. 10
          sasap = sasap*sasaij
          !> compute the gradient wrt the neighbor
          dsasaij = dsasaij/(sasaij*sqtj)
          nni = nni+1
          grdi(nni) = ia
          grds(:,nni) = dsasaij*tj(:)
        end if
      end if
    end do

  end subroutine compute_w_sp

!========================================================================================!
!========================================================================================!
end module pv_engrad

