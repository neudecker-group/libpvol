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
! along with xhcfflib. If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

program xhcfflib_main_tester
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use testmol
  use xhcff_interface
  use xhcff_surface_module
  use xhcff_engrad
  implicit none

  integer :: nat
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  integer :: chrg
  integer :: uhf
  integer :: i,j,k,l

!========================================================================================!
  real(wp) :: p,probe
  real(wp) :: energy
  real(wp),allocatable :: gradient(:,:)
  real(wp),allocatable :: gradDiff(:,:)
  real(wp) :: gnorm

  logical :: fail,pr
  integer :: io
  type(xhcff_calculator) :: xhcff
  type(surface_calculator) :: surf

  real(wp),parameter :: tolDiffD3 = 1e-6 !> tolerable difference for gradient units
  real(wp),parameter :: tolDiffBondi = 2e-4 !> tolerable difference for gradient units
!========================================================================================!
  fail = .false.
  pr = .true.

  nat = testnat
  allocate (at(nat),xyz(3,nat))
  at = testat
  xyz = testxyz
  p = testpressure
  probe = testproberad
  chrg = 0

  energy = 0.0_wp
  gnorm = 0.0_wp
  allocate (gradient(3,nat),source=0.0_wp)
  allocate (graddiff(3,nat),source=0.0_wp)
  write (*,*) nat
  write (*,*)
  do i = 1,nat
    write (*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
  end do
  call writetestcoord()

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  write (*,*)
  write (*,*) '================================================================'
  write (*,*) '===================== SASA CALCULATION ========================='
  write (*,*) '================================================================'

  !> small test
  call surf%deallocate()
  call surf%setup(nat,at,xyz,.true.,io,ngrid=lebedev%extreme,probe=0.0_wp)

  write (*,*)
  write (*,*) '========================== END ================================='
  write (*,*) '==================== SASA CALCULATION =========================='
  write (*,*) '========================== END ================================='

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  write (*,*)
  write (*,*) '================================================================'
  write (*,*) '==================== XHCFF SINGLEPOINT ========================='
  write (*,*) '================================================================'

  !> Xhcff with D3 radii
  call xhcff%init(nat,at,xyz,testpressure,proberad=testproberad,verbose=.false.)
  call xhcff%singlepoint(nat,at,xyz,energy,gradient)
  call xhcff%info()

  !> test difference to reference gradient
  gradDiff = gradient-testGradD3
  fail = ANY(abs(graddiff(:,:)) > tolDiffD3)

  if (fail) then
    write (*,*) 'UNITTEST FAILED for D3 radii!'
    write (*,*) 'difference between calculated gradient and reference:'

    do i = 1,nat
      write (*,'(2x,i3,3x,3f16.12)') i,graddiff(1:3,i)
    end do
  else
    write (*,*) 'Test passed!'
  end if

  !> Xhcff with Bondi radii
  call xhcff%reset
  call xhcff%init(nat,at,xyz,testpressure,gridpts=5294, &
  &    proberad=testproberad,vdwSet=1,verbose=.false.)
  call xhcff%singlepoint(nat,at,xyz,energy,gradient)
  call xhcff%info()

  !> test difference to reference gradient
  !> gradients are mirrored
  gradDiff = gradient+testGradBondi
  fail = any(abs(graddiff(:,:)) > tolDiffBondi)

  if (fail) then
    write (*,*) 'UNITTEST FAILED for Bondi radii!'
    write (*,*) 'difference between calculated gradient and reference:'

    do i = 1,nat
      write (*,'(2x,i3,3x,3f16.12)') i,graddiff(1:3,i)
    end do
  else
    write (*,*) 'Test passed!'
  end if

  write (*,*)
  write (*,*) '========================== END ================================='
  write (*,*) '==================== XHCFF SINGLEPOINT ========================='
  write (*,*) '========================== END ================================='

!=======================================================================================!
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program xhcfflib_main_tester
