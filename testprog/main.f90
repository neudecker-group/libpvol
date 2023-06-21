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
    use iso_fortran_env, only: wp=>real64,stdout=>output_unit
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

   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:)
   real(wp) :: gnorm

   logical :: fail,pr
   integer :: io
   type(xhcff_data) :: dat
   type(surface_calculator) :: surf

!========================================================================================!
    fail = .false.
    pr = .true.

    nat = testnat
    allocate(at(nat),xyz(3,nat))
    at = testat
    xyz = testxyz
    chrg = 0

    energy = 0.0_wp
    gnorm  = 0.0_wp
    allocate(gradient(3,nat),source=0.0_wp)

    write(*,*) nat
    write(*,*)
    do i=1,nat
      write(*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
    enddo
    call writetestcoord()

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

    write(*,*)
    write(*,*) '================================================================'
    write(*,*) '===================== SASA CALCULATION ========================='
    write(*,*) '================================================================'


    !> default probe size (1.5 Angstroem)
!    write(*,*)
!    call surf%setup(nat,at,xyz,.true.)

!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%tight)

!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%vtight)

!    write(*,*)
!   call surf%deallocate()
!   call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%extreme)



    !> smaller probe
!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,probe=1.0_wp)

!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%tight, probe=1.0_wp)

!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%vtight, probe=1.0_wp)

!    write(*,*)
!    call surf%deallocate()
!    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%extreme, probe=1.0_wp)

    !> small test
    call surf%deallocate()
    call surf%setup(nat,at,xyz,.true.,ngrid=lebedev%extreme, probe=0.0_wp)

    write(*,*)
    write(*,*) '========================== END ================================='
    write(*,*) '==================== SASA CALCULATION =========================='
    write(*,*) '========================== END ================================='


!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

    write(*,*)
    write(*,*) '================================================================'
    write(*,*) '==================== XHCFF SINGLEPOINT ========================='
    write(*,*) '================================================================'


    call xhcff_eg(nat,at,xyz,0.0_wp,0.00034_wp, energy, gradient)
    print *, 'xhcff gradient:'
    do i=1,nat
        write (*,'(2x,i3,3x,3f16.6)') , i, gradient(1:3,i)
    end do

    write(*,*)
    write(*,*) 'gnorm:',sqrt(sum(gradient**2))

    write(*,*)
    write(*,*) '========================== END ================================='
    write(*,*) '==================== XHCFF SINGLEPOINT ========================='
    write(*,*) '========================== END ================================='


!=======================================================================================!
   deallocate(gradient)
   deallocate(xyz,at)
!=======================================================================================!
end program xhcfflib_main_tester
