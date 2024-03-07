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

program xhcfflib_pv_tester
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use lih
  use pv_interface
  use xhcff_surface_module
  use pv_engrad
  use xhcff_type_timer
  implicit none

  integer :: nat
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  real(wp),allocatable :: stencil(:,:) !> stencil geom
  integer :: chrg
  integer :: uhf
  integer :: i,j,k,l
  character(len=40) :: atmp
!========================================================================================!
  real(wp) :: p, probe
  real(wp) :: fw, bw !> stencils
  real(wp) :: energy
  real(wp),allocatable :: gradient(:,:)
  real(wp),allocatable :: numGrad(:,:)
  real(wp),allocatable :: gradDiff(:,:)
  real(wp) :: gnorm

  logical :: fail,pr
  integer :: io
  type(pv_calculator) :: pv
  type(surface_calculator) :: surf

  integer :: ntimes
  type(xhcff_timer) :: timer

  real(wp) :: eDiff

  real(wp),parameter :: tolDiffD3 = 1e-6 !> tolerable difference for gradient units
  real(wp),parameter :: tolDiffBondiEnergy = 5e-4 !> tolerable difference for PV term

  integer :: threads
!========================================================================================!
  fail = .false.
  pr = .true.

  nat = testnat
  allocate (at(nat),xyz(3,nat),stencil(3,nat))
  at = testat
  xyz = testxyz
  p = testpressure
  probe = testproberad
  chrg = 0

  energy = 0.0_wp
  gnorm = 0.0_wp
  allocate (gradient(3,nat),source=0.0_wp)
  allocate (graddiff(3,nat),source=0.0_wp)
  allocate (numGrad(3,nat), source=0.0_wp)
  write (*,*) nat
  write (*,*)
  do i = 1,nat
    write (*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
  end do
  call writetestcoord()

!=======================================================================================!
!> parse optional command line args and their defaults
  threads = 0
  call ParseCommandLineArgs(threads)

!=======================================================================================!
!=======================================================================================!
!> Surface calculation USAGE
!=======================================================================================!
!=======================================================================================!

  write (*,*)
  write (*,*) '==========================BEGIN================================='
  write (*,*) '===================== SASA CALCULATION ========================='
  write (*,*) '==========================BEGIN================================='

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
  write (*,*) '==========================BEGIN================================='
  write (*,*) '====================  PV SINGLEPOINT  =========================='
  write (*,*) '==========================BEGIN================================='
  write (*,*)

  !=======================!
  !>  PV with D3 radii
  !=======================!
 ! call pv%init(nat,at,xyz,testpressure,proberad=testproberad,verbose=.false.,printlevel=2)
 ! call pv%singlepoint(nat,at,xyz,energy,gradient)
 ! call pv%info()

  !> test difference to reference volume
  

  !==========================!
  !>  PV with Bondi radii
  !==========================!
  !call pv%reset
  call pv%init(nat,at,xyz,testpressure,gridpts=5294, &
  &    proberad=testproberad,vdwSet=1,verbose=.false., printlevel=2)
  call pv%singlepoint(nat,at,xyz,energy,gradient)
  call pv%info()

  !> test energy difference
  eDiff = energy-testpv_bondi
  fail = abs(eDiff) > tolDiffBondiEnergy
  if (fail) then 
    write (*,*) '**** UNITTEST FAILED while Comparing Energies! *****'
    write (*,'(3A16)') 'Energy', 'Ref', 'Diff'
    write (*,'(3f16.12)') energy, testpv_bondi, eDiff
    stop
  end if

  stencil = xyz
  do i = 1,nat
    do j = 1,3
      stencil(j,i) = stencil(j,i) + 0.001
      call pv%singlepoint(nat,at,stencil,fw,gradient)
      stencil(j,i) = stencil(j,i) -0.002
      call pv%singlepoint(nat,at,stencil,bw,gradient)
      numGrad(j,i) = (fw - bw)/0.002
      stencil(j,i) = xyz(j,i)
    end do
  end do

  do i=1,nat
    write(*,'(3f16.12)') numGrad(1:3,i)
  end do

  write (*,*)
  write (*,*) '========================== END ================================='
  write (*,*) '=====================  PV SINGLEPOINT  ========================='
  write (*,*) '========================== END ================================='
!=======================================================================================!

  if (threads > 0) then
    write (*,*)
    write (*,*) '==========================BEGIN================================='
    write (*,*) '=====================  PV OpenMP TEST  ========================='
    write (*,*) '==========================BEGIN================================='

    ntimes = threads
    call timer%new(ntimes,.true.)
    do i = 1,ntimes
      call OMP_Set_Num_Threads(i)
#ifdef WITH_MKL
      call MKL_Set_Num_Threads(i)
#endif
      call ompprint_intern(atmp)

      call timer%measure(i,atmp)
      call pv%reset
      call pv%init(nat,at,xyz,testpressure,gridpts=5294, &
     &    proberad=testproberad,vdwSet=1,verbose=.false.)
      do j = 1,250 !> a few repetitions to actually see some CPU time...
        call pv%singlepoint(nat,at,xyz,energy,gradient)
      end do
      call timer%measure(i)
      write (*,*)
      call timer%write_timing(stdout,i)
    end do

    write (*,*)
    write (*,*) '========================== END ================================='
    write (*,*) '==================== XHCFF OpenMP TEST ========================='
    write (*,*) '========================== END ================================='
  end if
!=======================================================================================!
  deallocate (gradient,gradDiff,numGrad)
  deallocate (xyz,at,stencil)
!=======================================================================================!
end program xhcfflib_pv_tester

!=======================================================================================!
subroutine ompprint_intern(str)
 use omp_lib
  implicit none
  integer :: nproc,TID
  character(len=*) :: str
!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  IF (TID .EQ. 0) THEN
    nproc = OMP_GET_NUM_THREADS()
#ifdef WITH_MKL
    write (str,'(a,i3)') 'OMP/MKL threads = ',nproc
#else
    write (str,'(a,i3)') 'OMP threads = ',nproc
#endif
  END IF
!$OMP END PARALLEL
end subroutine ompprint_intern

!=======================================================================================!
subroutine ParseCommandLineArgs(threads)
  implicit none
  character(len=256) :: arg,arg2 ! Buffer to hold each argument
  integer :: numArgs,i     ! Variables to store argument count and loop index
  integer :: io,dumi
  !> IN/OUTPUTS
  integer,intent(inout) :: threads

  ! Get the number of command-line arguments
  numArgs = COMMAND_ARGUMENT_COUNT()

  ! Check if there are any arguments passed
  if (numArgs == 0) then
    return
  end if

  ! Loop through the command-line arguments
  do i = 1,numArgs
    ! Fetch each argument and store it in 'arg'
    call GET_COMMAND_ARGUMENT(i,arg)

    select case (trim(arg))
    case ('-T')
      call GET_COMMAND_ARGUMENT(i+1,arg2)
      read (arg2,*,iostat=io) dumi
      if (io == 0) threads = dumi
    end select
  end do

end subroutine ParseCommandLineArgs
