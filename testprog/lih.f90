module lih
!> Lithiumhydride
  use iso_fortran_env,only:wp => real64

!&<
  integer,parameter :: testnat = 2
  !> Atomtypes
  integer, parameter :: testat(testnat) = [3, 1]

  !> geometry in bohr
    real(wp),parameter :: testxyz(3,testnat) =    reshape(&
     [-6.69361629282015_wp,  1.98625390237610_wp,  0.122449150639047_wp, & 
    & -4.70384466526237_wp,  1.93397236145680_wp,  0.721668265764252_wp], &
    &  shape(testxyz))
  !> testpressure in GPa
  real(wp), parameter :: testpressure = 10.0_wp

  !> proberad in GPa
  real(wp), parameter :: testproberad = 0.0_wp

  !> internally generated Gradient as consistency check
  !> TODO: add 
!  real(wp),parameter :: testGradD3(3,testnat) =    reshape(&
!      &  shape(testGradD3))

  !> generated with Qchem and unscaled Radii as reference
!  real(wp),parameter :: testGradBondi(3,testnat) = reshape(&
!  &  shape(testGradBondi))
!&>

  ! Volume using unscaled Bondi Radii in Bohr ** 3
  real(wp), parameter :: testvol_bondi = 173.14791628632412_wp
  real(wp), parameter :: testpv_bondi = 0.0588517815670768_wp
  !public :: testnat
  !public :: testat
  !public :: testxyz
  !public :: testpressure
  !public :: testproberad
  !public :: testgrad
  !public :: writetestcoord

  character(len=2),parameter :: testelem(8) = ['h ','he','li','be', &
      &                        'b ','c ','n ','o ']
  private :: testelem

contains

  subroutine writetestcoord

    integer :: ich,i

    open (newunit=ich,file='coord')
    write (ich,'(a)') '$coord'
    do i = 1,testnat
      write (ich,'(3f20.14,2x,a2)') testxyz(1:3,i),testelem(testat(i))
    end do
    write (ich,'(a)') '$end'
    close (ich)

  end subroutine writetestcoord

end module lih
