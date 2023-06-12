! Created by felix on 6/12/23.


! try off a tesspoint class...
! TODO write correctly
module tesspoints
    use iso_fortran_env, only : wp => real64, stdout => output_unit
    implicit none
    private
    public :: tesspts !, fpt, printout



    !> tesselation points of one atom
    type tesspts
        !> number of tesselation points of this atom
    integer :: n

        !> coordinates of tesselation points
        real(wp), allocatable :: xyz

        !> areas of corresponding tesselation points
        real(wp), allocatable :: areas

!        contains
 !       procedure :: fill => fillpt
 !       procedure :: print => printout
    end type tesspts

!    contains
    !> fill information for tesspoint in tesspt object
 !   subroutine fillpt(this, ipt, coords, area)
 !   class(tesspts), intent(inout) :: this
 !   integer, intent(in) :: ipt
 !   real(wp), intent(in) :: coords(3)
 !   real(wp), intent(in) :: area
 !   integer :: i = 0

  !  do i = 1, 3
  !      this%xyz(ipt, i) = coords(i)
  !  end do
  !  this%areas(ipt) = area
  !  end subroutine fillpt

  !  subroutine printout(this)
  !      class(tesspts), intent(inout) :: this
  !      integer :: i = 1
  !      print *,"x, y, z, area"
  !      do i = 1,this%n
  !          print '(f4.2)', this%xyz(i, 1), " ", this%xyz(i, 2), " ", this%xyz(i, 3), " ", this%areas(i)
  !      end do
  !      print '(f4.2)', "total area: ", sum(this%areas)
  !  endsubroutine

end module tesspoints