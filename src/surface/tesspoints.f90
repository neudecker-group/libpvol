! Created by felix on 6/12/23.


! try off a tesspoint class...
! TODO write correctly
module tesspoints
    use iso_fortran_env, only : wp => real64, stdout => output_unit
    implicit none
    private
    public :: tesspts, fpt, printout



    !> tesselation points of one atom
    type tesspts(n)
        !> number of tesselation points of this atom
    integer, len :: n = 1

        !> coordinates of tesselation points
        real(wp), dimension(n,3) :: xyz

        !> areas of corresponding tesselation points
        real(wp), dimension(n)  :: areas

    end type tesspts

    contains
    !> fill information for tesspoint in tesspt object
    subroutine fpt(this, ipt, coords, area)
    type(tesspts), intent(in) :: this
    integer, intent(in) :: ipt
    real(wp), intent(in) :: coords(3)
    real(wp), intent(in) :: area
    integer :: i = 0

    do i = 1, 3
        this%xyz(ipt, i) = coords(i)
    end do
    this%areas(ipt) = area
    end subroutine fpt

    subroutine printout(this)
        type(tesspts), intent(in) :: this
        integer :: i = 1
        print *,"x, y, z, area"
        do i = 1,this%n
            print '(f4.2)', this%xyz(i, 1), " ", this%xyz(i, 2), " ", this%xyz(i, 3), " ", this%areas(i)
        end do
        print '(f4.2)', "total area: ", sum(this%areas)
    endsubroutine

end module tesspoints