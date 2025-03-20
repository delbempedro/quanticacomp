!------------------------------------------------------------
! File: L1-5255417-ex-1.f90
!
! Description:
!   Computes the area of a circle.
!
! Dependencies:
!   - None
!
! Since:
!   - 03/2025
!
! Authors:
!   - Pedro C. Delbem <pedrodelbem@usp.br>
!------------------------------------------------------------
program circle_area

    !deactivate implicit typing
    implicit none

    !define variables
    real radio, area
        
    !request radius value to the user
    write(*,*) 'Insert radius value:'

    !read user input
    read(*,*) radio

    !call calculate_area
    call calculate_area(radio, area)

    !print result
    write(*,*) 'Area of the circle is:', area

contains

    subroutine calculate_area(radio, area)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(in) :: radio
        real, intent(out) :: area

        !calculate area
        area = 4*atan(1.)*radio**2

    end subroutine calculate_area

end program circle_area