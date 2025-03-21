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
    real radius, area
        
    !request radius value to the user
    write(*,*) 'Insert radius value:'

    !read user input
    read(*,*) radius

    !call calculate_area
    call calculate_area(radius, area)

    !print result
    write(*,*) 'Area of the circle is:', area

contains

    subroutine calculate_area(radius, area)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(in) :: radius
        real, intent(out) :: area

        !calculate area (4*atan(1) = pi)
        area = 4*atan(1.)*radius**2

    end subroutine calculate_area

end program circle_area