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

    implicit none

    real radio, area
        
    !Request radius value to the user
    write(*,*) 'Insert radius value:'

    !Read user input
    read(*,*) radio

    !Call calculate_area
    call calculate_area(radio, area)

    !Print result
    write(*,*) 'Area of the circle is:', area

contains

    !Subroutine to calculate area
    subroutine calculate_area(radio, area)

        implicit none

        real, intent(in) :: radio !Input (raio)
        real, intent(out) :: area !Output (área)

        !Calculate area
        area = 4*atan(1.)*radio**2

    end subroutine calculate_area

end program circle_area