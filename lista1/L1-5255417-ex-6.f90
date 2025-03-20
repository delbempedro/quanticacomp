!------------------------------------------------------------
! File: L1-5255417-ex-6.f90
!
! Description:
!   Computes taylor series for exponential of -x with function of N
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
program exponential_taylor_series

    !deactivate implicit typing
    implicit none

    !define variables
    integer N(14), i, j
    real exponential_of_minus_x, exponential_of_x, x(5)

    !initialize variables
    N = [1,2,3,4,5,6,7,9,10**1,10**2,10**3,10**4,10**5,10**6]
    x = [0.1,1.0,10.0,100.0,1000.0]

    open(unit=1, file='exponential_taylor_series.txt', status='replace')

    do i=1,5

        !initialize variables
        exponential_of_minus_x = 0.0
        exponential_of_x = 0.0

        write(1,*) 'x:', x(i)

        do j=1,14
            !compute series
            call compute_series(exponential_of_minus_x, exponential_of_x, N(j), x(i))

            !compute exponential of -x
            exponential_of_x = 1/exponential_of_x

            !print result
            write(1,*) 'N:', N(j), exponential_of_minus_x, exponential_of_x, exp(-x(i))

        end do

    end do

    close(1)

contains

    subroutine compute_series(exponential_of_minus_x, exponential_of_x, N, x)

        !deactivate implicit typing
        implicit none

        !define variables
        integer i
        integer, intent(in) :: N
        real j, next_term_minus_x, next_term_x
        real, intent(in) :: x
        real, intent(inout) :: exponential_of_minus_x
        real, intent(inout) :: exponential_of_x

        !initialize next terms
        next_term_minus_x = 1.0
        next_term_x = 1.0
    
        !compute sum1
        do i=1,N

            !update j
            j = i

            !compute next terms
            next_term_minus_x = next_term_minus_x*(-x)/j
            next_term_x = next_term_x*x/j

            !compute sums
            exponential_of_minus_x = exponential_of_minus_x + next_term_minus_x
            exponential_of_x = exponential_of_x + next_term_x

        end do

    end subroutine compute_series

end program exponential_taylor_series