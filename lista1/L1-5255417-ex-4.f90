!------------------------------------------------------------
! File: L1-5255417-ex-4.f90
!
! Description:
!   Computes taylor series for exponential of -x
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
    integer i
    real*16 x(5), e_x, index, old_term, next_term, under_precision

    !initialize variables
    x = [0.1,1.0,10.0,100.0,1000.0]
    under_precision = 10**8

    do i=1,5

        !initialize variables
        e_x = 1.0
        index = 1.0
        next_term = 1.0
        old_term = 1.0

        !call compute_exponential
        call compute_exponential(x(i), e_x, index, old_term, next_term, under_precision)

        !print result
        write(*,*) "x = ", x(i)
        write(*,*) "    Computed e^(-", x(i), ") = ", e_x
        write(*,*) "    Real e^(-", x(i), ") = ", exp(-x(i))

    end do

contains

    subroutine compute_exponential(x, e_x, index, old_term, next_term, under_precision)

        !deactivate implicit typing
        implicit none

        !define variables
        real*16, intent(in) :: x
        real*16, intent(in) :: under_precision
        real*16, intent(inout) :: index
        real*16, intent(inout) :: next_term
        real*16, intent(inout) :: old_term
        real*16, intent(out) :: e_x

        !compute next term
        old_term = next_term
        next_term = next_term*(-x)/index

        !update e_x
        e_x = e_x + next_term

        !update index and factorial
        index = index + 1.0
        
        !update e_x while the absolute value of next_term is greater than precision
        do while (under_precision*abs(next_term-old_term) > abs(e_x))

            !compute next term
            old_term = next_term
            next_term = next_term*(-x)/index

            !update e_x
            e_x = e_x + next_term

            !update index and factorial
            index = index + 1.0

        end do

        write(*,*) next_term, old_term, e_x

    end subroutine compute_exponential

end program exponential_taylor_series