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
    real*16 x(5), e_x, index, next_term, precision

    !initialize variables
    x = [0.1,1.0,10.0,100.0,1000.0]
    precision = 1.0e-8

    do i=1,5

        !initialize variables
        e_x = 1.0
        index = 1.0
        next_term = 1.0

        !call compute_exponential
        call compute_exponential(x(i), e_x, index, next_term, precision)

        !print result
        write(*,*) "    x = ", x(i)
        write(*,*) "    Computed e^(-", x(i), ") = ", e_x
        write(*,*) "    Real e^(-", x(i), ") = ", exp(-x(i))

    end do

contains

    subroutine compute_exponential(x, e_x, index, next_term, precision)

        !deactivate implicit typing
        implicit none

        !define variables
        real*16, intent(in) :: x
        real*16, intent(in) :: precision
        real*16, intent(inout) :: index
        real*16, intent(inout) :: next_term
        real*16, intent(out) :: e_x
        
        !update e_x while the absolute value of next_term is greater than precision
        do while (abs(next_term) > precision)

            !compute next term
            next_term = next_term*(-x)/index

            !update e_x
            e_x = e_x + next_term

            !update index and factorial
            index = index + 1

        end do

    end subroutine compute_exponential

end program exponential_taylor_series