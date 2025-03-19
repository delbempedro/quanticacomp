!------------------------------------------------------------
! File: L1-5255417-ex-5.f90
!
! Description:
!   Computes sums of series
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
    integer N(14), i
    real sum1, sum2, sum3

    !initialize variables
    N = [1,2,3,4,5,6,7,9,10**1,10**2,10**3,10**4,10**5,10**6]

    do i=1,14

        !initialize variables
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0

        !compute series
        call compute_series(sum1, sum2, sum3, N(i))

        !print result
        write(*,*) N(i), sum1, sum2, sum3

    end do

contains

    subroutine compute_series(sum1, sum2, sum3, N)

        !deactivate implicit typing
        implicit none

        !define variables
        integer i
        integer, intent(in) :: N
        real x, sum2a, sum2b
        real, intent(inout) :: sum1
        real, intent(inout) :: sum2
        real, intent(inout) :: sum3
        
        !initialize partial sums
        sum2a = 0.0
        sum2b = 0.0

        !compute sum1
        do i=1,N

            !update x
            x = i

            !compute sums
            sum1 = sum1 + (-1)**(x)*x/(x+1)
            sum2a = sum2a + (2.0*x-1)/(2.0*x)
            sum2b = sum2b + (2.0*x)/(2.0*x+1)
            sum3 = sum3 + 1/(2.0*x*(2.0*x+1))

        end do

        !!! first serie is actually 1 to 2N?  !!!

        !!update sum1
        !do i=N,2*N

            !!update x
            !x = i

            !!compute sums
            !sum1 = sum1 + (-1)**(x)*x/(x+1)

        !end do

        !update sum2
        sum2 = -sum2a + sum2b

    end subroutine compute_series

end program exponential_taylor_series