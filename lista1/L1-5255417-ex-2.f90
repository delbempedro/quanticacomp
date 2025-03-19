!------------------------------------------------------------
! File: L1-5255417-ex-2.f90
!
! Description:
!   Computes overflow and underflow
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
program testing_overflow_and_underflow

    implicit none

    real test_overflow, real_overflow, test_underflow, real_underflow

    real*8 test_overflow8, real_overflow8, test_underflow8, real_underflow8

    test_overflow = 1.0
    test_overflow8 = 1.0
    call compute_overflow(test_overflow,real_overflow)
    call compute_overflow8(test_overflow8,real_overflow8)

    test_underflow = 1.0
    test_underflow8 = 1.0
    call compute_underflow(test_underflow,real_underflow)
    call compute_underflow8(test_underflow8,real_underflow8)
    
    write(*,*)'Overflow simple real order is: ', real_overflow
    write(*,*)'Overflow double real order is: ', real_overflow8

    write(*,*)'Underflow simple real order is: ', real_underflow
    write(*,*)'Underflow double real order is: ', real_underflow8

contains

    subroutine compute_overflow(test_overflow,real_overflow)

        implicit none

        real, intent(inout) :: test_overflow
        real, intent(inout) :: real_overflow
        
        infinity = huge(1.0_real(4))
        do while (test_overflow < infinity)
            write(*,*) test_overflow
            test_overflow = test_overflow*10
        end do

    end subroutine compute_overflow

    subroutine compute_overflow8(test_overflow8,real_overflow8)

        implicit none

        real*8, intent(inout) :: test_overflow8
        real*8, intent(inout) :: real_overflow8
        real*8 :: infinity

        infinity = huge(1.0_real(8))
        do while (test_overflow8 < infinity)
            real_overflow8 = test_overflow8
            write(*,*) test_overflow8
            test_overflow8 = test_overflow8*10
        end do

    end subroutine compute_overflow8

    subroutine compute_underflow(test_underflow,real_underflow)

        implicit none

        real, intent(inout) :: test_underflow
        real, intent(inout) :: real_underflow
        
        do while (test_underflow > 3.0E-38)
            real_underflow = test_underflow
            write(*,*) test_underflow
            test_underflow = test_underflow/10
        end do

    end subroutine compute_underflow

    subroutine compute_underflow8(test_underflow8,real_underflow8)

        implicit none

        real*8, intent(inout) :: test_underflow8
        real*8, intent(inout) :: real_underflow8
        
        do while (test_underflow8 > 3.0E-38)
            real_underflow8 = test_underflow8
            write(*,*) test_underflow8
            test_underflow8 = test_underflow8/10
        end do

    end subroutine compute_underflow8

end program testing_overflow_and_underflow