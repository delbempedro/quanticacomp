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

    !deactivate implicit typing
    implicit none

    !define variables
    real test_overflow, real_overflow, test_underflow, real_underflow
    real*8 test_overflow8, real_overflow8, test_underflow8, real_underflow8

    !initialize overflow variables
    test_overflow = 1.0
    test_overflow8 = 1.0

    !call overflows computations
    call compute_overflow(test_overflow,real_overflow)
    call compute_overflow8(test_overflow8,real_overflow8)

    !initialize underflow variables
    test_underflow = 1.0
    test_underflow8 = 1.0

    !call underflows computations
    call compute_underflow(test_underflow,real_underflow)
    call compute_underflow8(test_underflow8,real_underflow8)
    
    !print results
    write(*,*)'Overflow simple real is: ', real_overflow
    write(*,*)'Overflow double real is: ', real_overflow8
    write(*,*)'Underflow simple real is: ', real_underflow
    write(*,*)'Underflow double real is: ', real_underflow8

contains

    !computes single precision overflow
    subroutine compute_overflow(test_overflow,real_overflow)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(inout) :: test_overflow
        real, intent(inout) :: real_overflow
        
        !compute overflow
        do while (test_overflow < 2.0*test_overflow) !test_overflow => 2*test_overflow implies overflow

            !update real_overflow variable which saves the value before exceeding the overflow
            real_overflow = test_overflow

            !update test_overflow variable to test the next value
            test_overflow = test_overflow*10

        end do

    end subroutine compute_overflow

    !computes double precision overflow
    subroutine compute_overflow8(test_overflow8,real_overflow8)

        !deactivate implicit typing
        implicit none

        !define variables
        real*8, intent(inout) :: test_overflow8
        real*8, intent(inout) :: real_overflow8
        
        !compute overflow
        do while (test_overflow8 < 2.0*test_overflow8) !test_overflow8 => 2*test_overflow8 implies overflow

            !update real_overflow8 variable which saves the value before exceeding the overflow
            real_overflow8 = test_overflow8

            !update test_overflow8 variable to test the next value
            test_overflow8 = test_overflow8*10

        end do

    end subroutine compute_overflow8

    !computes single precision underflow
    subroutine compute_underflow(test_underflow,real_underflow)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(inout) :: test_underflow
        real, intent(inout) :: real_underflow
        
        !compute underflow
        do while (test_underflow > 0.5*test_underflow) !test_underflow > 0.5*test_underflow implies underflow

            !update real_underflow variable which saves the value before exceeding the underflow
            real_underflow = test_underflow

            !update test_underflow variable to test the next value
            test_underflow = test_underflow/10

        end do

    end subroutine compute_underflow

    !computes double precision underflow
    subroutine compute_underflow8(test_underflow8,real_underflow8)

        !deactivate implicit typing
        implicit none

        !define variables
        real*8, intent(inout) :: test_underflow8
        real*8, intent(inout) :: real_underflow8
        
        !compute underflow
        do while (test_underflow8 > 0.5*test_underflow) !test_underflow8 > 0.5*test_underflow8 implies underflow

            !update real_underflow8 variable which saves the value before exceeding the underflow
            real_underflow8 = test_underflow8

            !update test_underflow8 variable to test the next value
            test_underflow8 = test_underflow8/10
            
        end do

    end subroutine compute_underflow8

end program testing_overflow_and_underflow