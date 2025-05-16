!------------------------------------------------------------
! File: L1-5255417-ex-3.f90
!
! Description:
!   Finds computer precision
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
program computer_precision

    !deactivate implicit typing
    implicit none

    !define variables
    real*4 real4, sum4, aux4
    real*8 real8, sum8, aux8

    !initialize variables
    real4 = 1.0
    sum4 = 1.0
    aux4 = 1.0
    real8 = 1.0
    sum8 = 1.0
    aux8 = 1.0

    !find machine precision
    do while (sum4+real4 /= sum4) !while sum4+real4 is not equal to sum4

        !define aux4 to save the value before exceeding machine precision
        real4 = aux4

        !divide real4 by 2 and save the value in aux4
        aux4 = real4*0.5

    end do
    
    !find machine precision
    do while (sum8+real8 /= sum8) !while sum4+real4 is not equal to sum4
        
        !define aux8 to save the value before exceeding machine precision
        real8 = aux8

        !divide real4 by 2 and save the value in aux8
        aux8 = real8*0.5d0
        
    end do

    !print results
    write(*,*) "Machine precision for simple precision is:", real4
    write(*,*) "Machine precision for double precision is:", real8

end program computer_precision