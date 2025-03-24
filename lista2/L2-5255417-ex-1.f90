!------------------------------------------------------------
! File: L2-5255417-ex-1.f90
!
! Description:
!   Computes Newton-Raphson method convergence
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
program newton_raphson

    !deactivate implicit typing
    implicit none

    !define variables
    real x_k, x_kplus1, f_x, df_x
    integer iteration

    !initialize variables
    x_k = 10.0
    x_kplus1 = 0.0
    f_x = 0.0
    df_x = 0.0
    iteration = 1

    open(unit=1, file='output.txt', action='write')

    f_x = x_k**3 - 1
    df_x = 3*x_k**2
    x_kplus1 = x_k - f_x/df_x
    iteration = iteration + 1
    write(1,*) iteration, x_kplus1**3 - 1

    !Newton-Raphson method
    do while (abs(x_kplus1 - x_k) > 1e-6)
        x_k = x_kplus1
        f_x = x_k**3 - 1
        df_x = 3*x_k**2
        x_kplus1 = x_k - f_x/df_x
        iteration = iteration + 1
        write(1,*) iteration, x_kplus1**3 - 1
    end do

end program newton_raphson