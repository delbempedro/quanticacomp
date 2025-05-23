!------------------------------------------------------------
! File: L3-5255417-ex-2.f90
!
! Description:
!   Solve second order differential equation using Runge-Kutta method
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
program Runge_Kutta

    !deactivate implicit typing
    implicit none

    !declare parameters
    integer, parameter :: n = 100
    real, parameter :: lambda = -4.0*(2.0*atan(1.0))**2.0

    !declare variables
    integer :: i
    real y(0:n), v(0:n), k1, k2, k3, k4, l1, l2, l3, l4, h

    !define h
    write(*,*) "Insert h:"
    read(*,*) h

    !initialize y and v
    y(0) = 1.0
    v(0) = 0.0

    !open file for writing results
    open(1, file='results.txt', status='replace')
    
    !compute Runge-Kutta method
    do i = 0, n-1

        !print current step
        write(1,*) i, y(i), v(i)

        !compute coefficients
        k1 = v(i)
        l1 = h*y(i)
        k2 = v(i) + h*l1/2.0
        l2 = lambda*(y(i) + h*k1/2.0)
        k3 = v(i) + h*l2/2.0
        l3 = lambda*(y(i) + h*k2/2.0)
        k4 = v(i) + h*l3
        l4 = lambda*(y(i) + h*k3)

        !update y and v
        y(i+1) = y(i) + h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        v(i+1) = v(i) + h*(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0

    end do

    !print last step
    write(1,*) n, y(n), v(n)

    !close file
    close(1)

end program Runge_Kutta