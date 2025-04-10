!------------------------------------------------------------
! File: L3-5255417-ex-3.f90
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
program Numerov

    !deactivate implicit typing
    implicit none

    !declare parameters
    integer, parameter :: n = 100
    real, parameter :: lambda = 4.0*(2.0*atan(1.0))**2.0 !y''(x) + lambda*y(x) = 0

    !declare variables
    integer :: i
    real y(0:n), h

    !define h
    write(*,*) "Insert h:"
    read(*,*) h

    !initialize y and v
    y(0) = 1.0
    y(1) = 1.0 + lambda*(h**2.0)/2.0

    !open file for writing results
    open(1, file='results.txt', status='replace')
    
    !compute Numerov method
    do i = 1, n-1

        !print current step
        write(1,*) i, y(i)

        !update y
        y(i+1) = ( 2.0*y(i)*(1.0 - 5.0*h**2*lambda/12.0) - y(i-1)*(1.0 + h**2*lambda/12.0) ) / (1.0 + h**2*lambda/12.0)
        
    end do

    !print last step
    write(1,*) n, y(n)

    !close file
    close(1)

end program Numerov