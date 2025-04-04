!------------------------------------------------------------
! File: L2-5255417-ex-3.f90
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
program Runge-Kutta

    !deactivate implicit typing
    implicit none
    !declarations
    integer, parameter :: n = 1000
    real :: t0, t1, h
    real :: y0, y1
    real :: k1, k2, k3, k4
    real :: t(n), y(n)
    integer :: i
    !initial conditions
    t0 = 0.0
    t1 = 10.0
    y0 = 1.0
    y1 = 0.0
    !step size
    h = (t1 - t0) / n
    !initialize arrays
    t(1) = t0
    y(1) = y0
    !Runge-Kutta method
    do i = 1, n
        k1 = h * f(t(i), y(i))
        k2 = h * f(t(i) + h / 2, y(i) + k1 / 2)
        k3 = h * f(t(i) + h / 2, y(i) + k2 / 2)
        k4 = h * f(t(i) + h, y(i) + k3)
        y(i + 1) = y(i) + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t(i + 1) = t(i) + h
    end do
    !print results
    do i = 1, n
        print *, t(i), y(i)
    end do

end program Runge-Kutta