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
    real x_k, x_kplus1, f_x_k, df_x_k, f_x(5)
    integer iteration

    write(*,*) "Insert x^4 coefficient:"
    read(*,*) f_x(5)
    write(*,*) "Insert x^3 coefficient:"
    read(*,*) f_x(4)
    write(*,*) "Insert x^2 coefficient:"
    read(*,*) f_x(3)
    write(*,*) "Insert x coefficient:"
    read(*,*) f_x(2)
    write(*,*) "Insert constant coefficient:"
    read(*,*) f_x(1)

    !initialize variables
    x_k = 10.0
    x_kplus1 = 0.0
    f_x_k = f(x_k,f_x)
    df_x_k = df(x_k,f_x)
    iteration = 1

    open(unit=1, file='output.txt', action='write')

    f_x_k = f(x_k,f_x)
    df_x_k = df(x_k,f_x)
    x_kplus1 = x_k - f_x_k/df_x_k
    iteration = iteration + 1
    write(1,*) iteration, f(x_kplus1,f_x)

    !Newton-Raphson method
    do while (abs(x_kplus1 - x_k) > 1e-6)
        x_k = x_kplus1
        f_x_k = f(x_k,f_x)
        df_x_k = df(x_k,f_x)
        x_kplus1 = x_k - f_x_k/df_x_k
        iteration = iteration + 1
        write(1,*) iteration, f(x_kplus1,f_x)
    end do

    close(1)

contains

    function f(x,f_x) result(result)
        real, intent(in) :: x
        real, intent(in) :: f_x(5)
        real result
        
        result = f_x(5)*x**4. + f_x(4)*x**3. + f_x(3)*x**2. + f_x(2)*x + f_x(1)
    end function f

    function df(x,f_x) result(result)
        real, intent(in) :: x
        real, intent(in) :: f_x(5)
        real result
        
        result = 4.*f_x(5)*x**3. + 3.*f_x(4)*x**2. + 2.*f_x(3)*x + f_x(2)
    end function df

end program newton_raphson