!------------------------------------------------------------
! File: L2-5255417-ex-2.f90
!
! Description:
!
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
program find_roots

    !deactivate implicit typing
    implicit none

    !define variables
    real x_kminus1, x_k, x_kplus1, f_x_kminus1, f_x_k, df_x_k, initial_guess, pre_initial_guess
    integer iteration

    !request initial guess
    write(*,*) "Insert initial guess:"
    read(*,*) initial_guess

    !open first output file
    open(unit=1, file='newtonraphson1.txt', action='write')

    !initialize variables
    x_k = initial_guess
    x_kplus1 = 0.0
    f_x_k = f1(x_k)
    df_x_k = df1(x_k)
    iteration = 1

    !print header
    write(1,*) 'Newton-Raphson Method'
    write(1,*) 'Initial guess: ', x_k
    write(1,*) 'f(x) = x^2 - 5'

    !print first iteration
    write(1,*) iteration, f1(x_kplus1)

    !update f_x and df_x
    f_x_k = f1(x_k)
    df_x_k = df1(x_k)

    !update x_kplus1
    x_kplus1 = x_k - f_x_k/df_x_k

    !update iteration
    iteration = iteration + 1

    !print second iteration
    write(1,*) iteration, f1(x_kplus1)

    !Newton-Raphson method
    do while (abs(x_kplus1 - x_k) > 1e-6)

        !update x_k
        x_k = x_kplus1

        !update f_x and df_x
        f_x_k = f1(x_k)
        df_x_k = df1(x_k)

        !update x_kplus1
        x_kplus1 = x_k - f_x_k/df_x_k

        !update iteration
        iteration = iteration + 1

        !print iteration
        write(1,*) iteration, f1(x_kplus1)
        
    end do

    !print root
    write(1,*) 'Root:', x_kplus1

    !close first output file
    close(1)

    !open first output file
    open(unit=2, file='newtonraphson2.txt', action='write')

    !initialize variables
    x_k = initial_guess
    x_kplus1 = 0.0
    f_x_k = f2(x_k)
    df_x_k = df2(x_k)
    iteration = 1

    !print header
    write(2,*) 'Newton-Raphson Method'
    write(2,*) 'Initial guess: ', x_k
    write(2,*) 'f(x) = 5x^3 - 5x - 24'

    !print first iteration
    write(2,*) iteration, f2(x_kplus1)

    !update f_x and df_x
    f_x_k = f2(x_k)
    df_x_k = df2(x_k)

    !update x_kplus1
    x_kplus1 = x_k - f_x_k/df_x_k

    !update iteration
    iteration = iteration + 1

    !print second iteration
    write(2,*) iteration, f2(x_kplus1)

    !Newton-Raphson method
    do while (abs(x_kplus1 - x_k) > 1e-6)

        !update x_k
        x_k = x_kplus1

        !update f_x and df_x
        f_x_k = f2(x_k)
        df_x_k = df2(x_k)

        !update x_kplus1
        x_kplus1 = x_k - f_x_k/df_x_k

        !update iteration
        iteration = iteration + 1

        !print iteration
        write(2,*) iteration, f2(x_kplus1)
        
    end do

    !print root
    write(2,*) 'Root:', x_kplus1

    !close first output file
    close(2)

    !request initial guess
    write(*,*) "Insert pre-initial guess:"
    read(*,*) pre_initial_guess

    !open second output file
    open(unit=3, file='secant1.txt', action='write')

    !reinitialize variables
    x_kminus1 = pre_initial_guess
    x_k = initial_guess
    x_kplus1 = 0.0
    f_x_kminus1 = f1(x_kminus1)
    f_x_k = f1(x_k)
    iteration = 1

    !print header
    write(3,*) 'Secant Method'
    write(3,*) 'Initial guess: ', x_k
    write(3,*) 'Pre-initial guess: ', x_kminus1
    write(3,*) 'f(x) = x^2 - 5'

    !print first iteration
    write(3,*) iteration, f1(x_kplus1)

    !update f_x_k and f_x_k-1
    f_x_kminus1 = f1(x_kminus1)
    f_x_k = f1(x_k)

    !update x_kplus1 and x_kminus1
    x_kplus1 = x_k - f_x_k*(x_k - x_kminus1)/(f_x_k - f_x_kminus1)
    x_kminus1 = x_k

    !update iteration
    iteration = iteration + 1

    !print second iteration
    write(3,*) iteration, f1(x_kplus1)

    !Secant method
    do while (abs(x_kplus1 - x_k) > 1e-6)

        !update x_k
        x_k = x_kplus1

        !update f_x_k and f_x_k-1
        f_x_kminus1 = f1(x_kminus1)
        f_x_k = f1(x_k)

        !update x_kplus1 and x_kminus1
        x_kplus1 = x_k - f_x_k*(x_k - x_kminus1)/(f_x_k - f_x_kminus1)
        x_kminus1 = x_k

        !update iteration
        iteration = iteration + 1

        !print iteration
        write(3,*) iteration, f1(x_kplus1)

    end do

    !print root
    write(3,*) 'Root:', x_kplus1

    !close second output file
    close(3)

    !open second output file
    open(unit=4, file='secant2.txt', action='write')

    !reinitialize variables
    x_kminus1 = pre_initial_guess
    x_k = initial_guess
    x_kplus1 = 0.0
    f_x_kminus1 = f2(x_kminus1)
    f_x_k = f2(x_k)
    iteration = 1

    !print header
    write(4,*) 'Secant Method'
    write(4,*) 'Initial guess: ', x_k
    write(4,*) 'Pre-initial guess: ', x_kminus1
    write(4,*) 'f(x) = 5x^3 - 5x - 24'

    !print first iteration
    write(4,*) iteration, f2(x_kplus1)

    !update f_x_k and f_x_k-1
    f_x_kminus1 = f2(x_kminus1)
    f_x_k = f2(x_k)

    !update x_kplus1 and x_kminus1
    x_kplus1 = x_k - f_x_k*(x_k - x_kminus1)/(f_x_k - f_x_kminus1)
    x_kminus1 = x_k

    !update iteration
    iteration = iteration + 1

    !print second iteration
    write(4,*) iteration, f2(x_kplus1)

    !Secant method
    do while (abs(x_kplus1 - x_k) > 1e-6)

        !update x_k
        x_k = x_kplus1

        !update f_x_k and f_x_k-1
        f_x_kminus1 = f2(x_kminus1)
        f_x_k = f2(x_k)

        !update x_kplus1 and x_kminus1
        x_kplus1 = x_k - f_x_k*(x_k - x_kminus1)/(f_x_k - f_x_kminus1)
        x_kminus1 = x_k

        !update iteration
        iteration = iteration + 1

        !print iteration
        write(4,*) iteration, f2(x_kplus1)

    end do

    !print root
    write(4,*) 'Root:', x_kplus1

    !close second output file
    close(4)

contains

    function f1(x) result(result)
        real, intent(in) :: x
        real result
        
        result = x**2. - 5.
    end function f1

    function f2(x) result(result)
        real, intent(in) :: x
        real result
        
        result = 5.*x**3. - 5.*x - 24.
    end function f2

    function df1(x) result(result)
        real, intent(in) :: x
        real result
        
        result = 2.*x
    end function df1

    function df2(x) result(result)
        real, intent(in) :: x
        real result
        
        result = 15.*x**2. - 5.
    end function df2

end program find_roots