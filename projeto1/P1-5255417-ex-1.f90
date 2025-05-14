!------------------------------------------------------------
! File: P1-5255417-ex-1.f90
!
! Description:
!   Solve Poisson's equation using the Numerov algorithm
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

    !declare variables
    integer(8) :: i
    real(8) :: deltar, phi_i_minus_1, phi_i, phi_i_plus_1, r_i_minus_1, r_i, r_i_plus_1

    !initialize i
    i = 0

    !define delta r
    write(*,*) "Insert delta r:"
    read(*,*) deltar

    !initialize phi(0) and phi(delta r)
    phi_i_minus_1 = -1.0d0
    phi_i = -1.0d0 + deltar/2.0d0

    !open files for writing results
    open(1, file='results-analytic.txt', status='replace')
    open(2, file='results-zero-infinity.txt', status='replace')
    open(3, file='results-infinity-zero.txt', status='replace')
    
    !compute Numerov method
    do while (abs(phi_i) > deltar)
        
        !compute i
        i = i + 1

        !compute r_i-1, r_i, r_i+1
        r_i_minus_1 = (i-1)*deltar
        r_i         = i*deltar
        r_i_plus_1  = (i+1)*deltar

        !compute phi(i+1)
        phi_i_plus_1 = 2.0d0*phi_i - phi_i_minus_1 + (deltar**2.0d0)/12.0d0*( f(r_i_minus_1) + 10.0d0*f(r_i) + f(r_i_plus_1) )
        phi_i_minus_1 = phi_i
        phi_i = phi_i_plus_1

        !print analytic solution
        write(1,'(F12.8,1X,F12.8)') r_i, analytic_solution(r_i)

        !print zero-infinity solution
        write(2,'(F12.8,1X,F12.8)') r_i, phi_i

    end do

    !compute r_i
    r_i = (i+1)*deltar

    !print analytic last step
    write(1,'(F12.8,1X,F12.8)') r_i, analytic_solution(r_i)

    !print last step
    write(2,'(F12.8,1X,F12.8)') r_i, phi_i

    !close file
    close(1)
    close(2)

    !initialize phi(infinity) and phi(infinity-delta r)
    phi_i_plus_1 = 0.0d0
    phi_i = 0.0d0

    do while (i>0)
        
        !compute i
        i = i - 1

        !compute r_i-1, r_i, r_i+1
        r_i_minus_1 = (i-1)*deltar
        r_i         = i*deltar
        r_i_plus_1  = (i+1)*deltar

        !compute phi(i+1)
        phi_i_minus_1 = 2.0d0*phi_i - phi_i_plus_1 + (deltar**2.0d0)/12.0d0*( f(r_i_minus_1) + 10.0d0*f(r_i) + f(r_i_plus_1) )
        phi_i_plus_1 = phi_i
        phi_i = phi_i_minus_1

        !print infinity-zero solution
        write(3,'(F12.8,1X,F12.8)') r_i, phi_i

    end do

    !compute r_i
    r_i = i*deltar

    !print last step
    write(3,'(F12.8,1X,F12.8)') r_i, phi_i

    !close file
    close(3)

contains

    real(8) function f(r)

        !declare parameters
        real(8), intent(in) :: r

        !compute function
        f = -r*exp(-r)/2.0d0

    end function f

    real(8) function analytic_solution(r)

        !declare parameters
        real(8), intent(in) :: r

        !compute function
        analytic_solution = -exp(-r)*(r/2.0d0 + 1.0d0)

    end function analytic_solution

end program Numerov