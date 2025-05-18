!------------------------------------------------------------
! File: L4-5255417-ex-1.f90
!
! Description:
!   Finds the eigenvalues of the matrix discretization of the second derivative
!
! Dependencies:
!   - None
!
! Since:
!   - 05/2025
!
! Authors:
!   - Pedro C. Delbem <pedrodelbem@usp.br>
!------------------------------------------------------------
program Eigenvalues

    !deactivate implicit typing
    implicit none

    !define parameters
    real, parameter :: pi = 4.0*atan(1.0)

    !declare variables
    integer :: matrix_dimension, n
    real :: Pn, dPn, tolerance, real_lambda
    real, allocatable :: lambda(:)

    !define n
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !define tolerance
    write(*,*) "Insert tolerance"
    read(*,*) tolerance

    !allocate lambda
    allocate(lambda(0:matrix_dimension))

    !initialize lambda, Pn and dPn
    lambda(0) = 0.0
    Pn = 1.0
    dPn = 1.0/tolerance

    !open results file
    open(unit=1, file="L4-5255417-ex-1-results.txt", status="replace")

    !search eingenvalues
    do n = 1, matrix_dimension

        !reinitialize lambda
        lambda(n) = 1.1*lambda(n-1)

        GO TO 10

20      lambda(n) = lambda(n) + 4.0/matrix_dimension
        write(*,*) lambda

        !reinitialize Pn and dPn
10      Pn = 1.0
        dPn = 1.0/tolerance
        write(*,*) Pn, lambda(n), n

        !update Pn and dPn
        call computeP(lambda(n), matrix_dimension, Pn, dPn)
        write(*,*) Pn, lambda(n), n

        do while(abs(Pn) > tolerance)

            !update lambda
            lambda(n) = lambda(n) - Pn/dPn

            !update Pn and dPn
            call computeP(lambda(n), matrix_dimension, Pn, dPn)
            !write(*,*) Pn, n

        end do

        !write(*,*)isIn(lambda, n, matrix_dimension, tolerance)
        if (abs(lambda(n)-lambda(n-1))<tolerance) then
            write(*,*) "I'm Here"
            GO TO 20
        end if

        real_lambda = -4.0*( sin( n*pi/(2.0*(matrix_dimension+1.0 ) ) ) )**2.0
        write(1,*)"lambda:",lambda(n),"real lambda:",real_lambda,"diference:",abs(lambda(n)-real_lambda)

    end do

    !close file
    close(1)

contains

    subroutine computeP(lambda, matrix_dimension, Pn, dPn)
        
        !deactivate implicit typing
        implicit none

        !define parameters
        real, parameter :: A(0:1) = [-2.0,1.0]

        !declare variables
        integer, intent(in) :: matrix_dimension
        real, intent(in)    :: lambda
        real, intent(out)   :: Pn, dPn

        !declare local variables
        integer :: i
        real :: Pn_minus_1, Pn_minus_2, dPn_minus_1, dPn_minus_2

        !define P0 and P1
        Pn_minus_2 = A(1)
        Pn_minus_1 = A(0) - lambda

        !define dP0 and dP1
        dPn_minus_2 = 0.0
        dPn_minus_1 = -A(1)

        !if matrix dimension is zero or one
        if (matrix_dimension == 0) then
            Pn  = Pn_minus_2
            dPn = dPn_minus_2
            return
        else if (matrix_dimension == 1) then
            Pn  = Pn_minus_1
            dPn = dPn_minus_1
            return
        end if

        !if matrix dimension is two or bigger
        do i = 2, matrix_dimension

            !update Pn and dPn
            Pn = (A(0) - lambda)*Pn_minus_1 - A(1)*A(1)*Pn_minus_2
            dPn = (A(0) - lambda)*dPn_minus_1 - A(1)*A(1)*dPn_minus_2 - Pn_minus_1

            !update old ones
            Pn_minus_2 = Pn_minus_1
            Pn_minus_1 = Pn
            dPn_minus_2 = dPn_minus_1
            dPn_minus_1 = dPn

        end do

    end subroutine computeP

    logical function isIn(lambda, n, matrix_dimension, tolerance)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n, matrix_dimension
        real, intent(in) :: lambda(0:matrix_dimension), tolerance
        integer :: m

        isIn = .false.

        do m = 1, n-1
            write(*,*)abs(lambda(m)-lambda(n))
            if (abs(lambda(m)-lambda(n)) < tolerance) then
                isIn = .true.
                RETURN
            end if

        end do

        do m = n+1, matrix_dimension
            write(*,*)abs(lambda(m)-lambda(n))
            if (abs(lambda(m)-lambda(n)) < tolerance) then
                isIn = .true.
                RETURN
            end if

        end do

    end function isIn

end program Eigenvalues