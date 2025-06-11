!------------------------------------------------------------
! File: L5-5255417-ex-1.f90
!
! Description:
!   Finds the eigenvalues of the matrix discretization of the second derivative with power method
!
! Dependencies:
!   - None
!
! Since:
!   - 06/2025
!
! Authors:
!   - Pedro C. Delbem <pedrodelbem@usp.br>
!------------------------------------------------------------
program PowerMethod

    !deactivate implicit typing
    implicit none

    !define precision kind parameter
    integer, parameter :: dp = selected_real_kind(15, 307)

    !define parameters
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

    !declare variables
    integer :: matrix_dimension, n, control
    real(kind=dp) :: tolerance, old_lambda, lambda, old_under_lambda, under_lambda, old_lambdaB, lambdaB, capital_lambda_squared, real_lambda, A_main_diagonal, A_other_diagonals, difference, norm
    real(kind=dp), allocatable :: v(:), w(:)

    !define n
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !define tolerance
    write(*,*) "Insert Capital Lambda^2"
    read(*,*) capital_lambda_squared

    !allocate vectors
    allocate(v(matrix_dimension))
    allocate(w(matrix_dimension))

    !define diagonals
    A_main_diagonal = -2.0_dp
    A_other_diagonals = 1.0_dp

    !define tolerance
    tolerance = 1.0e-12

    !open file
    open(unit=1, file="eingenvalues.txt", status="replace")

    !A

    !initialize v vector with random numbers
    call random_number(v)

    !normalize the initial vector
    call vTw(norm, v, v, matrix_dimension)
    norm = sqrt(real(norm, dp))
    do n = 1, matrix_dimension
        v(n) = v(n) / norm
    end do

    !initialize old_lambda
    old_lambda = 0.0_dp

    !set control
    control = 1

    !compute difference
    call compute_difference(v, difference, lambda, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    do while (difference.gt.tolerance .or. abs(old_lambda-lambda).gt.tolerance)

        !compute next v
        call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call normalizew(v, w, matrix_dimension, tolerance)

        !save old_lambda
        old_lambda = lambda

        !compute difference
        call compute_difference(v, difference, lambda, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    end do

    !save eingenvector of max absolute eingenvalue
    write(1, *) "Eigenvector of max absolute eingenvalue:"
    do n = 1, matrix_dimension
        write(1, '(F8.4, 1X)') v(n)
    end do

    !compute real lambda
    real_lambda = -4.0_dp*( dsin( matrix_dimension*pi/(2.0_dp*(matrix_dimension+1.0_dp) ) ) )**2.0_dp

    write(*,*) "biggest absolut eingenvalue: ",lambda, "real eigenvalue: ",real_lambda

    !A^-1

    !initialize v vector with random numbers
    call random_number(v)

    !normalize the initial vector
    call vTw(norm, v, v, matrix_dimension)
    norm = sqrt(real(norm, dp))
    do n = 1, matrix_dimension
        v(n) = v(n) / norm
    end do

    !initialize old_under_lambda
    old_under_lambda = 0.0_dp

    !set control
    control = 2

    !compute difference
    call compute_difference(v, difference, under_lambda, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    do while (difference.gt.tolerance .or. abs(old_under_lambda-under_lambda).gt.tolerance)

        !compute next v
        call solve_Aw(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call normalizew(v, w, matrix_dimension, tolerance)

        !save old_under_lambda
        old_under_lambda = under_lambda

        !compute difference
        call compute_difference(v, difference, under_lambda, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    end do

    !save eingenvector of min absolute eingenvalue
    write(1, *) "Eigenvector of min absolute eingenvalue:"
    do n = 1, matrix_dimension
        write(1, '(F8.4, 1X)') v(n)
    end do

    !compute real under_lambda
    real_lambda = -4.0_dp*( dsin( pi/(2.0_dp*(matrix_dimension+1.0_dp) ) ) )**2.0_dp

    write(*,*) "smallest absolut eingenvalue: ",under_lambda**(-1.0_dp), "real eigenvalue: ",real_lambda

    !B = (I*capital_lambda^2 - A^2)

    !initialize v vector with random numbers
    call random_number(v)

    !normalize the initial vector
    call vTw(norm, v, v, matrix_dimension)
    norm = sqrt(real(norm, dp))
    do n = 1, matrix_dimension
        v(n) = v(n) / norm
    end do

    !initialize old_lambda
    old_lambdaB = 0.0_dp

    !set control
    control = 3

    !compute difference
    call compute_difference(v, difference, lambdaB, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    do while (difference.gt.tolerance .or. abs(old_lambdaB-lambdaB).gt.tolerance)

        !compute next v
        call Bv(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals, capital_lambda_squared)
        call normalizew(v, w, matrix_dimension, tolerance)

        !save old_lambda
        old_lambdaB = lambdaB

        !compute difference
        call compute_difference(v, difference, lambdaB, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

    end do

    !save eingenvector of last case
    write(1, *) "Eigenvector of last case:"
    do n = 1, matrix_dimension
        write(1, '(F8.4, 1X)') v(n)
    end do

    !compute lambda
    lambda = -sqrt(real(capital_lambda_squared-lambdaB, dp))

    write(*,*) "eingenvalue of last case: ",lambda

    !deallocate arrays
    deallocate(v)
    deallocate(w)

    !close file
    close(1)

contains

    subroutine Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp), intent(out) :: w(matrix_dimension)
        
        if (matrix_dimension.eq.1) then

            w(1) = v(1)*A_main_diagonal

        elseif (matrix_dimension.eq.2) then

            w(1) = v(1)*A_main_diagonal + v(2)*A_other_diagonals
            w(2) = v(1)*A_other_diagonals + v(2)*A_main_diagonal

        else

            !first element
            w(1) = v(1)*A_main_diagonal + v(2)*A_other_diagonals

            !middle elements
            do i = 2,matrix_dimension-1

                w(i) = v(i-1)*A_other_diagonals + v(i)*A_main_diagonal + v(i+1)*A_other_diagonals

            end do

            !last element
            w(matrix_dimension) = v(matrix_dimension-1)*A_other_diagonals + v(matrix_dimension)*A_main_diagonal

        end if

    end subroutine Av

    subroutine solve_Aw(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: A_main_diagonal, A_other_diagonals, v(matrix_dimension)
        real(kind=dp), intent(out) :: w(matrix_dimension)
        real(kind=dp) :: denominator, v_line(matrix_dimension), A_other_diagonals_line(matrix_dimension - 1) 

        if (matrix_dimension == 1) then
            w(1) = v(1) / A_main_diagonal
            return
        else

            !forward step
            A_other_diagonals_line(1) = A_other_diagonals / A_main_diagonal
            v_line(1) = v(1) / A_main_diagonal

            !compute c' and d'"
            do i = 2, matrix_dimension - 1

                denominator = A_main_diagonal - A_other_diagonals * A_other_diagonals_line(i-1)
                A_other_diagonals_line(i) = A_other_diagonals / denominator
                v_line(i) = (v(i) - A_other_diagonals * v_line(i-1)) / denominator

            end do

            !compute d' last element
            denominator = A_main_diagonal - A_other_diagonals * A_other_diagonals_line(matrix_dimension-1)
            v_line(matrix_dimension) = (v(matrix_dimension) - A_other_diagonals * v_line(matrix_dimension-1)) / denominator

            !back substitution
            w(matrix_dimension) = v_line(matrix_dimension)
            do i = matrix_dimension - 1, 1, -1
                w(i) = v_line(i) - A_other_diagonals_line(i) * w(i+1)
            end do

        end if

    end subroutine solve_Aw

    subroutine Bv(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals, capital_lambda_squared)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: v(matrix_dimension)
        real(kind=dp), intent(out) :: w(matrix_dimension)
        real(kind=dp), intent(in) :: A_main_diagonal, A_other_diagonals, capital_lambda_squared
        real(kind=dp) :: aux(matrix_dimension), A2v(matrix_dimension)
        

        !Av = A*v
        call Av(v, aux, matrix_dimension, A_main_diagonal, A_other_diagonals)

        !A2v = A * Av  (A^2*v)
        call Av(aux, A2v, matrix_dimension, A_main_diagonal, A_other_diagonals)

        !w = capital_lambda_squared * v - A2v_temp
        do i = 1, matrix_dimension
            w(i) = capital_lambda_squared*v(i) - A2v(i)
        end do

    end subroutine Bv

    subroutine vTw(result, v, w, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: v(matrix_dimension), w(matrix_dimension)
        real(kind=dp), intent(inout) :: result

        !initialize results
        result = 0.0_dp

        do i = 1,matrix_dimension

            result = result + v(i)*w(i)

        end do

    end subroutine vTw

    subroutine computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(inout) :: lambda
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp) :: w(matrix_dimension)

        call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call vTw(lambda, v, w, matrix_dimension)

    end subroutine computeLambda

    subroutine computeLambdaB(lambda_squared, v, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(inout) :: lambda_squared
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals, capital_lambda_squared
        real(kind=dp) :: w(matrix_dimension)

        call Bv(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals, capital_lambda_squared)
        call vTw(lambda_squared, v, w, matrix_dimension)

    end subroutine computeLambdaB

    subroutine compute1UnderLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(inout) :: lambda
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp) :: w(matrix_dimension)

        call solve_Aw(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call vTw(lambda, v, w, matrix_dimension)


    end subroutine compute1UnderLambda

    subroutine lambdav(v, w, lambda, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: lambda, v(matrix_dimension)
        real(kind=dp), intent(out) :: w(matrix_dimension)

        do i = 1,matrix_dimension

            w(i) = v(i)*lambda

        end do

    end subroutine lambdav

    subroutine normalizew(v, w, matrix_dimension, tolerance)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: w(matrix_dimension), tolerance
        real(kind=dp), intent(out) :: v(matrix_dimension)
        real(kind=dp) normalization

        call vTw(normalization, w, w, matrix_dimension)
        normalization = sqrt(real(normalization, dp))

        if (normalization.lt.tolerance) then

            write(*,*) "Warning!"
        else

            do i = 1,matrix_dimension

                v(i) = w(i)/normalization

            end do
        end if

    end subroutine normalizew

    subroutine compute_difference(v, difference, lambda, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared, control)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension, control
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp), intent(inout) :: difference, lambda, capital_lambda_squared
        real(kind=dp) :: aux1(matrix_dimension), aux2(matrix_dimension), residual(matrix_dimension)

        !compute difference
        if (control.eq.1) then
            call computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
            call Av(v, aux1, matrix_dimension, A_main_diagonal, A_other_diagonals)
        elseif (control.eq.2) then
            call compute1UnderLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
            call solve_Aw(v, aux1, matrix_dimension, A_main_diagonal, A_other_diagonals)
        elseif (control.eq.3) then
            call computeLambdaB(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension, capital_lambda_squared)
            call Bv(v, aux1, matrix_dimension, A_main_diagonal, A_other_diagonals, capital_lambda_squared)
        end if
        call lambdav(v, aux2, lambda, matrix_dimension)

        do i = 1,matrix_dimension
            residual(i) = aux1(i) - aux2(i)
        end do

        call vTw(difference, residual, residual, matrix_dimension)
        difference = sqrt(real(difference, dp))

    end subroutine compute_difference

end program PowerMethod