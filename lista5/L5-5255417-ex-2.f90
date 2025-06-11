!------------------------------------------------------------
! File: L5-5255417-ex-2.f90
!
! Description:
!   Determining the ground-state energy (E0) of an infinite potential well by finding the eigenvalue of the discretized second-derivative matrix (in the limit n->infinity) via the power method.
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
program InfinitePotentialWell

    !deactivate implicit typing
    implicit none

    !define precision kind parameter
    integer, parameter :: dp = selected_real_kind(15, 307)

    !define parameters
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

    !declare variables
    integer :: matrix_dimension, n
    real(kind=dp) :: tolerance, old_under_lambda, under_lambda, computed_E0, real_E0, A_main_diagonal, A_other_diagonals, difference, delta_x, norm
    real(kind=dp), allocatable :: v(:), w(:)

    !define n
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !define tolerance
    write(*,*) "Insert tolerance"
    read(*,*) tolerance

    !allocate vectors
    allocate(v(matrix_dimension))
    allocate(w(matrix_dimension))

    !compute delta_x
    delta_x = 1.0_dp/(matrix_dimension+1.0_dp)

    !define diagonals
    A_main_diagonal = -2.0_dp
    A_other_diagonals = 1.0_dp

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

    !compute difference
    call compute_difference(v, difference, under_lambda, A_main_diagonal, A_other_diagonals, matrix_dimension)

    do while (difference.gt.tolerance .or. abs(old_under_lambda-under_lambda).gt.tolerance)

        !compute next v
        call solve_Aw(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call normalizew(v, w, matrix_dimension, tolerance)

        !save old_under_lambda
        old_under_lambda = under_lambda

        !compute difference
        call compute_difference(v, difference, under_lambda, A_main_diagonal, A_other_diagonals, matrix_dimension)

    end do 

    !real E0
    real_E0 = pi**2.0_dp

    !computed E0
    computed_E0 = -under_lambda**(-1.0_dp)/(delta_x**2.0_dp)

    write(*,'(A, F12.7, A, F12.7)') "computed E0: ",computed_E0, " real E0: ",real_E0

    !open results file
    open(unit=1, file="L5-5255417-ex-2-results.txt", status="replace")

    !print phi(n=0)
    write(1,*)0.0_dp,0.0_dp

    do n = 1, matrix_dimension

        !print point phi(delta_x*n) =  v(n)/sqrt(delta_x)
        write(1,*)delta_x*n,v(n)/sqrt(real(delta_x, dp))

    end do

    !print phi(n=matrix_dimension+1)
    write(1,*)0.0_dp,0.0_dp

    !close file
    close(1)

    !deallocate arrays
    deallocate(v)
    deallocate(w)

contains

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

    subroutine computE0UnderLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(inout) :: lambda
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp) :: w(matrix_dimension)

        call solve_Aw(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call vTw(lambda, v, w, matrix_dimension)


    end subroutine computE0UnderLambda

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

    subroutine compute_difference(v, difference, lambda, A_main_diagonal, A_other_diagonals, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: i
        integer, intent(in) :: matrix_dimension
        real(kind=dp), intent(in) :: v(matrix_dimension), A_main_diagonal, A_other_diagonals
        real(kind=dp), intent(inout) :: difference, lambda
        real(kind=dp) :: aux1(matrix_dimension), aux2(matrix_dimension), residual(matrix_dimension)

        !compute difference
        call computE0UnderLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
        call solve_Aw(v, aux1, matrix_dimension, A_main_diagonal, A_other_diagonals)
        call lambdav(v, aux2, lambda, matrix_dimension)

        do i = 1,matrix_dimension
            residual(i) = aux1(i) - aux2(i)
        end do

        call vTw(difference, residual, residual, matrix_dimension)
        difference = sqrt(real(difference, dp))

    end subroutine compute_difference

end program InfinitePotentialWell