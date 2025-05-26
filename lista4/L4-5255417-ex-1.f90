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

    !define double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    !define parameters
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

    !declare variables
    integer :: matrix_dimension, n, number_of_intervals
    real(dp), allocatable :: sign_change_intervals(:,:)
    real(dp) :: Pn, dPn, tolerance, real_lambda, lambda_initial, mean_error, error, max_error
    real(dp), allocatable :: lambda(:)

    !define matrix dimension
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !define tolerance
    write(*,*) "Insert tolerance"
    read(*,*) tolerance

    !allocate lambda
    allocate(lambda(0:matrix_dimension))

    !find intervals where polynomial changes sign
    call find_sign_change_intervals(-4.0_dp, 0.0_dp, matrix_dimension, 1000*matrix_dimension, sign_change_intervals, number_of_intervals)

    !for each interval found, apply Newton-Raphson starting from midpoint of interval
    do n = 1, number_of_intervals

        !initialize lambda
        lambda_initial = (sign_change_intervals(1,n) + sign_change_intervals(2,n)) / 2.0_dp
        lambda(n) = lambda_initial

        !compute Pn and dPn
        call computeP(lambda(n), matrix_dimension, Pn, dPn)

        !Newton-Raphson
        do while (abs(Pn) > tolerance)

            !update lambda
            lambda(n) = lambda(n) - Pn / dPn

            !compute Pn and dPn
            call computeP(lambda(n), matrix_dimension, Pn, dPn)

        end do
    
    end do

    !sort lambda
    call sort_lambda(lambda, matrix_dimension)

    !open results file
    open(unit=1, file="L4-5255417-ex-1-results.txt", status="replace")

    !initialize mean_error and max_error
    mean_error = 0.0_dp
    max_error = 0.0_dp

    do n = 1, matrix_dimension
    
        !compute real lambda
        real_lambda = -4.0_dp*( sin( n*pi/(2.0_dp*(matrix_dimension+1.0_dp) ) ) )**2.0_dp

        !compute error
        error = abs(lambda(n)-real_lambda)

        !update max_error
        if (max_error < error) max_error = error

        !update mean_error
        mean_error = mean_error + error

        !write results
        write(1,'(A,F18.12,1X,A,F18.12,1X,A,F18.12,1X,A,ES10.2)') &
        "eigenvalue:", lambda(n), &
        "real eigenvalue:", real_lambda, &
        "difference:", error, &
        "+/-", tolerance
        
    end do

    !normalize
    mean_error = mean_error/matrix_dimension

    !write mean_error, max_error and tolerance
    write(*,'(A,F18.12,1X,A,F18.12,1X,A,ES10.2)') &
    "Mean Error: ",mean_error, &
    "Max Error:",max_error, &
    "+/-", tolerance

    !close file
    close(1)

contains

    subroutine computeP(lambda, matrix_dimension, Pn, dPn)
        
        !deactivate implicit typing
        implicit none

        !define parameters
        real(dp), parameter :: A(0:1) = [-2.0_dp,1.0_dp]

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in)    :: lambda
        real(dp), intent(out)   :: Pn, dPn

        !declare local variables
        integer :: i
        real(dp) :: Pn_minus_1, Pn_minus_2, dPn_minus_1, dPn_minus_2

        !define P0 and P1
        Pn_minus_2 = A(1)
        Pn_minus_1 = A(0) - lambda

        !define dP0 and dP1
        dPn_minus_2 = 0.0_dp
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

    subroutine find_sign_change_intervals(lambda_lower_bound, lambda_upper_bound, matrix_dimension, number_of_points, sign_change_intervals, number_of_intervals)

        !deactivate impliciting typing
        implicit none

        !define variables
        real(dp), intent(in) :: lambda_lower_bound
        real(dp), intent(in) :: lambda_upper_bound
        integer, intent(in) :: matrix_dimension
        integer, intent(in) :: number_of_points
        integer, intent(out) :: number_of_intervals
        real(dp), allocatable, intent(out) :: sign_change_intervals(:,:)

        !define local variables:
        real(dp) :: step_size
        real(dp) :: previous_lambda, current_lambda
        real(dp) :: previous_P, current_P
        integer :: index
        real(dp) :: dummy_dP  ! We don't need dP here, but computeP requires it

        !calculate step size to sample points uniformly between the bounds
        step_size = (lambda_upper_bound - lambda_lower_bound) / real(number_of_points - 1)
        number_of_intervals = 0

        !allocate maximum possible size for intervals (number_of_points - 1 intervals maximum)
        allocate(sign_change_intervals(2, number_of_points - 1))

        !evaluate P(lambda) at the first sample point
        previous_lambda = lambda_lower_bound
        call computeP(previous_lambda, matrix_dimension, previous_P, dummy_dP)

        !loop over the rest of the sample points to detect sign changes
        do index = 1, number_of_points - 1

            current_lambda = lambda_lower_bound + index * step_size
            call computeP(current_lambda, matrix_dimension, current_P, dummy_dP)

            !check if there's a sign change between two consecutive points
            if (previous_P * current_P < 0.0_dp) then

                !found an interval where P(lambda) changes sign
                number_of_intervals = number_of_intervals + 1
                sign_change_intervals(1, number_of_intervals) = previous_lambda
                sign_change_intervals(2, number_of_intervals) = current_lambda

            end if

            !update variables for next iteration
            previous_lambda = current_lambda
            previous_P = current_P

        end do

        !resize array to actual number of intervals found
        if (number_of_intervals > 0) then
            sign_change_intervals = sign_change_intervals(:, 1:number_of_intervals)
        else
            deallocate(sign_change_intervals)
        end if

    end subroutine find_sign_change_intervals

    subroutine sort_lambda(lambda, n)

        !deactivate impliciting typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(inout) :: lambda(0:n)
        integer :: i, j
        real(dp) :: temp

        do i = 0, n-1
            do j = 0, n-1-i
                if (lambda(j) < lambda(j+1)) then

                    !switch
                    temp = lambda(j)
                    lambda(j) = lambda(j+1)
                    lambda(j+1) = temp

                end if
            end do
        end do

    end subroutine sort_lambda


end program Eigenvalues