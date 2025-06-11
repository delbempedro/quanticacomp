!------------------------------------------------------------
! File: L6-5255417-ex-1.f90
!
! Description:
!   Use the Householder method to tridiagonalize a matrix, find its largest and lowest eigenvalue/eigenvector using the Power Method, and check the solution.
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
program Tridiagonalization

    !deactivate implicit typing
    implicit none

    !define precision kind parameter
    integer, parameter :: dp = selected_real_kind(15, 307)

    !declare variables
    integer :: matrix_dimension, i
    real(dp), allocatable :: A(:,:), At(:,:), O(:,:), d(:), e(:)
    real(dp), allocatable :: yt_min(:), yt_max(:), y_min(:), y_max(:)
    real(dp) :: lambda_min, lambda_max

    !define matrix dimension
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !allocate vectors
    allocate(A(matrix_dimension, matrix_dimension))
    allocate(At(matrix_dimension, matrix_dimension))
    allocate(O(matrix_dimension, matrix_dimension))
    allocate(d(matrix_dimension), e(matrix_dimension))
    allocate(yt_min(matrix_dimension), yt_max(matrix_dimension))
    allocate(y_min(matrix_dimension), y_max(matrix_dimension))

    write(*, '(/A, I0, A)') "-------> N = ", matrix_dimension, " <-------"

    !initialize matrix A and a copy At
    call initialize_A(A, matrix_dimension)
    call initialize_A(At, matrix_dimension)
    write(*,*) "Original Matrix A (first 5x5):"
    call print_matrix(A, matrix_dimension, min(matrix_dimension, 5))

    !tridiagonalize At with Householder, saving the transformation matrix in O
    call householder_reduction(At, O, matrix_dimension)
    write(*,*) "Tridiagonal Matrix T (first 5x5):"
    call print_matrix(At, matrix_dimension, min(matrix_dimension, 5))

    !verify transformation T = O' * A * O
    call verify_transformation(matrix_dimension, At, O, A)

    !extract the diagonals d (main) and e (off-diagonal) of At
    call get_tridiagonal_elements(At, d, e, matrix_dimension)

    !find the smallest and largest eigenvalues and the corresponding eigenvectors of At
    call PowerMethod(yt_max, lambda_max, yt_min, lambda_min, d, e, matrix_dimension)

    !save yt max and min
    open(unit=1, file="At_eingenvalues.txt", status="replace")
    write(1, *) "Eigenvector yt max:"
    do i = 1, matrix_dimension
        write(1, '(F8.4, 1X)') yt_max(i)
    end do
    write(1, *) "Eigenvector yt min:"
    do i = 1, matrix_dimension
        write(1, '(F8.4, 1X)') yt_min(i)
    end do
    close(1)

    !get eigenvectors of A, transform eigenvectors of At: y = O*yt
    y_max = matmul(O, yt_max)
    y_min = matmul(O, yt_min)

    !save y max and min
    open(unit=2, file="A_eingenvalues.txt", status="replace")
    write(2, *) "Eigenvector y max:"
    do i = 1, matrix_dimension
        write(2, '(F8.4, 1X)') y_max(i)
    end do
    write(2, *) "Eigenvector y min:"
    do i = 1, matrix_dimension
        write(2, '(F8.4, 1X)') y_min(i)
    end do
    close(2)

    !write eingenvalues
    write(*, '(/A, F12.8)') "Smallest eigenvalue (lambda_min): ", lambda_min
    write(*, '(A, F12.8)')  "Largest eigenvalue (lambda_max):  ", lambda_max
    
    !verify A*y = lambda*y for both cases
    call verify_solution(lambda_min, y_min, matrix_dimension, "smallest")
    call verify_solution(lambda_max, y_max, matrix_dimension, "biggest")

    !deallocate arrays
    deallocate(A, At, O, d, e, yt_min, yt_max, y_min, y_max)

contains

    subroutine initialize_A(A, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(out) :: A(matrix_dimension,matrix_dimension)
        integer :: i

        !initialize A with zeros
        A = 0.0_dp
        do i = 1, matrix_dimension

            !main diagonal
            A(i,i) = -5.0_dp / 2.0_dp

            !secondary diagonal
            if (i <= matrix_dimension-1) then
                A(i, i+1) = 4.0_dp / 3.0_dp
                A(i+1, i) = 4.0_dp / 3.0_dp
            end if

            !tertiary diagonal
            if (i <= matrix_dimension-2) then
                A(i, i+2) = -1.0_dp / 12.0_dp
                A(i+2, i) = -1.0_dp / 12.0_dp
            end if

        end do

    end subroutine initialize_A

    subroutine householder_reduction(A, O, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(inout) :: A(matrix_dimension,matrix_dimension)
        real(dp), intent(out) :: O(matrix_dimension,matrix_dimension)
        integer :: i, j

        real(dp) :: sigma, beta, mu
        real(dp), allocatable :: u(:), p(:), q(:), aux(:)

        !allocate vectors
        allocate(u(matrix_dimension), p(matrix_dimension), q(matrix_dimension), aux(matrix_dimension))

        !initialize O as the identity matrix
        O = 0.0_dp
        do i = 1, matrix_dimension
            O(i,i) = 1.0_dp
        end do

        !Householder reduction loop
        do i = 1, matrix_dimension - 2

            !calculate sigma
            sigma = norm2(A(i+1:matrix_dimension, i))

            if (sigma > 1.0e-15_dp) then
                if (A(i + 1, i) < 0.0_dp) then
                    sigma = -sigma
                end if

                !calculate Householder vector u
                u = 0.0_dp
                u(i + 1) = A(i + 1, i) + sigma
                do j = i + 2, matrix_dimension
                    u(j) = A(j, i)
                end do

                !calculate beta
                beta = 1.0_dp / (sigma * u(i + 1))

                !update matrix A using similarity transformation A = P*A*P
                p = matmul(A, u)
                p = beta * p 
                mu = beta * dot_product(p, u) / 2.0_dp
                q = p - mu * u
                A = A - outer_product(u, q) - outer_product(q, u)

                !update transformation matrix O = O*P
                O = O - beta * matmul(O, outer_product(u, u))

            end if

        end do

        !deallocate vectors
        deallocate(u, p, q)

    end subroutine householder_reduction

    function outer_product(v1, v2) result(res)

        !deactivate implicit typing
        implicit none

        !declare variables
        real(dp), intent(in) :: v1(:), v2(:)

        integer :: matrix_dimension
        real(dp), allocatable :: res(:,:)

        !define matrix_dimension
        matrix_dimension = size(v1)

        !allocate vector
        allocate(res(matrix_dimension, matrix_dimension))

        !do outer product
        res = spread(v1, dim=2, ncopies=matrix_dimension) * spread(v2, dim=1, ncopies=matrix_dimension)

    end function outer_product

    subroutine verify_transformation(matrix_dimension, A_t, O, A_original)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: A_t(matrix_dimension, matrix_dimension), O(matrix_dimension, matrix_dimension), A_original(matrix_dimension, matrix_dimension)

        real(dp) :: diff_norm
        real(dp), allocatable :: C(:,:), TEMP(:,:)
        
        !allocate vectors
        allocate(C(matrix_dimension,matrix_dimension), TEMP(matrix_dimension,matrix_dimension))
        
        !TEMP = O^T A
        TEMP = matmul(transpose(O), A_original)

        !C = O^T A O
        C = matmul(TEMP, O)
        
        !compute norm of the residual matrix
        diff_norm = norm2(A_t - C)

        !write results
        write(*,*)
        write(*,*) "--- Transformation Verification O^T A O = At ---"
        write(*, '(A, E12.5)') "Difference norm ||At - O^T A O||: ", diff_norm

        !deallocate vectors
        deallocate(C, TEMP)

    end subroutine verify_transformation

    subroutine get_tridiagonal_elements(At, d, e, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: At(matrix_dimension,matrix_dimension)
        real(dp), intent(out) :: d(matrix_dimension), e(matrix_dimension)
        integer :: i

        e = 0.0_dp
        do i = 1, matrix_dimension
            d(i) = At(i,i) !main diagonal
            if (i < matrix_dimension) then
                e(i) = At(i, i+1) !secondary diagonal
            end if
        end do

    end subroutine get_tridiagonal_elements
    
    subroutine compute_Ay(y, Ay, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: y(matrix_dimension)
        real(dp), intent(out) :: Ay(matrix_dimension)
        integer :: i
        real(dp), parameter :: d0 = -5.0_dp/2.0_dp, d1 = 4.0_dp/3.0_dp, d2 = -1.0_dp/12.0_dp
        
        !initialize vector
        Ay = 0.0_dp
        if (matrix_dimension <= 0) return

        !main loop for the central part of the matrix where all 5 diagonals are safe to access
        do i = 3, matrix_dimension - 2
            Ay(i) = d2*y(i-2) + d1*y(i-1) + d0*y(i) + d1*y(i+1) + d2*y(i+2)
        end do

        !handle the edges explicitly to avoid compiler warnings
        !row 1
        Ay(1) = d0*y(1)
        if (matrix_dimension >= 2) Ay(1) = Ay(1) + d1*y(2)
        if (matrix_dimension >= 3) Ay(1) = Ay(1) + d2*y(3)

        !row 2
        if (matrix_dimension >= 2) then
            Ay(2) = d1*y(1) + d0*y(2)
            if (matrix_dimension >= 3) Ay(2) = Ay(2) + d1*y(3)
            if (matrix_dimension >= 4) Ay(2) = Ay(2) + d2*y(4)
        end if

        !second to last row (only if matrix_dimension > 3, as matrix_dimension=1,2,3 cases are handled by above logic)
        if (matrix_dimension > 3) then
            Ay(matrix_dimension-1) = d2*y(matrix_dimension-3) + d1*y(matrix_dimension-2) + d0*y(matrix_dimension-1) + d1*y(matrix_dimension)
        end if

        !last row (only if matrix_dimension > 2)
        if (matrix_dimension > 2) then
            Ay(matrix_dimension) = d2*y(matrix_dimension-2) + d1*y(matrix_dimension-1) + d0*y(matrix_dimension)
        end if

    end subroutine compute_Ay

    subroutine verify_solution(lambda, y, matrix_dimension, label)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: lambda, y(matrix_dimension)
        character(len=*), intent(in) :: label

        real(dp) :: residual_norm
        real(dp), allocatable :: Ay(:), residual(:)

        !allocate vectors
        allocate(Ay(matrix_dimension), residual(matrix_dimension))
        
        !compute Ay = A * y for the original 5-diagonal matrix
        call compute_Ay(y, Ay, matrix_dimension)

        !compute residual vector, residual = Ay - lambda*y
        residual = Ay - lambda * y
        
        !compute residual norm ||residual||
        residual_norm = norm2(residual)

        !write results
        write(*, '(/A, A, A)') "--- Verification for Eigenvalue ", trim(label), " ---"
        write(*, '(A, F12.8)') "Eigenvalue:", lambda
        write(*, '(A, E12.5)') "Norm of residual ||Ay - lambday||: ", residual_norm

        !deallocate vectors
        deallocate(Ay, residual)

    end subroutine verify_solution

    subroutine print_matrix(matrix, matrix_dimension, n_print)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension, n_print
        real(dp), intent(in) :: matrix(matrix_dimension, matrix_dimension)
        integer :: i

        do i = 1, n_print
            write(*, '(10F12.6)') matrix(i, 1:n_print)
        end do

    end subroutine print_matrix

    !----------------------------- POWER METHOD SUBROUTINES -----------------------------!

    subroutine Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
    
        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: w(matrix_dimension)
        integer :: i
        
        !matrix_dimensions = 1
        if (matrix_dimension == 1) then
            w(1) = A_main_diagonal(1) * v(1)
            return
        end if
        
        !matrix_dimensions > 1
        w(1) = A_main_diagonal(1) * v(1) + A_other_diagonals(1) * v(2)
        do i = 2, matrix_dimension - 1
            w(i) = A_other_diagonals(i-1) * v(i-1) + A_main_diagonal(i) * v(i) + A_other_diagonals(i) * v(i+1)
        end do
        w(matrix_dimension) = A_other_diagonals(matrix_dimension-1) * v(matrix_dimension-1) + A_main_diagonal(matrix_dimension) * v(matrix_dimension)

    end subroutine Av

    subroutine solve_Aw(w, v, matrix_dimension, A_main_diagonal, A_other_diagonals)
        
        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: w(matrix_dimension)

        real(dp), allocatable :: c_prime(:), d_prime(:)
        integer :: i
        real(dp) :: m
        
        !matrix_dimension = 1
        if (matrix_dimension == 1) then
            w(1) = v(1) / A_main_diagonal(1)
            return
        end if
        
        !matrix_dimendion > 1

        !allocate vectors
        allocate(c_prime(matrix_dimension-1), d_prime(matrix_dimension))

        c_prime(1) = A_other_diagonals(1) / A_main_diagonal(1)
        d_prime(1) = v(1) / A_main_diagonal(1)

        do i = 2, matrix_dimension - 1
            m = A_main_diagonal(i) - A_other_diagonals(i-1) * c_prime(i-1)
            c_prime(i) = A_other_diagonals(i) / m
            d_prime(i) = (v(i) - A_other_diagonals(i-1) * d_prime(i-1)) / m
        end do
        
        m = A_main_diagonal(matrix_dimension) - A_other_diagonals(matrix_dimension-1) * c_prime(matrix_dimension-1)
        d_prime(matrix_dimension) = (v(matrix_dimension) - A_other_diagonals(matrix_dimension-1) * d_prime(matrix_dimension-1)) / m

        w(matrix_dimension) = d_prime(matrix_dimension)
        do i = matrix_dimension - 1, 1, -1
            w(i) = d_prime(i) - c_prime(i) * w(i+1)
        end do

        !deallocate vectors
        deallocate(c_prime, d_prime)

    end subroutine solve_Aw

    subroutine normalizew(v, w, matrix_dimension, tolerance)
        
        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: w(matrix_dimension), tolerance
        real(dp), intent(out) :: v(matrix_dimension)

        real(dp) :: normalization
        
        normalization = norm2(w)
        if (normalization < tolerance) then
            write(*,*) "Warning: vector norm is close to zero in normalizew."
            v = 0.0_dp
        else
            v = w / normalization
        end if

    end subroutine normalizew
    
    subroutine computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
        
        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: lambda

        real(dp), allocatable :: w(:)
        real(dp) :: v_dot_w, v_dot_v
        
        !allocate vectors
        allocate(w(matrix_dimension))

        !w = AV
        call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)

        !compute vw
        v_dot_w = dot_product(v, w)

        !compute vv
        v_dot_v = dot_product(v, v)

        !verify if v ~ 0 and compute lambda
        if (v_dot_v > 1.0e-20_dp) then
           lambda = v_dot_w / v_dot_v
        else
           lambda = 0.0_dp
        end if

        !deallocate vectors
        deallocate(w)

    end subroutine computeLambda

    subroutine PowerMethod(vector_max, lambda_max, vector_min, lambda_min, A_main_diagonal, A_other_diagonals, matrix_dimension)
        
        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: vector_max(matrix_dimension), lambda_max, vector_min(matrix_dimension), lambda_min
        
        real(dp), allocatable :: v(:), w(:)
        real(dp) :: tolerance, old_lambda, lambda, old_under_lambda, under_lambda
        integer :: iter
        
        !define tolerance
        tolerance = 1.0e-12_dp

        !allocate vectors
        allocate(v(matrix_dimension), w(matrix_dimension))

        ! --- Max Eigenvalue ---

            !initialize v vector with random numbers
        call random_number(v)

        !normalize the initial vector
        call vTw(difference, v, v, matrix_dimension)
        difference = sqrt(real(difference, dp))
        do n = 1, matrix_dimension
            v(n) = v(n) / difference
        end do

        old_lambda = 0.0_dp
        iter = 0
        do while(iter < 1000)
            iter = iter + 1
            call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
            call normalizew(v, w, matrix_dimension, tolerance)
            call computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
            if (abs(lambda - old_lambda) < tolerance) exit
            old_lambda = lambda
        end do

        !save lambda_max and vector_max
        lambda_max = lambda
        vector_max = v

        ! --- Min Eigenvalue ---

        !initialize v vector with random numbers
        call random_number(v)

        !normalize the initial vector
        call vTw(difference, v, v, matrix_dimension)
        difference = sqrt(real(difference, dp))
        do n = 1, matrix_dimension
            v(n) = v(n) / difference
        end do
        
        old_under_lambda = 0.0_dp
        iter = 0
        do while(iter < 1000)
           iter = iter + 1
           call solve_Aw(w, v, matrix_dimension, A_main_diagonal, A_other_diagonals)
           under_lambda = dot_product(v, w) ! under_lambda is now v' * A^-1 * v
           call normalizew(v, w, matrix_dimension, tolerance)
           if (iter > 1 .and. abs(1.0_dp/under_lambda - 1.0_dp/old_under_lambda) < tolerance) exit
           old_under_lambda = under_lambda
        end do

        !save lambda_min and vector_min
        lambda_min = 1.0_dp / under_lambda
        vector_min = v
        
        !deallocate vectors
        deallocate(v, w)

    end subroutine PowerMethod
    
end program Tridiagonalization
