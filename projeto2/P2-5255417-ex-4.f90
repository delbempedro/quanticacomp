!------------------------------------------------------------
! File: P2-5255417-ex-4.f90
!
! Description:
!   Use the Householder (with optimized memory) method to tridiagonalize a matrix, find its largest and lowest eigenvalue/eigenvector using the Power Method, and check the solution.
!
! Dependencies:
!   - None
!
! Since:
!   - 07/2025
!
! Authors:
!   - Pedro C. Delbem <pedrodelbem@usp.br>
!------------------------------------------------------------
program Tridiagonalization

    !deactivate implicit typing
    implicit none

    !define a kind parameter for double precision
    integer, parameter :: dp = selected_real_kind(15, 307)

    !declare variables
    integer :: matrix_dimension, i
    real(dp), allocatable :: At(:,:)
    real(dp), allocatable :: d(:), e(:)
    real(dp), allocatable :: yt_min(:), yt_max(:), y_min(:), y_max(:)
    real(dp) :: lambda_min, lambda_max
    real(dp), allocatable :: betas(:), u1_storage(:)

    
    !define matrix dimension from user input.
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !allocate vectors and matrix
    allocate(At(matrix_dimension, matrix_dimension))
    allocate(d(matrix_dimension), e(matrix_dimension))
    allocate(yt_min(matrix_dimension), yt_max(matrix_dimension))
    allocate(y_min(matrix_dimension), y_max(matrix_dimension))
    allocate(betas(matrix_dimension), u1_storage(matrix_dimension))

    write(*, '(/A, I0, A)') "-------> N = ", matrix_dimension, " <-------"


    !initialize At
    call initialize_A(At, matrix_dimension)

    !tridiagonalize At with Householder
    call householder_reduction(At, betas, u1_storage, matrix_dimension)

    !extract the main (d) and off-diagonal (e) elements from the modified At matrix.
    call get_tridiagonal_elements(At, d, e, matrix_dimension)

    !find the eigenvalues (lambda) and eigenvectors (yt)
    call PowerMethod(yt_max, lambda_max, yt_min, lambda_min, d, e, matrix_dimension)

    !transform the eigenvectors from the tridiagonal basis (yt) back to the original basis (y)
    call backward_transformation(At, betas, u1_storage, yt_max, y_max, matrix_dimension)
    call backward_transformation(At, betas, u1_storage, yt_min, y_min, matrix_dimension)

    !save yt max and min
    open(unit=1, file="At_eingenvector-4.txt", status="replace")
    write(1, *) "Eigenvector yt max:"
    do i = 1, matrix_dimension
        write(1, '(F8.4, 1X)') yt_max(i)
    end do
    write(1, *) "Eigenvector yt min:"
    do i = 1, matrix_dimension
        write(1, '(F8.4, 1X)') yt_min(i)
    end do
    close(1)

    !save y max and min
    open(unit=2, file="A_eingenvector-4.txt", status="replace")
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
    
    !verify the final solution by checking if A*y = lambda*y holds true.
    call verify_solution(lambda_min, y_min, matrix_dimension, "smallest")
    call verify_solution(lambda_max, y_max, matrix_dimension, "biggest")

    !deallocate vectors
    deallocate(At, d, e, yt_min, yt_max, y_min, y_max, betas, u1_storage)

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

    subroutine householder_reduction(A, betas, u1_storage, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(inout) :: A(matrix_dimension, matrix_dimension)
        real(dp), intent(out) :: betas(matrix_dimension), u1_storage(matrix_dimension)
        integer :: i,j,k
        real(dp) :: sigma, mu
        real(dp), allocatable :: u(:), p(:)

        !allocate variables
        allocate(u(matrix_dimension), p(matrix_dimension))

        do i = 1, matrix_dimension - 2

            !initialize sigma, betas and u1
            sigma = norm2(A(i+1:matrix_dimension, i))
            betas(i) = 0.0_dp
            u1_storage(i) = 0.0_dp

            if (sigma > 1.0e-15_dp) then
                if (A(i + 1, i) < 0.0_dp) then
                    sigma = -sigma
                end if


                !the Householder vector that defines the reflection plane.
                u = 0.0_dp
                u(i + 1) = A(i + 1, i) + sigma
                u(i + 2:matrix_dimension) = A(i + 2:matrix_dimension, i)

        
                !the scalar coefficient for the reflection matrix P = I - beta*u*u^T
                betas(i) = 1.0_dp / (sigma * u(i + 1))

                !store the first element of u, as it's lost after the update.
                u1_storage(i) = u(i+1)


                p = matmul(A, u)
                p = betas(i) * p 
                mu = betas(i) * dot_product(p, u) / 2.0_dp
                p = p - mu * u

                do j = 1, matrix_dimension
                    do k = 1, matrix_dimension
                        A(j, k) = A(j, k) - u(j) * p(k) - p(j) * u(k)
                    end do
                end do


                !store the "tail" of the u vector in the lower part of matrix A, 
                A(i+2:matrix_dimension, i) = u(i+2:matrix_dimension)

            end if

        end do

        !deallocate vectors
        deallocate(u, p)

    end subroutine householder_reduction

    subroutine verify_transformation_min_memory(A_original_parameters, tridiagonal_matrix_stored, beta_values, u1_vector_storage, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: tridiagonal_matrix_stored(matrix_dimension, matrix_dimension)
        real(dp), intent(in) :: beta_values(matrix_dimension), u1_vector_storage(matrix_dimension)
        character(len=1), intent(in) :: A_original_parameters
        integer :: j

        !temporary variables
        real(dp), allocatable :: previous_O_column(:), current_O_column(:), next_O_column(:)
        real(dp), allocatable :: left_hand_side_column(:), right_hand_side_column(:)
        real(dp), allocatable :: basis_vector(:)
        real(dp) :: total_difference_squared, difference_norm

        !allocate vectors
        allocate(previous_O_column(matrix_dimension), current_O_column(matrix_dimension), next_O_column(matrix_dimension))
        allocate(left_hand_side_column(matrix_dimension), right_hand_side_column(matrix_dimension))
        allocate(basis_vector(matrix_dimension))

        total_difference_squared = 0.0_dp

        !pre-calculate the first and second columns of O to start the "sliding window"
        basis_vector = 0.0_dp
        basis_vector(1) = 1.0_dp
        call backward_transformation(tridiagonal_matrix_stored, beta_values, u1_vector_storage, basis_vector, current_O_column, matrix_dimension)

        if (matrix_dimension > 1) then
            basis_vector = 0.0_dp
            basis_vector(2) = 1.0_dp
            call backward_transformation(tridiagonal_matrix_stored, beta_values, u1_vector_storage, basis_vector, next_O_column, matrix_dimension)
        else
            next_O_column = 0.0_dp
        end if
        
        previous_O_column = 0.0_dp

        !loop through each column j to check if (A*O)_j == (O*At)_j
        do j = 1, matrix_dimension

            !compute the j-th column of the right-hand side (RHS = A * O)
            call compute_Ay(current_O_column, right_hand_side_column, matrix_dimension)

            !compute the j-th column of the left-hand side (LHS = O * At) since At is tridiagonal, (At)_j has at most 3 non-zero elements.
            left_hand_side_column = tridiagonal_matrix_stored(j,j) * current_O_column
            if (j > 1) then
                left_hand_side_column = left_hand_side_column + tridiagonal_matrix_stored(j-1,j) * previous_O_column
            end if
            if (j < matrix_dimension) then
                left_hand_side_column = left_hand_side_column + tridiagonal_matrix_stored(j+1,j) * next_O_column
            end if

            !accumulate the squared difference between the columns
            total_difference_squared = total_difference_squared + norm2(left_hand_side_column - right_hand_side_column)**2

            !slide the window of O columns for the next iteration
            previous_O_column = current_O_column
            current_O_column = next_O_column
            if (j < matrix_dimension - 1) then
                basis_vector = 0.0_dp
                basis_vector(j + 2) = 1.0_dp
                call backward_transformation(tridiagonal_matrix_stored, beta_values, u1_vector_storage, basis_vector, next_O_column, matrix_dimension)
            else
                next_O_column = 0.0_dp
            end if
        end do

        !sqrt value
        difference_norm = sqrt(total_difference_squared)
        
        write(*,*)
        write(*,*) "--- Transformation Verification (min memory) ||A*O - O*At|| ---"
        write(*, '(A, E12.5)') "Difference norm: ", difference_norm

        !deallocate all temporary memory
        deallocate(previous_O_column, current_O_column, next_O_column, left_hand_side_column, right_hand_side_column, basis_vector)

    end subroutine verify_transformation_min_memory

    subroutine backward_transformation(A_stored, betas, u1_storage, yt, y, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: A_stored(matrix_dimension, matrix_dimension)
        real(dp), intent(in) :: betas(matrix_dimension), u1_storage(matrix_dimension)
        real(dp), intent(in) :: yt(matrix_dimension)
        real(dp), intent(out) :: y(matrix_dimension)
        integer :: i
        real(dp) :: s
        real(dp), allocatable :: u(:)

        !allocate variables
        allocate(u(matrix_dimension))

        !start with y = yt
        y = yt
        
        !apply the transformations P_{n-2}, ..., P_1 in reverse order
        do i = matrix_dimension - 2, 1, -1

            !if beta is zero, no transformation was applied
            if (abs(betas(i)) > 1.0e-15_dp) then

                !reconstruct the Householder vector u_i
                u = 0.0_dp
                u(i+1) = u1_storage(i)
                u(i+2:matrix_dimension) = A_stored(i+2:matrix_dimension, i)

                !apply the reflection P*y = y - beta * u * (u' * y)
                s = dot_product(u(i+1:matrix_dimension), y(i+1:matrix_dimension))
                y(i+1:matrix_dimension) = y(i+1:matrix_dimension) - betas(i) * s * u(i+1:matrix_dimension)

            end if

        end do

        !deallocate vector u
        deallocate(u)

    end subroutine backward_transformation

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

    !======================================================================
    !                       POWER METHOD SUBROUTINES
    !======================================================================
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
        real(dp) :: tolerance, old_lambda, lambda, old_under_lambda, under_lambda, norm
        integer :: iter, n
        
        !define tolerance
        tolerance = 1.0e-12_dp

        !allocate vectors
        allocate(v(matrix_dimension), w(matrix_dimension))

        ! --- Max Eigenvalue ---

            !initialize v vector with random numbers
        call random_number(v)

        !normalize the initial vector
        norm = dot_product(v,v)
        do n = 1, matrix_dimension
            v(n) = v(n) / norm
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
        norm = dot_product(v,v)
        do n = 1, matrix_dimension
            v(n) = v(n) / norm
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