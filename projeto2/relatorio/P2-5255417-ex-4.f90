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
    real(dp) :: lambda_min, lambda_maxr.
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
    write(*,*) "Original Matrix (first 5x5):"
    call print_matrix(At, matrix_dimension, min(matrix_dimension, 5))

    !tridiagonalize At with Householder
    call householder_reduction_optimized(At, betas, u1_storage, matrix_dimension)
    write(*,*) "Tridiagonal Matrix T (first 5x5):"
    call print_matrix(At, matrix_dimension, min(matrix_dimension, 5))

    !extract the main (d) and off-diagonal (e) elements from the modified At matrix.
    call get_tridiagonal_elements(At, d, e, matrix_dimension)

    !find the eigenvalues (lambda) and eigenvectors (yt)
    call PowerMethod(yt_max, lambda_max, yt_min, lambda_min, d, e, matrix_dimension)

    !transform the eigenvectors from the tridiagonal basis (yt) back to the original basis (y)
    call backward_transformation(At, betas, u1_storage, yt_max, y_max, matrix_dimension)
    call backward_transformation(At, betas, u1_storage, yt_min, y_min, matrix_dimension)

    !save y max and min
    open(unit=2, file="A_eingenvector-b.txt", status="replace")
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
        integer :: i, j
        real(dp) :: sigma, mu
        real(dp), allocatable :: u(:), p(:), q(:)

        !allocate variables
        allocate(u(matrix_dimension), p(matrix_dimension), q(matrix_dimension))

        do i = 1, matrix_dimension - 2

            !compute sigma
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
                q = p - mu * u
                A = A - outer_product(u, q) - outer_product(q, u)


                !store the "tail" of the u vector in the lower part of matrix A, 
                A(i+2:matrix_dimension, i) = u(i+2:matrix_dimension)

            end if

        end do

        !deallocate vectors
        deallocate(u, p, q)

    end subroutine householder_reduction


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


    function outer_product(v1, v2) result(res)

        !deactivate implicit typing
        implicit none

        !declare variables
        real(dp), intent(in) :: v1(:), v2(:)
        integer :: matrix_dimension
        real(dp), allocatable :: res(:,:)

        !compute matrix dimension
        matrix_dimension = size(v1)

        !allocate vectors
        allocate(res(matrix_dimension, matrix_dimension))

        !do outer product
        res = spread(v1, dim=2, ncopies=matrix_dimension) * spread(v2, dim=1, ncopies=matrix_dimension)

    end function outer_product

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
            d(i) = At(i,i)
            if (i < matrix_dimension) then
                e(i) = At(i, i+1)
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
        
        Ay = 0.0_dp
        if (matrix_dimension <= 0) return

        ! This specialized matrix-vector product is memory-efficient as it does not 
        ! require storing the full 5-diagonal matrix. It calculates the product
        ! on-the-fly based on the known structure of the original matrix.
        do i = 1, matrix_dimension
            Ay(i) = d0 * y(i)
            if (i > 1)   Ay(i) = Ay(i) + d1 * y(i-1)
            if (i > 2)   Ay(i) = Ay(i) + d2 * y(i-2)
            if (i < matrix_dimension)     Ay(i) = Ay(i) + d1 * y(i+1)
            if (i < matrix_dimension - 1) Ay(i) = Ay(i) + d2 * y(i+2)
        end do
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

        allocate(Ay(matrix_dimension), residual(matrix_dimension))
        ! Calculate the left side of the equation: Ay = A*y
        call compute_Ay(y, Ay, matrix_dimension)
        ! Calculate the residual vector: r = Ay - lambda*y
        residual = Ay - lambda * y
        ! The norm of the residual should be close to zero if y and lambda are correct.
        residual_norm = norm2(residual)

        write(*, '(/A, A, A)') "--- Verification for Eigenvalue ", trim(label), " ---"
        write(*, '(A, F12.8)') "Eigenvalue:", lambda
        write(*, '(A, E12.5)') "Norm of residual ||Ay - lambday||: ", residual_norm
        deallocate(Ay, residual)
    end subroutine verify_solution

    subroutine print_matrix(matrix, matrix_dimension, n_print)
        implicit none
        integer, intent(in) :: matrix_dimension, n_print
        real(dp), intent(in) :: matrix(matrix_dimension, matrix_dimension)
        integer :: i
        do i = 1, n_print
            write(*, '(10F12.6)') matrix(i, 1:n_print)
        end do
    end subroutine print_matrix

    !======================================================================
    !                       POWER METHOD SUBROUTINES
    !======================================================================
    ! Calculates w = A*v for a tridiagonal matrix A.
    subroutine Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        implicit none
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: w(matrix_dimension)
        integer :: i
        
        if (matrix_dimension == 1) then
            w(1) = A_main_diagonal(1) * v(1)
            return
        end if
        
        w(1) = A_main_diagonal(1) * v(1) + A_other_diagonals(1) * v(2)
        do i = 2, matrix_dimension - 1
            w(i) = A_other_diagonals(i-1) * v(i-1) + A_main_diagonal(i) * v(i) + A_other_diagonals(i) * v(i+1)
        end do
        w(matrix_dimension) = A_other_diagonals(matrix_dimension-1) * v(matrix_dimension-1) + A_main_diagonal(matrix_dimension) * v(matrix_dimension)
    end subroutine Av

    ! Solves the system A*w = v for a tridiagonal matrix A, which is equivalent
    ! to calculating w = A_inverse * v. This uses the Thomas Algorithm.
    subroutine solve_Aw(w, v, matrix_dimension, A_main_diagonal, A_other_diagonals)
        implicit none
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: w(matrix_dimension)
        real(dp), allocatable :: c_prime(:), d_prime(:)
        integer :: i
        real(dp) :: m
        
        if (matrix_dimension == 1) then
            w(1) = v(1) / A_main_diagonal(1)
            return
        end if
        
        allocate(c_prime(matrix_dimension-1), d_prime(matrix_dimension))
        ! Forward elimination step
        c_prime(1) = A_other_diagonals(1) / A_main_diagonal(1)
        d_prime(1) = v(1) / A_main_diagonal(1)
        do i = 2, matrix_dimension - 1
            m = A_main_diagonal(i) - A_other_diagonals(i-1) * c_prime(i-1)
            c_prime(i) = A_other_diagonals(i) / m
            d_prime(i) = (v(i) - A_other_diagonals(i-1) * d_prime(i-1)) / m
        end do
        m = A_main_diagonal(matrix_dimension) - A_other_diagonals(matrix_dimension-1) * c_prime(matrix_dimension-1)
        d_prime(matrix_dimension) = (v(matrix_dimension) - A_other_diagonals(matrix_dimension-1) * d_prime(matrix_dimension-1)) / m
        
        ! Backward substitution step
        w(matrix_dimension) = d_prime(matrix_dimension)
        do i = matrix_dimension - 1, 1, -1
            w(i) = d_prime(i) - c_prime(i) * w(i+1)
        end do
        deallocate(c_prime, d_prime)
    end subroutine solve_Aw

    ! Normalizes a vector w into v.
    subroutine normalizew(v, w, matrix_dimension, tolerance)
        implicit none
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
    
    ! Computes the Rayleigh quotient lambda = (v' * A * v) / (v' * v)
    subroutine computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
        implicit none
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: v(matrix_dimension), A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: lambda
        real(dp), allocatable :: w(:)
        real(dp) :: v_dot_w, v_dot_v
        
        allocate(w(matrix_dimension))
        call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals)
        v_dot_w = dot_product(v, w)
        v_dot_v = dot_product(v, v)
        if (v_dot_v > 1.0e-20_dp) then
           lambda = v_dot_w / v_dot_v
        else
           lambda = 0.0_dp
        end if
        deallocate(w)
    end subroutine computeLambda

    ! Main routine to find the largest and smallest eigenvalues and eigenvectors
    subroutine PowerMethod(vector_max, lambda_max, vector_min, lambda_min, A_main_diagonal, A_other_diagonals, matrix_dimension)
        implicit none
        integer, intent(in) :: matrix_dimension
        real(dp), intent(in) :: A_main_diagonal(matrix_dimension), A_other_diagonals(matrix_dimension)
        real(dp), intent(out) :: vector_max(matrix_dimension), lambda_max, vector_min(matrix_dimension), lambda_min
        real(dp), allocatable :: v(:), w(:)
        real(dp) :: tolerance, old_lambda, lambda, old_under_lambda, under_lambda
        integer :: iter
        
        tolerance = 1.0e-12_dp
        allocate(v(matrix_dimension), w(matrix_dimension))

        ! --- Max Eigenvalue (Standard Power Method) ---
        call random_number(v)
        v = v / norm2(v)
        old_lambda = 0.0_dp
        iter = 0
        do while(iter < 1000)
            iter = iter + 1
            call Av(v, w, matrix_dimension, A_main_diagonal, A_other_diagonals) ! w = A*v
            call normalizew(v, w, matrix_dimension, tolerance)
            call computeLambda(lambda, v, A_main_diagonal, A_other_diagonals, matrix_dimension)
            if (abs(lambda - old_lambda) < tolerance) exit
            old_lambda = lambda
        end do
        lambda_max = lambda
        vector_max = v

        ! --- Min Eigenvalue (Inverse Power Method) ---
        call random_number(v)
        v = v / norm2(v)
        old_under_lambda = 0.0_dp
        iter = 0
        do while(iter < 1000)
           iter = iter + 1
           call solve_Aw(w, v, matrix_dimension, A_main_diagonal, A_other_diagonals) ! w = A_inverse * v
           under_lambda = dot_product(v, w) ! This is now the eigenvalue of A_inverse
           call normalizew(v, w, matrix_dimension, tolerance)
           if (iter > 1 .and. abs(1.0_dp/under_lambda - 1.0_dp/old_under_lambda) < tolerance) exit
           old_under_lambda = under_lambda
        end do
        lambda_min = 1.0_dp / under_lambda ! Eigenvalue of A is 1 / eigenvalue of A_inverse
        vector_min = v
        
        deallocate(v, w)
    end subroutine PowerMethod
    
end program Tridiagonalization