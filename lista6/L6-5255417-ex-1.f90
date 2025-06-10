!------------------------------------------------------------
! File: L6-5255417-ex-1.f90
!
! Description:
!   Use the Householder method to tridagonize uma matrix, found out its largest and lowest self-value/autovetor and check the solution
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
    integer :: matrix_dimension
    real(dp), allocatable :: A(:,:), At(:,:), O(:,:), d(:), e(:)
    real(dp), allocatable :: y_min(:), y_max(:)
    real(dp) :: lambda_min, lambda_max
    integer :: min_index, max_index

    !define matrix dimension
    write(*,*) "Inset matrix dimension:"
    read(*,*) matrix_dimension

    !alocate vectors
    allocate(A(matrix_dimension, matrix_dimension))
    allocate(At(matrix_dimension, matrix_dimension))
    allocate(O(matrix_dimension, matrix_dimension))
    allocate(d(matrix_dimension), e(matrix_dimension))
    allocate(y_min(matrix_dimension), y_max(matrix_dimension))

    write(*, '(/A, I0, A)') "-------> N = ", matrix_dimension, " <-------"

    !initialize matrix and the 5 diagonals
    call initialize_A(A, matrix_dimension)
    call initialize_A(At, matrix_dimension)
    write(*,*) "A (first 5x5):"
    call print_matrix(A, matrix_dimension, min(matrix_dimension, 5))

    !tridiagonalize A with Householder savind in At
    call householder_reduction(At, O, matrix_dimension)
    write(*,*) "Tridiagonal Matrix (T) (first 5x5):"
    call print_matrix(At, matrix_dimension, min(matrix_dimension, 5))

    !verify transformation
    call verify_transformation(matrix_dimension, At, O, A)

    !extract the diagonals d (main) e e (secondary) of T
    call get_tridiagonal_elements(At, d, e, matrix_dimension)

    !compute eingenvalues (d) and eigenvectors (O) of T
    !call tqli_eigen(d, e, O, matrix_dimension)

    !found the smallest and the biggest eigenvalues and the correspondent eigenvectors
    !lambda_min = d(1)
    !min_index = 1
    !lambda_max = d(1)
    !max_index = 1
    !call find_min_max_eigenpair(d, O, lambda_min, y_min, lambda_max, y_max, matrix_dimension)
    call PowerMethod(y_max, lambda_max, y_min, lambda_min, d, e, matrix_dimension)

    write(*, '(/A, F12.8)') "smallest eigenvalues: ", lambda_min
    write(*, '(A, F12.8)')  "biggest eigenvalues: ", lambda_max
    
    !verify Ay = lambda*y for both cases
    call verify_solution(lambda_min, y_min, matrix_dimension, "smallest")
    call verify_solution(lambda_max, y_max, matrix_dimension, "biggest")

    !deallocate arrays
    deallocate(A, O, d, e, y_min, y_max)

contains

    subroutine initialize_A(A, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(out) :: A(n,n)
        integer :: i

        !initialize with zeros
        A = 0.0_dp

        !initialize the 5 diagonals
        do i = 1, n

            !main diagonal
            A(i,i) = -5.0_dp / 2.0_dp

            !secondary diagonals
            if (i <= n-1) then
                A(i, i+1) = 4.0_dp / 3.0_dp
                A(i+1, i) = 4.0_dp / 3.0_dp
            end if

            !tertiary diagonals
            if (i <= n-2) then
                A(i, i+2) = -1.0_dp / 12.0_dp
                A(i+2, i) = -1.0_dp / 12.0_dp
            end if

        end do

    end subroutine initialize_A

    subroutine householder_reduction(A, O, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out) :: O(n,n)
        integer :: i, j, k
        real(dp) :: scale, h, g, f
        real(dp), allocatable :: e(:), p(:)

        !alocate vectors
        allocate(e(n), p(n))

        !main Householder reduction loop
        do i = n, 2, -1

            j = i-1
            h = 0.0_dp
            scale = 0.0_dp

            if (j > 1) then
                do k = 1, j
                    scale = scale + abs(A(i,k))
                end do
                if (scale == 0.0_dp) then
                    e(i) = A(i,j)
                else

                    !scale the sub-column to avoid numerical overflow/underflow
                    do k = 1, j
                        A(i,k) = A(i,k)/scale
                        h = h + A(i,k)**2
                    end do
                    f = A(i,j)
                    g = -sign(sqrt(h), f)
                    e(i) = scale*g
                    h = h - f*g
                    A(i,j) = f - g
                    f = 0.0_dp

                    !compute the Householder vector
                    do k = 1, j
                        p(k) = 0.0_dp
                        do l = 1, k
                            p(k) = p(k) + A(k,l)*A(i,l)
                        end do
                        do l = k+1, j
                            p(k) = p(k) + A(l,k)*A(i,l)
                        end do
                    end do

                    !compute the Householder matrix
                    do k = 1, j
                        p(k) = p(k)/h
                        f = f + p(k)*A(i,k)
                    end do
                    f = f/(h+h)
                    do k = 1, j
                        p(k) = p(k) - f*A(i,k)
                    end do
                    do k = 1, j
                        do l = 1, k
                            A(k,l) = A(k,l) - p(k)*A(i,l) - p(l)*A(i,k)
                        end do
                    end do
                end if

            !if j = 1, there is no sub-column
            else
                e(i) = A(i,j)

            end if

        end do

        !accumulates the transformations in O

        !initialize O with identity
        O = 0.0_dp
        do i = 1, n
            O(i,i) = 1.0_dp
        end do

        !accumulates the transformations
        do i = n, 2, -1
            j = i-1
            if (e(i) /= 0.0_dp) then
                h = e(i)*A(i,j)
                do k = 1, j
                    g = 0.0_dp
                    do l = 1, j
                        g = g + O(l,k)*A(i,l)
                    end do
                    g = g/h
                    do l = 1, j
                        O(l,k) = O(l,k) - g*A(i,l)
                    end do
                end do
            end if
        end do
        
        !restores matrix A to be tridiagonal
        do i = 1, n
            do j = 1, n
                if (i /= j) then
                    A(i,j) = 0.0_dp
                end if
            end do
            A(i,i) = e(i+1)
            if (i /= 1) A(i,i-1) = e(i)
            if (i /= n) A(i,i+1) = e(i+1)
        end do

        !deallocates temporary arrays
        deallocate(e, p)

    end subroutine householder_reduction

    subroutine verify_transformation(n, A_t, O, A_original)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: A_t(n, n), O(n, n), A_original(n, n)
        real(dp) :: diff_norm, sum_sq
        real(dp), allocatable :: OT(:,:), TEMP(:,:), C(:,:), DIFF(:,:)
        integer :: i, j, k

        !allocate vectors
        allocate(OT(n,n), TEMP(n,n), C(n,n), DIFF(n,n))

        !O -> O^T.
        do i = 1, n
            do j = 1, n
                OT(j, i) = O(i, j)
            end do
        end do

        !TEMP = A_original * O.
        do i = 1, n
            do j = 1, n

                !initialize TEMP with zeros
                TEMP(i, j) = 0.0_dp

                !do the product
                do k = 1, n
                    TEMP(i, j) = TEMP(i, j) + A_original(i, k) * O(k, j)
                end do

            end do
        end do

        !C = O^T * TEMP (O^T * A_original * O).
        do i = 1, n
            do j = 1, n

                !initialize C with zeros
                C(i, j) = 0.0_dp

                !do the product
                do k = 1, n
                    C(i, j) = C(i, j) + OT(i, k) * TEMP(k, j)
                end do

            end do
        end do

        !compute difference (DIFF = A_t - C)
        do i = 1, n
            do j = 1, n
                DIFF(i, j) = A_t(i, j) - C(i, j)
            end do
        end do

        !compute squared norm difference
        sum_sq = 0.0_dp
        do i = 1, n
            do j = 1, n
                sum_sq = sum_sq + DIFF(i, j)**2
            end do
        end do

        !take root
        diff_norm = sqrt(sum_sq)

        !print result
        write(*,*)
        write(*,*) "--- Transformation Verification O^T * A * O ---"
        write(*, '(A, E12.5)') "Difference norm ||A_t - O^T*A*O||: ", diff_norm

        !deallocate vectors
        deallocate(OT, TEMP, C, DIFF)

    end subroutine verify_transformation

    subroutine get_tridiagonal_elements(T, d, e, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: T(n,n)
        real(dp), intent(out) :: d(n), e(n)
        integer :: i

        !computes the tridiagonal elements
        do i = 1, n

            !main diagonal
            d(i) = T(i,i)

            !secondary diagonal
            if (i < n) then
                e(i) = T(i, i+1)
            else
                e(i) = 0.0_dp
            end if
            
        end do

    end subroutine get_tridiagonal_elements
    
    subroutine tqli_eigen(d, e, O, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(inout) :: d(n), e(n), O(n,n)
        integer :: i, k, l, m, iter
        real(dp) :: s, r, p, g, f, dd, c, b

        do i = 2, n
            e(i-1) = e(i)
        end do
        e(n) = 0.0_dp

        do l = 1, n
            iter = 0
            do
                do m = l, n-1
                    dd = abs(d(m)) + abs(d(m+1))
                    if (abs(e(m)) + dd == dd) exit
                end do
                if (m == l) exit
                if (iter == 30) then
                    write(*,*) "Error"
                    return
                end if
                iter = iter + 1
                g = (d(l+1) - d(l)) / (2.0_dp * e(l))
                r = sqrt(g*g + 1.0_dp)
                g = d(m) - d(l) + e(l) / (g + sign(r, g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp
                do i = m-1, l, -1
                    f = s * e(i)
                    b = c * e(i)
                    r = sqrt(f*f + g*g)
                    e(i+1) = r
                    if (r == 0.0_dp) then
                        d(i+1) = d(i+1) - p
                        e(m) = 0.0_dp
                        exit
                    end if
                    s = f/r
                    c = g/r
                    g = d(i+1) - p
                    r = (d(i) - g)*s + 2.0_dp*c*b
                    p = s*r
                    d(i+1) = g + p
                    g = c*r - b
                    do k = 1, n
                        f = O(k,i+1)
                        O(k,i+1) = s*O(k,i) + c*f
                        O(k,i) = c*O(k,i) - s*f
                    end do
                end do
                if (r == 0.0_dp .and. i >= l) cycle
                d(l) = d(l) - p
                e(l) = g
                e(m) = 0.0_dp
            end do
        end do

    end subroutine tqli_eigen

    subroutine find_min_max_eigenpair(d, O, lambda_min, y_min, lambda_max, y_max, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: d(n), O(n,n)
        real(dp), intent(out) :: lambda_min, lambda_max, y_min(n), y_max(n)
        integer :: i, min_index, max_index
        
        lambda_min = d(1)
        min_index = 1
        lambda_max = d(1)
        max_index = 1
        
        do i = 2, n
            if (d(i) < lambda_min) then
                lambda_min = d(i)
                min_index = i
            end if
            if (d(i) > lambda_max) then
                lambda_max = d(i)
                max_index = i
            end if
        end do
        
        y_min = O(:, min_index)
        y_max = O(:, max_index)

    end subroutine find_min_max_eigenpair

    subroutine matvec_A(y, Ay, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: y(n)
        real(dp), intent(out) :: Ay(n)
        integer :: i
        real(dp), parameter :: d0 = -5.0_dp/2.0_dp, d1 = 4.0_dp/3.0_dp, d2 = -1.0_dp/12.0_dp

        do i = 1, n
            Ay(i) = d0 * y(i)
            if (i > 1)   Ay(i) = Ay(i) + d1 * y(i-1)
            if (i > 2)   Ay(i) = Ay(i) + d2 * y(i-2)
            if (i < n)   Ay(i) = Ay(i) + d1 * y(i+1)
            if (i < n-1) Ay(i) = Ay(i) + d2 * y(i+2)
        end do

    end subroutine matvec_A

    subroutine verify_solution(lambda, y, n, label)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda, y(n)
        character(len=*), intent(in) :: label
        real(dp) :: residual_norm
        real(dp), allocatable :: Ay(:), residual(:)

        allocate(Ay(n), residual(n))
        
        !compute Ay = A * y
        call matvec_A(y, Ay, n)

        !compute residue vector, residual = Ay - lambda*y
        do i = 1, n
            residual(i) = Ay(i) - lambda * y(i)
        end do

        !compute residual norm ||residual||
        call dot_product(residual_norm, residual, residual, n)
        residual_norm = sqrt(residual_norm)

        write(*, '(/A, A, A)') "--- Verificacao para Autovalor ", trim(label), " ---"
        write(*, '(A, F12.8)') "Autovetor (primeiros 5 elementos):", (y(i), i=1,min(n,5))
        write(*, '(A, E12.5)') "Norma do residuo ||Ay - lambda*y||: ", residual_norm

        !deallocate arrays
        deallocate(Ay, residual)

    end subroutine verify_solution

    subroutine dot_product(result, v1, v2, n)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: n
        real(dp), intent(in) :: v1(n), v2(n)
        real(dp), intent(out) :: result
        integer :: i

        result = 0.0_dp
        do i = 1, n
            result = result + v1(i) * v2(i)
        end do

    end subroutine dot_product

    subroutine print_matrix(matrix, matrix_dimension, number_of_lines_to_print)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer, intent(in) :: matrix_dimension, number_of_lines_to_print
        real(dp), intent(in) :: matrix(matrix_dimension, matrix_dimension)
        integer :: i


        do i = 1, number_of_lines_to_print
            write(*, '(10F12.6)') matrix(i, 1:number_of_lines_to_print)
        end do
        
    end subroutine print_matrix

    !!!----------------------------- POWER METHOD SUBROUTINES -----------------------------!!!

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
        end if
        call lambdav(v, aux2, lambda, matrix_dimension)

        do i = 1,matrix_dimension
            residual(i) = aux1(i) - aux2(i)
        end do

        call vTw(difference, residual, residual, matrix_dimension)
        difference = sqrt(real(difference, dp))

    end subroutine compute_difference

    subroutine PowerMethod(vector_max, lambda_max, vector_min, lambda_min, A_main_diagonal, A_other_diagonals, matrix_dimension)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer :: matrix_dimension, n, control
        real(kind=dp) :: vector_max(matrix_dimension), lambda_max, vector_min(matrix_dimension), lambda_min, A_main_diagonal, A_other_diagonals
        real(kind=dp), allocatable :: v(matrix_dimension), w(matrix_dimension), tolerance, old_lambda, lambda, old_under_lambda, under_lambda, difference

        !define tolerance
        tolerance=1.0e-10

        !A

        !initialize v vector (v = 1/sqrt(dimension) * (1 ... 1)^T)
        do n = 1,matrix_dimension

            v(n) = 1.0_dp / sqrt(real(matrix_dimension, dp))

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

        lambda_max = lambda
        vector_max = v

        !A^-1

        !initialize v vector (v = 1/sqrt(dimension) * (1 ... 1)^T)
        do n = 1,matrix_dimension

            v(n) = 1.0_dp / sqrt(real(matrix_dimension, dp))

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

        lambda_min=under_lambda**(-1.0_dp)
        do n = 1,matrix_dimension
            vector_min(n) = v(n)**(-1.0_dp)
        end do

    end subroutine PowerMethod


end program Tridiagonalization