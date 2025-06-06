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
    real(dp), allocatable :: A(:,:), O(:,:), d(:), e(:)
    real(dp), allocatable :: y_min(:), y_max(:)
    real(dp) :: lambda_min, lambda_max
    integer :: min_idx, max_idx

    !define matrix dimension
    write(*,*) "Inset matrix dimension:"
    read(*,*) matrix_dimension

    !alocate vectors
    allocate(A(matrix_dimension, matrix_dimension))
    allocate(O(matrix_dimension, matrix_dimension))
    allocate(d(matrix_dimension), e(matrix_dimension))
    allocate(y_min(matrix_dimension), y_max(matrix_dimension))

    write(*, '(/A, I0, A)') "-------> N = ", matrix_dimension, " <-------"

    !initialize matrix and the 5 diagonals
    call initialize_A(A, matrix_dimension)
    write(*,*) "A (first 5x5):"
    call print_matrix(A, matrix_dimension, min(matrix_dimension, 5))

    !tridiagonalize A with Householder
    call householder_reduction(A, O, matrix_dimension)
    write(*,*) "Tridiagonal Matrix (T) (first 5x5):"
    call print_matrix(A, matrix_dimension, min(matrix_dimension, 5))

    !extract the diagonals d (main) e e (secondary) of T
    call get_tridiagonal_elements(A, d, e, matrix_dimension)

    !compute eingenvalues (d) and eigenvectors (O) of T
    call tqli_eigen(d, e, O, matrix_dimension)

    !found the smallest and the biggest eigenvalues and the correspondent eigenvectors
    lambda_min = d(1)
    min_idx = 1
    lambda_max = d(1)
    max_idx = 1
    call find_min_max_eigenpair(d, O, lambda_min, y_min, lambda_max, y_max, matrix_dimension)

    write(*, '(/A, F12.8)') "smallest eigenvalues: ", lambda_min
    write(*, '(A, F12.8)')  "biggest eigenvalues: ", lambda_max
    
    !verify Ay = lambda*y for both cases
    call verify_solution(lambda_min, y_min, matrix_dimension, "smallest")
    call verify_solution(lambda_max, y_max, matrix_dimension, "biggest")

    !deallocate arrays
    deallocate(A, O, d, e, y_min, y_max)

contains

    subroutine initialize_A(A, n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(out) :: A(n,n)
        integer :: i

        A = 0.0_dp
        do i = 1, n
            A(i,i) = -5.0_dp / 2.0_dp
            if (i <= n-1) then
                A(i, i+1) = 4.0_dp / 3.0_dp
                A(i+1, i) = 4.0_dp / 3.0_dp
            end if
            if (i <= n-2) then
                A(i, i+2) = -1.0_dp / 12.0_DP
                A(i+2, i) = -1.0_dp / 12.0_DP
            end if
        end do
    end subroutine initialize_A

    subroutine householder_reduction(A, O, n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out) :: O(n,n)
        integer :: i, j, k
        real(dp) :: scale, h, g, f
        real(dp), allocatable :: e(:), p(:)

        allocate(e(n), p(n))

        ! Realiza a redução de Householder para uma matriz simétrica A
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
                    do k = 1, j
                        p(k) = 0.0_dp
                        do l = 1, k
                            p(k) = p(k) + A(k,l)*A(i,l)
                        end do
                        do l = k+1, j
                            p(k) = p(k) + A(l,k)*A(i,l)
                        end do
                    end do
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
            else
                e(i) = A(i,j)
            end if
        end do

        ! Acumula as transformações em O
        O = 0.0_dp
        do i = 1, n
            O(i,i) = 1.0_dp
        end do
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
        
        ! Restaura a matriz A para ser tridiagonal
        do i = 1, n
            do j = 1, n
                if (i /= j) then
                    A(i,j) = 0.0_dp
                end if
            end do
            A(i,i) = e(i+1) ! Usa 'e' que virará diagonal principal
            if (i /= 1) A(i,i-1) = e(i)
            if (i /= n) A(i,i+1) = e(i+1)
        end do

        deallocate(e, p)
    end subroutine householder_reduction

    subroutine get_tridiagonal_elements(T, d, e, n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: T(n,n)
        real(dp), intent(out) :: d(n), e(n)
        integer :: i
        do i = 1, n
            d(i) = T(i,i)
            if (i < n) then
                e(i) = T(i, i+1)
            else
                e(i) = 0.0_dp
            end if
        end do
    end subroutine get_tridiagonal_elements
    
    subroutine tqli_eigen(d, e, O, n)
        implicit none
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
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: d(n), O(n,n)
        real(dp), intent(out) :: lambda_min, lambda_max, y_min(n), y_max(n)
        integer :: i, min_idx, max_idx
        
        lambda_min = d(1)
        min_idx = 1
        lambda_max = d(1)
        max_idx = 1
        
        do i = 2, n
            if (d(i) < lambda_min) then
                lambda_min = d(i)
                min_idx = i
            end if
            if (d(i) > lambda_max) then
                lambda_max = d(i)
                max_idx = i
            end if
        end do
        
        y_min = O(:, min_idx)
        y_max = O(:, max_idx)
    end subroutine find_min_max_eigenpair

    subroutine matvec_A(y, Ay, n)
        implicit none
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
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda, y(n)
        character(len=*), intent(in) :: label
        real(dp) :: residual_norm
        real(dp), allocatable :: Ay(:), residual(:)

        allocate(Ay(n), residual(n))
        
        ! Calcula Ay = A * y
        call matvec_A(y, Ay, n)

        ! Calcula o resíduo: residual = Ay - lambda*y
        do i = 1, n
            residual(i) = Ay(i) - lambda * y(i)
        end do

        ! Calcula a norma do resíduo: ||residual||
        call dot_product(residual_norm, residual, residual, n)
        residual_norm = sqrt(residual_norm)

        write(*, '(/A, A, A)') "--- Verificacao para Autovalor ", trim(label), " ---"
        write(*, '(A, F12.8)') "Autovetor (primeiros 5 elementos):", (y(i), i=1,min(n,5))
        write(*, '(A, E12.5)') "Norma do residuo ||Ay - lambda*y||: ", residual_norm

        deallocate(Ay, residual)
    end subroutine verify_solution

    subroutine dot_product(result, v1, v2, n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v1(n), v2(n)
        real(dp), intent(out) :: result
        integer :: i
        result = 0.0_dp
        do i = 1, n
            result = result + v1(i) * v2(i)
        end do
    end subroutine dot_product

    subroutine print_matrix(mat, n_total, n_print)
        implicit none
        integer, intent(in) :: n_total, n_print
        real(dp), intent(in) :: mat(n_total, n_total)
        integer :: i
        do i = 1, n_print
            write(*, '(10F12.6)') mat(i, 1:n_print)
        end do
    end subroutine print_matrix

end program Tridiagonalization