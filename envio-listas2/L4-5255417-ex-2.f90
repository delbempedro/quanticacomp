!------------------------------------------------------------
! File: L4-5255417-ex-2.f90
!
! Description:
!   Finds the eigenvector of the matrix discretization of the second derivative
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
program Eingenvector

    !deactivate implicit typing
    implicit none

    !define double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    !declare variables
    integer :: matrix_dimension, j, i
    real(dp), allocatable :: d(:), c_line(:), d_line(:), x(:), r(:)
    real(dp) :: a, b, c, denominator, random_variable

    !define matrix dimension
    write(*,*) "Insert matrix dimension:"
    read(*,*) matrix_dimension

    !define j randomly
    call random_seed()
    call random_number(random_variable)
    j = int(random_variable * matrix_dimension) + 1

    !allocate vectors
    allocate(d(1:matrix_dimension), c_line(1:matrix_dimension), d_line(1:matrix_dimension), x(1:matrix_dimension), r(1:matrix_dimension))

    !initialize diagonals
    a = 1.0_dp
    b = -2.0_dp
    c = 1.0_dp

    !initialize d vector
    d = 0.0_dp
    d(j) = 1.0_dp

    !foward step
    c_line(1) = c / b
    d_line(1) = d(1) / b

    !compute c' and d'
    do i = 2, matrix_dimension-1

        denominator = b - a * c_line(i-1)
        c_line(i) = c / denominator
        d_line(i) = (d(i) - a * d_line(i-1)) / denominator

    end do

    !compute d' last element
    denominator = b - a * c_line(matrix_dimension-1)
    d_line(matrix_dimension) = (d(matrix_dimension) - a * d_line(matrix_dimension-1)) / denominator

    !back substitution
    x(matrix_dimension) = d_line(matrix_dimension)
    do i = matrix_dimension-1, 1, -1
        x(i) = d_line(i) - c_line(i) * x(i+1)
    end do

    !compute residue vector
    r(1) = b*x(1) + c*x(2) - d(1) !first term
    do i = 2, matrix_dimension-1
        r(i) = a*x(i-1) + b*x(i) + c*x(i+1) - d(i)
    end do
    r(matrix_dimension) = a*x(matrix_dimension-1) + b*x(matrix_dimension) - d(matrix_dimension)

    !print results
    write(*,'(A,I0,A)') "The solution x of Ax = d (with j = ", j, ") and the residue vector are:"
    do i = 1, matrix_dimension
        write(*,'(I3,2X,F12.6,2X,F12.6)') i, clean_zero(x(i)), clean_zero(r(i))
    end do

    deallocate(d, c_line, d_line, x, r)

contains

    pure function clean_zero(x) result(y)

        !deactivate impliciting typing
        implicit none

        !define precisions
        real(selected_real_kind(15, 307)), intent(in) :: x
        real(selected_real_kind(15, 307)) :: y

        !eliminate -0
        if (abs(x) < 1.0e-12) then
            y = 0.0
        else
            y = x
        end if

    end function clean_zero

end program Eingenvector