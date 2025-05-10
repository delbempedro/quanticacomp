!------------------------------------------------------------
! File: P1-5255417-ex-3.f90
!
! Description:
!   Solve the Schrodinger (with Lennard-Jones) equation using the Numerov algorithm in conjunction with the Matching Method
!
! Dependencies:
!   - None
!
! Since:
!   - 03/2025
!
! Authors:
!   - Pedro C. Delbem <pedrodelbem@usp.br>
!------------------------------------------------------------
program MatchingMethod

    !deactivate implicit typing
    implicit none

    !declare variables
    integer*8, parameter :: max_points = 100000
    integer*8 i, j, infinity, k, l, n, i_min!, number_of_energies
    real*8 deltar, k2, initial_k2, k_line2, deltak2, diff_phi, diff_dphi, dphi_plus, dphi_minus, r_min, r_0, r_plus, r_minus
    real*8, dimension(0:max_points) :: phi_plus, phi_minus

    ! define r_0
    r_0 = 2.0d0**(1.0d0/6.0d0)

    !initialize i
    i = 0

    !initialize j
    write(*,*) "Where is infinity? (Insert j)"
    read(*,*) infinity
    j = infinity

    !define delta r
    write(*,*) "Insert delta r:"
    read(*,*) deltar
    
    !initialize k2
    write(*,*) "Insert first k2:"
    read(*,*) initial_k2

    !initialize deltak2
    write(*,*) "Insert delta k2:"
    read(*,*) deltak2

    !ask for the number of energies
    !write(*,*) "Insert the number of energies:"
    !read(*,*) number_of_energies

    !define r min
    r_min = 0.5d0

    !compute first i minimum
    i_min = 0
    do while (i_min*deltar.lt.r_min)
        i_min = i_min + 1
    end do

    !open file for writing results
    open(1, file='matching-method.txt', status='replace')
    open(2, file='eigenfunctions.txt', status='replace')

    do k = 1, 2

        !initialize k2
        k2 = initial_k2

        !initialize k_line2
        k_line2 = 1.0d0 * k

        write(1,'(A,I1)') "k'² = ",k


        do n = 1, 4!number_of_energies

            write(1,'(A,I13)') "    Energy Number = ",n

            !reinitialize i and j before each numerov execution
            i = i_min
            j = infinity

            !reinitialize variables before each numerov execution
            phi_plus(i-1) = 0.0d0
            phi_plus(i) = 1.0d0
            phi_minus(j+1) = 0.0d0
            phi_minus(j) = 1.0d0
            diff_phi  = 1.0d0
            diff_dphi = 1.0d0

            call numerov(i, j, deltar, k2, k_line2, phi_plus, r_plus, phi_minus, r_minus, r_0, r_min)

            ! Compute difference of phi+ and phi- in r_0 
            diff_phi = abs(phi_plus(j) - phi_minus(i))

            ! Compute difference of dphi+ and dphi- in r_0 
            dphi_plus = (phi_plus(j+1) - phi_plus(j-1)) / (2.0d0 * deltar)
            dphi_minus = (phi_minus(i+1) - phi_minus(i-1)) / (2.0d0 * deltar)
            diff_dphi = abs(dphi_plus - dphi_minus)

            !write(*,'(A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6)') &
            !    " φ₊ = ", phi_plus(j), " | φ₋ = ", phi_minus(i), " | Δφ = ", diff_phi, &
            !    " | Δφ' = ", diff_dphi, " | k² = ", k2

            do while (diff_phi > deltar*1e-3 .and. diff_dphi > deltar*1e-3)

                !write(*,'(A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6)') &
                !" φ₊ = ", phi_plus(j), " | φ₋ = ", phi_minus(i), " | Δφ = ", diff_phi, &
                !" | Δφ' = ", diff_dphi, " | k² = ", k2

                ! Update k2
                k2 = k2 - deltak2


                !reinitialize i and j before each numerov execution
                i = i_min
                j = infinity

                !reinitialize variables before each numerov execution
                phi_plus(i-1) = 0.0d0
                phi_plus(i) = 1.0d0
                phi_minus(j+1) = 0.0d0
                phi_minus(j) = 1.0d0
                diff_phi  = 1.0d0
                diff_dphi = 1.0d0
                
                call numerov(i, j, deltar, k2, k_line2, phi_plus, r_plus, phi_minus, r_minus, r_0, r_min)

                ! Compute difference of phi+ and phi- in r_0 
                diff_phi = abs(phi_plus(j) - phi_minus(i))

                ! Compute difference of dphi+ and dphi- in r_0 
                dphi_plus = (phi_plus(j+1) - phi_plus(j-1)) / (2.0d0 * deltar)
                dphi_minus = (phi_minus(i+1) - phi_minus(i-1)) / (2.0d0 * deltar)
                diff_dphi = abs(dphi_plus - dphi_minus)

            end do

            ! Write results to file
            write(1,'(A,F12.6)') "  k² = ", k2

            k2 = k2 - deltak2

            !write eingenfunction
            l = i_min
            do while (l*deltar.lt.r_0)

                write(2,*) phi_plus(l), l*deltar
                l = l + 1

            end do
            do while (l.lt.infinity)

                write(2,*) phi_plus(l), l*deltar
                l = l + 1

            end do
            write(2,*) "----------------------------------------------"

        end do

    end do

    close(1)
    close(2)

contains

    real function f(x, k2, k_line2)

        implicit none
        real*8, intent(in) :: x
        real*8, intent(in) :: k2, k_line2

        if (x <= 0.0d0) then
            f = 1.0d30  ! valor grande para forçar decaimento
        else
            f = k_line2 * (x**(-12.0d0) - x**(-6.0d0)) - k2
        end if

    end function f

    subroutine numerov(i, j, deltar, k2, k_line2, phi_plus, r_plus, phi_minus, r_minus, r_0, r_min)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer*8 :: i, j
        real*8 :: deltar, k2, k_line2, r_min, r_0, r_plus, r_minus
        real*8, dimension(0:max_points) :: phi_plus, phi_minus

        do while (i*deltar.lt.r_0)

            !write(*,*) " φ₊ = ", phi_plus(i), " | r₊ = ", r_plus, " | i", i

            ! Advance i
            i = i + 1

            ! Compute r+
            r_plus = i*deltar

            ! Numerov for forward
            phi_plus(i+1) = (2.0d0*phi_plus(i)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_plus, k2, k_line2)) - &
                                phi_plus(i-1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_plus-deltar, k2, k_line2))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_plus+deltar, k2, k_line2))

        end do

        do while(r_0.lt.j*deltar)

            !write(*,*) " φ₋ = ", phi_minus(i), " | r₋ = ", r_minus, " | j", j

            ! Advance j
            j = j - 1

            ! Compute r-
            r_minus = j*deltar

            ! Numerov for backward
            phi_minus(j-1) = (2.0d0*phi_minus(j)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_minus, k2, k_line2)) - &
                                phi_minus(j+1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_minus+deltar, k2, k_line2))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_minus-deltar, k2, k_line2))

        end do

    end subroutine numerov

end program MatchingMethod
