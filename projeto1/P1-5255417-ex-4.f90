!------------------------------------------------------------
! File: P1-5255417-ex-4.f90
!
! Description:
!   Solve the Schrodinger (with Lennard-Jones approximating by Harmonic Oscilator) equation using the Numerov algorithm in conjunction with the Matching Method
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
    integer(8), parameter :: max_points = 100000
    integer(8) i, j, infinity, e, l, n, i_min
    real(8) deltar, k, k_line, deltak, diff_phi, diff_dphi, dphi_plus, dphi_minus, r_min, r_0, r_plus, r_minus, r, profundidade
    real(8), dimension(0:max_points) :: phi_plus, phi_minus

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

    !initialize deltak
    write(*,*) "Insert delta k:"
    read(*,*) deltak

    !define r min
    r_min = 0.95d0

    !compute first i minimum
    i_min = 0
    do while (i_min*deltar.lt.r_min)
        i_min = i_min + 1
    end do

    !open file for writing results
    open(1, file='eigenfunctions1.txt', status='replace')
    open(2, file='eigenfunctions2.txt', status='replace')
    open(3, file='matching-method.txt', status='replace')

    profundidade = 1.0d0

    do e = 1, 2

        !define k_line
        k_line = profundidade * e

        !intialize k
        k = 0.0d0
        write(*,*) k

        write(3,'(A,I1)') "k' = ",e

        do n = 1, 4

            write(3,'(A,I13)') "    Energy Number = ",n

            !reinitialize i and j before each numerov execution
            i = i_min
            j = infinity

            !reinitialize variables before each numerov execution
            phi_plus(i-1) = 0.0d0
            phi_plus(i)   = 1.0d-6
            phi_minus(j+1) = 0.0d0
            phi_minus(j)   = 1.0d-6
            diff_phi  = 1.0d0
            diff_dphi = 1.0d0

            call numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

            ! Compute difference of phi+ and phi- in r_0 
            diff_phi = abs(phi_plus(i) - phi_minus(j))

            ! Compute difference of dphi+ and dphi- in r_0 
            dphi_plus = (phi_plus(i+1) - phi_plus(i-1)) / (2.0d0 * deltar)
            dphi_minus = (phi_minus(j+1) - phi_minus(j-1)) / (2.0d0 * deltar)
            diff_dphi = abs(dphi_plus - dphi_minus)

            do while (diff_phi > deltar .or. diff_dphi > deltar)

                ! Update k
                k = k + deltak

                !reinitialize i and j before each numerov execution
                i = i_min
                j = infinity

                !reinitialize variables before each numerov execution
                phi_plus(i-1) = 0.0d0
                phi_plus(i)   = 1.0d-6
                phi_minus(j+1) = 0.0d0
                phi_minus(j)   = 1.0d-6
                diff_phi  = 1.0d0
                diff_dphi = 1.0d0
                
                call numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

                ! Compute difference of phi+ and phi- in r_0 
                diff_phi = abs(phi_plus(i) - phi_minus(j))

                ! Compute difference of dphi+ and dphi- in r_0 
                dphi_plus = (phi_plus(i+1) - phi_plus(i-1)) / (2.0d0 * deltar)
                dphi_minus = (phi_minus(j+1) - phi_minus(j-1)) / (2.0d0 * deltar)
                diff_dphi = abs(dphi_plus - dphi_minus)

                write(*,*) k
                write(*,*) phi_plus(i), phi_minus(j)
                write(*,*) diff_phi, diff_dphi

            end do

            ! Write results to file
            write(3,'(A,F12.6)') "  k = ", k

            k = k + deltak

            !write eingenfunction
            l = i_min
            do while (l*deltar.lt.r_0)

                write(e,*) phi_plus(l), l*deltar
                l = l + 1

            end do
            do while (l.lt.infinity)

                write(e,*) phi_minus(l), l*deltar
                l = l + 1

            end do
            write(e,*) "----------------------------------------------"

        end do

        ! plot V(r)
        r = r_min
        do while(r.lt.infinity*deltar)

            write(e,*) V_harmonic(r, k_line, r_0), r
            r = r + deltar

        end do

    end do

    close(3)
    close(2)
    close(1)

contains

    real(8) function f(x, k, k_line)

        !deactivate implicit typing
        implicit none

        !declare variables
        real(8), intent(in) :: x
        real(8), intent(in) :: k, k_line

        f = V_harmonic(x, k_line, r_0) - k

    end function 

    real(8) function V_harmonic(x, k_line, r_0)

        !deactivate implicit typing
        implicit none

        !declare variables
        real(8), intent(in) :: x, k_line, r_0
        real(8) :: omega2

        omega2 = k_line * (156.0d0 * r_0**(-14) - 42.0d0 * r_0**(-8)) ! V''(r_0)
        V_harmonic = 0.5d0 * omega2 * (x - r_0)**2

    end function

    subroutine numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

        !deactivate implicit typing
        implicit none

        !declare variables
        integer(8) :: i, j
        real(8) :: deltar, k, k_line, r_0, r_plus, r_minus
        real(8), dimension(0:max_points) :: phi_plus, phi_minus

        do while (i*deltar.lt.r_0)

            ! Compute r+
            r_plus = i*deltar

            ! Numerov for forward
            phi_plus(i+1) = (2.0d0*phi_plus(i)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_plus, k, k_line)) - &
                                phi_plus(i-1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_plus-deltar, k, k_line))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_plus+deltar, k, k_line))

            ! Advance i
            i = i + 1

        end do

        do while(r_0.lt.j*deltar)

            ! Compute r-
            r_minus = j*deltar

            ! Numerov for backward
            phi_minus(j-1) = (2.0d0*phi_minus(j)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_minus, k, k_line)) - &
                                phi_minus(j+1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_minus+deltar, k, k_line))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_minus-deltar, k, k_line))

            ! Advance j
            j = j - 1

        end do

    end subroutine numerov

end program MatchingMethod
