!------------------------------------------------------------
! File: P1-5255417-ex-2.f90
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

    implicit none

    !declare variables
    integer*8 :: i, j, k, n, number_of_energies
    real*8 :: deltar, phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1
    real*8 :: r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1
    real*8 :: phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1
    real*8 :: r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1
    real*8 :: k2, k_line2, deltak2
    real*8 :: diff_phi, diff_dphi
    real*8 :: dphi_plus, dphi_minus

    !initialize i
    i = 0

    !initialize j
    write(*,*) "Where is infinity? (Insert j)"
    read(*,*) j

    !define delta r
    write(*,*) "Insert delta r:"
    read(*,*) deltar
    
    !initialize k2
    write(*,*) "Insert first k2:"
    read(*,*) k2

    !initialize deltak2
    write(*,*) "Insert delta k2:"
    read(*,*) deltak2

    !ask for the number of energies
    write(*,*) "Insert the number of energies:"
    read(*,*) number_of_energies

    do k = 1, 2

        !initialize k_line2
        k_line2 = 1.0d0 * k

        do n = 1, number_of_energies

            !initialize phi+(0) and phi+(delta r)
            phi_plus_i_minus_1 = 0.0d0
            phi_plus_i = 0.0d0

            !initialize phi-(0) and phi-(delta r)
            phi_minus_i_plus_1 = 0.0d0
            phi_minus_i = 0.0d0

            !initialize diff_phi and diff_dphi
            diff_phi = 1.0d0
            diff_dphi = 1.0d0

            !open file for writing results
            open(k, file='matching-method.txt', status='replace')

            call numerov_step(i, j, deltar, k2, k_line2, &
                                phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1, &
                                r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1, &
                                phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1, &
                                r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1)

            ! Compute differences
            diff_phi = abs(phi_plus_i - phi_minus_i)

            dphi_plus = (phi_plus_i_plus_1 - phi_plus_i_minus_1) / (2.0d0 * deltar)
            dphi_minus = (phi_minus_i_plus_1 - phi_minus_i_minus_1) / (2.0d0 * deltar)
            diff_dphi = abs(dphi_plus - dphi_minus)

            do while (diff_phi > deltar .and. diff_dphi > deltar)

                ! Update k2
                k2 = k2 + deltak2
                
                call numerov_step(i, j, deltar, k2, k_line2, &
                                phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1, &
                                r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1, &
                                phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1, &
                                r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1)

                ! Compute differences
                diff_phi = abs(phi_plus_i - phi_minus_i)

                dphi_plus = (phi_plus_i_plus_1 - phi_plus_i_minus_1) / (2.0d0 * deltar)
                dphi_minus = (phi_minus_i_plus_1 - phi_minus_i_minus_1) / (2.0d0 * deltar)
                diff_dphi = abs(dphi_plus - dphi_minus)

            end do

            ! Write results to file
            write(k, '(F12.8,1X,F12.8)') k_line2, k2

            close(k)

        end do

    end do

    

contains

    real function f(x, k2, k_line2)
        real*8, intent(in) :: x
        real*8, intent(in) :: k2, k_line2

        f = k_line2 * (x**(-12.0d0) - x**(-6.0d0)) - k2
    end function f

    subroutine numerov_step(i, j, deltar, k2, k_line2, &
                            phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1, &
                            r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1, &
                            phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1, &
                            r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1)

        implicit none

        ! Inputs/Outputs
        integer*8, intent(inout) :: i, j
        real*8, intent(in) :: deltar, k2, k_line2
        real*8, intent(inout) :: phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1
        real*8, intent(inout) :: r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1
        real*8, intent(inout) :: phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1
        real*8, intent(inout) :: r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1

        ! Advance i and j
        i = i + 1
        j = j - 1

        ! Compute positions
        r_plus_i_minus_1 = (i-1)*deltar
        r_plus_i = i*deltar
        r_plus_i_plus_1 = (i+1)*deltar

        r_minus_i_minus_1 = (j-1)*deltar
        r_minus_i = j*deltar
        r_minus_i_plus_1 = (j+1)*deltar

        ! Numerov for forward
        phi_plus_i_plus_1 = (2.0d0*phi_plus_i*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_plus_i, k2, k_line2)) - &
                             phi_plus_i_minus_1*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_plus_i_minus_1, k2, k_line2))) / &
                            (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_plus_i_plus_1, k2, k_line2))

        phi_plus_i_minus_1 = phi_plus_i
        phi_plus_i = phi_plus_i_plus_1

        ! Numerov for backward
        phi_minus_i_minus_1 = (2.0d0*phi_minus_i*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_minus_i, k2, k_line2)) - &
                               phi_minus_i_plus_1*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_minus_i_plus_1, k2, k_line2))) / &
                              (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_minus_i_minus_1, k2, k_line2))

        phi_minus_i_plus_1 = phi_minus_i
        phi_minus_i = phi_minus_i_minus_1

        ! Save outputs (optionally)
        write(n,'(F12.8,1X,F12.8)') r_plus_i, phi_plus_i
        write(n,'(F12.8,1X,F12.8)') r_minus_i, phi_minus_i

    end subroutine numerov_step

end program MatchingMethod
