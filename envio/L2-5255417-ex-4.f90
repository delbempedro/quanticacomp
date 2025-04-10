!------------------------------------------------------------
! File: L2-5255417-ex-4.f90
!
! Description:
!   Find first three energy levels of quantum 1D infity square well using shooting method
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
program shooting_method

    !deactivate implicit typing
    implicit none

    !define variables
    real*8 deltax, deltak, phi_deltax, dphi_0, k, phi_xminus1, phi_x, phi_xplus1, x, energy_values(3)
    integer i, j, number_of_iterations

    !define constants
    real*8, parameter :: phi_0 = 0.0d0

    !initialize x
    x = 0.0d0

    write(*,*) "Insert the number of iterations:"
    read(*,*) number_of_iterations

    write(*,*) "Insert k:"
    read(*,*) k

    write(*,*) "Insert deltak:"
    read(*,*) deltak

    write(*,*) "Insert phi_deltax (non zero):"
    read(*,*) phi_deltax
    do while (phi_deltax == 0.0)
        write(*,*) "phi_deltax cannot be zero, input again"
        read(*,*) phi_deltax
    end do

    !define deltax
    deltax = 1.0d0/number_of_iterations

    do i = 1,3

        !initialize phi_x and phi_xplus1
        phi_xminus1 = 0.1d0
        phi_x = 1.0d0
        phi_xplus1 = 2.0d0

        do while (phi_x*phi_xplus1 > 0.0d0) !update k until phi(1) and phi(1+deltax) have different signs

            !do the first iteration
            phi_xminus1 = phi_0
            phi_x = phi_deltax
            phi_xplus1 = 2.0*phi_x - phi_xminus1 - (k**2.0)*(deltax**2.0)*phi_x
            x = deltax

            !update phi
            call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations, x)

            !update k
            k = k + deltak

        end do

        energy_values(i) = k**2.0d0

        !update k
        k = k + 2.0d0*deltak

    end do
    
    !print the energy levels
    write(*,*) "First energy level: ", energy_values(1)
    write(*,*) "Second energy level: ", energy_values(2)
    write(*,*) "Third energy level: ", energy_values(3)

contains

    subroutine update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations, x)

        !deactivate implicit typing
        implicit none

        !define variables
        real*8, intent(inout) :: phi_xminus1, phi_x, phi_xplus1, k, deltax, x
        integer, intent(in) :: number_of_iterations

        do j = 2, number_of_iterations
            phi_xminus1 = phi_x
            phi_x = phi_xplus1
            phi_xplus1 = 2.0d0*phi_x - phi_xminus1 - (k**2.0d0)*(deltax**2.0d0)*phi_x
            x = x + deltax
        end do

    end subroutine update_phi

end program shooting_method