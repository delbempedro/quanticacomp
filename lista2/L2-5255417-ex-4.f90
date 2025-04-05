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
    real deltax, deltak, phi_deltax, dphi_0, k, phi_xminus1, phi_x, phi_xplus1, x
    integer i, number_of_iterations

    !define constants
    real phi_0

    !define phi(0)
    phi_0 = 0.0

    !initialize x
    x = 0.0

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
    deltax = 1.0/number_of_iterations
    write(*,*) "deltax: ", deltax

    !initialize phi
    phi_x = 1.0
    
    !update k until phi(1) >= deltak
    do while (phi_x >= deltak)

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

    write(*,*) "First energy level: ", k

    !update k
    k = k + 2.0*deltak

    !second level    
    do while (phi_x >= deltak) !update k until phi(1) >= deltak

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

    write(*,*) "Second energy level: ", k

    !update k
    k = k + 2.0*deltak

    !third level
    do while (phi_x >= deltak) !update k until phi(1) >= deltak

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

    write(*,*) "Third energy level: ", k

    write(*,*) deltax*number_of_iterations

contains

    subroutine update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations, x)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(inout) :: phi_xminus1, phi_x, phi_xplus1, k, deltax, x
        integer, intent(in) :: number_of_iterations

        do i = 2, number_of_iterations
            phi_xminus1 = phi_x
            phi_x = phi_xplus1
            phi_xplus1 = 2.0*phi_x - phi_xminus1 - (k**2.0)*(deltax**2.0)*phi_x
            x = x + deltax
            write(*,*) "x: ", x, "i", i
        end do

    end subroutine update_phi

end program shooting_method