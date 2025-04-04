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
    real deltax, deltak, phi_deltax, dphi_0, k, phi_xminus1, phi_x, phi_xplus1
    integer i, number_of_iterations

    !define constants
    real phi_0

    !define phi(0)
    phi_0 = 0.0

    write(*,*) "Insert k:"
    read(*,*) k

    write(*,*) "Insert deltax:"
    read(*,*) deltax
    
    write(*,*) "Insert deltak:"
    read(*,*) deltak

    write(*,*) "Insert phi_deltax (non zero):"
    read(*,*) phi_deltax
    do while (phi_deltax == 0.0)
        write(*,*) "phi_deltax cannot be zero, input again"
        read(*,*) phi_deltax
    end do

    write(*,*) "Insert dphi_0 (non zero):"
    read(*,*) dphi_0
    do while (dphi_0 == 0.0)
        write(*,*) "dphi_0 cannot be zero, input again"
        read(*,*) dphi_0
    end do

    !define number of iterations
    number_of_iterations = int(1.0/deltax)

    !update phi
    phi_xminus1 = phi_0
    phi_x = phi_deltax
    phi_xplus1 = 2.0*phi_x - phi_xminus1 - k**2*deltax**2*phi_x

    !do first try
    !call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)
    
    !update k until phi(1) >= deltak
    do while (phi_xplus1 >= deltak)

        !update phi
        call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)
        write(*,*) "phi(1) = ", phi_xplus1
        !update k
        k = k + deltak

    end do

    write(*,*) "First energy level: ", k

    k = k + 2.0*deltak
    phi_xminus1 = phi_0
    phi_x = phi_deltax
    phi_xplus1 = 2.0*phi_x - phi_xminus1 - k**2*deltax**2*phi_x

    !second level    
    do while (phi_xplus1 >= deltak) !update k until phi(1) >= deltak

        !update k
        k = k + deltak

        !update phi
        call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)

    end do

    write(*,*) "Second energy level: ", k

    k = k + 2.0*deltak
    phi_xminus1 = phi_0
    phi_x = phi_deltax
    phi_xplus1 = 2.0*phi_x - phi_xminus1 - k**2*deltax**2*phi_x

    !third level
    do while (phi_xplus1 >= deltak) !update k until phi(1) >= deltak

        !update k
        k = k + deltak

        !update phi
        call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)

    end do

    write(*,*) "Third energy level: ", k

contains

    subroutine update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)

        !deactivate implicit typing
        implicit none

        !define variables
        real, intent(inout) :: phi_xminus1, phi_x, phi_xplus1, k, deltax
        integer, intent(in) :: number_of_iterations

        do i = 1, number_of_iterations
            phi_xminus1 = phi_x
            phi_x = phi_xplus1
            phi_xplus1 = 2.0*phi_x - phi_xminus1 - k**2*deltax**2*phi_x
        end do

    end subroutine update_phi

end program shooting_method