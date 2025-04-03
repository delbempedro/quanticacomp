!------------------------------------------------------------
! File: L2-5255417-ex-3.f90
!
! Description:
!   Find first three energy levels of quantum 1D infinite square well using shooting method
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

    ! deactivate implicit typing
    implicit none

    ! define variables
    real :: deltax, deltak, phi_deltax, dphi_0, k, phi_xminus1, phi_x, phi_xplus1
    integer :: i, number_of_iterations

    ! define constants
    real, parameter :: phi_0 = 0.0

    ! user input
    write(*,*) "Insert initial k:"
    read(*,*) k

    write(*,*) "Insert deltax:"
    read(*,*) deltax

    write(*,*) "Insert deltak:"
    read(*,*) deltak

    write(*,*) "Insert phi_deltax (non-zero):"
    read(*,*) phi_deltax
    do while (phi_deltax == 0.0)
        write(*,*) "phi_deltax cannot be zero, input again"
        read(*,*) phi_deltax
    end do

    write(*,*) "Insert dphi_0 (non-zero):"
    read(*,*) dphi_0
    do while (dphi_0 == 0.0)
        write(*,*) "dphi_0 cannot be zero, input again"
        read(*,*) dphi_0
    end do

    ! define number of iterations
    number_of_iterations = int(1.0/deltax)

    ! find energy levels
    call find_energy_level(k, deltak, deltax, number_of_iterations, phi_deltax, "First")
    call find_energy_level(k, deltak, deltax, number_of_iterations, phi_deltax, "Second")
    call find_energy_level(k, deltak, deltax, number_of_iterations, phi_deltax, "Third")

contains

    subroutine update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)
        implicit none
        real, intent(inout) :: phi_xminus1, phi_x, phi_xplus1
        real, intent(in) :: k, deltax
        integer, intent(in) :: number_of_iterations
        
        do i = 1, number_of_iterations
            phi_xminus1 = phi_x
            phi_x = phi_xplus1
            phi_xplus1 = 2.0 * phi_x - phi_xminus1 - k**2 * deltax**2 * phi_x
        end do
    end subroutine update_phi

    subroutine find_energy_level(k, deltak, deltax, number_of_iterations, phi_deltax, label)
        implicit none
        real, intent(inout) :: k
        real, intent(in) :: deltak, deltax, phi_deltax
        integer, intent(in) :: number_of_iterations
        character(len=*), intent(in) :: label
        real :: phi_xminus1, phi_x, phi_xplus1

        ! initialize phi values
        phi_xminus1 = 0.0
        phi_x = phi_deltax
        phi_xplus1 = 2.0 * phi_x - phi_xminus1 - k**2 * deltax**2 * phi_x

        do while (phi_xplus1 > deltak)
            call update_phi(phi_xminus1, phi_x, phi_xplus1, k, deltax, number_of_iterations)
            
            k = k + deltak
        end do

        write(*,*) trim(label), " energy level: ", k

        ! Ensure next search starts from a higher value
        k = k + 2.0 * deltak
    end subroutine find_energy_level

end program shooting_method
