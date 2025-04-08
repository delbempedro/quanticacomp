program secant_method

    !deactivate implicit typing
    implicit none

    !declare variables
    real*8 :: deltax, phi_deltax, f_kminus1, f_k, k, kminus1, kplus1
    real*8 :: energy_values(3)
    integer :: i, number_of_iterations, iteration
    real*8, parameter :: tol = 1.0d-6

    !user inputs
    write(*,*) "Insert the number of iterations:"
    read(*,*) number_of_iterations

    write(*,*) "Insert k:"
    read(*,*) k

    write(*,*) "Insert kminus1:"
    read(*,*) kminus1

    write(*,*) "Insert phi_deltax (non zero):"
    read(*,*) phi_deltax
    do while (phi_deltax == 0.0d0)
        write(*,*) "phi_deltax cannot be zero, input again"
        read(*,*) phi_deltax
    end do

    !defines deltax
    deltax = 1.0d0 / number_of_iterations

    do i = 1, 3

        !initializes iteration
        iteration = 0

        !initializes k and kminus1
        kminus1 = kminus1 + 1.0d0*i
        k = k + 1.0d0*i
        
        !first iteration
        f_kminus1 = final_phi(kminus1, deltax, phi_deltax, number_of_iterations)
        f_k = final_phi(k, deltax, phi_deltax, number_of_iterations)

        do
            kplus1 = k - f_k * (k - kminus1) / (f_k - f_kminus1)
            iteration = iteration + 1

            if (abs(kplus1 - k) < tol) exit

            !update k and kminus1
            kminus1 = k
            k = kplus1

            !update f_k and f_kminus1
            f_kminus1 = f_k
            f_k = final_phi(k, deltax, phi_deltax, number_of_iterations)

        end do

        !save the energy value
        energy_values(i) = kplus1**2.0d0

    end do

    !print the results
    write(*,*) "First energy level: ", energy_values(1)
    write(*,*) "Second energy level: ", energy_values(2)
    write(*,*) "Third energy level: ", energy_values(3)

contains

    function final_phi(k, deltax, phi_deltax, number_of_iterations) result(phi)

        !deactivate implicit typing
        implicit none

        !declare variables
        real*8 :: k, deltax, phi_deltax
        integer :: number_of_iterations
        real*8 :: phi, phi_xminus1, phi_x, phi_xplus1
        integer :: j

        !initialize phi
        phi_xminus1 = 0.0d0
        phi_x = phi_deltax

        !calculate phi using the finite difference method
        do j = 1, number_of_iterations
            phi_xplus1 = 2.0d0*phi_x - phi_xminus1 - (k**2.0d0)*(deltax**2.0d0)*phi_x
            phi_xminus1 = phi_x
            phi_x = phi_xplus1
        end do

        !calculate the final phi value
        phi = phi_x

    end function final_phi

end program secant_method
