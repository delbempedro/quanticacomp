!------------------------------------------------------------
! File: P1-5255417-ex-3.f90
!
! Description:
!   Solve the radial Schrödinger equation with the Lennard-Jones potential
!   using the Numerov integration method combined with the Matching Method.
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

    ! Disable implicit typing of variables
    implicit none

    ! Define double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Declare variables
    integer :: i, j
    integer :: number_of_points, number_of_states, index_match
    real(dp) :: r_min, r_max, delta_r, delta_k, tolerance
    real(dp) :: k_line, r0, k_trial, delta_k_trial
    real(dp) :: diff_dpsi, diff_dpsi_prev
    real(dp) :: dpsi_left, dpsi_right, scaling_factor
    logical :: sign_changed

    ! Arrays for wavefunctions and eigenvalues
    real(dp), allocatable :: psi_left(:), psi_right(:), psi_total(:)
    real(dp), allocatable :: ks(:), potential(:)

    real(dp) :: r  ! Radial coordinate

    !---------------------------------------
    ! Input section: Read parameters from user
    !---------------------------------------
    print *, "Insert k_line:"
    read(*,*) k_line
    print *, "Insert r_min, r_max:"
    read(*,*) r_min, r_max
    print *, "Insert delta_k:"
    read(*,*) delta_k
    print *, "Insert delta_r:"
    read(*,*) delta_r
    print *, "Insert tolerance:"
    read(*,*) tolerance
    print *, "Insert number of states to find:"
    read(*,*) number_of_states

    ! Compute number of integration points and matching index
    number_of_points = int((r_max - r_min)/delta_r) + 1
    r0 = (2.0_dp)**(1.0_dp/6.0_dp)                       ! Equilibrium point of Lennard-Jones potential
    index_match = int((r0 - r_min)/delta_r) + 1          ! Index corresponding to r0

    ! Allocate arrays
    allocate(psi_left(0:number_of_points+1))
    allocate(psi_right(0:number_of_points+1))
    allocate(psi_total(0:number_of_points+1))
    allocate(ks(number_of_states))
    allocate(potential(0:number_of_points+1))

    ! Open output file for writing eigenfunctions
    open(unit=1, file="eigenfunctions_total.txt", status='replace', action='write')

    !---------------------------------------
    ! Begin eigenvalue search loop
    !---------------------------------------
    k_trial = -k_line         ! Initial guess for energy (negative)
    delta_k_trial = delta_k   ! Initial step in energy

    i = 1
    do while (i <= number_of_states)

        ! Compute left and right wavefunctions using current trial energy
        call compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)

        ! Scale left wavefunction to match right one at the matching point
        scaling_factor = psi_right(index_match) / psi_left(index_match)
        psi_left = psi_left * scaling_factor

        ! Compute derivatives at the matching point (central difference)
        dpsi_left  = (psi_left(index_match+1)  - psi_left(index_match-1))  / (2.0_dp * delta_r)
        dpsi_right = (psi_right(index_match+1) - psi_right(index_match-1)) / (2.0_dp * delta_r)

        diff_dpsi = dpsi_left - dpsi_right
        diff_dpsi_prev = diff_dpsi

        ! Iterative refinement of k_trial using the matching condition
        do while (abs(diff_dpsi) > tolerance)

            k_trial = k_trial + delta_k_trial
            call compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)
            scaling_factor = psi_right(index_match) / psi_left(index_match)
            psi_left = psi_left * scaling_factor

            dpsi_left  = (psi_left(index_match+1)  - psi_left(index_match-1))  / (2.0_dp * delta_r)
            dpsi_right = (psi_right(index_match+1) - psi_right(index_match-1)) / (2.0_dp * delta_r)

            diff_dpsi_prev = diff_dpsi
            diff_dpsi = dpsi_left - dpsi_right

            ! If the sign of the derivative difference changed, halve and invert the step
            sign_changed = (diff_dpsi * diff_dpsi_prev < 0.0_dp)
            if (sign_changed) delta_k_trial = -delta_k_trial / 2.0_dp

            ! Stop if the step is too small
            if (abs(delta_k_trial) < 1.0e-15_dp) exit
        end do

        ! Save eigenvalue found
        ks(i) = k_trial

        ! Combine the two wavefunctions into a single total solution
        do j = 1, index_match
            psi_total(j) = psi_left(j)
        end do
        do j = index_match+1, number_of_points
            psi_total(j) = psi_right(j)
        end do

        ! Normalize to better visualization
        call normalize_wavefunction(psi_total, number_of_points, delta_r)

        ! Write state number, energy and wavefunction to file
        r = r_min
        write(1,*) "State ", i, " Energy k = ", k_trial
        do j = 1, number_of_points
            potential(j) = V(r, k_line)
            write(1,*) r, potential(j), psi_total(j)
            r = r + delta_r
        end do
        write(1,*)

        ! Prepare next trial energy
        i = i + 1
        k_trial = k_trial + delta_k
        delta_k_trial = delta_k
    end do

    close(1)

contains

    ! Lennard-Jones potential function V(r)
    real(dp) function V(r, k_line)
        implicit none
        real(dp), intent(in) :: r, k_line
        if (r <= 1.0e-10_dp) then
            V = 1.0e10_dp  ! Prevent divergence at r = 0
        else
            V = k_line * ((1.0_dp/r)**12 - (1.0_dp/r)**6)
        end if
    end function V

    ! Subroutine that computes wavefunctions using the Numerov method
    subroutine compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)
        implicit none
        integer, intent(in) :: number_of_points
        real(dp), intent(in) :: delta_r, k_trial, r_min, k_line
        real(dp), intent(out) :: psi_left(0:number_of_points+1), psi_right(0:number_of_points+1)

        real(dp) :: r, k1, k2, k3, c
        integer :: i, j

        c = delta_r**2 / 12.0_dp  ! Coefficient for Numerov scheme

        ! Initialize left wavefunction boundary conditions
        psi_left = 0.0_dp
        psi_left(0) = -delta_r**4
        psi_left(1) = delta_r

        ! Initialize right wavefunction boundary conditions
        psi_right = 0.0_dp
        psi_right(number_of_points+1) = -delta_r**4
        psi_right(number_of_points) = delta_r

        ! Integrate left wavefunction forward (from r_min)
        r = r_min
        do i = 2, number_of_points
            k1 = k_trial - V(r, k_line)
            k2 = k_trial - V(r + delta_r, k_line)
            k3 = k_trial - V(r - delta_r, k_line)
            psi_left(i) = (2.0_dp * (1.0_dp - 5.0_dp * c * k1) * psi_left(i-1) &
                          - (1.0_dp + c * k3) * psi_left(i-2)) / (1.0_dp + c * k2)
            r = r + delta_r
        end do

        ! Integrate right wavefunction backward (from r_max)
        r = r_min + delta_r * real(number_of_points, dp)
        do j = number_of_points - 1, 1, -1
            k1 = k_trial - V(r, k_line)
            k2 = k_trial - V(r - delta_r, k_line)
            k3 = k_trial - V(r + delta_r, k_line)
            psi_right(j) = (2.0_dp * (1.0_dp - 5.0_dp * c * k1) * psi_right(j+1) &
                          - (1.0_dp + c * k3) * psi_right(j+2)) / (1.0_dp + c * k2)
            r = r - delta_r
        end do

    end subroutine compute_wavefunctions

        subroutine normalize_wavefunction(psi, npoints, delta_r)
        implicit none
        integer, intent(in) :: npoints
        real(dp), intent(in) :: delta_r
        real(dp), intent(inout) :: psi(0:npoints+1)
        integer :: i
        real(dp) :: norm

        ! Sum all values
        norm = 0.0_dp
        do i = 1, npoints
            norm = norm + psi(i)**2 * delta_r
        end do
        norm = sqrt(norm)

        ! Normalize
        if (norm > 1.0e-15_dp) then
            do i = 1, npoints
                psi(i) = psi(i) / norm
            end do
        end if
    end subroutine normalize_wavefunction

end program MatchingMethod
