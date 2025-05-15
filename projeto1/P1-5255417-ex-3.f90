program MatchingMethod

    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    integer :: i, j
    integer :: number_of_points, number_of_states, index_match
    real(dp) :: r_min, r_max, delta_r, delta_k, tolerance
    real(dp) :: k_line, r0, k_trial, delta_k_trial
    real(dp) :: diff_dpsi, diff_dpsi_prev
    real(dp) :: dpsi_left, dpsi_right, scaling_factor
    logical :: sign_changed

    real(dp), allocatable :: psi_left(:), psi_right(:), psi_total(:)
    real(dp), allocatable :: ks(:), potential(:)

    real(dp) :: r

    !--- Leitura dos parâmetros ---
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

    number_of_points = int((r_max - r_min)/delta_r) + 1
    r0 = (2.0_dp)**(1.0_dp/6.0_dp)
    index_match = int((r0 - r_min)/delta_r) + 1

    allocate(psi_left(0:number_of_points+1))
    allocate(psi_right(0:number_of_points+1))
    allocate(psi_total(0:number_of_points+1))
    allocate(ks(number_of_states))
    allocate(potential(0:number_of_points+1))

    open(unit=1, file="eigenfunctions_total.txt", status='replace', action='write')

    k_trial = -k_line
    delta_k_trial = delta_k

    i = 1
    do while (i <= number_of_states)

        call compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)

        ! Escala psi_left para casar com psi_right no ponto de matching
        scaling_factor = psi_right(index_match) / psi_left(index_match)
        psi_left = psi_left * scaling_factor

        ! Calcula derivadas no ponto de matching
        dpsi_left = (psi_left(index_match+1) - psi_left(index_match-1)) / (2.0_dp * delta_r)
        dpsi_right = (psi_right(index_match+1) - psi_right(index_match-1)) / (2.0_dp * delta_r)

        diff_dpsi = dpsi_left - dpsi_right
        diff_dpsi_prev = diff_dpsi

        do while (abs(diff_dpsi) > tolerance)
            k_trial = k_trial + delta_k_trial
            call compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)
            scaling_factor = psi_right(index_match) / psi_left(index_match)
            psi_left = psi_left * scaling_factor
            dpsi_left = (psi_left(index_match+1) - psi_left(index_match-1)) / (2.0_dp * delta_r)
            dpsi_right = (psi_right(index_match+1) - psi_right(index_match-1)) / (2.0_dp * delta_r)
            diff_dpsi_prev = diff_dpsi
            diff_dpsi = dpsi_left - dpsi_right
            sign_changed = (diff_dpsi * diff_dpsi_prev < 0.0_dp)
            if (sign_changed) delta_k_trial = -delta_k_trial / 2.0_dp
            if (abs(delta_k_trial) < 1.0e-15_dp) exit
        end do

        ks(i) = k_trial

        ! Monta a função de onda total combinando psi_left e psi_right
        do j = 1, index_match
            psi_total(j) = psi_left(j)
        end do
        do j = index_match+1, number_of_points
            psi_total(j) = psi_right(j)
        end do

        ! Normalização simples (opcional, só para melhor visualização)
        call normalize_wavefunction(psi_total, number_of_points, delta_r)

        ! Escreve resultados no arquivo: r, V(r), psi_total
        r = r_min
        write(1,*) "State ", i, " Energy k = ", k_trial
        do j = 1, number_of_points
            potential(j) = V(r, k_line)
            write(1,'(F12.6,1X,F12.6,1X,F12.6)') r, potential(j), psi_total(j)
            r = r + delta_r
        end do
        write(1,*)
        i = i + 1
        k_trial = k_trial + delta_k
        delta_k_trial = delta_k
    end do

    close(1)

contains

    real(dp) function V(r, k_line)
        implicit none
        real(dp), intent(in) :: r, k_line
        if (r <= 1.0e-10_dp) then
            V = 1.0e10_dp
        else
            V = k_line * ((1.0_dp/r)**12 - (1.0_dp/r)**6)
        end if
    end function V

    subroutine compute_wavefunctions(psi_left, psi_right, number_of_points, delta_r, k_trial, r_min, k_line)
        implicit none
        integer, intent(in) :: number_of_points
        real(dp), intent(in) :: delta_r, k_trial, r_min, k_line
        real(dp), intent(out) :: psi_left(0:number_of_points+1), psi_right(0:number_of_points+1)

        real(dp) :: r, k1, k2, k3, c
        integer :: i, j

        c = delta_r**2 / 12.0_dp

        ! Inicializa psi_left
        psi_left = 0.0_dp
        psi_left(0) = -delta_r**4
        psi_left(1) = delta_r

        ! Inicializa psi_right
        psi_right = 0.0_dp
        psi_right(number_of_points+1) = -delta_r**4
        psi_right(number_of_points) = delta_r

        ! Integra psi_left (r_min -> r_max)
        r = r_min
        do i = 2, number_of_points
            k1 = k_trial - V(r, k_line)
            k2 = k_trial - V(r + delta_r, k_line)
            k3 = k_trial - V(r - delta_r, k_line)
            psi_left(i) = (2.0_dp * (1.0_dp - 5.0_dp * c * k1) * psi_left(i-1) &
                          - (1.0_dp + c * k3) * psi_left(i-2)) / (1.0_dp + c * k2)
            r = r + delta_r
        end do

        ! Integra psi_right (r_max -> r_min)
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

        norm = 0.0_dp
        do i = 1, npoints
            norm = norm + psi(i)**2 * delta_r
        end do
        norm = sqrt(norm)

        if (norm > 1.0e-15_dp) then
            do i = 1, npoints
                psi(i) = psi(i) / norm
            end do
        end if
    end subroutine normalize_wavefunction

end program MatchingMethod
