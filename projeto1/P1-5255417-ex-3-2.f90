program MatchingMethod

    implicit none

    integer(8), parameter :: max_points = 100000
    integer(8) :: i, j, infinity, e, l, n, i_min
    real(8) :: deltar, k, k_line, deltak, diff_phi, diff_dphi, dphi_plus, dphi_minus, r_min, r_0, r_plus, r_minus, r, profundidade
    real(8), dimension(0:max_points) :: phi_plus, phi_minus
    real(8) :: diff_dphi_ant

    r_0 = 2.0d0**(1.0d0/6.0d0)

    write(*,*) "Where is infinity? (Insert j)"
    read(*,*) infinity
    j = infinity

    write(*,*) "Insert delta r:"
    read(*,*) deltar

    write(*,*) "Insert delta k:"
    read(*,*) deltak

    r_min = 1.05d0

    i_min = 0
    do while (i_min*deltar.lt.r_min)
        i_min = i_min + 1
    end do

    open(1, file='eigenfunctions1.txt', status='replace')
    open(2, file='eigenfunctions2.txt', status='replace')
    open(3, file='matching-method.txt', status='replace')

    profundidade = 100.0d0

    do e = 1, 2

        k_line = profundidade * e

        k = -k_line/2.0d0
        write(*,*) k

        write(3,'(A,I1)') "k' = ",e

        do n = 1, 4

            write(3,'(A,I13)') "    Energy Number = ",n

            i = i_min
            j = infinity

            phi_plus(i-1) = 0.0d0
            phi_plus(i)   = deltar
            phi_minus(j+1) = 0.0d0
            phi_minus(j)   = deltar

            call numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

            diff_phi = abs(phi_plus(i) - phi_minus(j))

            dphi_plus = (phi_plus(i+1) - phi_plus(i-1)) / (2.0d0 * deltar)
            dphi_minus = (phi_minus(j+1) - phi_minus(j-1)) / (2.0d0 * deltar)
            diff_dphi = dphi_plus - dphi_minus

            deltak = abs(deltak)
            diff_dphi_ant = diff_dphi

            !write(*,*) diff_phi, diff_dphi

            do while ( (diff_phi > deltar .or. abs(diff_dphi) > deltar) .and. k < 0.0d0 )

                !write(*,*)k

                if (diff_dphi * diff_dphi_ant < 0.0d0) then
                    deltak = -deltak / 2.0d0
                end if

                k = k + deltak

                i = i_min
                j = infinity

                phi_plus(i-1) = 0.0d0
                phi_plus(i)   = deltar
                phi_minus(j+1) = 0.0d0
                phi_minus(j)   = deltar

                call numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

                diff_phi = abs(phi_plus(i) - phi_minus(j))

                dphi_plus = (phi_plus(i+1) - phi_plus(i-1)) / (2.0d0 * deltar)
                dphi_minus = (phi_minus(j+1) - phi_minus(j-1)) / (2.0d0 * deltar)
                diff_dphi_ant = diff_dphi
                diff_dphi = dphi_plus - dphi_minus

                !write(*,*) k, diff_phi, abs(diff_dphi)

            end do

            write(3,'(A,F12.6)') "  k = ", k

            k = k + deltak

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

        r = r_min
        do while(r.lt.infinity*deltar)

            write(e,*) V(r, k_line), r
            r = r + deltar

        end do

    end do

    close(3)
    close(2)
    close(1)

contains

    real(8) function f(x, k, k_line)

        implicit none

        real(8), intent(in) :: x
        real(8), intent(in) :: k, k_line

        f = k - V(x, k_line)

    end function 

    real(8) function V(x, k_line)

        implicit none

        real(8), intent(in) :: x, k_line

        V = k_line * (1.0d0/x**(12.0d0) - 1.0d0/x**(6.0d0))

    end function

    subroutine numerov(i, j, deltar, k, k_line, phi_plus, r_plus, phi_minus, r_minus, r_0)

        implicit none

        integer(8) :: i, j
        real(8) :: deltar, k, k_line, r_0, r_plus, r_minus
        real(8), dimension(0:max_points) :: phi_plus, phi_minus

        do while (i*deltar.lt.r_0)

            i = i + 1
            r_plus = i*deltar

            phi_plus(i+1) = (2.0d0*phi_plus(i)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_plus, k, k_line)) - &
                                phi_plus(i-1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_plus-deltar, k, k_line))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_plus+deltar, k, k_line))

        end do

        do while(r_0.lt.j*deltar)

            j = j - 1
            r_minus = j*deltar

            phi_minus(j-1) = (2.0d0*phi_minus(j)*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_minus, k, k_line)) - &
                                phi_minus(j+1)*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_minus+deltar, k, k_line))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_minus-deltar, k, k_line))

        end do

    end subroutine numerov

end program MatchingMethod
