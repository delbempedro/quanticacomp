!------------------------------------------------------------
! File: P1-5255417-ex-3.f90
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


contains

    subroutine Numerov

        implicit none

        ! Inputs/Outputs
        integer*8, intent(inout) :: i, j
        real*8, intent(in) :: deltar, k2, k_line2
        real*8, intent(inout) :: phi_plus_i_minus_1, phi_plus_i, phi_plus_i_plus_1
        real*8, intent(inout) :: r_plus_i_minus_1, r_plus_i, r_plus_i_plus_1
        real*8, intent(inout) :: phi_minus_i_minus_1, phi_minus_i, phi_minus_i_plus_1
        real*8, intent(inout) :: r_minus_i_minus_1, r_minus_i, r_minus_i_plus_1
        real*8, intent(in) :: r_0

        do while (i*deltar.ne.r_0)

            write(*,*) " φ₊ = ", phi_plus_i, " | r₊ = ", r_plus_i, " | i", i

            ! Advance i
            i = i + 1

            ! Compute positions
            r_plus_i_minus_1 = (i-1)*deltar
            r_plus_i = i*deltar
            r_plus_i_plus_1 = (i+1)*deltar

            ! Numerov for forward
            phi_plus_i_plus_1 = (2.0d0*phi_plus_i*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_plus_i, k2, k_line2)) - &
                                phi_plus_i_minus_1*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_plus_i_minus_1, k2, k_line2))) / &
                                (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_plus_i_plus_1, k2, k_line2))

            phi_plus_i_minus_1 = phi_plus_i
            phi_plus_i = phi_plus_i_plus_1

    end do

    do while(r_0.ne.j*deltar)

        write(*,*) " φ₋ = ", phi_minus_i, " | r₋ = ", r_minus_i, " | j", j

        ! Advance j
        j = j - 1

        ! Compute positions
        r_minus_i_minus_1 = (j-1)*deltar
        r_minus_i = j*deltar
        r_minus_i_plus_1 = (j+1)*deltar

        ! Numerov for backward
        phi_minus_i_minus_1 = (2.0d0*phi_minus_i*(1.0d0 - 5.0d0*(deltar**2.0d0)/6.0d0*f(r_minus_i, k2, k_line2)) - &
                            phi_minus_i_plus_1*(1.0d0 - (deltar**2.0d0)/12.0d0*f(r_minus_i_plus_1, k2, k_line2))) / &
                            (1.0d0 + (deltar**2.0d0)/12.0d0*f(r_minus_i_minus_1, k2, k_line2))

        phi_minus_i_plus_1 = phi_minus_i
        phi_minus_i = phi_minus_i_minus_1

    end do

    end subroutine Numerov

end program MatchingMethod