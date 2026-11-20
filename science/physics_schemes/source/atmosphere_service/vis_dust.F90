!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate contribution of dust to extinction for visibility.
!>        Based on the UM's pws_vis2_diag, with pws_dustmmr1_em set to zero.
!>
!>        Note: The UM blends TWO dust mmr variables;
!>        the ground emission and the general dust tracking.
!>        It blends them according to the pws_dustmmr1_em ratio.

module vis_dust_mod

implicit none

contains


  !> @brief Calculate aviation diagnostics.
  !> @param[in]     base_vis              The base visibility on which to include dust visibility.
  !> @param[in]     t1p5m                 1.5m temperature.
  !> @param[in]     p_star                Surface pressure.
  !> @param[in]     acc_ins_du            Accumulated dust MMR.
  !> @param[in]     cor_ins_du            Coarse dust MMR.
  !> @param[out]    vis_with_dust         Result, visibility with dust.
  !> @param[out]    vis_overall           Overall visibility, to add dust.
  subroutine vis_dust(base_vis, t1p5m, p_star, &
                      acc_ins_du, cor_ins_du, vis_with_dust, vis_overall)

    use um_types, only: real_umphys
    use visbty_constants_mod, only: lnliminalcontrast, recipvisair
    use planet_constants_mod, only: r  ! Gas constant for dry air
    use parkind1, only: jprb, jpim
    use constants_mod, only : r_def, r_um

    implicit none


    ! Arguments
    real(kind=real_umphys), intent(in) :: base_vis(1)
    real(kind=r_def), intent(in) :: t1p5m(1)
    real(r_um), intent(in) :: p_star(1)
    real(kind=r_def), intent(in) :: acc_ins_du(1)
    real(kind=r_def), intent(in) :: cor_ins_du(1)
    real(kind=real_umphys), intent(out) :: vis_with_dust(1)
    real(kind=real_umphys), intent(out) :: vis_overall(1)


    ! Locals

    ! extinction in clean air
    real(kind=real_umphys) :: beta_air

    ! specific extinction coeffs at 550nm for each bin
    real(kind=real_umphys), parameter :: k_ext(2) = [ 700.367_real_umphys, 141.453_real_umphys ]

    ! air density
    real(kind=real_umphys) :: rho

    ! extinction due to dust
    real(kind=real_umphys) :: beta_dust

    ! running total of the extinction
    real(kind=real_umphys) :: beta_tot


    ! Calculate the extinction in clean air from recipvisair.
    beta_air = -lnliminalcontrast * recipvisair

    ! Calculate the extinction due to dust.
    rho = p_star(1) / (t1p5m(1) * r)
    beta_dust = rho * ( acc_ins_du(1)*k_ext(1) + cor_ins_du(1)*k_ext(2) )


    ! Invert visibility to calculate total extinction.
    beta_tot = -lnliminalcontrast / base_vis(1)

    ! Add the extinction from dust.
    beta_tot = beta_tot + beta_dust

    ! Invert back to get visibilities.
    vis_overall(1)  = -lnliminalcontrast / beta_tot


    ! Include small contribution from air to limit vis model's max value
    vis_with_dust(1) = -lnliminalcontrast / (beta_dust + beta_air)


  end subroutine vis_dust

end module vis_dust_mod
