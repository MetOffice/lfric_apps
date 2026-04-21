! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate contribution of dust to extinction for visibility.

! Description:
!   tbd

! Documentation:
!   tbd

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_visibility

! Code description:
!   This code is written to UMDP3 standards.
module vis_dust_mod

implicit none

! Note: The UM blends TWO dust mmr variables;
! the ground emission and the general dust tracking.
! It blends them according to the pws_dustmmr1_em ratio.

character(len=*), parameter, private :: ModuleName = 'VIS_DUST_MOD'
contains

! Based on the UM's pws_vis2_diag, with pws_dustmmr1_em set to zero.
subroutine vis_dust(                                                   &
            vis_no_precip,                                             &
            t1p5m,                                                     &
            p_star,                                                    &
            acc_ins_du,                                                &
            cor_ins_du,                                                &
            vis_with_dust,                                             &
            vis_overall)

  !--------
  ! Modules
  !--------
  use um_types, only: real_umphys
  use visbty_constants_mod, only: lnliminalcontrast, recipvisair
  use planet_constants_mod, only: planet_radius
  use yomhook, only: lhook, dr_hook
  use parkind1, only: jprb, jpim
  use constants_mod, only : r_def, r_um


  implicit none

  !-------
  ! Inputs
  !-------

  ! Todo: the calling code has the number of points hard coded to 1, but we should consider making this more general.

  ! Visibility without precip, to which we will add extinction from dust.
  real(kind=real_umphys), intent(in) :: vis_no_precip(1)

  ! 1.5m temperature
  real(kind=r_def), intent(in) :: t1p5m(1)

  ! Surface Pressure
  real(r_um), intent(in) :: p_star(1)

  ! Accumulated dust, ~= UM's dust_div1.
  real(kind=r_def), intent(in) :: acc_ins_du(1)

  ! Coarse dust, ~= UM's dust_div2.
  real(kind=r_def), intent(in) :: cor_ins_du(1)

  ! Result, visibility with dust.
  real(kind=real_umphys), intent(inout) :: vis_with_dust(1)

  ! Overall visibility.
  real(kind=real_umphys), intent(out) :: vis_overall(1)


  !-------
  ! Locals
  !-------

  ! extinction in clean air
  REAL(KIND=real_umphys) :: beta_air

  ! specific extinction coeffs at 550nm for each bin
  REAL(KIND=real_umphys) :: k_ext(2)

  ! air density
  REAL(KIND=real_umphys) :: rho

  ! extinction due to dust
  REAL(KIND=real_umphys) :: beta_dust

  ! running total of the extinction
  REAL(KIND=real_umphys) :: beta_tot



  DATA k_ext(1), k_ext(2) / 700.367, 141.453 /

  integer(kind=jpim), parameter :: zhook_in  = 0
  integer(kind=jpim), parameter :: zhook_out = 1
  real(kind=jprb)               :: zhook_handle

  character(len=*), parameter :: RoutineName='VIS_DUST'

  ! tracing
  if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! no bounds check necessary?


  ! Calculate the extinction in clean air from recipvisair
  beta_air = -lnliminalcontrast * recipvisair


  ! Calculate the extinction due to dust
  rho = p_star(1) / (t1p5m(1) * planet_radius)

  beta_dust = rho * ( acc_ins_du(1)*k_ext(1) + cor_ins_du(1)*k_ext(2) )

  ! Calculate visibility fields

  !   Invert visibility to calculate total extinction
  beta_tot = -lnliminalcontrast / vis_no_precip(1)

  !   Add the extinction from dust to the total
  beta_tot = beta_tot + beta_dust


  !   Invert back to get visibilities
  vis_overall(1)  = -lnliminalcontrast / beta_tot


  !   Include small contribution from air to limit vis model's max value
  vis_with_dust(1) = -lnliminalcontrast / (beta_dust + beta_air)


  if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return

end subroutine vis_dust
end module vis_dust_mod
