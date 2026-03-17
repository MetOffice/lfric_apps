! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate contribution of precipitation to extinction for visibility.

! Description:
!   Process fields of precipitation intensity to give scattering coefft
!   in 1/metres. This is added to an input visibility to give a total.
!   Calculated at a single model level or level within surface layer
!   e.g. screen height (1.5m)

! Documentation:
!   Forecasting Research Scientific Paper NO.4
!   Diagnosis of visibility in the UK Met Office Mesoscale Model
!   and the use of a visibility analysis to constrain initial
!   conditions.  SP Ballard, BJ Wright, BW Golding    1992
!     NIMROD diagnostic:
!   Wright, B. J., 1997: Improvements to the Nimrod Visibility
!     Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!   Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!     for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_visibility

! Code description:
!   This code is written to UMDP3 standards.
module vis_precip_mod

implicit none




! Note: The UM blends TWO dust mmr variables;
! the ground emission and the general dust tracking.
! It blends them according to the pws_dustmmr1_em ratio.


! todo: do we want to add the visibility with or without precip?


! todo: discuss these omimtted variables from vis_precip, which are always 1:
! p_field (Number of points in field)
! points (Number of gridpoints being processed)
! k1stpt (First gridpoint processed within complete field)


! Based on the UM's pws_vis2_diag, with pws_dustmmr1_em set to zero.
character(len=*), parameter, private :: ModuleName = 'VIS__MOD'
contains
subroutine vis_dust(                                                   &
            t1p5m,                                                     &
            vis_no_precip,                                             &
            acc_ins_du,                                                &
            cor_ins_du,                                                &
            vis_with_dust)

  !--------
  ! Modules
  !--------
  use um_types, only: real_umphys

  implicit none

  !-------
  ! Inputs
  !-------

  ! Todo: the calling code has the number of points hard coded to 1, but we should consider making this more general.

  ! Visibility without precip, to which we will add extinction from dust.
  real(kind=real_umphys), intent(in) :: vis_no_precip(1)

  ! 1.5m temperature
  real(kind=real_umphys), intent(in) :: t1p5m(1)

  ! Accumulated dust, ~= UM's dust_div1.
  real(kind=real_umphys), intent(in) :: acc_ins_du(1)

  ! Coarse dust, ~= UM's dust_div2.
  real(kind=real_umphys), intent(in) :: cor_ins_du(1)

  ! Result, visibility with dust.
  real(kind=real_umphys), intent(inout) :: vis_with_dust(1)
                                        
  !-------
  ! Locals
  !-------
  
  REAL(KIND=real_umphys) :: k_ext(2)
                        ! specific extinction coeffs at 550nm for each bin
  DATA k_ext(1), k_ext(2) / 700.367, 141.453 /


  ! tracing
  if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
  
  ! no bounds check necessary?
  

  ! Calculate the extinction in clean air from RecipVisAir
  beta_air = -LnLiminalContrast * RecipVisAir

  ! Calculate the extinction due to dust
  rho = pstar(1) / (t1p5m(i,j) * r)

  beta_dust(i,j)= rho * ( acc_ins_du(i,j)*k_ext(1) +                          &
                          cor_ins_du(i,j)*k_ext(2) )

  ! Calculate visibility fields

  !   Invert visibility to calculate total extinction
  beta_tot(i,j) = -LnLiminalContrast / vis_no_precip(i,j)

  !   Add the extinction from dust to the total
  beta_tot(i,j) = beta_tot(i,j) + beta_dust(i,j)



  !   Invert back to get visibilities
  pws_1p5m_vis_tot(i,j)  = -LnLiminalContrast / beta_tot(i,j)


  !   Include small contribution from air to limit vis model's max value
  pws_1p5m_vis_dust(i,j) = -LnLiminalContrast / (beta_dust(i,j) + beta_air)


  if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return

end subroutine vis_dust
end module vis_precip_mod
