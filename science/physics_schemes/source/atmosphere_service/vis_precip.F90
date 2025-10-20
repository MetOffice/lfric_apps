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

character(len=*), parameter, private :: ModuleName = 'VIS_PRECIP_MOD'
contains
subroutine vis_precip                                                          &
           (vis_no_precip,                                                     &
                                                      !INPUT
            lca,cca,pct,                                                       &
                                                      !INPUT
            beta_ls_rain, beta_ls_snow,                                        &
                                                      !INPUT
            beta_c_rain, beta_c_snow,                                          &
                                                      !INPUT
            p_field,points,k1stpt,                                             &
                                                      !INPUT
            vis_overall,vis_lsp,vis_cp,                                        &
                                                      !OUTPUT
            i_error)                                  !OUTPUT

! Modules
use um_types, only: real_umphys
use visbty_constants_mod, only: lnliminalcontrast
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
integer, intent(in) ::                                                         &
 p_field,                                                                      &
                                      ! in NO. points in field.
 points,                                                                       &
              ! in Number of gridpoints being processed.
 k1stpt
              ! in First gridpoint processed within complete field.

real(kind=real_umphys), intent(in) ::                                          &
 vis_no_precip(p_field),                                                       &
                                      ! in Vis outside precip.
 lca(p_field),                                                                 &
                                      ! in Total Layer Cloud.
 cca(p_field),                                                                 &
                                      ! in Convective Cloud.
 beta_ls_rain(p_field),                                                        &
                                      ! in Scattering in LS Rain.
 beta_ls_snow(p_field),                                                        &
                                      ! in Scattering in LS Snow.
 beta_c_rain(p_field),                                                         &
                                      ! in Scattering in Conv Rain
 beta_c_snow(p_field)             ! in Scattering in Conv Snow

logical, intent(in) ::                                                         &
 pct                              ! in T:Cloud amounts are in %
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
real(kind=real_umphys), intent(out) ::                                         &
 vis_overall(p_field),                                                         &
                                     ! out Visibility overall
 vis_lsp(p_field),                                                             &
                                     ! out Visibility in LS Precip.
 vis_cp(p_field)                     ! out Visibility in Conv Precip.
integer, intent(out) :: i_error      ! out Error code
!----------------------------------------------------------------------
! Local variables:-----------------------------------------------------
!----------------------------------------------------------------------
integer :: i       ! Loop counters: I - horizontal field index;

real(kind=real_umphys) ::                                                      &
 beta_no_precip,                                                               &
 p_lsp(p_field),                                                               &
 p_cp(p_field)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='VIS_PRECIP'

!---------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
i_error=0
if ((k1stpt+points-1) >  p_field) then
  i_error=1
  return
end if
if (pct) then
  do i = k1stpt, k1stpt+points-1
    p_cp(i)=cca(i)/100.0
    p_lsp(i)=(1.0-p_cp(i))*lca(i)/100.0
  end do
else
  do i = k1stpt, k1stpt+points-1
    p_cp(i)=cca(i)
    p_lsp(i)=(1.0-p_cp(i))*lca(i)
  end do
end if

do i = k1stpt, k1stpt+points-1

  beta_no_precip=-lnliminalcontrast/vis_no_precip(i)

  if (p_lsp(i)  >   0.0) then
    vis_lsp(i) = -lnliminalcontrast /                                          &
      (beta_no_precip +                                                        &
       beta_ls_rain(i) + beta_ls_snow(i))
  else
    vis_lsp(i)=vis_no_precip(i)
  end if

  if (p_cp(i)  >   0.0) then
    vis_cp(i) = -lnliminalcontrast /                                           &
      (beta_no_precip +                                                        &
       beta_c_rain(i) + beta_c_snow(i))
  else
    vis_cp(i)=vis_no_precip(i)
  end if

  ! Ensure no rounding problems lead to vis > vis in clear air
  vis_overall(i) = min((1.0-p_cp(i)-p_lsp(i))*vis_no_precip(i) +               &
                   p_lsp(i)*vis_lsp(i) +p_cp(i)*vis_cp(i),                     &
                   vis_no_precip(i))

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine vis_precip
end module vis_precip_mod
