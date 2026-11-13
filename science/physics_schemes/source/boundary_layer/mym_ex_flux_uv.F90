! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To calculate momentum fluxes in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_ex_flux_uv_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_EX_FLUX_UV_MOD'
contains

subroutine mym_ex_flux_uv(                                                     &
        dimsi, dimsi_s, dimso, bl_levels,                                      &
        rdz_u_v, rhokm_u_v, rhogamuv_uv, u_v, tau_xy_fd_uv,                    &
        tau_x_y, tau_grad, tau_count_grad)

use atm_fields_bounds_mod, only: array_dims
use jules_surface_mod, only: formdrag, explicit_stress
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Intent IN Variables
type(array_dims), intent(in) ::                                                &
   dimsi,      & ! Array dimensions for the inputs
   dimsi_s,    & ! Array dimensions for input u or v (has haloes).
   dimso         ! Array dimensions for the outputs and work variables

integer, intent(in) :: bl_levels
                 ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   rdz_u_v (dimsi%i_start:dimsi%i_end,                                         &
            dimsi%j_start:dimsi%j_end, 2:bl_levels),                           &
                 ! Reciprocal of the vertical
                 ! distance from level K-1 to
                 ! level K. (K > 1) on wind levels
   rhokm_u_v (dimsi%i_start:dimsi%i_end,                                       &
              dimsi%j_start:dimsi%j_end, bl_levels),                           &
                 ! Exchange coefficients for
                 ! momentum, on UV-grid with
                 ! first and last j_end ignored.
                 ! for K>=2, between rho level K and K-1.
                 ! i.e. assigned at theta level K-1
   rhogamuv_uv(dimsi%i_start:dimsi%i_end,                                      &
               dimsi%j_start:dimsi%j_end, 2:bl_levels),                        &
                 ! Counter Gradient Term for U or V
                 ! defined on UV-grid
   u_v(dimsi_s%i_start:dimsi_s%i_end,                                          &
       dimsi_s%j_start:dimsi_s%j_end,bl_levels),                               &
                 ! Westerly_Southerly component of wind.
   tau_xy_fd_uv(dimsi%i_start:dimsi%i_end,                                     &
                dimsi%j_start:dimsi%j_end, bl_levels)
                 ! X/Y-component of form-drag stress
                 !    at a UV point

! Intent INOUT Variables
real(kind=r_bl), intent(in out) ::                                      &
   tau_x_y (dimso%i_start:dimso%i_end,                                         &
            dimso%j_start:dimso%j_end, bl_levels)
                 ! explicit x_y-component of
                 ! turbulent stress at levels
                 ! k-1/2; eg. TAUX(,1) is surface
                 ! stress. UV-grid, 1st and last j_end
                 ! set to "missing data". (N/sq m)

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   tau_grad(dimso%i_start:dimso%i_end,                                         &
            dimso%j_start:dimso%j_end,bl_levels),                              &
                 ! k*du/dz grad stress (kg/m/s2)
   tau_count_grad(dimso%i_start:dimso%i_end,                                   &
                  dimso%j_start:dimso%j_end,bl_levels)
                 ! Counter gradient stress (kg/m/s2)

! LOCAL VARIABLES.
integer ::                                                                     &
   i, j, k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_EX_FLUX_UV'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k=1
do j = dimso%j_start, dimso%j_end
  do i = dimso%i_start, dimso%i_end
    tau_grad(i,j,k) = 0.0
    tau_count_grad(i,j,k) = 0.0
  end do
end do

do k = 2, bl_levels
  do j = dimso%j_start, dimso%j_end
    do i = dimso%i_start, dimso%i_end

      tau_grad(i,j,k) = rhokm_u_v(i,j,k) *                                     &
                     ( u_v(i,j,k) - u_v(i,j,k-1) ) *rdz_u_v(i,j,k)
      tau_count_grad(i,j,k) = rhogamuv_uv(i, j, k)
      tau_x_y(i,j,k) = tau_grad(i,j,k) + tau_count_grad(i,j,k)

    end do
  end do
end do

! Add explicit orographic stress, noting that the surface stress
! is to be added later

if (formdrag  ==  explicit_stress) then
  do k = 2, bl_levels
    do j = dimso%j_start, dimso%j_end
      do i = dimso%i_start, dimso%i_end
        tau_x_y(i,j,k) = tau_x_y(i,j,k) + tau_xy_fd_uv(i,j,k)
      end do
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_ex_flux_uv
end module mym_ex_flux_uv_mod
