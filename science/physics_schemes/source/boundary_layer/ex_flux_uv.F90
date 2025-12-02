! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate explicit flux of momentum in u or v direction

!  Programming standard: UMDP 3

!  Documentation: UM Documentation Paper No 24.

!  subroutine EX_FLUX_UV

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module ex_flux_uv_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'EX_FLUX_UV_MOD'
contains

subroutine ex_flux_uv (                                                        &
  dimsi, dimsi_s, dimso, bl_levels,                                            &
  u_v, zhnl, rdz_u_v, rhokm_u_v, f_ngstress_uv, tau_xy_fd_uv,                  &
  tau_x_y, tau_grad,tau_non_grad                                               &
  )

use atm_fields_bounds_mod, only: array_dims
use bl_option_mod, only: max_stress_grad, ng_stress,                           &
     BrownGrant97_limited, BrownGrant97_original, zero, one
use jules_surface_mod, only: formdrag, explicit_stress
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! ARGUMENTS WITH intent in. IE: INPUT VARIABLES.

type(array_dims), intent(in) ::                                                &
   dimsi,      & ! Array dimensions for the inputs
   dimsi_s,    & ! Array dimensions for input u or v (has halos).
   dimso         ! Array dimensions for the outputs and work variables

integer, intent(in) :: bl_levels
                             ! in No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.

real(kind=r_bl), intent(in) ::                                                 &
  rdz_u_v (dimsi%i_start:dimsi%i_end,                                          &
           dimsi%j_start:dimsi%j_end, 2:bl_levels),                            &
!                                in Reciprocal of the vertical
!                                   distance from level K-1 to
!                                   level K. (K > 1) on wind levels
    rhokm_u_v (dimsi%i_start:dimsi%i_end,                                      &
               dimsi%j_start:dimsi%j_end, bl_levels),                          &
!                                in Exchange coefficients for
!                                   momentum, on UV-grid with
!                                   first and last j_end ignored.
!                                   for K>=2 (from KMKH).
    f_ngstress_uv(dimsi%i_start:dimsi%i_end,                                   &
                  dimsi%j_start:dimsi%j_end,2:bl_levels),                      &
!                                in dimensionless function for
                               !    non-gradient wind stress,
                               !    either U or V depending on call
    u_v (dimsi_s%i_start:dimsi_s%i_end,                                        &
         dimsi_s%j_start:dimsi_s%j_end,bl_levels),                             &
!                                in Westerly_Southerly component of
!                                   wind.
    tau_xy_fd_uv(dimsi%i_start:dimsi%i_end,                                    &
                 dimsi%j_start:dimsi%j_end, bl_levels),                        &
                               ! in X/Y-component of form-drag stress
                               !    at a UV point
    zhnl(dimso%i_start:dimso%i_end,                                            &
         dimso%j_start:dimso%j_end)    ! in non-local BL depth
    ! Note: this is a work variable passed in from bdy_expl3, so has
    ! the same dimensions as the outputs.

! INOUT variables
real(kind=r_bl), intent(in out) ::                                             &
  tau_x_y (dimso%i_start:dimso%i_end,                                          &
           dimso%j_start:dimso%j_end, bl_levels)
!                                out explicit x_y-component of
!                                    turbulent stress at levels
!                                    k-1/2; eg. TAUX(,1) is surface
!                                    stress. UV-grid, 1st and last j_end
!                                    set to "missing data". (N/sq m)

! ARGUMENTS WITH intent out. IE: INPUT VARIABLES CHANGED ON OUTPUT.
real(kind=r_bl), intent(out) ::                                                &
  tau_grad(dimso%i_start:dimso%i_end,                                          &
           dimso%j_start:dimso%j_end,bl_levels),                               &
!                                out k*du/dz grad stress (kg/m/s2)
    tau_non_grad(dimso%i_start:dimso%i_end,                                    &
                 dimso%j_start:dimso%j_end,bl_levels)
!                                out Non-grad stress (kg/m/s2)

integer ::                                                                     &
  i,                                                                           &
  j,                                                                           &
  k

real(kind=r_bl) ::                                                             &
  tau_surf(dimso%i_start:dimso%i_end,                                          &
           dimso%j_start:dimso%j_end),                                         &
                             ! Explicit surface stress
  sign_tau,                                                                    &
                             ! Sign of surface stress
  bl_stress_grad

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EX_FLUX_UV'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!  1.  Calculate "explicit" surface fluxes of momentum
!-----------------------------------------------------------------------

! Set diagnostics to zero in level 1

!$OMP PARALLEL DEFAULT(SHARED) private(i, j, k, sign_tau, bl_stress_grad)
if ( ng_stress == BrownGrant97_limited .or.                                    &
     ng_stress == BrownGrant97_original ) then
      ! Limit the explicitly calculated surface stress used to scale
      ! the non-gradient stress parametrization such that the implied
      ! stress gradient across the BL is less than MAX_STRESS_GRAD.
      ! This has been found by experimentation to be sufficient to
      ! stop large T increments being generated in the dynamics
      ! solver, via large increments to W

!$OMP do SCHEDULE(STATIC)
  do j = dimso%j_start, dimso%j_end
    do i = dimso%i_start, dimso%i_end
      sign_tau = sign(one, tau_x_y(i,j,1) )
      bl_stress_grad = abs( tau_x_y(i,j,1) )/zhnl(i,j)
      bl_stress_grad = min( max_stress_grad, bl_stress_grad )
      tau_surf(i,j) = sign_tau * zhnl(i,j) * bl_stress_grad
    end do
  end do
!$OMP end do
else
!$OMP do SCHEDULE(STATIC)
  do j = dimso%j_start, dimso%j_end
    do i = dimso%i_start, dimso%i_end
      tau_surf(i,j) = tau_x_y(i,j,1)
    end do
  end do
!$OMP end do
end if

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  if ( k>1 ) then
    do j = dimso%j_start, dimso%j_end
      do i = dimso%i_start, dimso%i_end

        tau_grad(i,j,k) = rhokm_u_v(i,j,k) *                                   &
                        ( u_v(i,j,k) - u_v(i,j,k-1) ) *rdz_u_v(i,j,k)
        tau_non_grad(i,j,k) = f_ngstress_uv(i,j,k) * tau_surf(i,j)
        tau_x_y(i,j,k) = tau_grad(i,j,k) + tau_non_grad(i,j,k)

        ! Add explicit orographic stress, noting that the surface stress
        ! is to be added later
        if (formdrag  ==  explicit_stress) then
          tau_x_y(i,j,k) = tau_x_y(i,j,k) + tau_xy_fd_uv(i,j,k)
        end if

      end do
    end do
  else
    do j = dimso%j_start, dimso%j_end
      do i = dimso%i_start, dimso%i_end
        tau_grad(i,j,k) = tau_x_y(i,j,1)
        tau_non_grad(i,j,k) = zero
      end do
    end do
  end if
end do
!$OMP end do
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ex_flux_uv
end module ex_flux_uv_mod
