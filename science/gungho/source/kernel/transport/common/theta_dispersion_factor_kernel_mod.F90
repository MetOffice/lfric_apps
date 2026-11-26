!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates the "theta factor" used in the theta dispersion correction
!> @details Calculates the "theta_factor" used in to correct the horizontal
!!          fluxes used in consistent theta transport, to capture the correct
!!          dispersion relation. This factor is
!!              0.25 * dz * dtheta/dz
!!          The kernel calculates the gradient of a Wtheta field at its native
!!          points, either by fitting a quadratic through neighbouring points
!!          or, if specified, exp(a0 + a1*z + a2*z**2).
!!          Only implemented for the lowest-order elements.

module theta_dispersion_factor_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,     &
                                    CELL_COLUMN, GH_REAL, GH_WRITE,  &
                                    GH_SCALAR, GH_LOGICAL
use fs_continuity_mod,       only :  Wtheta
use constants_mod,           only : r_tran, r_def, i_def, EPS_R_TRAN, l_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: theta_dispersion_factor_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                     &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta), &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)           &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: theta_dispersion_factor_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public theta_dispersion_factor_code

contains

!> @brief Calculates the "theta factor" used in the theta dispersion correction
!> @param[in]     nlayers      Number of layers in the shifted mesh
!> @param[in,out] dtheta       The vertical gradient factor of a Wtheta field
!!                             at its own points
!> @param[in]     theta        The transported Wtheta variable, in Wtheta
!!                             on the prime mesh
!> @param[in]     height_wt    Heights of Wtheta points on the prime mesh
!> @param[in]     logspace     Whether to do interpolation of log(theta)
!> @param[in]     ndf_wt       Num of DoFs per cell for Wtheta
!> @param[in]     undf_wt      Num of DoFs per partition for Wtheta
!> @param[in]     map_wt       Base cell DoF-map for Wtheta
subroutine theta_dispersion_factor_code( nlayers,    &
                                         dtheta,     &
                                         theta,      &
                                         height_wt,  &
                                         logspace,   &
                                         ndf_wt,     &
                                         undf_wt,    &
                                         map_wt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_wt, ndf_wt
  integer(kind=i_def), intent(in)    :: map_wt(ndf_wt)
  logical(kind=l_def), intent(in)    :: logspace
  real(kind=r_tran),   intent(inout) :: dtheta(undf_wt)
  real(kind=r_tran),   intent(in)    :: theta(undf_wt)
  real(kind=r_def),    intent(in)    :: height_wt(undf_wt)

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)   :: z_k, z_km1, z_kp1, dz, z_w2, a1, a2
  real(kind=r_tran)   :: theta_km1, theta_k, theta_kp1, theta_factor

  ! -------------------------------------------------------------------------- !
  ! Compute dz/4 * dtheta/dz
  ! -------------------------------------------------------------------------- !

  ! Set bottom value: 0.25* dz / dtheta/dz. dz factor will cancel
  ! Unclear whether to assume dz is that of the shifted cell (so 0.5*dz)
  ! Just use linear fit at bottom, as theta should not have large values here
  dtheta(map_wt(1)) = 0.25_r_tran*(theta(map_wt(1)+1) - theta(map_wt(1)))

  ! Loop through internal layers
  do k = 1, nlayers - 1
    ! Simplify code by extracting heights
    z_k = real(height_wt(map_wt(1)+k), r_tran)
    z_kp1 = real(height_wt(map_wt(1)+k+1), r_tran)
    z_km1 = real(height_wt(map_wt(1)+k-1), r_tran)
    theta_kp1 = theta(map_wt(1)+k+1)
    theta_k = theta(map_wt(1)+k)
    theta_km1 = theta(map_wt(1)+k-1)
    ! dz needs to correspond to dz of layer on shifted mesh
    dz = 0.5_r_tran*(z_kp1 - z_km1)

    if (.not. logspace .or. MIN(theta_kp1, theta_k, theta_km1) < EPS_R_TRAN) then
      ! ---- Method for fitting a quadratic through neighbouring points ------ !
      ! The point is the centre of the shifted mesh
      ! ---------------------------------------------------------------------- !
      z_w2 = 0.25_r_tran*(2.0_r_tran*z_k + z_kp1 + z_km1)
      ! Fit quadratic in z to neighbouring theta values
      theta_factor = 0.25_r_tran * dz *                                        &
        ( theta_km1 * (2.0_r_tran*z_w2 - z_k - z_kp1)                          &
                       / ((z_km1 - z_k) * (z_km1 - z_kp1))                     &
        + theta_k * (2.0_r_tran*z_w2 - z_km1 - z_kp1)                          &
                     / ((z_k - z_km1) * (z_k - z_kp1))                         &
        + theta_kp1 * (2.0_r_tran*z_w2 - z_k - z_km1)                          &
                       / ((z_kp1 - z_k) * (z_kp1 - z_km1)) )

    else
      ! ---- Method for fitting to an exponential of a quadratic ------------- !
      ! The point is the theta point on the original mesh
      ! ---------------------------------------------------------------------- !
      ! Fit quadratic in z to neighbouring theta values
      ! Expand theta = exp(a0 + a1*z + a2*z**2)
      a1 =                                                                     &
      ((log(theta_kp1)-log(theta_k))*(z_k**2.0_r_tran - z_km1**2.0_r_tran)     &
       - (log(theta_k)-log(theta_km1))*(z_kp1**2.0_r_tran - z_k**2.0_r_tran))  &
          / ((z_kp1 - z_k)*(z_k - z_km1)*(z_km1 - z_kp1))
      a2 = - ((log(theta_kp1) - log(theta_k))*(z_k - z_km1)                    &
              - (log(theta_k) - log(theta_km1))*(z_kp1 - z_k))                 &
              / ((z_kp1 - z_k)*(z_k - z_km1)*(z_km1 - z_kp1))

      theta_factor = 0.25_r_tran * dz * (a1 + 2.0_r_tran*a2*z_k) * theta_k
    end if

    dtheta(map_wt(1)+k) = theta_factor
  end do

  ! Set top value
  k = nlayers
  z_k = real(height_wt(map_wt(1)+k), r_tran)
  z_km1 = real(height_wt(map_wt(1)+k-1), r_tran)
  theta_k = theta(map_wt(1)+k)
  theta_km1 = theta(map_wt(1)+k-1)
  ! dz needs to correspond to dz of layer on shifted mesh, but unclear if this
  ! should take into account the half-level
  dz = z_k - z_km1
  if (.not. logspace .or. MIN(theta_k, theta_km1) < EPS_R_TRAN) then
    ! Gradient of a linear field
    dtheta(map_wt(1)+k) = 0.25_r_tran*(theta_k - theta_km1)
  else
    ! Fit in log space
    a1 = (log(theta_k) - log(theta_km1)) / (z_k - z_km1)
    dtheta(map_wt(1)+k) = 0.25_r_tran*dz*a1*theta_k
  end if

end subroutine theta_dispersion_factor_code

end module theta_dispersion_factor_kernel_mod