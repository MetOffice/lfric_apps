!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! TODO #3706: when r_tran operator types are supported by PSyclone, this
! kernel should use that instead of an r_def div operator
!-------------------------------------------------------------------------------

!> @brief Iterative update of a mass flux to ensure positivity of a density
!> @details Modifies a mass flux field through an iterative process, to ensure
!!          positivity of the density that will be transported using this flux,
!!          Let the scaling be s, which must be chosen so that
!!          f(n+1) = f(n) - dt*div(flux_in) - dt*div(s*flux_out) >= min_value
!!          Rearranging this for s gives:
!!          s = [f(n) - dt*div(flux_in) - min_value] / [- dt*div(flux_out)]

module iterate_min_flux_kernel_mod
use argument_mod,            only : arg_type,                           &
                                    GH_FIELD, GH_OPERATOR, GH_READ,     &
                                    CELL_COLUMN, GH_REAL, GH_READINC,   &
                                    GH_SCALAR
use fs_continuity_mod,       only : W3, W2
use constants_mod,           only : r_def, r_tran, i_def, EPS_R_TRAN
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: iterate_min_flux_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                      &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,    W3),     & ! field
       arg_type(GH_FIELD,    GH_REAL, GH_READINC, W2),     & ! flux
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,    W3, W2), & ! div
       arg_type(GH_SCALAR,   GH_REAL, GH_READ),            & ! min_value
       arg_type(GH_SCALAR,   GH_REAL, GH_READ)             & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: iterate_min_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public iterate_min_flux_code

contains

!> @brief Iterative update of a mass flux to ensure positivity of a density
!! @param[in]     cell        Horizontal cell index
!! @param[in]     nlayers     Number of layers
!! @param[in]     field       Density field before application of div(flux)
!! @param[in,out] flux        Flux field to be updated
!! @param[in]     ncell_3d    Total number of cells related to div
!! @param[in]     div         Divergence operator used in the update
!! @param[in]     min_value   The minimum value to enforce for the density field
!! @param[in]     dts         The time-step used in the update
!! @param[in]     ndf_w3      Num of DoFs per cell for the density field
!! @param[in]     undf_w3     Num of DoFs per partition for the density field
!! @param[in]     map_w3      Base cell DoF-map for the density field
!! @param[in]     ndf_w2      Num of DoFs per cell for the flux field
!! @param[in]     undf_w2     Num of DoFs per partition for the flux field
!! @param[in]     map_w2      Base cell DoF-map for the flux field
subroutine iterate_min_flux_code(cell,                    &
                                 nlayers,                 &
                                 field, flux,             &
                                 ncell_3d,                &
                                 div,                     &
                                 min_value,               &
                                 dts,                     &
                                 ndf_w3, undf_w3, map_w3, &
                                 ndf_w2, undf_w2, map_w2  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: cell, nlayers
  integer(kind=i_def), intent(in)    :: ncell_3d
  integer(kind=i_def), intent(in)    :: undf_w3, ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w2, ndf_w2
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  real(kind=r_tran),   intent(inout) :: flux(undf_w2)
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_def),    intent(in)    :: div(ndf_w3,ndf_w2,ncell_3d)
  real(kind=r_tran),   intent(in)    :: min_value
  real(kind=r_tran),   intent(in)    :: dts

  ! Internal variables
  integer(kind=i_def) :: k, ik, df_w3, df_w2
  real(kind=r_tran)   :: a, b, inc, inc_out, inc_in
  real(kind=r_tran)   :: field_new, flux_scaler
  integer(kind=i_def) :: flux_change_id(ndf_w2)

  do k = 0, nlayers-1
    ik = (cell-1)*nlayers + k + 1

    do df_w3 = 1, ndf_w3

      ! Compute positive and negative increments from the flux
      inc_out = 0.0_r_tran
      inc_in  = 0.0_r_tran
      flux_change_id = 0

      do df_w2 = 1, ndf_w2
        inc = - dts*real(div(df_w3,df_w2,ik),r_tran)*flux(map_w2(df_w2)+k)

        ! Outgoing fluxes are indicated by flux_change_id
        if ( inc < 0.0_r_tran ) then
          inc_out = inc_out - inc
          flux_change_id(df_w2) = 1
        else
          inc_in = inc_in + inc
        end if
      end do

      ! Estimate of new field value
      field_new = field(map_w3(df_w3)+k) + inc_in - inc_out

      if ( field_new < 0.0_r_tran) then
        ! Compute scaling for outgoing fluxes
        ! Scaling should be 1 if there is no danger of going below min_value
        ! s = [f(n) - dt*div(flux_in) - min_value] / [- dt*div(flux_out)]
        !   = a / inc_out

        a = field(map_w3(df_w3)+k) + inc_in - (min_value + EPS_R_TRAN)
        b = a / max(inc_out, EPS_R_TRAN)
        ! Adjust scaling to ensure it is between 0 and 1
        flux_scaler = min(max(0.0_r_tran,b),1.0_r_tran)

        ! Scale outgoing fluxes (indicated by flux_change_id)
        do df_w2 = 1, ndf_w2
          if ( flux_change_id(df_w2) == 1 ) then
            flux(map_w2(df_w2)+k) = flux(map_w2(df_w2)+k) * flux_scaler
          end if
        end do

      end if

    end do

  end do

end subroutine iterate_min_flux_code

end module iterate_min_flux_kernel_mod
