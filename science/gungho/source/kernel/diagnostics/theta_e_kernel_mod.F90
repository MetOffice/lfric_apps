!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Computes theta_e, the wet equivalent potential temperature.
!>
!> @details Computes theta_e, the wet equivalent potential temperature, at
!>          Wtheta points from the potential temperature, the Exner pressure
!>          and the mixing ratio of water vapour.
!>
module theta_e_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, &
                                         GH_FIELD, GH_REAL,   &
                                         GH_WRITE, GH_READ,   &
                                         CELL_COLUMN
  use constants_mod,              only : r_def, i_def
  use driver_water_constants_mod, only : latent_heat_h2o_condensation
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: theta_e_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: theta_e_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: theta_e_code

contains

!> @brief   Computes theta_e, the wet equivalent potential temperature.
!! @param[in]     nlayers      Number of layers
!! @param[in,out] theta_e      Output theta_e field (wet equiv. pot. temperature)
!! @param[in]     theta        Input potential temperature field
!! @param[in]     exner_at_wt  Exner pressure at Wtheta points
!! @param[in]     mr_v         Mixing ratio of water vapour
!! @param[in]     cp           Heat capacity of dry air at constant pressure
!! @param[in]     ndf_wt       Number of degrees of freedom per cell for Wtheta
!! @param[in]     undf_wt      Total number of degrees of freedom for Wtheta
!! @param[in]     map_wt       Dofmap for the cell at the base of the column for Wtheta
subroutine theta_e_code( nlayers,     &
                         theta_e,     &
                         theta,       &
                         exner_at_wt, &
                         mr_v,        &
                         cp,          &
                         ndf_wt,      &
                         undf_wt,     &
                         map_wt       )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt

  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_def),    dimension(undf_wt), intent(inout) :: theta_e
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: exner_at_wt
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: mr_v
  real(kind=r_def),                        intent(in)    :: cp

  ! Internal variables
  integer(kind=i_def) :: df, min_col_index, max_col_index
  real(kind=r_def)    :: temperature

  ! Find minimum and maximum DoF numberings for this column
  min_col_index = minval(map_wt)
  max_col_index = maxval(map_wt) + nlayers - 1

  ! Directly loop over all DoFs in the column
  do df = min_col_index, max_col_index

    temperature = theta(df) * exner_at_wt(df)

    theta_e(df) = theta(df) * &
          exp( latent_heat_h2o_condensation * mr_v(df) / (cp * temperature) )

  end do

end subroutine theta_e_code

end module theta_e_kernel_mod
