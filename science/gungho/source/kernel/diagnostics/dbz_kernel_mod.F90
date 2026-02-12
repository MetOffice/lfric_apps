!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Computes the relative humidity from prognostic variables.
!>
!> @details Computes the relative humidity field, at Wtheta points from the
!>          potential temperature, the Exner pressure and the mixing ratio of
!>          water vapour.
!>
module dbz_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, &
                                         GH_FIELD, GH_REAL,   &
                                         GH_WRITE, GH_READ,   &
                                         CELL_COLUMN
  use constants_mod,              only : r_def, i_def
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type
  use physics_common_mod,         only : qsaturation

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: dbz_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                 &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: dbz_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: dbz_code

contains

subroutine dbz_code( nlayers,       &
                     dbz,           &
                     exner_at_wt,   &
                     theta,         &
                     rho_at_wt,     &
                     mr_v,          &
                     mr_cl,         &
                     mr_r,          &
                     recip_epsilon, &
                     kappa,         &
                     p_zero,        &
                     ndf_wt,        &
                     undf_wt,       &
                     map_wt         )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt

  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_def),    dimension(undf_wt), intent(inout) :: dbz
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: exner_at_wt
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: rho_at_wt
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: mr_v, mr_cl, mr_r
  real(kind=r_def),                        intent(in)    :: recip_epsilon
  real(kind=r_def),                        intent(in)    :: kappa
  real(kind=r_def),                        intent(in)    :: p_zero

  ! Internal variables
  integer(kind=i_def) :: df, min_col_index, max_col_index
  real(kind=r_def)    :: temperature, pressure, mr_sat, rho_sqrt, rel_hum
  real(kind=r_def)    :: rho_wet, vel_r, z_e

  real(kind=r_def), parameter :: vconr = 2503.23638966667_r_def
  real(kind=r_def), parameter :: normr = 25132741228.7183_r_def

  ! Find minimum and maximum DoF numberings for this column
  min_col_index = minval(map_wt)
  max_col_index = maxval(map_wt) + nlayers - 1

  ! Directly loop over all DoFs in the column
  do df = min_col_index, max_col_index

    temperature = theta(df) * exner_at_wt(df)
    pressure = p_zero * exner_at_wt(df) ** (1.0_r_def/kappa)
    ! Pressure for saturation curve is needed in mbar
    mr_sat = qsaturation(temperature, 0.01_r_def*pressure)
    rel_hum = mr_v(df) / mr_sat * ( 1.0_r_def + mr_sat * recip_epsilon ) / &
                                  ( 1.0_r_def + mr_v(df) * recip_epsilon )
    rho_sqrt = SQRT(MIN(10.0_r_def, rho_at_wt(min_col_index) / rho_at_wt(df)))
    rho_wet = rho_at_wt(df) * MAX(1.0e-12_r_def, mr_r(df) + DIM(mr_cl(df), 1.0e-3_r_def))
    vel_r = MAX(1.0e-3_r_def, vconr * rho_sqrt * EXP(0.2_r_def*LOG(rho_wet / normr)))
    z_e = 200.0_r_def * EXP(1.6_r_def * LOG(3.6e6_r_def * (rho_wet / 1.0e3_r_def) * vel_r))
    dbz(df) = 10.0_r_def * LOG10(MAX(z_e, 0.01_r_def))
  end do

end subroutine dbz_code

end module dbz_kernel_mod
