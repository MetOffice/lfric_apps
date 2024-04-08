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
module relative_humidity_kernel_mod

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
  type, public, extends(kernel_type) :: relative_humidity_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: relative_humidity_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: relative_humidity_code

contains

!> @brief   Computes the relative humidity from prognostic variables.
!! @param[in]     nlayers       Number of layers
!! @param[in,out] rel_hum       Output relative humidity field
!! @param[in]     mr_v          Mixing ratio of water vapour
!! @param[in]     theta         Input potential temperature field
!! @param[in]     exner_at_wt   Exner pressure at Wtheta points
!! @param[in]     recip_epsilon Ratio of gas constants for water vapour and dry
!!                              air. 1/epsilon = Rv/Rd
!! @param[in]     kappa         Power relating atmospheric pressure and Exner
!!                              pressure. kappa = Rd / cp
!! @param[in]     p_zero        Reference atmospheric pressure
!! @param[in]     ndf_wt        Number of degrees of freedom per cell for Wtheta
!! @param[in]     undf_wt       Total number of degrees of freedom for Wtheta
!! @param[in]     map_wt        Dofmap for the cell at the base of the column for Wtheta
subroutine relative_humidity_code( nlayers,       &
                                   rel_hum,       &
                                   mr_v,          &
                                   theta,         &
                                   exner_at_wt,   &
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

  real(kind=r_def),    dimension(undf_wt), intent(inout) :: rel_hum
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: exner_at_wt
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: mr_v
  real(kind=r_def),                        intent(in)    :: recip_epsilon
  real(kind=r_def),                        intent(in)    :: kappa
  real(kind=r_def),                        intent(in)    :: p_zero

  ! Internal variables
  integer(kind=i_def) :: df, min_col_index, max_col_index
  real(kind=r_def)    :: temperature, pressure, mr_sat

  ! Find minimum and maximum DoF numberings for this column
  min_col_index = minval(map_wt)
  max_col_index = maxval(map_wt) + nlayers - 1

  ! Directly loop over all DoFs in the column
  do df = min_col_index, max_col_index

    temperature = theta(df) * exner_at_wt(df)
    pressure = p_zero * exner_at_wt(df) ** (1.0_r_def/kappa)
    ! Pressure for saturation curve is needed in mbar
    mr_sat = qsaturation(temperature, 0.01_r_def*pressure)
    rel_hum(df) = mr_v(df) / mr_sat * ( 1.0_r_def + mr_sat * recip_epsilon ) / &
                                      ( 1.0_r_def + mr_v(df) * recip_epsilon )
  end do

end subroutine relative_humidity_code

end module relative_humidity_kernel_mod
