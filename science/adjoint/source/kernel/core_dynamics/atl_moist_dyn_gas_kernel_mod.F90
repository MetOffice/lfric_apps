!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the adjoint for moisture-dependent factors
module atl_moist_dyn_gas_kernel_mod

  use argument_mod,                  only: arg_type,          &
                                           GH_FIELD, GH_REAL, &
                                           GH_SCALAR,         &
                                           GH_READWRITE,      &
                                           CELL_COLUMN
  use constants_mod,                 only: r_def, i_def
  use fs_continuity_mod,             only: Wtheta
  use kernel_mod,                    only: kernel_type
  use planet_config_mod,             only: recip_epsilon

  implicit none

  private

  !--------------------------------------------------------------------------
  ! Public types
  !--------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: atl_moist_dyn_gas_kernel_type
      private
      type(arg_type) :: meta_args(2) = (/                     &
           arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
           arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta)  &
           /)
      integer :: operates_on = CELL_COLUMN
  contains
      procedure, nopass :: atl_moist_dyn_gas_code
  end type

  !--------------------------------------------------------------------------
  ! Contained functions/subroutines
  !--------------------------------------------------------------------------
  public :: atl_moist_dyn_gas_code

  contains

  !> @brief Compute the tangent linear moist dynamical factors
  !! @param[in]      nlayers        Integer the number of layers
  !! @param[in,out]  moist_dyn_gas  Change in Gas factor (m_v / epsilon)
  !! @param[in,out]  mr_v           Change in Water vapour mixing ratio
  !! @param[in]      ndf_wtheta     The number of degrees of freedom per cell
  !!                                for wtheta
  !! @param[in]      udf_wtheta     The number of total degrees of freedom for
  !!                                wtheta
  !! @param[in]      map_wtheta     Integer array holding the dofmap for the
  !!                                cell at the base of the column
  subroutine atl_moist_dyn_gas_code( nlayers, moist_dyn_gas, mr_v, &
                                     ndf_wtheta, undf_wtheta, map_wtheta )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta

    real(kind=r_def), dimension(undf_wtheta),   intent(inout) :: moist_dyn_gas
    real(kind=r_def), dimension(undf_wtheta),   intent(inout) :: mr_v
    integer(kind=i_def), dimension(ndf_wtheta), intent(in)    :: map_wtheta

    ! Internal variables
    integer(kind=i_def)                 :: k, df

    do k = nlayers - 1, 0, -1
      do df = ndf_wtheta, 1, -1
        mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k) + moist_dyn_gas(map_wtheta(df) + k) * recip_epsilon
        moist_dyn_gas(map_wtheta(df) + k) = 0.0_r_def
      end do
    end do

  end subroutine atl_moist_dyn_gas_code

end module atl_moist_dyn_gas_kernel_mod
