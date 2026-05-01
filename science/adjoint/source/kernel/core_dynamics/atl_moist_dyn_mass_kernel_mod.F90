!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the tangent linear for moisture-dependent factors
module atl_moist_dyn_mass_kernel_mod

  use argument_mod,                  only: arg_type,          &
                                           GH_FIELD, GH_REAL, &
                                           GH_SCALAR,         &
                                           GH_READWRITE,      &
                                           CELL_COLUMN
  use constants_mod,                 only: r_def, i_def
  use fs_continuity_mod,             only: Wtheta
  use kernel_mod,                    only: kernel_type

  implicit none

  private

  !--------------------------------------------------------------------------
  ! Public types
  !--------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: atl_moist_dyn_mass_kernel_type
      private
      type(arg_type) :: meta_args(2) = (/                       &
           arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta), &
           arg_type(GH_FIELD*6, GH_REAL, GH_READWRITE, Wtheta)  &
           /)
      integer :: operates_on = CELL_COLUMN
  contains
      procedure, nopass :: atl_moist_dyn_mass_code
  end type

  !--------------------------------------------------------------------------
  ! Contained functions/subroutines
  !--------------------------------------------------------------------------
  public :: atl_moist_dyn_mass_code

  contains

  !> @brief Compute the tangent linear moist dynamical factors
  !! @param[in]      nlayers        Integer the number of layers
  !! @param[in,out]  moist_dyn_tot  Change in Total mass factor (sum m_x)
  !! @param[in,out]  mr_v           Change in Water vapour mixing ratio
  !! @param[in]      mr_cl          Change in Liquid cloud mixing ratio
  !! @param[in]      mr_r           Change in Rain mixing ratio
  !! @param[in]      mr_s           Change in Snow mixing ratio
  !! @param[in]      mr_g           Change in Graupel mixing ratio
  !! @param[in]      mr_ci          Change in Ice cloud mixing ratio
  !! @param[in]      ndf_wtheta     The number of degrees of freedom per cell
  !!                                for wtheta
  !! @param[in]      udf_wtheta     The number of total degrees of freedom for
  !!                                wtheta
  !! @param[in]      map_wtheta     Integer array holding the dofmap for the
  !!                                cell at the base of the column
  subroutine atl_moist_dyn_mass_code(nlayers, moist_dyn_tot,               &
                                     mr_v, mr_cl, mr_r, mr_s, mr_g, mr_ci, &
                                     ndf_wtheta, undf_wtheta, map_wtheta)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta

    real(kind=r_def), dimension(undf_wtheta),   intent(inout) :: moist_dyn_tot
    real(kind=r_def), dimension(undf_wtheta),   intent(in)    :: mr_v, mr_cl, mr_r
    real(kind=r_def), dimension(undf_wtheta),   intent(in)    :: mr_ci, mr_s, mr_g
    integer(kind=i_def), dimension(ndf_wtheta), intent(in)    :: map_wtheta

    ! Internal variables
    integer(kind=i_def)                 :: k, df

    do k = nlayers - 1, 0, -1
      do df = ndf_wtheta, 1, -1
        mr_g(map_wtheta(df) + k) = mr_g(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        mr_s(map_wtheta(df) + k) = mr_s(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        mr_ci(map_wtheta(df) + k) = mr_ci(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        mr_r(map_wtheta(df) + k) = mr_r(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        mr_cl(map_wtheta(df) + k) = mr_cl(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k) + moist_dyn_tot(map_wtheta(df) + k)
        moist_dyn_tot(map_wtheta(df) + k) = 0.0_r_def
      end do
    end do

  end subroutine atl_moist_dyn_mass_code

end module atl_moist_dyn_mass_kernel_mod
