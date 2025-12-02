!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the initial buoyancy field.
!!
!> @details Computes initial buoyancy perturbation field from a given
!!          analytic expression.
module initial_swe_buoyancy_kernel_mod

  use argument_mod,         only: arg_type, func_type,                &
                                  GH_FIELD, GH_INC, GH_READ,          &
                                  ANY_SPACE_9, ANY_SPACE_1, GH_BASIS, &
                                  GH_DIFF_BASIS, GH_REAL, GH_SCALAR,  &
                                  CELL_COLUMN, GH_EVALUATOR,          &
                                  ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use shallow_water_settings_config_mod, &
                            only: swe_test

  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: initial_swe_buoyancy_kernel_type
      private
      type(arg_type) :: meta_args(4) = (/                                     &
          arg_type(GH_FIELD,   GH_REAL, GH_INC,   ANY_SPACE_1),               &
          arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
          arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
          /)
      type(func_type) :: meta_funcs(1) = (/                    &
          func_type(ANY_SPACE_9, GH_BASIS)                     &
          /)
      integer :: operates_on = CELL_COLUMN
      integer :: gh_shape = GH_EVALUATOR
  contains
      procedure, nopass :: initial_swe_buoyancy_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public initial_swe_buoyancy_code

contains

  !> @brief Compute the initial buoyancy field.
  !> @param[in]     nlayers  The number of model layers
  !> @param[in,out] buoyancy The field to compute
  !> @param[in]     chi_1    Real array, the 1st component of the coordinate field
  !> @param[in]     chi_2    Real array, the 2nd component of the coordinate field
  !> @param[in]     chi_3    Real array, the 3rd component of the coordinate field
  !> @param[in]     panel_id Field containing the ID of the mesh panel
  !> @param[in]     domain_x Domain size in x direction
  !> @param[in]     ndf_wb   The number of degrees of freedom per cell for wb
  !> @param[in]     undf_wb  The number of total degrees of freedom for wb
  !> @param[in]     map_wb   Integer array holding the dofmap for the cell at the base of the column
  !> @param[in]     ndf_wx   The number of degrees of freedom per cell
  !> @param[in]     undf_wx  The total number of degrees of freedom
  !> @param[in]     map_wx   Integer array holding the dofmap for the cell at the base of the column
  !> @param[in]     wx_basis Real 3-dim array holding basis functions evaluated at quadrature points
  !> @param[in]     ndf_pid  The number of DoFs per cell for the panel ID
  !> @param[in]     undf_pid The number of DoFs for this partition for the panel ID
  !> @param[in]     map_pid  DoF-map for the panel ID
  subroutine initial_swe_buoyancy_code(nlayers,                           &
                                       buoyancy,                          &
                                       chi_1, chi_2, chi_3,               &
                                       panel_id,                          &
                                       domain_x,                          &
                                       ndf_wb, undf_wb, map_wb,           &
                                       ndf_wx, undf_wx, map_wx, wx_basis, &
                                       ndf_pid, undf_pid, map_pid)

    use analytic_swe_buoyancy_profiles_mod, only : analytic_swe_buoyancy
    use sci_chi_transform_mod,              only : chi2xyz

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wb, ndf_wx, undf_wb, undf_wx
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), dimension(ndf_wb),  intent(in) :: map_wb
    integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
    integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
    real(kind=r_def), dimension(undf_wb),          intent(inout) :: buoyancy
    real(kind=r_def),                              intent(in)    :: domain_x
    real(kind=r_def), dimension(undf_wx),          intent(in)    :: chi_1, chi_2, chi_3
    real(kind=r_def), dimension(undf_pid),         intent(in)    :: panel_id
    real(kind=r_def), dimension(1,ndf_wx,ndf_wb),  intent(in)    :: wx_basis

    ! Internal variables
    integer(kind=i_def)                 :: df, df0, ipanel
    real(kind=r_def), dimension(ndf_wx) :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                    :: coord(3), xyz(3)

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Avoid compile warning for unused variables
    df = nlayers

    ! compute the pointwise buoyancy profile
    do df0 = 1, ndf_wx
      chi_1_e(df0) = chi_1( map_wx(df0) )
      chi_2_e(df0) = chi_2( map_wx(df0) )
      chi_3_e(df0) = chi_3( map_wx(df0) )
    end do

    do df = 1, ndf_wb
      coord(:) = 0.0_r_def
      do df0 = 1, ndf_wx
        coord(1) = coord(1) + chi_1_e(df0)*wx_basis(1,df0,df)
        coord(2) = coord(2) + chi_2_e(df0)*wx_basis(1,df0,df)
        coord(3) = coord(3) + chi_3_e(df0)*wx_basis(1,df0,df)
      end do

      call chi2xyz(coord(1), coord(2), coord(3), ipanel, xyz(1), xyz(2), xyz(3))
      buoyancy(map_wb(df)) = analytic_swe_buoyancy(xyz, swe_test, domain_x)
    end do

  end subroutine initial_swe_buoyancy_code

end module initial_swe_buoyancy_kernel_mod
