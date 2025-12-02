!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the initial buoyancy field

!> @details Computes initial buoyancy perturbation field from a given
!!          analytic expression
module initial_buoyancy_kernel_mod

    use argument_mod,  only: arg_type, func_type,                &
                             GH_FIELD, GH_REAL, GH_INC, GH_READ, &
                             ANY_SPACE_9, ANY_SPACE_1, GH_BASIS, &
                             CELL_COLUMN, GH_EVALUATOR,          &
                             ANY_DISCONTINUOUS_SPACE_3
    use constants_mod, only: r_def, i_def
    use kernel_mod,    only: kernel_type

    implicit none

    private

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_buoyancy_kernel_type
        private
        type(arg_type) :: meta_args(3) = (/                                    &
             arg_type(GH_FIELD,   GH_REAL, GH_INC,  ANY_SPACE_1),              &
             arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
             arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
             /)
        type(func_type) :: meta_funcs(1) = (/                                  &
             func_type(ANY_SPACE_9, GH_BASIS)                                  &
             /)
        integer :: operates_on = CELL_COLUMN
        integer :: gh_shape = GH_EVALUATOR
    contains
        procedure, nopass :: initial_buoyancy_code
    end type

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public :: initial_buoyancy_code
contains

    !> @brief Compute the initial buoyancy field
    !! @param[in] nlayers The number of model layers
    !! @param[in,out] buoyancy The field to compute
    !! @param[in] chi_1 Real array, the x component of the coordinate field
    !! @param[in] chi_2 Real array, the y component of the coordinate field
    !! @param[in] chi_3 Real array, the z component of the coordinate field
    !! @param[in] panel_id Field with IDs of mesh panels
    !! @param[in] ndf_wt The number of degrees of freedom per cell for wt
    !! @param[in] undf_wt The number of total degrees of freedom for wt
    !! @param[in] map_wt Integer array holding the dofmap for the cell at the base of the column
    !! @param[in] ndf_wx The number of degrees of freedom per cell
    !! @param[in] undf_wx The total number of degrees of freedom
    !! @param[in] map_wx Integer array holding the dofmap for the cell at the base of the column
    !! @param[in] wx_basis Real 5-dim array holding basis functions evaluated at quadrature points,
    !! @param[in] ndf_pid The number of degrees of freedom per cell for panel id
    !! @param[in] ndf_pid The total number of degrees of freedom for panel id
    !! @param[in] map_pid Integer array with the dofmap for panel_id
    subroutine initial_buoyancy_code(nlayers,                           &
                                     buoyancy,                          &
                                     chi_1, chi_2, chi_3,               &
                                     panel_id,                          &
                                     ndf_wt, undf_wt, map_wt,           &
                                     ndf_wx, undf_wx, map_wx, wx_basis, &
                                     ndf_pid, undf_pid, map_pid         )

        use analytic_buoyancy_profiles_mod, only : analytic_buoyancy
        use sci_chi_transform_mod,          only : chi2xyz

        implicit none

        ! Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wt, ndf_wx, ndf_pid
        integer(kind=i_def), intent(in) :: undf_wt, undf_wx, undf_pid
        integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
        integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
        integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
        real(kind=r_def), dimension(undf_wt),          intent(inout) :: buoyancy
        real(kind=r_def), dimension(undf_wx),          intent(in)    :: chi_1, chi_2, chi_3
        real(kind=r_def), dimension(undf_pid),         intent(in)    :: panel_id
        real(kind=r_def), dimension(1,ndf_wx,ndf_wt),  intent(in)    :: wx_basis

        ! Internal variables
        integer(kind=i_def)                 :: df, df0, k, ipanel
        real(kind=r_def), dimension(ndf_wx) :: chi_1_e, chi_2_e, chi_3_e
        real(kind=r_def)                    :: chi_wt(3), xyz(3)

        ipanel = int(panel_id(map_pid(1)), i_def)

        ! compute the pointwise buoyancy profile
        do k = 0, nlayers-1
          do df0 = 1, ndf_wx
            chi_1_e(df0) = chi_1( map_wx(df0) + k )
            chi_2_e(df0) = chi_2( map_wx(df0) + k )
            chi_3_e(df0) = chi_3( map_wx(df0) + k )
          end do

          do df = 1, ndf_wt
            chi_wt(:) = 0.0_r_def
            do df0 = 1, ndf_wx
              chi_wt(1) = chi_wt(1) + chi_1_e(df0)*wx_basis(1,df0,df)
              chi_wt(2) = chi_wt(2) + chi_2_e(df0)*wx_basis(1,df0,df)
              chi_wt(3) = chi_wt(3) + chi_3_e(df0)*wx_basis(1,df0,df)
            end do

            call chi2xyz(chi_wt(1), chi_wt(2), chi_wt(3), &
                         ipanel, xyz(1), xyz(2), xyz(3))

            buoyancy(map_wt(df) + k) = analytic_buoyancy(xyz)
          end do
        end do

    end subroutine initial_buoyancy_code

end module initial_buoyancy_kernel_mod
