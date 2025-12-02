!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute sum of dry and moisture mixing ratios
module moist_dyn_mass_kernel_mod

    use argument_mod,                  only: arg_type,          &
                                             GH_FIELD, GH_REAL, &
                                             GH_SCALAR,         &
                                             GH_WRITE, GH_READ, &
                                             CELL_COLUMN
    use constants_mod,                 only: r_def, i_def
    use fs_continuity_mod,             only: Wtheta
    use kernel_mod,                    only: kernel_type

    implicit none

    private

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: moist_dyn_mass_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                  &
             arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),  &
             arg_type(GH_FIELD*6, GH_REAL, GH_READ,  Wtheta) &
             /)
        integer :: operates_on = CELL_COLUMN
    contains
        procedure, nopass :: moist_dyn_mass_code
    end type

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public :: moist_dyn_mass_code

contains

    !> @brief Compute the total of all the mass mixing ratios
    !! @param[in] nlayers Integer the number of layers
    !! @param[in,out] moist_dyn_tot Total mass factor (1 + sum m_x)
    !! @param[in] mr_v Water vapour mixing ratio
    !! @param[in] mr_cl Liquid cloud mixing ratio
    !! @param[in] mr_r Rain mixing ratio
    !! @param[in] mr_s Snow mixing ratio
    !! @param[in] mr_g Graupel mixing ratio
    !! @param[in] mr_ci Cloud ice mixing ratio
    !! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
    !! @param[in] udf_wtheta The number of total degrees of freedom for wtheta
    !! @param[in] map_wtheta Integer array holding the dofmap for the cell at the base of the column

    subroutine moist_dyn_mass_code(nlayers, moist_dyn_tot,               &
                                   mr_v, mr_cl, mr_r, mr_s, mr_g, mr_ci, &
                                   ndf_wtheta, undf_wtheta, map_wtheta )

        implicit none

        ! Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta

        real(kind=r_def), dimension(undf_wtheta),   intent(inout) :: moist_dyn_tot
        real(kind=r_def), dimension(undf_wtheta),   intent(in)    :: mr_v, mr_cl, mr_r
        real(kind=r_def), dimension(undf_wtheta),   intent(in)    :: mr_ci, mr_s, mr_g
        integer(kind=i_def), dimension(ndf_wtheta), intent(in)    :: map_wtheta

        ! Internal variables
        integer(kind=i_def)                 :: k, df

        real(kind=r_def)                    :: mr_v_at_dof, mr_cl_at_dof, mr_r_at_dof
        real(kind=r_def)                    :: mr_ci_at_dof, mr_s_at_dof, mr_g_at_dof

        ! compute the pointwise mr profile
        do k = 0, nlayers-1

          do df = 1, ndf_wtheta
            mr_v_at_dof  = mr_v(map_wtheta(df) + k)
            mr_cl_at_dof = mr_cl(map_wtheta(df) + k)
            mr_r_at_dof  = mr_r(map_wtheta(df) + k)
            mr_ci_at_dof = mr_ci(map_wtheta(df) + k)
            mr_s_at_dof  = mr_s(map_wtheta(df) + k)
            mr_g_at_dof  = mr_g(map_wtheta(df) + k)
            moist_dyn_tot(map_wtheta(df) + k) = 1.0_r_def + mr_v_at_dof  +  &
                                                            mr_cl_at_dof +  &
                                                            mr_r_at_dof  +  &
                                                            mr_s_at_dof  +  &
                                                            mr_g_at_dof  +  &
                                                            mr_ci_at_dof
          end do

        end do

    end subroutine moist_dyn_mass_code

end module moist_dyn_mass_kernel_mod
