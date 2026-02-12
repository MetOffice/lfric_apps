!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial mixing ratio fields

!> @details The kernel computes initial mixing ratio fields for mr in the same
!>          space as that of theta

module initial_rel_humidity_kernel_mod

    use argument_mod,                  only: arg_type, func_type,       &
                                             GH_FIELD, GH_REAL,         &
                                             GH_SCALAR, GH_BASIS,       &
                                             GH_WRITE, GH_READ,         &
                                             ANY_SPACE_9,               &
                                             ANY_DISCONTINUOUS_SPACE_3, &
                                             GH_EVALUATOR, CELL_COLUMN
    use fs_continuity_mod,             only: W3, Wtheta
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use idealised_config_mod,          only: test

    implicit none
    private

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_rel_humidity_kernel_type
        private
        type(arg_type) :: meta_args(3) = (/                                     &
             arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                   &
             arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),              &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
             /)
        type(func_type) :: meta_funcs(1) = (/                                   &
             func_type(ANY_SPACE_9, GH_BASIS)                                   &
             /)
        integer :: operates_on = CELL_COLUMN
        integer :: gh_shape = GH_EVALUATOR
    contains
        procedure, nopass :: initial_rel_humidity_code
    end type

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public :: initial_rel_humidity_code
contains

  subroutine initial_rel_humidity_code(nlayers,                          &
                                  rel_hum,                               &
                                  chi_1, chi_2, chi_3,                   &
                                  panel_id,                              &
                                  ndf_wtheta, undf_wtheta, map_wtheta,   &
                                  ndf_chi, undf_chi, map_chi, chi_basis, &
                                  ndf_pid, undf_pid, map_pid             )

    use analytic_moisture_profiles_mod, only : analytic_rel_humidity
    use sci_chi_transform_mod,          only : chi2xyz

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), dimension(ndf_wtheta),  intent(in) :: map_wtheta
    integer(kind=i_def), dimension(ndf_chi),     intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: rel_hum
    real(kind=r_def), dimension(undf_chi),    intent(in)    :: chi_1, chi_2, chi_3
    real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id

    real(kind=r_def), dimension(1,ndf_chi,ndf_wtheta), intent(in) :: chi_basis

    ! Internal variables
    integer(kind=i_def) :: k, df, dfc, ipanel
    real(kind=r_def)    :: coords(3), xyz(3)
    real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e

    ipanel = int(panel_id(map_pid(1)), i_def)


    do k = 0, nlayers-1
      do dfc = 1, ndf_chi
        chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
        chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
        chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
      end do

      do df = 1, ndf_wtheta
        coords(:) = 0.0_r_def
        do dfc = 1, ndf_chi
          coords(1) = coords(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
          coords(2) = coords(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
          coords(3) = coords(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
        end do

        call chi2xyz(coords(1), coords(2), coords(3), &
                      ipanel, xyz(1), xyz(2), xyz(3))
        rel_hum(map_wtheta(df) + k) = analytic_rel_humidity(xyz, test)

      end do
    end do

  end subroutine initial_rel_humidity_code

end module initial_rel_humidity_kernel_mod
