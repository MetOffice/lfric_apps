!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes initial exner field.
!>
module initial_exner_sample_kernel_mod

  use argument_mod,         only : arg_type, func_type,       &
                                   GH_FIELD, GH_SCALAR,       &
                                   GH_READ, GH_WRITE,         &
                                   GH_REAL, ANY_SPACE_9,      &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   GH_BASIS, CELL_COLUMN,     &
                                   GH_EVALUATOR
  use constants_mod,        only : r_def, i_def
  use fs_continuity_mod,    only : W3
  use idealised_config_mod, only : test
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: initial_exner_sample_kernel_type
      private
      type(arg_type) :: meta_args(4) = (/                                      &
           arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
           arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
           arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
           arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
           /)
      type(func_type) :: meta_funcs(1) = (/                                    &
           func_type(ANY_SPACE_9, GH_BASIS)                                    &
           /)
      integer :: operates_on = CELL_COLUMN
      integer :: gh_shape = GH_EVALUATOR
  contains
      procedure, nopass :: initial_exner_sample_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: initial_exner_sample_code

contains

  !> @brief Computes the initial Exner field
  !> @param[in] nlayers Number of layers
  !> @param[in,out] exner Pressure field
  !> @param[in] chi_1 First component of the chi coordinate field
  !> @param[in] chi_2 Second component of the chi coordinate field
  !> @param[in] chi_3 Third component of the chi coordinate field
  !> @param[in] panel_id A field giving the ID for mesh panels
  !> @param[in] current_time Current time of the model run
  !> @param[in] ndf_w3 Number of degrees of freedom per cell for W3
  !> @param[in] undf_w3 Total number of degrees of freedom for W3
  !> @param[in] map_w3 Dofmap for the cell at the base of the column for W3
  !> @param[in] ndf_chi Number of degrees of freedom per cell for Wchi
  !> @param[in] undf_chi Number of degrees of freedom for Wchi
  !> @param[in] map_chi Dofmap for the cell at the base of the column for Wchi
  !> @param[in] basis_chi_on_w3 Basis functions for Wchi evaluated at
  !!                            W3 nodal points
  !> @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !> @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !> @param[in] map_pid  Dofmap for the cell at the base of the column
  !!                     for panel_id
  subroutine initial_exner_sample_code(nlayers,                    &
                                       exner,                      &
                                       chi_1, chi_2, chi_3,        &
                                       panel_id,                   &
                                       current_time,               &
                                       ndf_w3, undf_w3, map_w3,    &
                                       ndf_chi, undf_chi, map_chi, &
                                       basis_chi_on_w3,            &
                                       ndf_pid, undf_pid, map_pid  )

    use analytic_pressure_profiles_mod, only : analytic_pressure
    use sci_chi_transform_mod,          only : chi2xyz

    implicit none

    ! Arguments
    integer(kind=i_def),                               intent(in)    :: nlayers
    integer(kind=i_def),                               intent(in)    :: ndf_w3
    integer(kind=i_def),                               intent(in)    :: ndf_chi
    integer(kind=i_def),                               intent(in)    :: ndf_pid
    integer(kind=i_def),                               intent(in)    :: undf_w3
    integer(kind=i_def),                               intent(in)    :: undf_chi
    integer(kind=i_def),                               intent(in)    :: undf_pid
    integer(kind=i_def), dimension(ndf_w3),            intent(in)    :: map_w3
    integer(kind=i_def), dimension(ndf_chi),           intent(in)    :: map_chi
    integer(kind=i_def), dimension(ndf_chi),           intent(in)    :: map_pid
    real(kind=r_def),    dimension(undf_w3),           intent(inout) :: exner
    real(kind=r_def),    dimension(undf_chi),          intent(in)    :: chi_1
    real(kind=r_def),    dimension(undf_chi),          intent(in)    :: chi_2
    real(kind=r_def),    dimension(undf_chi),          intent(in)    :: chi_3
    real(kind=r_def),    dimension(undf_pid),          intent(in)    :: panel_id
    real(kind=r_def),    dimension(1,ndf_chi, ndf_w3), intent(in)    :: basis_chi_on_w3
    real(kind=r_def),                                  intent(in)    :: current_time

    ! Internal variables
    real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                     :: coords(3), xyz(3)
    integer(kind=i_def)                  :: df1, df, k, ipanel

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Compute the RHS & LHS integrated over one cell and solve
    do k = 0, nlayers-1
      do df1 = 1, ndf_chi
        chi_1_e(df1) = chi_1(map_chi(df1) + k)
        chi_2_e(df1) = chi_2(map_chi(df1) + k)
        chi_3_e(df1) = chi_3(map_chi(df1) + k)
      end do

      do df = 1, ndf_w3
        coords = 0.0_r_def
        do df1 = 1, ndf_chi
          coords(1) = coords(1) + chi_1_e(df1)*basis_chi_on_w3(1,df1,df)
          coords(2) = coords(2) + chi_2_e(df1)*basis_chi_on_w3(1,df1,df)
          coords(3) = coords(3) + chi_3_e(df1)*basis_chi_on_w3(1,df1,df)
        end do

        call chi2xyz(coords(1), coords(2), coords(3), &
                     ipanel, xyz(1), xyz(2), xyz(3))

        exner(map_w3(df) + k) = analytic_pressure(xyz, test, current_time)

      end do
    end do

  end subroutine initial_exner_sample_code

end module initial_exner_sample_kernel_mod
