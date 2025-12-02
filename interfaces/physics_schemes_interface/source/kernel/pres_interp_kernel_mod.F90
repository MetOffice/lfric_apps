!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to pressure interpolation code.

module pres_interp_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_INTEGER,                &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: pres_interp_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ),                          &
!        arg_type(GH_SCALAR_ARRAY,GH_REAL, GH_READ, 1), see PSyclone issue #1312
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pres_interp_code
  end type pres_interp_kernel_type

  public :: pres_interp_code

contains

  !> @details Interpolates a field onto pressure levels, with appropriate
  !!          extrapolation below the surface and above the model top.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     data_in       Input data to interpolate
  !> @param[in]     exner         Exner pressure at levels of input data
  !> @param[in]     nplev         Number of pressure levels
  !> @param[in]     plevs         Pressure level values
  !> @param[in,out] data_out      Interpolated data
  !> @param[in]     p_zero        Reference pressure
  !> @param[in]     kappa         Rd / cp
  !> @param[in]     ndf_in        Number of degrees of freedom per cell for in fields
  !> @param[in]     undf_in       Number of total degrees of freedom for in fields
  !> @param[in]     map_in        Dofmap for the cell at the base of the column for in fields
  !> @param[in]     ndf_out       Number of degrees of freedom per cell for out fields
  !> @param[in]     undf_out      Number of total degrees of freedom for out fields
  !> @param[in]     map_out       Dofmap for the cell at the base of the column for out fields
  subroutine pres_interp_code(nlayers,                         &
                              data_in,                         &
                              ex_at_data,                      &
                              nplev,                           &
                              plevs,                           &
                              data_out,                        &
                              p_zero,                          &
                              kappa,                           &
                              ndf_in, undf_in, map_in,         &
                              ndf_out, undf_out, map_out)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in), dimension(ndf_out)  :: map_out
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in), dimension(ndf_in)  :: map_in

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_in) :: data_in
    real(kind=r_def),    intent(in), dimension(undf_in) :: ex_at_data
    real(kind=r_def),    intent(inout), dimension(undf_out) :: data_out

    ! Constants passed explicitly from algorithm
    real(kind=r_def),    intent(in) :: p_zero, kappa
    real(kind=r_def),    intent(in), dimension(nplev) :: plevs

    ! Internal variables
    integer(kind=i_def) :: k, level_above, top_df, kp
    real(kind=r_def) :: desired_ex

    do kp = 1, nplev

      ! Level we want
      desired_ex = (plevs(kp)/p_zero)**kappa

      ! For Wtheta ndf = 2, loop k = 0, nlayers
      ! For W3 ndf = 1, loop k = 0, nlayers - 1
      top_df = nlayers + ndf_in - 2

      ! Find first level where real pressure is below desired pressure
      level_above = -1_i_def
      do k = 0, top_df
        if ( (ex_at_data(map_in(1)+k) < desired_ex) ) then
          level_above = k
          exit
        end if
      end do

      if (level_above == -1_i_def) then
        ! Desired level is above model top, use highest level value
        data_out(map_out(1)+kp-1) = data_in(map_in(1)+top_df)
      else if (level_above == 0_i_def) then
        ! Desired level is below surface, use lowest level value
        data_out(map_out(1)+kp-1) = data_in(map_in(1)+level_above)
      else
        ! Linear interpolation
        data_out(map_out(1)+kp-1) = ( (desired_ex - &
                                       ex_at_data(map_in(1)+level_above-1)) * &
                                      data_in(map_in(1)+level_above) &
                                    - (desired_ex - &
                                       ex_at_data(map_in(1)+level_above)) * &
                                      data_in(map_in(1)+level_above-1) ) / &
                                    (ex_at_data(map_in(1)+level_above) - &
                                     ex_at_data(map_in(1)+level_above-1))
      end if

    end do

  end subroutine pres_interp_code

end module pres_interp_kernel_mod
