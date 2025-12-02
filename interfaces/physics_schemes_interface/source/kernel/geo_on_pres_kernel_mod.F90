!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to geopotential on pressure levels code.

module geo_on_pres_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_INTEGER,                &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  use fs_continuity_mod,    only: WTHETA
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: geo_on_pres_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ),                          &
!        arg_type(GH_SCALAR_ARRAY,GH_REAL, GH_READ, 1), see PSyclone issue #1312
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: geo_on_pres_code
  end type geo_on_pres_kernel_type

  public :: geo_on_pres_code

contains

  !> @details Computes the geopotential height on pressure levels, with
  !!          appropriate extrapolation below the surface and above the model
  !!          top.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     data_in       Input data to interpolate
  !> @param[in]     ex_at_data    Exner pressure at levels of input data
  !> @param[in]     theta         Potential temperature
  !> @param[in]     height_wth    Height at theta levels
  !> @param[in]     exner_wth     Exner at theta levels
  !> @param[in]     nplev         Number of pressure levels
  !> @param[in]     plevs         Pressure level values
  !> @param[in,out] data_out      Interpolated data
  !> @param[in]     p_zero        Reference pressure
  !> @param[in]     kappa         Rd / cp
  !> @param[in]     cp            Specific heat at constant pressure
  !> @param[in]     gravity       Acceleration due to gravity
  !> @param[in]     ex_power      (cp * lapse_rate ) / g
  !> @param[in]     ndf_in        Number of degrees of freedom per cell for in fields
  !> @param[in]     undf_in       Number of total degrees of freedom for in fields
  !> @param[in]     map_in        Dofmap for the cell at the base of the column for in fields
  !> @param[in]     ndf_out       Number of degrees of freedom per cell for out fields
  !> @param[in]     undf_out      Number of total degrees of freedom for out fields
  !> @param[in]     map_out       Dofmap for the cell at the base of the column for out fields
  subroutine geo_on_pres_code(nlayers,                         &
                              data_in,                         &
                              ex_at_data,                      &
                              theta,                           &
                              height_wth,                      &
                              exner_wth,                       &
                              nplev,                           &
                              plevs,                           &
                              data_out,                        &
                              p_zero,                          &
                              kappa,                           &
                              cp,                              &
                              gravity,                         &
                              ex_power,                        &
                              ndf_in, undf_in, map_in,         &
                              ndf_wth, undf_wth, map_wth,      &
                              ndf_out, undf_out, map_out)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in), dimension(ndf_out)  :: map_out
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in), dimension(ndf_in)  :: map_in
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_in) :: data_in
    real(kind=r_def),    intent(in), dimension(undf_in) :: ex_at_data
    real(kind=r_def),    intent(in), dimension(undf_wth) :: theta
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def),    intent(inout), dimension(undf_out) :: data_out

    ! Constants passed explicitly from algorithm
    real(kind=r_def),    intent(in) :: p_zero, kappa, cp, gravity, ex_power
    real(kind=r_def),    intent(in), dimension(nplev) :: plevs

    ! Internal variables
    integer(kind=i_def) :: k, level_above, top_df, kp, level_extrap
    real(kind=r_def) :: desired_ex, h_ref_lev
    real(kind=r_def), parameter:: extrap_height = 2000.0_r_def

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
        ! Desired level is above model top, extrapolate up
        data_out(map_out(1)+kp-1) = data_in(map_in(1)+top_df) - &
                                   (desired_ex - ex_at_data(map_in(1)+top_df)) &
                                    * cp * theta(map_wth(1)+nlayers) / gravity
      else if (level_above == 0_i_def) then
        ! Desired level is below surface, extrapolate down
        do k = 1, nlayers
          if ((height_wth(map_wth(1)+k) - height_wth(map_wth(1))) &
               > extrap_height) then
            level_extrap = k-1
            exit
          end if
        end do
        h_ref_lev = (height_wth(map_wth(1)+level_extrap) - data_in(map_in(1))) &
                   / (1.0_r_def - (exner_wth(map_wth(1)+level_extrap) / &
                                   ex_at_data(map_in(1)))**ex_power )
        data_out(map_out(1)+kp-1) = data_in(map_in(1)+level_above) &
                                   + h_ref_lev * (1.0_r_def - &
                                  (desired_ex/ex_at_data(map_in(1)))**ex_power)
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

  end subroutine geo_on_pres_code

end module geo_on_pres_kernel_mod
