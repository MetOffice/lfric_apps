! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!> @brief Calculate geopotential thickness and snow probabbility.
!>
module aviation_thickness_snow_kernel_mod

  use argument_mod,         only: arg_type,                        &
                                  gh_field, gh_scalar, gh_logical, &
                                  gh_read, gh_write, gh_integer,   &
                                  gh_real, cell_column,            &
                                  any_discontinuous_space_1,       &
                                  any_discontinuous_space_2
  use kernel_mod,           only: kernel_type
  use constants_mod,        only: r_def, i_def, l_def

  implicit none

  ! The aviation diagnostics kernel type.
  type, extends(kernel_type) :: aviation_thickness_snow_kernel_type
    type(arg_type), dimension(10) :: meta_args = (/ &

      ! Output fields: thickness 850 & 500, and snow probability.
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &

      ! Source field plev_geopot.
      arg_type(gh_field, gh_real, gh_read, any_discontinuous_space_2), &

      ! Field request flags.
      arg_type(gh_scalar, gh_logical, gh_read), &
      arg_type(gh_scalar, gh_logical, gh_read), &
      arg_type(gh_scalar, gh_logical, gh_read), &

      ! Level indices: 1000, 850 & 500.
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read) &
      /)

    integer :: operates_on = cell_column

    contains
      procedure, nopass :: code => aviation_thickness_snow_kernel_code
    end type aviation_thickness_snow_kernel_type

contains


  !> @brief Calculate geopotential thickness and snow probabbility
  !>        from geopotential height at pressure levels.
  !> @details Assumes lowest order W3 data, where ndf is always 1.
  !>          Thickness:
  !>            Subtract geopotential heights at 850 and 500hPa from 1000hPa.
  !>          Snow probability:
  !>            Implement Boyden (1964), using 850 and 1000 hPa.
  !>            See https://github.com/MetOffice/Section20/issues/21
  !> @param[in]     nlayers                 The number of layers in a column.
  !> @param[out]    thickness_850           Output geopot thickness between 1000 and 850 hPa in m.
  !> @param[out]    thickness_500           Output geopot thickness between 1000 and 500 hPa in m.
  !> @param[out]    snow_probability        Output snow probability in percentage.
  !> @param[in]     plev_geopot             Geopotential height at pressure levels in m.
  !> @param[in]     thickness_850_flag      thickness_850 request flag.
  !> @param[in]     thickness_500_flag      thickness_500 request flag.
  !> @param[in]     snow_probability_flag   snow_probability request flag.
  !> @param[in]     i1000                   Index of 1000 hPa level in plev_geopot.
  !> @param[in]     i850                    Index of 850 hPa level in plev_geopot.
  !> @param[in]     i500                    Index of 500 hPa level in plev_geopot.
  !> @param[in]     result_ndf              Number of DOFs in the result cell.
  !> @param[in]     result_undf             Number of DOFs in the result field.
  !> @param[in]     result_map              Dofmap to the result's bottom cell.
  !> @param[in]     source_ndf              Number of DOFs in the source cell.
  !> @param[in]     source_undf             Number of DOFs in the source field.
  !> @param[in]     source_map              Dofmap to the source's bottom cell.
  subroutine aviation_thickness_snow_kernel_code(nlayers, &
             ! output fields.
             thickness_850, thickness_500, snow_probability, &

             ! source field.
             plev_geopot, &

             ! request flags.
             thickness_850_flag, thickness_500_flag, snow_probability_flag, &

             ! level incides.
             i1000, i850, i500, &

             ! kernel stuff.
             result_ndf, result_undf, result_map, &
             source_ndf, source_undf, source_map)

    implicit none

    ! Arguments (kernel).
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: result_ndf, source_ndf
    integer(kind=i_def), intent(in) :: result_undf, source_undf
    integer(kind=i_def), intent(in), dimension(result_ndf) :: result_map
    integer(kind=i_def), intent(in), dimension(source_ndf) :: source_map

    ! Arguments (algorithm).
    real(kind=r_def), intent(out), dimension(result_undf) :: thickness_850
    real(kind=r_def), intent(out), dimension(result_undf) :: thickness_500
    real(kind=r_def), intent(out), dimension(result_undf) :: snow_probability
    real(kind=r_def), intent(in), dimension(source_undf) :: plev_geopot

    ! Request flags.
    logical(kind=l_def), intent(in) :: thickness_850_flag
    logical(kind=l_def), intent(in) :: thickness_500_flag
    logical(kind=l_def), intent(in) :: snow_probability_flag

    ! Level indices.
    integer(kind=i_def), intent(in) :: i1000, i850, i500


    ! Local variables.
    real(kind=r_def) :: gph_1000, gph_850


    ! Process the single DOF in this cell.
    gph_1000 = plev_geopot(source_map(1) + i1000-1)

    if (thickness_850_flag) then
      thickness_850(result_map(1)) = &
        plev_geopot(source_map(1)+i850-1) - gph_1000
    end if

    if (thickness_500_flag) then
      thickness_500(result_map(1)) = &
        plev_geopot(source_map(1)+i500-1) - gph_1000
    end if

    if(snow_probability_flag) then
      gph_850 = plev_geopot(source_map(1) + i850-1)
      snow_probability(result_map(1)) = &
        5220.0_r_def + 3.86666_r_def*gph_1000 - 4.0_r_def*gph_850

      if (gph_1000 > 1.0e8) then
          snow_probability(result_map(1)) = 0.0_r_def
      end if

      ! Limit to percentage.
      if (snow_probability(result_map(1)) < 0.0_r_def) then
          snow_probability(result_map(1)) = 0.0_r_def
      else if (snow_probability(result_map(1)) > 100.0_r_def) then
          snow_probability(result_map(1)) = 100.0_r_def
      end if

    end if

  end subroutine aviation_thickness_snow_kernel_code

end module aviation_thickness_snow_kernel_mod
