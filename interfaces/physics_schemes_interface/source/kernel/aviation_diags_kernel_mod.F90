!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
! Section 20 aviation diagnostics kernel.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file currently belongs in section: physics_schemes_interface
! whilst discussions are ongoing about its final location.
!
module aviation_diags_kernel_mod

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
  type, extends(kernel_type) :: aviation_diags_kernel_type
    type(arg_type), dimension(8) :: meta_args = (/ &

      ! Output fields.
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &

      ! Source field.
      arg_type(gh_field, gh_real, gh_read, any_discontinuous_space_2), &

      ! Request flags.
      arg_type(gh_scalar, gh_logical, gh_read), &
      arg_type(gh_scalar, gh_logical, gh_read), &

      ! Level indices.
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read) &
      /)

    integer :: operates_on = cell_column

    contains
      procedure, nopass :: code => aviation_diags_kernel_code
    end type aviation_diags_kernel_type

contains

  subroutine aviation_diags_kernel_code(nlayers,     &
             ! Output fields.
             thickness_850, thickness_500,           &

             ! Source field.
             plev_geopot,                            &

             ! Request flags.
             thickness_850_flag, thickness_500_flag, &

             ! Level incides.
             i1000, i850, i500,                      &

             ! Kernel stuff.
             result_ndf, result_undf, result_map,    &
             source_ndf, source_undf, source_map)

    ! Subtract geopotential heights at 850 and 500 hPa from that at 1000 hPa
    ! to calculate thickness fields.

    implicit none

    ! Arguments (kernel)

    ! The number of layers in a column.
    integer(kind=i_def), intent(in) :: nlayers

    ! Number of degrees of freedom (columns) in the cell we're processing.
    integer(kind=i_def), intent(in) :: result_ndf, source_ndf

    ! Number of unique degrees of freedom in the fields.
    integer(kind=i_def), intent(in) :: result_undf, source_undf

    ! Degrees of freedom maps. Offsets to the bottom of each column.
    integer(kind=i_def), intent(in), dimension(result_ndf) :: result_map
    integer(kind=i_def), intent(in), dimension(source_ndf) :: source_map


    ! Arguments (algorithm)

    ! Output thickness fields.
    real(kind=r_def), intent(out), dimension(result_undf) :: thickness_850
    real(kind=r_def), intent(out), dimension(result_undf) :: thickness_500

    ! Geopotential height at pressure levels.
    real(kind=r_def), intent(in), dimension(source_undf) :: plev_geopot

    ! Request flags.
    logical(kind=l_def), intent(in) :: thickness_850_flag, thickness_500_flag

    ! Level indices.
    integer(kind=i_def), intent(in) :: i1000, i850, i500


    ! Local variables
    integer(kind=i_def) :: df
    real(kind=r_def) :: gph_1000


    ! Process every DOF in this cell.
    do df = 1, result_ndf

      gph_1000 = plev_geopot(source_map(df) + i1000-1)

      if (thickness_850_flag .and. i850 /= -1) then
        thickness_850(result_map(df)) = &
          plev_geopot(source_map(df)+i850-1) - gph_1000
      end if

      if (thickness_500_flag .and. i500 /= -1) then
        thickness_500(result_map(df)) = &
          plev_geopot(source_map(df)+i500-1) - gph_1000
      end if

    end do

  end subroutine aviation_diags_kernel_code

end module aviation_diags_kernel_mod
