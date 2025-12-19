!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Aviation diagnostics

module aviation_diags_kernel_mod

  use argument_mod,         only: arg_type,            &
                                  GH_FIELD, GH_SCALAR, &
                                  GH_READ, GH_WRITE, GH_INTEGER, &
                                  GH_REAL, CELL_COLUMN, &
                                  ANY_DISCONTINUOUS_SPACE_1, ANY_DISCONTINUOUS_SPACE_2
!  use fs_continuity_mod,    only: WTHETA, W3
  use kernel_mod,           only: kernel_type
  use constants_mod,        only: r_def, i_def

  implicit none

  type, extends(kernel_type) :: aviation_diags_kernel_type
    type(arg_type), dimension(6) :: meta_args = (/ &
            ! output
            arg_type(gh_field, gh_real, gh_write, ANY_DISCONTINUOUS_SPACE_1), &
            arg_type(gh_field, gh_real, gh_write, ANY_DISCONTINUOUS_SPACE_1), &
            ! source field
            arg_type(gh_field, gh_real, gh_read, ANY_DISCONTINUOUS_SPACE_2), &
            ! level indices
            arg_type(gh_scalar, gh_integer, gh_read), &
            arg_type(gh_scalar, gh_integer, gh_read), &
            arg_type(gh_scalar, gh_integer, gh_read) &
            /)
    integer :: operates_on = cell_column
  contains
    procedure, nopass :: code => aviation_diags_kernel_code
  end type aviation_diags_kernel_type

contains

  SUBROUTINE aviation_diags_kernel_code(nlayers, &
          thickness_850, thickness_500, &
          plev_geopot, i1000, i850, i500, &
          result_ndf, result_undf, result_map, &
          source_ndf, source_undf, source_map)
    USE constants_mod

    IMPLICIT NONE


    ! kernel params

    ! the number of layers in a column
    INTEGER(KIND=i_def), intent(in) :: nlayers

    ! number of degrees of freedom for the particular column
    INTEGER(KIND=i_def), intent(in) :: result_ndf, source_ndf

    ! number of unique degrees of freedom
    INTEGER(KIND=i_def), intent(in) :: result_undf, source_undf

    ! degrees of freedom map (dofmap) which indicates the location of the required values in the field array
    INTEGER(KIND=i_def), intent(in), dimension(result_ndf) :: result_map
    INTEGER(KIND=i_def), intent(in), dimension(source_ndf) :: source_map


    ! algorithm params

    ! Level indices. For i850 and i500, -1 = not requested.
    INTEGER(KIND=i_def), intent(in) :: i1000, i850, i500


    ! input and output fields
    REAL(KIND=r_def), intent(inout), dimension(result_undf) :: thickness_850, thickness_500
    REAL(KIND=r_def), intent(in), dimension(source_undf) :: plev_geopot  ! geopotential height at pressure levels


    ! local variables
    INTEGER(KIND=i_def) :: df


    ! process every dof in this cell
    do df = 1, result_ndf

      if (i850 /= -1) then
        thickness_850(result_map(df)) = plev_geopot(source_map(df)+i850-1) - plev_geopot(source_map(df)+i1000-1)
      end if

      if (i500 /= -1) then
        thickness_500(result_map(df)) = plev_geopot(source_map(df)+i500-1) - plev_geopot(source_map(df)+i1000-1)
      end if

    end do

  END SUBROUTINE aviation_diags_kernel_code

end module aviation_diags_kernel_mod
