!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
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
MODULE aviation_diags_kernel_mod

  USE argument_mod,         ONLY: arg_type,            &
                                  GH_FIELD, GH_SCALAR, GH_LOGICAL, &
                                  GH_READ, GH_WRITE, GH_INTEGER, &
                                  GH_REAL, CELL_COLUMN, &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  USE kernel_mod,           ONLY: kernel_type
  USE constants_mod,        ONLY: r_def, i_def, l_def

  IMPLICIT NONE

  ! The aviation diagnostics kernel type.
  TYPE, EXTENDS(kernel_type) :: aviation_diags_kernel_type
    TYPE(arg_type), DIMENSION(8) :: meta_args = (/ &

      ! Output fields.
      arg_type(gh_field, gh_real, gh_write, ANY_DISCONTINUOUS_SPACE_1), &
      arg_type(gh_field, gh_real, gh_write, ANY_DISCONTINUOUS_SPACE_1), &

      ! Source field.
      arg_type(gh_field, gh_real, gh_read, ANY_DISCONTINUOUS_SPACE_2), &

      ! Request flags.
      arg_type(gh_scalar, gh_logical, gh_read), &
      arg_type(gh_scalar, gh_logical, gh_read), &

      ! Level indices.
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read), &
      arg_type(gh_scalar, gh_integer, gh_read) &
      /)

    INTEGER :: operates_on = cell_column

    CONTAINS
      PROCEDURE, NOPASS :: code => aviation_diags_kernel_code
    END TYPE aviation_diags_kernel_type

CONTAINS

  SUBROUTINE aviation_diags_kernel_code(nlayers, &
             ! Output fields.
             thickness_850, thickness_500, &

             ! Source field.
             plev_geopot, &

             ! Request flags.
             thickness_850_flag, thickness_500_flag, &

             ! Level incides.
             i1000, i850, i500, &

             ! Kernel stuff.
             result_ndf, result_undf, result_map, &
             source_ndf, source_undf, source_map)

    ! Subtract geopotential heights at 850 and 500 hPa from that at 1000 hPa
    ! to calculate thickness fields.

    IMPLICIT NONE

    ! Arguments (kernel)

    ! The number of layers in a column.
    INTEGER(KIND=i_def), INTENT(IN) :: nlayers

    ! Number of degrees of freedom (columns) in the cell we're processing.
    INTEGER(KIND=i_def), INTENT(IN) :: result_ndf, source_ndf

    ! Number of unique degrees of freedom in the fields.
    INTEGER(KIND=i_def), INTENT(IN) :: result_undf, source_undf

    ! Degrees of freedom maps. Offsets to the bottom of each column.
    INTEGER(KIND=i_def), INTENT(IN), DIMENSION(result_ndf) :: result_map
    INTEGER(KIND=i_def), INTENT(IN), DIMENSION(source_ndf) :: source_map


    ! Arguments (algorithm)

    ! Output thickness fields.
    REAL(KIND=r_def), INTENT(OUT), DIMENSION(result_undf) :: thickness_850
    REAL(KIND=r_def), INTENT(OUT), DIMENSION(result_undf) :: thickness_500

    ! Geopotential height at pressure levels.
    REAL(KIND=r_def), INTENT(IN), DIMENSION(source_undf) :: plev_geopot

    ! Request flags.
    LOGICAL(KIND=l_def), INTENT(IN) :: thickness_850_flag, thickness_500_flag

    ! Level indices. For i850 and i500, -1 means "not requested".
    INTEGER(KIND=i_def), INTENT(IN) :: i1000, i850, i500


    ! Local variables
    INTEGER(KIND=i_def) :: df
    REAL(KIND=r_def) :: gph_1000


    ! Process every DOF in this cell.
    DO df = 1, result_ndf

      gph_1000 = plev_geopot(source_map(df) + i1000-1)

      IF (i850 /= -1) THEN
        thickness_850(result_map(df)) = &
          plev_geopot(source_map(df)+i850-1) - gph_1000
      END IF

      IF (i500 /= -1) THEN
        thickness_500(result_map(df)) = &
          plev_geopot(source_map(df)+i500-1) - gph_1000
      END IF

    END DO

  END SUBROUTINE aviation_diags_kernel_code

END MODULE aviation_diags_kernel_mod
