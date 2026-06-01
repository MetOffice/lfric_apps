!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Logs values within a single column
!>
module log_column_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, &
                               CELL_COLUMN, GH_SCALAR, GH_WRITE,     &
                               GH_INTEGER,                           &
                               ANY_DISCONTINUOUS_SPACE_1,            &
                               ANY_DISCONTINUOUS_SPACE_2,            &
                               ANY_DISCONTINUOUS_SPACE_3,            &
                               ANY_DISCONTINUOUS_SPACE_4
  use constants_mod,     only: r_double, r_single, r_def, i_def
  use kernel_mod,        only: kernel_type
  use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: log_column_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                               &
    /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: log_column_code

  interface log_column_code
    module procedure &
      log_column_code_r_single, &
      log_column_code_r_double
    end interface

contains

!> @details       Logs all of the field values within a single column
!> @param[in]     nlayers    Number of layers in the mesh
!> @param[in]     dummy      Dummy field to allow psyclone to compile
!> @param[in]     field      Field to be logged
!> @param[in]     panel_id   2D W3 field, whose DoFmap is the column indices
!> @param[in]     column     The index of the column to log
!> @param[in]     ndf_3d     Number of dofs per cell for 3d field
!> @param[in]     undf_3d    Total number of Dofs in this partition for 3d field
!> @param[in]     map_3d     Dofmap for 3d field
!> @param[in]     ndf_2d     Number of dofs per cell for 2d field
!> @param[in]     undf_2d    Total number of Dofs in this partition for 2d field
!> @param[in]     map_2d     Dofmap for 2d field
subroutine log_column_code_r_single( nlayers,  &
                                     dummy,    &
                                     field,    &
                                     panel_id, &
                                     column,   &
                                     ndf_3d,   &
                                     undf_3d,  &
                                     map_3d,   &
                                     ndf_2d,   &
                                     undf_2d,  &
                                     map_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: column
  integer(kind=i_def), intent(in)    :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in)    :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in)    :: map_3d(ndf_3d)
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  real(kind=r_single), intent(inout) :: dummy(undf_3d)
  real(kind=r_single), intent(in)    :: field(undf_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_2d)

  ! Internal variables
  integer(kind=i_def) :: k, df

  if (map_2d(1) == column) then
    select case (ndf_3d)

    case (1)
      do k = 1, nlayers
        write(log_scratch_space, *) 'Column ', column, ' layer ', k, ':', field(map_3d(1)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case (2)
      do k = 1, nlayers + 1
        write(log_scratch_space, *) 'Column ', column, ' layer ', k, ':', field(map_3d(1)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case (4)
      do df = 1, 4
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do

    case (6)
      do df = 1, 4
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do
      df = 5
      do k = 1, nlayers + 1
        write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                    ' layer ', k, ':', field(map_3d(df)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case default
      do df = 1, ndf_3d
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do

    end select
  end if

end subroutine log_column_code_r_single

subroutine log_column_code_r_double( nlayers,  &
                                     dummy,    &
                                     field,    &
                                     panel_id, &
                                     column,   &
                                     ndf_3d,   &
                                     undf_3d,  &
                                     map_3d,   &
                                     ndf_2d,   &
                                     undf_2d,  &
                                     map_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: column
  integer(kind=i_def), intent(in)    :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in)    :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in)    :: map_3d(ndf_3d)
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  real(kind=r_double), intent(inout) :: dummy(undf_3d)
  real(kind=r_double), intent(in)    :: field(undf_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_2d)

  ! Internal variables
  integer(kind=i_def) :: k, df

  if (map_2d(1) == column) then
    select case (ndf_3d)

    case (1)
      do k = 1, nlayers
        write(log_scratch_space, *) 'Column ', column, ' layer ', k, ':', field(map_3d(1)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case (2)
      do k = 1, nlayers + 1
        write(log_scratch_space, *) 'Column ', column, ' layer ', k, ':', field(map_3d(1)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case (4)
      do df = 1, 4
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do

    case (6)
      do df = 1, 4
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do
      df = 5
      do k = 1, nlayers + 1
        write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                    ' layer ', k, ':', field(map_3d(df)+k-1)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end do

    case default
      do df = 1, ndf_3d
        do k = 1, nlayers
          write(log_scratch_space, *) 'Column ', column, 'dof ', df, &
                                      ' layer ', k, ':', field(map_3d(df)+k-1)
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
        end do
      end do

    end select
  end if

end subroutine log_column_code_r_double

end module log_column_kernel_mod
