!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Logs the coordinates of each column
!>
module log_location_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, &
                               GH_WRITE, GH_INTEGER,                 &
                               ANY_DISCONTINUOUS_SPACE_1,            &
                               ANY_DISCONTINUOUS_SPACE_3,            &
                               OWNED_AND_HALO_CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type
  use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: log_location_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3)   &
    /)
    integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
  contains
    procedure, nopass :: log_location_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: log_location_code

contains

!> @details       Logs all of the field values within a single column
!> @param[in]     nlayers    Number of layers in the mesh
!> @param[in,out] dummy      Dummy field to allow psyclone to compile
!> @param[in]     longitude  Longitude field to log
!> @param[in]     latitude   Latitude field to log
!> @param[in]     panel_id   2D W3 field, whose DoFmap is the column indices
!> @param[in]     ndf        Number of dofs per cell for long/lat
!> @param[in]     undf       Total number of Dofs in this partition for long/lat
!> @param[in]     map        Dofmap for long/lat
!> @param[in]     ndf_pid    Number of dofs per cell for panel ID
!> @param[in]     undf_pid   Total number of Dofs in this partition for panel ID
!> @param[in]     map_pid    Dofmap for panel ID
subroutine log_location_code( nlayers,   &
                              dummy,     &
                              longitude, &
                              latitude,  &
                              panel_id,  &
                              ndf,       &
                              undf,      &
                              map,       &
                              ndf_pid,   &
                              undf_pid,  &
                              map_pid )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf, undf
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: map(ndf)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  real(kind=r_def),    intent(inout) :: dummy(undf)
  real(kind=r_def),    intent(in)    :: longitude(undf)
  real(kind=r_def),    intent(in)    :: latitude(undf)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: df

  select case (ndf)
  case (1)
    write(log_scratch_space, '(A,I8,A,2E16.8,I8)')                             &
        'Coordinates, column ', map_pid(1), ' : ',                             &
        longitude(map(1)), latitude(map(1)), INT(panel_id(map_pid(1)), i_def)
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  case default
    do df = 1, ndf
      write(log_scratch_space, '(A,I8,A,I8,A,2E16.8,I8)')                      &
          'Coordinates, column ', map_pid(1), 'dof ', df, ' : ',               &
          longitude(map(df)), latitude(map(df)), INT(panel_id(map_pid(1)), i_def)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do
  end select

end subroutine log_location_code

end module log_location_kernel_mod
