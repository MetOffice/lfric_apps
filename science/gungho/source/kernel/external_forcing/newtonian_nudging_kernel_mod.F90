!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs Newtonian nudging.
!> @details This kernel computes increment for performing Newtonian nudging of
!!          a field towards a reference field.
!!          Only implemented for the lowest-order elements
module newtonian_nudging_kernel_mod

  use argument_mod,              only: arg_type,                               &
                                       GH_FIELD, GH_SCALAR, GH_SCALAR_ARRAY,   &
                                       GH_REAL, GH_INTEGER,                    &
                                       GH_READ, GH_WRITE,                      &
                                       ANY_DISCONTINUOUS_SPACE_1,              &
                                       ANY_DISCONTINUOUS_SPACE_2,              &
                                       CELL_COLUMN
  use constants_mod,             only: r_def, i_def
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: newtonian_nudging_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_SCALAR_ARRAY,                                              &
                            GH_REAL,    GH_READ, 1),                           &
        arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),  &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_REAL,    GH_READ)                               &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: newtonian_nudging_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: newtonian_nudging_code

contains

!> @brief Performs Newtonian nudging
!> @param[in]     nlayers             Number of layers
!> @param[in,out] increment           The increment to compute
!> @param[in]     field_in            Field to be nudged
!> @param[in]     field_ref           Reference field to nudge towards
!> @param[in]     nvert               Number of vertical levels
!> @param[in]     vert_weights        Vertical weights for nudging
!> @param[in]     tropopause_level    Tropopause level for each column
!> @param[in]     min_level           Minimum vertical level for nudging
!> @param[in]     max_level           Maximum vertical level for nudging
!> @param[in]     time_weight         Time weight for nudging
!> @param[in]     ndf_3d              Num of DoFs per cell in 3D field
!> @param[in]     undf_3d             Num of DoFs in this partition in 3D field
!> @param[in]     map_3d              DoF map for the 3D field
!> @param[in]     ndf_2d              Num of DoFs per cell in 2D field
!> @param[in]     undf_2d             Num of DoFs in this partition in 2D field
!> @param[in]     map_2d              DoF map for the 2D field
subroutine newtonian_nudging_code(nlayers,                                     &
                                  increment,                                   &
                                  field_in,                                    &
                                  field_ref,                                   &
                                  nvert,                                       &
                                  vert_weights,                                &
                                  tropopause_level,                            &
                                  min_level,                                   &
                                  max_level,                                   &
                                  time_weight,                                 &
                                  ndf_3d,                                      &
                                  undf_3d,                                     &
                                  map_3d,                                      &
                                  ndf_2d,                                      &
                                  undf_2d,                                     &
                                  map_2d)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_2d, ndf_3d
  integer(kind=i_def), intent(in)    :: undf_2d, undf_3d
  integer(kind=i_def), intent(in)    :: nvert(1)
  integer(kind=i_def), intent(in)    :: min_level, max_level
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  integer(kind=i_def), intent(in)    :: map_3d(ndf_3d)
  real(kind=r_def),    intent(in)    :: vert_weights(nvert(1))
  real(kind=r_def),    intent(in)    :: time_weight
  integer(kind=i_def), intent(in)    :: tropopause_level(undf_2d)
  real(kind=r_def),    intent(in)    :: field_ref(undf_3d)
  real(kind=r_def),    intent(in)    :: field_in(undf_3d)
  real(kind=r_def),    intent(inout) :: increment(undf_3d)

  ! Local variables
  integer(kind=i_def) :: nl
  integer(kind=i_def) :: idx_fb, idx_ft, idx_vb, idx_vt
  integer(kind=i_def) :: min_level_loc, max_level_loc

  nl = nlayers - 2 + ndf_3d  ! nlayers for Wtheta, nlayers - 1 for W3

  ! Wtheta: Level 0 corresponds to the surface
  ! At min_level, Wtheta weights will be 0 unless taper_bottom == level_bottom
  ! At max_level, Wtheta weights will be 0 unless taper_top == level_top
  ! The tropopause level is included in the calculation if it is below max_level
  ! W3: (k_wt < k_w3 < k_wt + 1)
  ! Below min_level, W3 weights are set to 0, so we can ignore these levels
  ! Above max_level, W3 weights are set to 0, so we can ignore these levels
  ! Tropopause level is at a Wtheta level, so include W3 levels below this

  ! The following indices take the above into account:
  min_level_loc = MAX(min_level, 0)
  ! max_level, not above tropopause. Subtract 1 for W3
  max_level_loc =                                                              &
      MAX(MIN(max_level, tropopause_level(map_2d(1)) - 2 + ndf_3d, nl), 0)

  idx_fb = map_3d(1) + min_level_loc
  idx_ft = map_3d(1) + max_level_loc

  ! Index of vertical weights. These start from 1 (whether W3 or Wtheta)
  idx_vb = 1 + min_level_loc
  idx_vt = 1 + max_level_loc

  increment(idx_fb:idx_ft) = (                                                 &
    (field_ref(idx_fb:idx_ft) - field_in(idx_fb:idx_ft))                       &
    * time_weight * vert_weights(idx_vb:idx_vt)                                &
  )

end subroutine newtonian_nudging_code

end module newtonian_nudging_kernel_mod
