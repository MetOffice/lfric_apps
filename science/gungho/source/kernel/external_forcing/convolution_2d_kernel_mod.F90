!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes weights for convolution for spectral nudging.
!> @details This kernel computes weights for performing a convolution that is
!!          used for nudging the spectrum of a field towards the spectrum of
!!          a reference field.
!!          Only implemented for the lowest-order elements
module convolution_2d_kernel_mod

  use argument_mod,              only: arg_type,                               &
                                       GH_FIELD, GH_SCALAR_ARRAY, GH_SCALAR,   &
                                       GH_REAL, GH_INTEGER, GH_READ, GH_WRITE, &
                                       CELL_COLUMN, STENCIL, REGION,           &
                                       ANY_DISCONTINUOUS_SPACE_1,              &
                                       ANY_DISCONTINUOUS_SPACE_2,              &
                                       ANY_DISCONTINUOUS_SPACE_9

  use constants_mod,             only: r_def, i_def
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: convolution_2d_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,   &
                                                            STENCIL(REGION)),  &
        arg_type(GH_SCALAR_ARRAY,                                              &
                            GH_REAL,    GH_READ, 1),                           &
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
        arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),  &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_REAL,    GH_READ)                               &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: convolution_2d_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: convolution_2d_code

contains

!> @brief Performs convolution for spectral nudging.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] filtered_field   The filtered field to compute
!> @param[in]     field_in         Field to be filtered
!> @param[in]     stencil_size     Size of the stencil used in the convolution
!> @param[in]     stencil_map      Map of the stencil used in the convolution
!> @param[in]     nvert            Number of vertical weights
!> @param[in]     vert_weights     1D array of vertical weights for nudging
!> @param[in]     conv_weights     2D weights for performing the convolution
!> @param[in]     tropopause_level Tropopause level for each column
!> @param[in]     min_level        Minimum vertical level for nudging
!> @param[in]     max_level        Maximum vertical level for nudging
!> @param[in]     time_weight      Time weight for nudging
!> @param[in]     ndf_3d           Number of DoFs per cell in 3D field
!> @param[in]     undf_3d          Number of DoFs in this partition in 3D field
!> @param[in]     map_3d           DoF map for the 3D field
!> @param[in]     ndf_2d_cw        Num of DoFs per cell in weights field
!> @param[in]     undf_2d_cw       Num DoFs in this partition in weights field
!> @param[in]     map_2d_cw        DoF map for the weights field
!> @param[in]     ndf_2d           Num of DoFs per cell in 2D trop field
!> @param[in]     undf_2d          Num DoFs in this partition in 2D trop field
!> @param[in]     map_2d           DoF map for the 2D trop field
subroutine convolution_2d_code(nlayers,                                        &
                               filtered_field,                                 &
                               field_in,                                       &
                               stencil_size,                                   &
                               stencil_map,                                    &
                               nvert,                                          &
                               vert_weights,                                   &
                               conv_weights,                                   &
                               tropopause_level,                               &
                               min_level,                                      &
                               max_level,                                      &
                               time_weight,                                    &
                               ndf_3d,                                         &
                               undf_3d,                                        &
                               map_3d,                                         &
                               ndf_2d_cw,                                      &
                               undf_2d_cw,                                     &
                               map_2d_cw,                                      &
                               ndf_2d,                                         &
                               undf_2d,                                        &
                               map_2d)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_2d, ndf_2d_cw, ndf_3d
  integer(kind=i_def), intent(in)    :: undf_2d, undf_2d_cw, undf_3d
  integer(kind=i_def), intent(in)    :: nvert(1)
  integer(kind=i_def), intent(in)    :: stencil_size
  integer(kind=i_def), intent(in)    :: min_level, max_level
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  integer(kind=i_def), intent(in)    :: map_2d_cw(ndf_2d_cw)
  integer(kind=i_def), intent(in)    :: map_3d(ndf_3d)
  integer(kind=i_def), intent(in)    :: stencil_map(ndf_3d, stencil_size)
  real(kind=r_def),    intent(in)    :: vert_weights(nvert(1))
  real(kind=r_def),    intent(in)    :: time_weight
  real(kind=r_def),    intent(in)    :: conv_weights(undf_2d_cw)
  integer(kind=i_def), intent(in)    :: tropopause_level(undf_2d)
  real(kind=r_def),    intent(in)    :: field_in(undf_3d)
  real(kind=r_def),    intent(inout) :: filtered_field(undf_3d)

  ! Local variables
  integer(kind=i_def) :: i, nl
  integer(kind=i_def) :: idx_w, idx_sb, idx_st, idx_fb, idx_ft, idx_vb, idx_vt
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

  do i = 1, stencil_size
    idx_w = map_2d_cw(1) + i - 1                        ! Multidata index
    idx_sb = stencil_map(1, i) + min_level_loc          ! Stencil index bottom
    idx_st = stencil_map(1, i) + max_level_loc          ! Stencil index top

    filtered_field(idx_fb:idx_ft) = (                                          &
        filtered_field(idx_fb:idx_ft)                                          &
        + field_in(idx_sb:idx_st) * conv_weights(idx_w)                        &
    )
  end do
  ! Apply vertical weights
  filtered_field(idx_fb:idx_ft) = (                                            &
    filtered_field(idx_fb:idx_ft) * time_weight * vert_weights(idx_vb:idx_vt)  &
  )

end subroutine convolution_2d_code

end module convolution_2d_kernel_mod
