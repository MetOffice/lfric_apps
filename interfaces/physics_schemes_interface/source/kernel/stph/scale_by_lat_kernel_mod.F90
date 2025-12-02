!-------------------------------------------------------------------------------
!(c) Crown copyright 2022 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Invert field in Northern hemisphere
module scale_by_lat_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,          &
                               GH_REAL, GH_WRITE, GH_READ,  &
                               ANY_DISCONTINUOUS_SPACE_1,   &
                               ANY_DISCONTINUOUS_SPACE_2,   &
                               CELL_COLUMN
  use fs_continuity_mod, only: W3, Wtheta
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: scale_by_lat_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field_out
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! field_in
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2)  & ! latitude
    /)
    integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: scale_by_lat_code
  end type scale_by_lat_kernel_type

  public scale_by_lat_code
contains

  !> @brief Invert field in northern hemisphere
  !> @param[in]    nlayers     The number of layers
  !> @param[inout] field_out   field_out
  !> @param[in]    field_in    field_in
  !> @param[in]    latitude    2D field of latitudes having the same ndf as field
  !> @param[in]    ndf_f       Number of degrees of freedom per cell for field
  !> @param[in]    undf_f      Number of unique degrees of freedom for field
  !> @param[in]    map_f       Dofmap for the cell at the base of the column for field
  !> @param[in]    ndf_2d      Number of degrees of freedom per cell for 2D space
  !> @param[in]    undf_2d     Number of unique degrees of freedom for 2D space
  !> @param[in]    map_2d      Dofmap for the cell at the base of the column for 2D space

  subroutine scale_by_lat_code(nlayers,   &
                               field_out, &
                               field_in,  &
                               latitude,  &
                               ndf_f,     &
                               undf_f,    &
                               map_f,     &
                               ndf_2d,    &
                               undf_2d,   &
                               map_2d     &
                               )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_f, ndf_2d
    integer(kind=i_def), intent(in) :: undf_f, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_f)  ::  map_f
    integer(kind=i_def), intent(in), dimension(ndf_2d) ::  map_2d

    !fields
    real(kind=r_def), intent(inout), dimension(undf_f)  :: field_out
    real(kind=r_def), intent(in),    dimension(undf_f)  :: field_in
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude

    integer(kind=i_def) :: k, df

    ! Invert field where latitude > 0
    do df = 1, ndf_f
      do k = 0, nlayers-1
        field_out(map_f(df) + k) = field_in(map_f(df) + k ) * &
             sign(1.0_r_def, latitude(map_2d(df))) * (-1.0_r_def)
      end do
    end do

  end subroutine scale_by_lat_code

end module scale_by_lat_kernel_mod
