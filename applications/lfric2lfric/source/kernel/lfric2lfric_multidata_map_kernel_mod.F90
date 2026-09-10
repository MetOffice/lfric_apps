!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Map multidata fields with different ndata numbers.
!> @details This is used for lfric2lfric reconfiguration, most likely
!!        for regridding 2D fields on surface tiles. The multidata_map
!!        is read from configuration and is an integer array (of size of
!!        the destination ndata) containing a list of the array positions
!!        from the source data.
module lfric2lfric_multidata_map_kernel_mod

  use argument_mod,      only: arg_type,        &
                               GH_FIELD,        &
                               GH_SCALAR,       &
                               GH_REAL,         &
                               GH_INTEGER,      &
                               GH_READ,         &
                               GH_WRITE,        &
                               GH_LOGICAL,      &
                               CELL_COLUMN,     &
                               ANY_SPACE_1,     &
                               ANY_SPACE_2
  use constants_mod,     only: i_def, r_def, l_def
  use kernel_mod,        only: kernel_type
  use dst_jules_surface_types_config_mod, &
                         only: multidata_map

  implicit none

  private

  type, public, extends(kernel_type) :: lfric2lfric_multidata_map_type
    private
    type(arg_type), dimension(5) :: meta_args = (/                &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,  ANY_SPACE_1), & ! field_dst
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! ndata_dst
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_SPACE_2),   & ! field_src
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! ndata_src
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                 & ! ndata_first 
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: lfric2lfric_multidata_map_code
  end type lfric2lfric_multidata_map_type

  public :: lfric2lfric_multidata_map_code

contains
  
!> @brief Map a multidata field to another multidata field
!! @param[in]     nlayers     Number of layers
!! @param[in,out] field_dst   Output (destination) multidata field
!! @param[in]     ndata_dst   Number of data points in field_dst
!! @param[in]     field_src   Input (source) multidata field
!! @param[in]     ndata_src   Number of data points in field_src
!! @param[in]     ndata_first Data layout (data- or layer-first)
!! @param[in]     ndf_dst     Number of dofs per cell for field_dst
!! @param[in]     undf_dst    Total number of dofs per cell for field_dst
!! @param[in]     map_dst     Cell dofmap for field_dst
!! @param[in]     ndf_src     Number of dofs per cell for field_src
!! @param[in]     undf_src    Total number of dofs per cell for field_src
!! @param[in]     map_src     Cell dofmap for field_src
subroutine lfric2lfric_multidata_map_code(nlayers,     &
                                          field_dst,   &
                                          ndata_dst,   &
                                          field_src,   &
                                          ndata_src,   &
                                          ndata_first, &
                                          ndf_dst,     &
                                          undf_dst,    &
                                          map_dst,     &
                                          ndf_src,     &
                                          undf_src,    &
                                          map_src)

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndata_src
  integer(kind=i_def), intent(in) :: ndf_src
  integer(kind=i_def), intent(in) :: undf_src
  integer(kind=i_def), intent(in) :: map_src(ndf_src)

  integer(kind=i_def), intent(in) :: ndata_dst
  integer(kind=i_def), intent(in) :: ndf_dst
  integer(kind=i_def), intent(in) :: undf_dst
  integer(kind=i_def), intent(in) :: map_dst(ndf_dst)

  logical(kind=l_def), intent(in) :: ndata_first

  real(kind=r_def),    intent(in)  :: field_src(undf_src)
  real(kind=r_def),    intent(out) :: field_dst(undf_dst)

  integer(kind=i_def) :: df, i_dst, i_src, k

  if (ndata_first) then
    do k = 0, nlayers - 1
      i_dst = map_dst(1) + k * ndata_dst
      i_src = map_src(1) + k * ndata_src
      do df = 0, ndata_dst -1
        field_dst( i_dst + df ) = field_src( i_src + multidata_map(df+1) -1 )
      end do
    end do
  else
    do df = 0, ndata_dst - 1
      i_dst = map_dst(1) + df * nlayers
      i_src = map_src(1) + (multidata_map(df+1) - 1) * nlayers
      do k =0, nlayers-1
         field_dst( i_dst + k ) = field_src( i_src + k )
      end do
    end do
  end if

end subroutine lfric2lfric_multidata_map_code

end module lfric2lfric_multidata_map_kernel_mod
