!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Map multidata fields with different ndata numbers.
!> @details This is used for lfric2lfric reconfiguration, most likely
!!          for regridding 2D fields on snow layers and tiles, which
!!          are of size number_snow_layers * number_land_tiles. But can also
!!          be used for fields on tiles (where the number of snow_layers is
!!          set to 1).
module lfric2lfric_multidata_kernel_mod

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

  type, public, extends(kernel_type) :: lfric2lfric_multidata_type
    private
    type(arg_type), dimension(7) :: meta_args = (/                &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,  ANY_SPACE_1), & ! field_dst
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,   ANY_SPACE_2), & ! field_src
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! n_snow_layers_dst
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! n_land_tiles_dst
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! n_snow_layers_src
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! n_land_tiles_src
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                 & ! ntiles_first 
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: lfric2lfric_multidata_code
  end type lfric2lfric_multidata_type

  public :: lfric2lfric_multidata_code

contains
  
!> @brief Map a multidata (snow) field to another multidata (snow) field
!! @param[in]     nlayers      Number of layers
!! @param[in,out] field_dst    Output (destination) multidata field
!! @param[in]     field_src    Input (source) multidata field
!! @param[in]     n_snow_layers_dst Number of snow layers in destination field
!! @param[in]     n_tiles_dst  Number of land tiles in destination fields
!! @param[in]     n_snow_layers_src Number of snow layers in source field
!! @param[in]     n_tiles_src  Number of land tiles in source field
!! @param[in]     ntiles_first Data layout (tiles- or layer-first)
!! @param[in]     ndf_dst      Number of dofs per cell for field_dst
!! @param[in]     undf_dst     Total number of dofs per cell for field_dst
!! @param[in]     map_dst      Cell dofmap for field_dst
!! @param[in]     ndf_src      Number of dofs per cell for field_src
!! @param[in]     undf_src     Total number of dofs per cell for field_src
!! @param[in]     map_src      Cell dofmap for field_src
subroutine lfric2lfric_multidata_code(nlayers,           &
                                      field_dst,         &
                                      field_src,         &
                                      n_snow_layers_dst, &
                                      n_tiles_dst,       &
                                      n_snow_layers_src, &
                                      n_tiles_src,       &
                                      ntiles_first,      &
                                      ndf_dst,           &
                                      undf_dst,          &
                                      map_dst,           &
                                      ndf_src,           &
                                      undf_src,          &
                                      map_src)

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_src
  integer(kind=i_def), intent(in) :: undf_src
  integer(kind=i_def), intent(in) :: map_src(ndf_src)

  integer(kind=i_def), intent(in) :: ndf_dst
  integer(kind=i_def), intent(in) :: undf_dst
  integer(kind=i_def), intent(in) :: map_dst(ndf_dst)

  integer(kind=i_def), intent(in) :: n_snow_layers_dst, n_tiles_dst
  integer(kind=i_def), intent(in) :: n_snow_layers_src, n_tiles_src

  logical(kind=l_def), intent(in) :: ntiles_first

  real(kind=r_def),    intent(in) :: field_src(undf_src)
  real(kind=r_def),    intent(out):: field_dst(undf_dst)

  integer(kind=i_def) :: df, i_dst, i_src, k

  if ( ntiles_first ) then
    ! Data layout is tiles T first (before layers L):
    ! L1T1, L1T2, L1T3, ..., L2T1, L2T2, L2T3, ... 
    do k = 0, n_snow_layers_dst - 1
      i_dst = map_dst(1) + k * n_tiles_dst
      i_src = map_src(1) + k * n_tiles_src
      do df = 0, n_tiles_dst -1
        field_dst( i_dst + df ) = field_src( i_src + multidata_map(df + 1) -1 )
      end do
    end do
   
  else
    ! Data layout is layers L first (before tiles T):
    ! L1T1, L2T1, ..., L1T2, L2T2, ..., L1T3, L2T3, ... 
    do df = 0, n_tiles_dst - 1
      i_dst = map_dst(1) + df * n_snow_layers_dst
      i_src = map_src(1) + ( multidata_map(df + 1) - 1 ) * n_snow_layers_src
      do k = 0, n_snow_layers_dst - 1
         field_dst( i_dst + k ) = field_src( i_src + k )
      end do
    end do
  end if

end subroutine lfric2lfric_multidata_code

end module lfric2lfric_multidata_kernel_mod
