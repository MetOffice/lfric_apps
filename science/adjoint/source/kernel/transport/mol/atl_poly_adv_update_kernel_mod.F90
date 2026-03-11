!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Line-by-line adjoint of poly_adv_update kernel, with different active variables
!>        to adj_poly_adv_update. The reconstruction field used is the linearisation state.

! @todo Kernel does not work correctly when run in parallel; fix is being worked on

module atl_poly_adv_update_kernel_mod

use argument_mod,      only : arg_type, func_type,   &
                              GH_FIELD,              &
                              GH_READ, GH_READWRITE, &
                              GH_REAL, GH_INC,       &
                              STENCIL, CROSS2D,      &
                              CELL_COLUMN,           &
                              ANY_DISCONTINUOUS_SPACE_1
use constants_mod,     only : i_def, l_def, r_tran
use fs_continuity_mod, only : W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: atl_poly_adv_update_kernel_type
  type(arg_type) :: META_ARGS(4) = (/ &
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                                 & ! advective
    arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), & ! ls_reconstruction
    arg_type(GH_FIELD, GH_REAL, GH_INC, W2),                                           & ! wind
    arg_type(GH_FIELD, GH_REAL, GH_READ, W2, STENCIL(CROSS2D))                         & ! wind_dir
  /)
  integer :: OPERATES_ON = cell_column
  contains
    procedure, nopass :: atl_poly_adv_update_code
end type atl_poly_adv_update_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: atl_poly_adv_update_code

contains

!> @brief Computes the adjoint of tangent linear horizontal advective update of a field
!> @param[in]     nlayers           Number of layers
!> @param[in,out] advective         Field containing the horizontal advective
!!                                  update: (u,v).grad_h(reconstruction)
!> @param[in]     ls_reconstruction Passive linearisation state multidata field
!!                                  containing the edge reconstruction for each cell
!!                                  if it is the upwind cell for that edge
!> @param[in]     smap_md_size      Size of the multidata stencil map in
!!                                  each direction
!> @param[in]     smap_md_max       Maximum size of the multidata stencil map
!> @param[in]     smap_md           Stencil map for the multidata fields
!> @param[in,out] wind              Active wind field
!> @param[in]     wind_dir          Wind field used to determine direction,
!!                                  equal to wind when used in gungho
!!                                  but ls_wind when used in the linear model
!> @param[in]     smap_w2_size      Size of the w2stencil map in each direction
!> @param[in]     smap_w2_max       Maximum size of the w2 stencil map
!> @param[in]     smap_w2           Stencil map for the w2 fields
!> @param[in]     ndf_wt            Number of degrees of freedom per cell
!> @param[in]     undf_wt           Number of unique degrees of freedom for the advective field
!> @param[in]     map_wt            Dofmap for the cell at the base of the column
!> @param[in]     ndf_md            Number of degrees of freedom per cell
!> @param[in]     undf_md           Number of unique degrees of freedom for the
!!                                  reconstructed field
!> @param[in]     map_md            Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2            Number of degrees of freedom per cell
!> @param[in]     undf_w2           Number of unique degrees of freedom for the wind field
!> @param[in]     map_w2            Dofmap for the cell at the base of the column
subroutine atl_poly_adv_update_code( nlayers,           &
                                     advective,         &
                                     ls_reconstruction, &
                                     smap_md_size,      &
                                     smap_md_max,       &
                                     smap_md,           &
                                     wind,              &
                                     wind_dir,          &
                                     smap_w2_size,      &
                                     smap_w2_max,       &
                                     smap_w2,           &
                                     ndf_wt,            &
                                     undf_wt,           &
                                     map_wt,            &
                                     ndf_md,            &
                                     undf_md,           &
                                     map_md,            &
                                     ndf_w2,            &
                                     undf_w2,           &
                                     map_w2 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: ndf_md
  integer(kind=i_def), intent(in)                    :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in) :: map_md
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  integer(kind=i_def), intent(in)                                  :: smap_md_max
  integer(kind=i_def), dimension(4), intent(in)                    :: smap_md_size
  integer(kind=i_def), dimension(ndf_md,smap_md_max,4), intent(in) :: smap_md
  integer(kind=i_def), intent(in)                                  :: smap_w2_max
  integer(kind=i_def), dimension(4), intent(in)                    :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2,smap_w2_max,4), intent(in) :: smap_w2

  real(kind=r_tran), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_tran), dimension(undf_w2), intent(inout) :: wind
  real(kind=r_tran), dimension(undf_md), intent(in)    :: ls_reconstruction
  real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind_dir

  ! Internal variables
  integer(kind=i_def), parameter :: nfaces = 4
  integer(kind=i_def), parameter :: w = 1
  integer(kind=i_def), parameter :: s = 2
  integer(kind=i_def), parameter :: e = 3
  integer(kind=i_def), parameter :: n = 4

  integer(kind=i_def) :: k
  integer(kind=i_def) :: df
  integer(kind=i_def) :: ijp
  integer(kind=i_def) :: df1
  integer(kind=i_def) :: df2

  integer(kind=i_def), dimension(4)      :: direction_dofs
  real(kind=r_tran)                      :: direction
  real(kind=r_tran), dimension(nfaces)   :: v_dot_n
  integer(kind=i_def), dimension(nfaces) :: opposite
  logical(kind=l_def), dimension(nfaces) :: missing_neighbour

  real(kind=r_tran), dimension(2,0:nlayers) :: uv_dir
  real(kind=r_tran), dimension(4,0:nlayers) :: ls_tracer
  real(kind=r_tran), dimension(2,0:nlayers) :: uv
  real(kind=r_tran)                         :: ls_dtdx
  real(kind=r_tran)                         :: ls_dtdy


  v_dot_n = (/ -1.0_r_tran, 1.0_r_tran, 1.0_r_tran, -1.0_r_tran /)
  opposite(:) = -1
  missing_neighbour(:) = .false.

  do df = 1, nfaces, 1
    df1 = map_w2(df)
    if (smap_w2_size(df) > 1) then
      do df2 = 1, nfaces, 1
        if (smap_w2(df2,2,df) == df1) then
          opposite(df) = df2
        end if
      enddo
    else
      opposite(df) = df
      missing_neighbour(df) = .true.
    end if
  enddo

  k = 0
  uv_dir(1,k) = 0.25_r_tran * wind_dir(map_w2(1)) + 0.25_r_tran * wind_dir(map_w2(3))
  uv_dir(2,k) = 0.25_r_tran * wind_dir(map_w2(2)) + 0.25_r_tran * wind_dir(map_w2(4))
  do k = 1, nlayers - 1, 1
    uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2(1)) + 0.25_r_tran * wind_dir(k + map_w2(3)) + &
                  0.25_r_tran * wind_dir(k + map_w2(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2(3) - 1)
    uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2(2)) + 0.25_r_tran * wind_dir(k + map_w2(4)) + &
                  0.25_r_tran * wind_dir(k + map_w2(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2(4) - 1)
  enddo
  k = nlayers
  uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2(3) - 1)
  uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2(4) - 1)

  direction_dofs(:) = (/ 1, 2, 1, 2 /)

  do df = 1, nfaces, 1
    do k = 0, nlayers, 1
      direction = uv_dir(direction_dofs(df),k) * v_dot_n(df)
      if (direction > 0.0_r_tran .OR. missing_neighbour(df)) then
        ijp = df * nlayers + df - nlayers + map_md(1) - 1
        ls_tracer(df,k) = ls_reconstruction(ijp + k)
      else
        ijp = nlayers * opposite(df) - nlayers + opposite(df) + smap_md(1,2,df) - 1
        ls_tracer(df,k) = ls_reconstruction(ijp + k)
      end if
    enddo
  enddo

  do k = nlayers, 0, -1
    ls_dtdx = ls_tracer(e,k) - ls_tracer(w,k)
    ls_dtdy = ls_tracer(n,k) - ls_tracer(s,k)
    uv(1,k) = ls_dtdx * advective(map_wt(1) + k)
    uv(2,k) = -ls_dtdy * advective(map_wt(1) + k)
    advective(map_wt(1) + k) = 0.0_r_tran
  enddo

  k = nlayers
  wind(k + map_w2(2) - 1) = wind(k + map_w2(2) - 1) + 0.25_r_tran * uv(2,k)
  wind(k + map_w2(4) - 1) = wind(k + map_w2(4) - 1) + 0.25_r_tran * uv(2,k)
  wind(k + map_w2(1) - 1) = wind(k + map_w2(1) - 1) + 0.25_r_tran * uv(1,k)
  wind(k + map_w2(3) - 1) = wind(k + map_w2(3) - 1) + 0.25_r_tran * uv(1,k)

  do k = nlayers - 1, 1, -1
    wind(k + map_w2(2)) = wind(k + map_w2(2)) + 0.25_r_tran * uv(2,k)
    wind(k + map_w2(4)) = wind(k + map_w2(4)) + 0.25_r_tran * uv(2,k)
    wind(k + map_w2(2) - 1) = wind(k + map_w2(2) - 1) + 0.25_r_tran * uv(2,k)
    wind(k + map_w2(4) - 1) = wind(k + map_w2(4) - 1) + 0.25_r_tran * uv(2,k)
    wind(k + map_w2(1)) = wind(k + map_w2(1)) + 0.25_r_tran * uv(1,k)
    wind(k + map_w2(3)) = wind(k + map_w2(3)) + 0.25_r_tran * uv(1,k)
    wind(k + map_w2(1) - 1) = wind(k + map_w2(1) - 1) + 0.25_r_tran * uv(1,k)
    wind(k + map_w2(3) - 1) = wind(k + map_w2(3) - 1) + 0.25_r_tran * uv(1,k)
  enddo

  k = 0
  wind(map_w2(2)) = wind(map_w2(2)) + 0.25_r_tran * uv(2,k)
  wind(map_w2(4)) = wind(map_w2(4)) + 0.25_r_tran * uv(2,k)
  wind(map_w2(1)) = wind(map_w2(1)) + 0.25_r_tran * uv(1,k)
  wind(map_w2(3)) = wind(map_w2(3)) + 0.25_r_tran * uv(1,k)

end subroutine atl_poly_adv_update_code

end module atl_poly_adv_update_kernel_mod
