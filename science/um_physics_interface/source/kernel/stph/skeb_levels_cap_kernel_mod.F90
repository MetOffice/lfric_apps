!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Apply vertical smoothing (linear ramp/down) for SPT increments
module skeb_levels_cap_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, &
                               GH_WRITE,  GH_READ,  &
                               GH_REAL, CELL_COLUMN, &
                               GH_SCALAR, GH_INTEGER

  use fs_continuity_mod, only: W2
  use constants_mod,     only: r_def, i_def, r_um
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: skeb_levels_cap_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/              &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),  & ! field_out
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),   & ! field_in
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),  & ! skeb_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),  & ! skeb_level_top
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ )   & ! level2km
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: skeb_levels_cap_code
  end type skeb_levels_cap_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public skeb_levels_cap_code
contains

  !> @brief Apply vertical ramp up/down on SPT buffer verfical levels
  !> @details Apply the vertical ramp up (down) over SPT increments in the
  !>          SPT vertical levels from spt_level_bottom to spt_level_begin_tapering_bottom
  !>          (from spt_level_begin_tapering_top to spt_level_top)
  !> @param[in]     nlayers       The number of layers
  !> @param[in,out] field_out     Field with scaled du wind perturbations
  !> @param[in]     field_in      Field with du wind perturbations
  !> @param[in]     skeb_level_bottom Bottom SKEB level
  !> @param[in]     skeb_level_top    Top SKEB level
  !> @param[in]     level2km      Lowest level above 2km
  !> @param[in]     ndf_w2        Number of DOFs per cell for W2 space
  !> @param[in]     undf_w2       Number of unique DOFs for W2 space
  !> @param[in]     map_w2        dofmap for the cell at the base of the column for W2 space

  subroutine skeb_levels_cap_code(nlayers,           &
                                  field_out,         &
                                  field_in,          &
                                  skeb_level_bottom, &
                                  skeb_level_top,    &
                                  level2km,          &
                                  ndf_w2,            &
                                  undf_w2,           &
                                  map_w2             &
                                  )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w2
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in), dimension(ndf_w2)  :: map_w2

    ! field with perturbation
    real(kind=r_def), intent(in), dimension(undf_w2) :: field_in
    real(kind=r_def), intent(inout), dimension(undf_w2) :: field_out

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top, level2km

    real(kind=r_def) :: lev_fraction,  logscale
    integer(kind=i_def) :: k, df

    ! --------------------------------------------------------------------
    !  Reduce wind increments to 0 at level=1 from level~2km using LOG10
    !  This avoids creating large values near the top of the boundary lyr
    ! --------------------------------------------------------------------
    logscale = 10.0_r_def/level2km

    do df = 1, ndf_w2

      ! Zero any increments below lowest skeb level
      do k = 1, skeb_level_bottom-1
        field_out(map_w2(df)+k-1) =  0.0_r_def
      end do

      ! Scale levels below 2km
      if (skeb_level_bottom < level2km) then
        do k = skeb_level_bottom, level2km
          lev_fraction = max(0.0_r_def , log10(logscale*k))   ! log decrease to zero
          field_out(map_w2(df)+k-1) = field_in(map_w2(df)+k-1)* lev_fraction
        end do
      end if

      ! Ramp off backscatter (linear) at top 3 levels of SKEB2 range
      do k = skeb_level_top-2, skeb_level_top
        lev_fraction = 0.25_r_def * (skeb_level_top - k + 1)
        field_out(map_w2(df)+k-1) =  field_in(map_w2(df)+k-1)* lev_fraction
      end do

      ! Zero any above top level
      do k = skeb_level_top+1, nlayers
        field_out(map_w2(df)+k-1) =  0.0_r_def
      end do

    end do ! loop over dofs

  end subroutine skeb_levels_cap_code

end module skeb_levels_cap_kernel_mod
