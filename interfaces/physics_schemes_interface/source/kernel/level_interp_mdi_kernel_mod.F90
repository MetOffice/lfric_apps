!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!> @brief Interpolate a field onto surfaces of constant value of another field.

module level_interp_mdi_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_INTEGER,                &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,        only: r_def, i_def, rmdi
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Interpolation orders, matching the UM's interpor_mod
  integer(kind=i_def), public, parameter :: interp_order_linear = 1_i_def
  integer(kind=i_def), public, parameter :: interp_order_cubic  = 3_i_def

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: level_interp_mdi_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ),                          &
!        arg_type(GH_SCALAR_ARRAY,GH_REAL, GH_READ, 1), see PSyclone issue #1312
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ)                           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: level_interp_mdi_code
  end type level_interp_mdi_kernel_type

  public :: level_interp_mdi_code

contains

  !> @details Interpolates a field onto surfaces of constant value of a second
  !!          "coordinate" field, for example potential vorticity onto theta
  !!          surfaces, or potential temperature onto a potential vorticity
  !!          surface. This is a port of the UM's vert_interp_mdi2 and shares
  !!          its conventions:
  !!
  !!          * The bracketing level is the *highest* level at which the
  !!            coordinate field is less than or equal to the target value. For
  !!            a non-monotonic coordinate such as |PV| this selects the
  !!            uppermost crossing, which is what makes the theta on PV=+/-2
  !!            diagnostic a dynamical tropopause rather than a boundary layer
  !!            artefact.
  !!          * Columns which do not span the target value are set to missing
  !!            data rather than extrapolated. This differs from
  !!            pres_interp_kernel_mod, which clamps to the end values.
  !!          * Cubic interpolation degrades to linear where the four point
  !!            stencil would run off either end of the column.
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     data_in       Input data to interpolate
  !> @param[in]     coord_in      Coordinate field at levels of the input data
  !> @param[in]     nlev_out      Number of output surfaces
  !> @param[in]     target_levs   Coordinate values of the output surfaces
  !> @param[in,out] data_out      Interpolated data
  !> @param[in]     interp_order  Interpolation order, linear or cubic
  !> @param[in]     ndf_in        Number of degrees of freedom per cell for in fields
  !> @param[in]     undf_in       Number of total degrees of freedom for in fields
  !> @param[in]     map_in        Dofmap for the cell at the base of the column for in fields
  !> @param[in]     ndf_out       Number of degrees of freedom per cell for out fields
  !> @param[in]     undf_out      Number of total degrees of freedom for out fields
  !> @param[in]     map_out       Dofmap for the cell at the base of the column for out fields
  subroutine level_interp_mdi_code(nlayers,                     &
                                   data_in,                     &
                                   coord_in,                    &
                                   nlev_out,                    &
                                   target_levs,                 &
                                   data_out,                    &
                                   interp_order,                &
                                   ndf_in, undf_in, map_in,     &
                                   ndf_out, undf_out, map_out)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nlev_out
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in), dimension(ndf_out) :: map_out
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in), dimension(ndf_in)  :: map_in

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in),    dimension(undf_in)  :: data_in
    real(kind=r_def),    intent(in),    dimension(undf_in)  :: coord_in
    real(kind=r_def),    intent(inout), dimension(undf_out) :: data_out

    ! Constants passed explicitly from algorithm
    integer(kind=i_def), intent(in) :: interp_order
    real(kind=r_def),    intent(in), dimension(nlev_out) :: target_levs

    ! Internal variables
    integer(kind=i_def) :: k, kl, level_below, top_df, base_in, base_out
    real(kind=r_def) :: target_lev
    real(kind=r_def) :: c_minus, c_here, c_plus, c_plus2

    base_in  = map_in(1)
    base_out = map_out(1)

    ! For Wtheta ndf = 2, loop k = 0, nlayers
    ! For W3 ndf = 1, loop k = 0, nlayers - 1
    top_df = nlayers + ndf_in - 2

    do kl = 1, nlev_out

      target_lev = target_levs(kl)

      ! Highest level whose coordinate value is at or below the target, leaving
      ! room for the level above. level_below = -1 flags missing data.
      level_below = -1_i_def
      if (target_lev <= coord_in(base_in+top_df) ) then
        do k = 0, top_df - 1
          if ( coord_in(base_in+k) <= target_lev ) then
            level_below = k
          end if
        end do
      end if

      if ( level_below == -1_i_def ) then

        ! Column does not span the target surface
        data_out(base_out+kl-1) = rmdi

      else if ( interp_order == interp_order_linear .or.  &
                level_below == 0_i_def .or.               &
                level_below == top_df - 1_i_def ) then

        ! Linear interpolation. Also used at both ends of the column, where the
        ! cubic stencil would be out of range.
        data_out(base_out+kl-1) =                                    &
             ( (target_lev - coord_in(base_in+level_below))          &
                 * data_in(base_in+level_below+1)                    &
             - (target_lev - coord_in(base_in+level_below+1))        &
                 * data_in(base_in+level_below) )                    &
             / ( coord_in(base_in+level_below+1)                     &
               - coord_in(base_in+level_below) )

      else

        ! Cubic interpolation over levels level_below-1 to level_below+2
        c_minus = coord_in(base_in+level_below-1)
        c_here  = coord_in(base_in+level_below)
        c_plus  = coord_in(base_in+level_below+1)
        c_plus2 = coord_in(base_in+level_below+2)

        data_out(base_out+kl-1) =                                    &
             ( (target_lev - c_here)                                 &
             * (target_lev - c_plus)                                 &
             * (target_lev - c_plus2) )                              &
             / ( (c_minus - c_here)                                  &
               * (c_minus - c_plus)                                  &
               * (c_minus - c_plus2) )                               &
             * data_in(base_in+level_below-1)                        &
             + ( (target_lev - c_minus)                              &
               * (target_lev - c_plus)                               &
               * (target_lev - c_plus2) )                            &
             / ( (c_here - c_minus)                                  &
               * (c_here - c_plus)                                   &
               * (c_here - c_plus2) )                                &
             * data_in(base_in+level_below)                          &
             + ( (target_lev - c_minus)                              &
               * (target_lev - c_here)                               &
               * (target_lev - c_plus2) )                            &
             / ( (c_plus - c_minus)                                  &
               * (c_plus - c_here)                                   &
               * (c_plus - c_plus2) )                                &
             * data_in(base_in+level_below+1)                        &
             + ( (target_lev - c_minus)                              &
               * (target_lev - c_here)                               &
               * (target_lev - c_plus) )                             &
             / ( (c_plus2 - c_minus)                                 &
               * (c_plus2 - c_here)                                  &
               * (c_plus2 - c_plus) )                                &
             * data_in(base_in+level_below+2)

      end if

    end do

  end subroutine level_interp_mdi_code

end module level_interp_mdi_kernel_mod
