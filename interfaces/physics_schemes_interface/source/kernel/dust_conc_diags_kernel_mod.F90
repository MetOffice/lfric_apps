!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Some of the content of this file has been produced with the assistance of
!> Met Office Claude Code Enterprise.
!> @brief Interface to the dust concentration diagnostics.

module dust_conc_diags_kernel_mod

  use argument_mod,       only : arg_type,                                  &
                                 GH_FIELD, GH_REAL,                         &
                                 GH_READ, GH_WRITE,                        &
                                 CELL_COLUMN,                               &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_def, i_def, rmdi
  use empty_data_mod,     only : empty_real_data
  use fs_continuity_mod,  only : Wtheta, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: dust_conc_diags_kernel_type
    private
    ! acc_ins_du, cor_ins_du, air_density, height_w3, height_wth,
    ! dust_conc_surface, dust_conc_2000_5000ft
    type(arg_type) :: meta_args(7) = (/                              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: dust_conc_diags_code
  end type

  public :: dust_conc_diags_code

contains

  !> @brief Calculate total dust concentration at the surface and averaged
  !>        over the 2000-5000 feet layer.
  !> @details Total dust concentration is the summed accumulation- and
  !>          coarse-mode insoluble dust mass mixing ratio times air
  !>          density. The surface value is that quantity on the first
  !>          theta level above the ground. The 2000-5000 feet value is a
  !>          thickness-weighted mean of the same quantity over the model
  !>          layers that reach into the band, excluding the lowest and
  !>          the topmost as the UM does. Each is counted at its full
  !>          thickness, so a layer straddling an edge is not clipped.
  !>          Layer bounds are the rho-level heights, measured from the
  !>          surface. Where no layer reaches the band, the value is
  !>          missing data.
  !>
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     acc_ins_du             Accumulation-mode insoluble dust mmr
  !> @param[in]     cor_ins_du             Coarse-mode insoluble dust mmr
  !> @param[in]     air_density            Air density in wth space
  !> @param[in]     height_w3              Height above sea level in w3
  !> @param[in]     height_wth             Height above sea level in wth
  !> @param[in,out] dust_conc_surface      Dust concentration at the surface
  !> @param[in,out] dust_conc_2000_5000ft  Dust concentration in the
  !>                                       2000-5000 feet layer
  !> @param[in]     ndf_wth                Number of degrees of freedom per
  !>                                       cell for potential temperature space
  !> @param[in]     undf_wth               Number of unique degrees of
  !>                                       freedom for potential temperature
  !>                                       space
  !> @param[in]     map_wth                Dofmap for the cell at the base
  !>                                       of the column for potential
  !>                                       temperature space
  !> @param[in]     ndf_w3                 Number of degrees of freedom per
  !>                                       cell for density space
  !> @param[in]     undf_w3                Number of unique degrees of
  !>                                       freedom for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base
  !>                                       of the column for density space
  !> @param[in]     ndf_2d                 Number of degrees of freedom per
  !>                                       cell for 2D fields
  !> @param[in]     undf_2d                Number of unique degrees of
  !>                                       freedom for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base
  !>                                       of the column for 2D fields
  subroutine dust_conc_diags_code( nlayers,               &
                                   acc_ins_du,            &
                                   cor_ins_du,            &
                                   air_density,           &
                                   height_w3,             &
                                   height_wth,            &
                                   dust_conc_surface,     &
                                   dust_conc_2000_5000ft, &
                                   ndf_wth,               &
                                   undf_wth,              &
                                   map_wth,               &
                                   ndf_w3,                &
                                   undf_w3,               &
                                   map_w3,                &
                                   ndf_2d,                &
                                   undf_2d,               &
                                   map_2d                 )

    use science_conversions_mod, only: feet_to_metres

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    real(kind=r_def), intent(in), dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(in), dimension(undf_wth) :: cor_ins_du
    real(kind=r_def), intent(in), dimension(undf_wth) :: air_density
    real(kind=r_def), intent(in), dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth

    real(kind=r_def), intent(inout), pointer :: dust_conc_surface(:)
    real(kind=r_def), intent(inout), pointer :: dust_conc_2000_5000ft(:)

    ! The 2000-5000 feet band.
    real(kind=r_def), parameter :: bottom_of_layer =                   &
                                     2000.0_r_def * feet_to_metres
    real(kind=r_def), parameter :: top_of_layer    =                   &
                                     5000.0_r_def * feet_to_metres

    real(kind=r_def) :: conc_tot
    real(kind=r_def) :: ref_height, upper, lower, thickness, conc
    real(kind=r_def) :: band_depth

    integer(kind=i_def) :: k

    if ( .not. associated(dust_conc_surface, empty_real_data) ) then
      conc_tot = total_dust_conc( acc_ins_du(map_wth(1)+1),  &
                                  cor_ins_du(map_wth(1)+1),  &
                                  air_density(map_wth(1)+1) )
      dust_conc_surface(map_2d(1)) = conc_tot
    end if

    if ( .not. associated(dust_conc_2000_5000ft, empty_real_data) ) then
      ! Measured from the surface, the bottom Wtheta level. The UM measures
      ! from its first rho level instead, which lifts its band a few metres
      ! off the ground; that is treated as a UM error, not reproduced.
      ref_height = height_wth(map_wth(1))
      conc = 0.0_r_def
      band_depth = 0.0_r_def

      ! Layer bounds rise monotonically with k, so the levels below the
      ! band can be skipped and the loop left for good once past its top.
      ! A layer straddling either edge still counts at its full thickness,
      ! as the UM counts it.
      do k = 2, nlayers-1
        upper = height_w3(map_w3(1)+k)   - ref_height
        lower = height_w3(map_w3(1)+k-1) - ref_height

        if (upper < bottom_of_layer) cycle
        if (lower > top_of_layer) exit

        conc_tot = total_dust_conc( acc_ins_du(map_wth(1)+k),  &
                                    cor_ins_du(map_wth(1)+k),  &
                                    air_density(map_wth(1)+k) )

        ! Differenced without the datum, as the UM differences it: the
        ! two forms are not bit-identical once ref_height is large.
        thickness = height_w3(map_w3(1)+k) - height_w3(map_w3(1)+k-1)

        conc = conc + thickness * conc_tot
        band_depth = band_depth + thickness
      end do

      ! Divided by the summed depth of the layers counted, making this
      ! their thickness-weighted mean. The UM divides by the nominal band
      ! depth instead, which overstates the mean whenever a layer straddles
      ! an edge.
      if (band_depth > 0.0_r_def) then
        dust_conc_2000_5000ft(map_2d(1)) = conc / band_depth
      else
        dust_conc_2000_5000ft(map_2d(1)) = rmdi
      end if
    end if

  end subroutine dust_conc_diags_code

  !> @brief Total dust concentration on a single level, in grams per
  !>        cubic metre.
  !> @details Mirrors the UM's total dust concentration calculation: the
  !>          summed insoluble dust mass mixing ratios times air density.
  !>          The UM derives that density as pressure over the dry gas
  !>          constant times temperature; here it is read from a density
  !>          field instead, so that which density is meant is explicit at
  !>          the call site.
  !>
  !>          The UM can instead weight the surface value by this
  !>          timestep's dust emission, when its emission-weighting
  !>          namelist value is set. That branch is not ported; this is
  !>          the UM's default, unweighted result.
  !>
  !> @param[in] acc_ins  Accumulation-mode insoluble dust mmr
  !> @param[in] cor_ins  Coarse-mode insoluble dust mmr
  !> @param[in] rho_air  Air density
  !> @return    conc     Total dust concentration in grams per m3
  pure function total_dust_conc( acc_ins, cor_ins, rho_air ) result(conc)

    implicit none

    real(kind=r_def), intent(in) :: acc_ins, cor_ins
    real(kind=r_def), intent(in) :: rho_air

    real(kind=r_def) :: conc

    ! Conversion from mass mixing ratio times air density, in kilograms
    ! per cubic metre, to the grams per cubic metre both diagnostics are
    ! output in. The UM reaches the same units in two steps, scaling to
    ! micrograms where it builds the concentration field and back again
    ! where it writes the diagnostic.
    real(kind=r_def), parameter :: kg_to_g = 1.0e3_r_def

    conc = (acc_ins + cor_ins) * rho_air * kg_to_g

  end function total_dust_conc

end module dust_conc_diags_kernel_mod
