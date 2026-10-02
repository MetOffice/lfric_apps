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
                                 GH_FIELD, GH_SCALAR, GH_REAL,              &
                                 GH_READ, GH_WRITE,                        &
                                 CELL_COLUMN,                               &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_def, i_def
  use empty_data_mod,     only : empty_real_data
  use fs_continuity_mod,  only : Wtheta, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: dust_conc_diags_kernel_type
    private
    ! acc_ins_du, cor_ins_du, theta, exner_in_wth, height_w3,
    ! dust_conc_surface, dust_conc_2000_5000ft, p_zero, kappa, rd
    type(arg_type) :: meta_args(10) = (/                             &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),              &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                  &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                      &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                      &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                       &
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
  !>          density, with density taken as pressure over the dry gas
  !>          constant times temperature. The surface value is that
  !>          quantity on the lowest theta level. The 2000-5000 feet value
  !>          is a thickness-weighted sum of the same quantity over the
  !>          model layers whose rho-level bounds (measured from the first
  !>          rho level) fall within the 609.6-1524.0 m band, divided by
  !>          the nominal 914.4 m band depth.
  !>
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     acc_ins_du             Accumulation-mode insoluble dust mmr
  !> @param[in]     cor_ins_du             Coarse-mode insoluble dust mmr
  !> @param[in]     theta                  Potential temperature
  !> @param[in]     exner_in_wth           Exner pressure in wth space
  !> @param[in]     height_w3              Height above sea level in w3
  !> @param[in,out] dust_conc_surface      Dust concentration at the surface
  !> @param[in,out] dust_conc_2000_5000ft  Dust concentration in the
  !>                                       2000-5000 feet layer
  !> @param[in]     p_zero                 Reference surface pressure
  !> @param[in]     kappa                  Ratio of rd to cp
  !> @param[in]     rd                     Gas constant for dry air
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
                                   theta,                 &
                                   exner_in_wth,          &
                                   height_w3,             &
                                   dust_conc_surface,     &
                                   dust_conc_2000_5000ft, &
                                   p_zero,                &
                                   kappa,                 &
                                   rd,                    &
                                   ndf_wth,               &
                                   undf_wth,              &
                                   map_wth,               &
                                   ndf_w3,                &
                                   undf_w3,               &
                                   map_w3,                &
                                   ndf_2d,                &
                                   undf_2d,               &
                                   map_2d                 )

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
    real(kind=r_def), intent(in), dimension(undf_wth) :: theta
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in), dimension(undf_w3)  :: height_w3

    real(kind=r_def), intent(inout), pointer :: dust_conc_surface(:)
    real(kind=r_def), intent(inout), pointer :: dust_conc_2000_5000ft(:)

    real(kind=r_def), intent(in) :: p_zero, kappa, rd

    ! Conversion from micrograms to grams per cubic metre. The kilograms
    ! to micrograms step is applied in total_dust_conc, keeping the two
    ! conversions as separate explicit steps as the UM source does.
    real(kind=r_def), parameter :: micg_to_g = 1.0e-6_r_def

    ! 2000 and 5000 feet in metres, and the nominal depth of that band.
    real(kind=r_def), parameter :: bottom_of_layer = 609.6_r_def
    real(kind=r_def), parameter :: top_of_layer    = 1524.0_r_def
    real(kind=r_def), parameter :: layer_thickness = 914.4_r_def

    real(kind=r_def) :: conc_tot
    real(kind=r_def) :: ref_height, upper, lower, thickness, conc

    integer(kind=i_def) :: k

    if ( .not. associated(dust_conc_surface, empty_real_data) ) then
      conc_tot = total_dust_conc( acc_ins_du(map_wth(1)+1),   &
                                  cor_ins_du(map_wth(1)+1),   &
                                  theta(map_wth(1)+1),        &
                                  exner_in_wth(map_wth(1)+1), &
                                  p_zero, kappa, rd )
      dust_conc_surface(map_2d(1)) = conc_tot * micg_to_g
    end if

    if ( .not. associated(dust_conc_2000_5000ft, empty_real_data) ) then
      ref_height = height_w3(map_w3(1))
      conc = 0.0_r_def

      do k = 2, nlayers-1
        conc_tot = total_dust_conc( acc_ins_du(map_wth(1)+k),   &
                                    cor_ins_du(map_wth(1)+k),   &
                                    theta(map_wth(1)+k),        &
                                    exner_in_wth(map_wth(1)+k), &
                                    p_zero, kappa, rd )

        upper = height_w3(map_w3(1)+k)   - ref_height
        lower = height_w3(map_w3(1)+k-1) - ref_height

        if (upper < bottom_of_layer) then
          thickness = 0.0_r_def
        else if (lower > top_of_layer) then
          thickness = 0.0_r_def
        else
          thickness = height_w3(map_w3(1)+k) - height_w3(map_w3(1)+k-1)
        end if

        conc = conc + thickness * conc_tot
      end do

      dust_conc_2000_5000ft(map_2d(1)) = conc * micg_to_g / layer_thickness
    end if

  end subroutine dust_conc_diags_code

  !> @brief Total dust concentration on a single level, in micrograms per
  !>        cubic metre.
  !> @details Mirrors the UM's total dust concentration calculation: the
  !>          summed insoluble dust mass mixing ratios times air density,
  !>          scaled from kilograms to micrograms. Air density is taken as
  !>          pressure over the dry gas constant times temperature, as the
  !>          UM does, rather than read from a density field. Both
  !>          diagnostics call this, so the single UM formula keeps a
  !>          single definition here and the two cannot drift apart.
  !>
  !> @param[in] acc_ins  Accumulation-mode insoluble dust mmr
  !> @param[in] cor_ins  Coarse-mode insoluble dust mmr
  !> @param[in] theta    Potential temperature
  !> @param[in] exner    Exner pressure
  !> @param[in] p_zero   Reference surface pressure
  !> @param[in] kappa    Ratio of rd to cp
  !> @param[in] rd       Gas constant for dry air
  !> @return    conc     Total dust concentration in micrograms per m3
  pure function total_dust_conc( acc_ins, cor_ins, theta, exner, &
                                 p_zero, kappa, rd ) result(conc)

    implicit none

    real(kind=r_def), intent(in) :: acc_ins, cor_ins
    real(kind=r_def), intent(in) :: theta, exner
    real(kind=r_def), intent(in) :: p_zero, kappa, rd

    real(kind=r_def) :: conc

    ! Conversion from mass mixing ratio times air density, in kilograms
    ! per cubic metre, to micrograms per cubic metre.
    real(kind=r_def), parameter :: kg_to_micg = 1.0e9_r_def

    real(kind=r_def) :: pressure, temperature, rho_air

    pressure    = p_zero * exner**(1.0_r_def/kappa)
    temperature = theta * exner
    rho_air     = pressure / (rd * temperature)
    conc        = (acc_ins + cor_ins) * rho_air * kg_to_micg

  end function total_dust_conc

end module dust_conc_diags_kernel_mod
