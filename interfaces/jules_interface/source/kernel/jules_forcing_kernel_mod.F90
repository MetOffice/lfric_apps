!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies JULES site forcing data to the fields used to drive JULES
!>
module jules_forcing_kernel_mod

  use argument_mod,       only : arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 GH_READ, GH_WRITE,         &
                                 GH_READWRITE,              &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 CELL_COLUMN
  use constants_mod,      only : i_def, r_def
  use fs_continuity_mod,  only : W3, Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_forcing_kernel_type
    private
    type(arg_type) :: meta_args(19) = (/                                      &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                    & ! theta
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                    & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                        & ! u_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                        & ! v_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                        & ! wetrho_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_wind
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_tair
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_qair
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_psurf
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_lwdown
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_swdown
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! forcing_rainf
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! ls_rain
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1)  & ! ls_snow
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: jules_forcing_code
  end type

  public :: jules_forcing_code

  !> Air temperature at or below which precipitation falls as snow (K),
  !> as the default t_for_snow in standalone JULES
  real(kind=r_def), parameter :: t_for_snow = 274.0_r_def

  !> Fraction of the downward SW flux that is photosynthetically active,
  !> as used by standalone JULES
  real(kind=r_def), parameter :: par_fraction = 0.5_r_def

contains

  !> @brief Applies JULES site forcing data to the fields used to drive JULES
  !> @details The near-surface state seen by JULES is set from the forcing:
  !>          exner at the surface gives the forcing surface pressure and
  !>          exner above it is hydrostatic at the forcing temperature and
  !>          humidity height, theta and m_v on the first theta level give the forcing
  !>          temperature and humidity, the wind on the first W3 level is the
  !>          forcing wind speed (v = 0) and the density on the first W3 level
  !>          is consistent with these. Downward radiation and precipitation
  !>          are set from the forcing, with the precipitation falling as snow
  !>          if the forcing temperature is at or below t_for_snow. All SW
  !>          is treated as direct and half of it as photosynthetically active.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] theta               Potential temperature (K)
  !> @param[in,out] m_v                 Vapour mixing ratio (kg kg-1)
  !> @param[in,out] exner_in_wth        Exner pressure in Wtheta
  !> @param[in,out] u_in_w3             Eastward wind in W3 (m s-1)
  !> @param[in,out] v_in_w3             Northward wind in W3 (m s-1)
  !> @param[in,out] wetrho_in_w3        Wet density in W3 (kg m-3)
  !> @param[in]     forcing_wind        Forcing wind speed (m s-1)
  !> @param[in]     forcing_tair        Forcing air temperature (K)
  !> @param[in]     forcing_qair        Forcing specific humidity (kg kg-1)
  !> @param[in]     forcing_psurf       Forcing surface pressure (Pa)
  !> @param[in]     forcing_lwdown      Forcing downward LW flux (W m-2)
  !> @param[in]     forcing_swdown      Forcing downward SW flux (W m-2)
  !> @param[in]     forcing_rainf       Forcing precipitation rate (kg m-2 s-1)
  !> @param[out]    sw_down_surf        Downward SW flux at surface (W m-2)
  !> @param[out]    lw_down_surf        Downward LW flux at surface (W m-2)
  !> @param[out]    sw_down_blue_surf   Downward PAR flux at surface (W m-2)
  !> @param[out]    sw_direct_blue_surf Direct downward PAR flux (W m-2)
  !> @param[out]    ls_rain             Large-scale rainfall rate (kg m-2 s-1)
  !> @param[out]    ls_snow             Large-scale snowfall rate (kg m-2 s-1)
  !> @param[in]     ndf_wth             Number of DOFs per cell for Wtheta
  !> @param[in]     undf_wth            Number of unique DOFs for Wtheta
  !> @param[in]     map_wth             Dofmap for the cell at the base of the column for Wtheta
  !> @param[in]     ndf_w3              Number of DOFs per cell for W3
  !> @param[in]     undf_w3             Number of unique DOFs for W3
  !> @param[in]     map_w3              Dofmap for the cell at the base of the column for W3
  !> @param[in]     ndf_2d              Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d             Number of unique DOFs for 2D fields
  !> @param[in]     map_2d              Dofmap for the cell at the base of the column for 2D fields
  subroutine jules_forcing_code( nlayers,                                &
                                 theta, m_v, exner_in_wth,               &
                                 u_in_w3, v_in_w3, wetrho_in_w3,         &
                                 forcing_wind, forcing_tair,             &
                                 forcing_qair, forcing_psurf,            &
                                 forcing_lwdown, forcing_swdown,         &
                                 forcing_rainf,                          &
                                 sw_down_surf, lw_down_surf,             &
                                 sw_down_blue_surf, sw_direct_blue_surf, &
                                 ls_rain, ls_snow,                       &
                                 ndf_wth, undf_wth, map_wth,             &
                                 ndf_w3, undf_w3, map_w3,                &
                                 ndf_2d, undf_2d, map_2d )

    use planet_constants_mod,      only: p_zero, kappa, r, c_virtual, g
    use initialization_config_mod, only: jules_forcing_z1_tq

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

    real(kind=r_def), intent(inout) :: theta(undf_wth)
    real(kind=r_def), intent(inout) :: m_v(undf_wth)
    real(kind=r_def), intent(inout) :: exner_in_wth(undf_wth)
    real(kind=r_def), intent(inout) :: u_in_w3(undf_w3)
    real(kind=r_def), intent(inout) :: v_in_w3(undf_w3)
    real(kind=r_def), intent(inout) :: wetrho_in_w3(undf_w3)
    real(kind=r_def), intent(in)    :: forcing_wind(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_tair(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_qair(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_psurf(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_lwdown(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_swdown(undf_2d)
    real(kind=r_def), intent(in)    :: forcing_rainf(undf_2d)
    real(kind=r_def), intent(inout) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(inout) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(inout) :: sw_down_blue_surf(undf_2d)
    real(kind=r_def), intent(inout) :: sw_direct_blue_surf(undf_2d)
    real(kind=r_def), intent(inout) :: ls_rain(undf_2d)
    real(kind=r_def), intent(inout) :: ls_snow(undf_2d)

    ! Local variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: tair, qair, tv, pressure_1

    tair = forcing_tair(map_2d(1))
    qair = max(forcing_qair(map_2d(1)), 0.0_r_def)
    tv   = tair * (1.0_r_def + c_virtual * qair)

    ! Exner at the surface from the forcing surface pressure, and above it
    ! from the hydrostatic pressure at the forcing temperature and humidity
    ! height. Only the surface and first theta levels are used by JULES.
    pressure_1 = forcing_psurf(map_2d(1)) * &
                 exp( -g * jules_forcing_z1_tq / (r * tv) )
    exner_in_wth(map_wth(1)) = (forcing_psurf(map_2d(1)) / p_zero)**kappa
    do k = 1, nlayers
      exner_in_wth(map_wth(1) + k) = (pressure_1 / p_zero)**kappa
    end do

    ! Temperature and humidity on the first theta level. The JULES kernels
    ! use m_v directly as the specific humidity of the first level.
    theta(map_wth(1) + 1) = tair / exner_in_wth(map_wth(1) + 1)
    m_v(map_wth(1) + 1)   = qair

    ! Wind on the first W3 level, with the forcing wind speed as the u component
    u_in_w3(map_w3(1)) = forcing_wind(map_2d(1))
    v_in_w3(map_w3(1)) = 0.0_r_def

    ! Density on the first W3 level consistent with the forcing state
    wetrho_in_w3(map_w3(1)) = pressure_1 / (r * tv)

    ! Downward radiation, with all SW taken as direct
    sw_down_surf(map_2d(1))        = forcing_swdown(map_2d(1))
    lw_down_surf(map_2d(1))        = forcing_lwdown(map_2d(1))
    sw_down_blue_surf(map_2d(1))   = par_fraction * forcing_swdown(map_2d(1))
    sw_direct_blue_surf(map_2d(1)) = sw_down_blue_surf(map_2d(1))

    ! Precipitation phase from the forcing temperature
    if (tair <= t_for_snow) then
      ls_rain(map_2d(1)) = 0.0_r_def
      ls_snow(map_2d(1)) = forcing_rainf(map_2d(1))
    else
      ls_rain(map_2d(1)) = forcing_rainf(map_2d(1))
      ls_snow(map_2d(1)) = 0.0_r_def
    end if

  end subroutine jules_forcing_code

end module jules_forcing_kernel_mod
