!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code.
!-----------------------------------------------------------------------------
!> @brief CAPE, CIN and equilibrium level of a near-surface mixed-layer parcel.
module cape_cin_mixed_layer_kernel_mod

  use argument_mod,      only: arg_type,                                      &
                               GH_FIELD, GH_SCALAR, GH_REAL,                  &
                               GH_READ, GH_WRITE, CELL_COLUMN,                &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def, l_def, rmdi
  use empty_data_mod,    only: empty_real_data
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata type.
  type, public, extends(kernel_type) :: cape_cin_mixed_layer_kernel_type
    private
    type(arg_type) :: meta_args(18) = (/                                      &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! theta
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! m_v
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! height_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                          & ! exner
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                          & ! rho
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                          & ! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                          & ! height_w3
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! cape
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! cin
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! equilibrium_height
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! gravity
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! rd
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! cp
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! kappa
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! one_over_kappa
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! p_zero
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                & ! epsilon
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cape_cin_mixed_layer_code
  end type

  public :: cape_cin_mixed_layer_code

contains

  !> @brief CAPE, CIN and equilibrium level of a near-surface mixed-layer
  !>        parcel.
  !> @details The parcel takes the mass-weighted mean temperature and water
  !>          vapour of the lowest 500 m and the pressure at its middle, and
  !>          is lifted from there without mixing. The equilibrium height is
  !>          the height of the level of neutral buoyancy when the CAPE is
  !>          significant, and missing data otherwise. Outputs that were not
  !>          requested are not written.
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     theta              Potential temperature (K)
  !> @param[in]     m_v                Water vapour mixing ratio
  !> @param[in]     exner_wth          Exner pressure on Wtheta
  !> @param[in]     height_wth         Height of Wtheta levels above sea
  !>                                   level (m)
  !> @param[in]     exner_w3           Exner pressure on W3
  !> @param[in]     rho_dry            Dry density on W3 (kg m-3)
  !> @param[in]     rho_wet            Wet density on W3 (kg m-3)
  !> @param[in]     height_w3          Height of W3 levels above sea level (m)
  !> @param[in,out] cape               Convective available potential energy
  !>                                   (J/kg)
  !> @param[in,out] cin                Convective inhibition (J/kg)
  !> @param[in,out] equilibrium_height Height of the equilibrium level above
  !>                                   sea level (m)
  !> @param[in]     gravity            Gravitational acceleration (m s-2)
  !> @param[in]     rd                 Gas constant for dry air (J kg-1 K-1)
  !> @param[in]     cp                 Specific heat of dry air (J kg-1 K-1)
  !> @param[in]     kappa              Ratio of rd to cp
  !> @param[in]     one_over_kappa     Ratio of cp to rd
  !> @param[in]     p_zero             Reference pressure (Pa)
  !> @param[in]     epsilon            Ratio of gas constants of dry air and
  !>                                   water vapour
  !> @param[in]     ndf_wth            Number of DOFs per cell for Wtheta
  !> @param[in]     undf_wth           Number of unique DOFs for Wtheta
  !> @param[in]     map_wth            Dofmap for the base of the column on
  !>                                   Wtheta
  !> @param[in]     ndf_w3             Number of DOFs per cell for W3
  !> @param[in]     undf_w3            Number of unique DOFs for W3
  !> @param[in]     map_w3             Dofmap for the base of the column on W3
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d            Number of unique DOFs for 2D fields
  !> @param[in]     map_2d             Dofmap for the 2D fields
  subroutine cape_cin_mixed_layer_code( nlayers,                              &
                                        theta, m_v, exner_wth, height_wth,    &
                                        exner_w3, rho_dry, rho_wet,           &
                                        height_w3,                            &
                                        cape, cin, equilibrium_height,        &
                                        gravity, rd, cp, kappa,               &
                                        one_over_kappa, p_zero, epsilon,      &
                                        ndf_wth, undf_wth, map_wth,           &
                                        ndf_w3, undf_w3, map_w3,              &
                                        ndf_2d, undf_2d, map_2d )

    use parcel_ascent_cape_cin_mod,     only: mixed_layer_parcel_start,       &
                                              parcel_ascent_cape_cin
    use science_aviation_constants_mod, only: parcel_mixed_layer_depth,       &
                                              parcel_min_cape

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    real(kind=r_def), intent(in), dimension(undf_wth) :: theta, m_v
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in), dimension(undf_w3)  :: exner_w3
    real(kind=r_def), intent(in), dimension(undf_w3)  :: rho_dry, rho_wet
    real(kind=r_def), intent(in), dimension(undf_w3)  :: height_w3

    real(kind=r_def), pointer, intent(inout) :: cape(:)
    real(kind=r_def), pointer, intent(inout) :: cin(:)
    real(kind=r_def), pointer, intent(inout) :: equilibrium_height(:)

    real(kind=r_def), intent(in) :: gravity, rd, cp, kappa
    real(kind=r_def), intent(in) :: one_over_kappa, p_zero, epsilon

    ! Local variables
    integer(kind=i_def) :: k_start
    real(kind=r_def)    :: p_wth(0:nlayers), p_w3(0:nlayers-1)
    real(kind=r_def)    :: t_parcel, mv_parcel, p_parcel, partial_depth
    real(kind=r_def)    :: column_cape, column_cin, lnb_height

    p_wth = p_zero * exner_wth(map_wth(1):map_wth(1)+nlayers)**one_over_kappa
    p_w3  = p_zero * exner_w3(map_w3(1):map_w3(1)+nlayers-1)**one_over_kappa

    call mixed_layer_parcel_start( nlayers,                                   &
                                   theta(map_wth(1):map_wth(1)+nlayers),      &
                                   m_v(map_wth(1):map_wth(1)+nlayers),        &
                                   height_wth(map_wth(1):map_wth(1)+nlayers), &
                                   exner_w3(map_w3(1):map_w3(1)+nlayers-1),   &
                                   p_w3, p_wth(0),                            &
                                   rho_dry(map_w3(1):map_w3(1)+nlayers-1),    &
                                   rho_wet(map_w3(1):map_w3(1)+nlayers-1),    &
                                   height_w3(map_w3(1):map_w3(1)+nlayers-1),  &
                                   t_parcel, mv_parcel, p_parcel, k_start,    &
                                   partial_depth )

    call parcel_ascent_cape_cin( nlayers, k_start, .true._l_def,              &
                                 parcel_mixed_layer_depth, partial_depth,     &
                                 t_parcel,                                    &
                                 t_parcel / ( p_parcel / p_zero )**kappa,     &
                                 mv_parcel, p_parcel,                         &
                                 theta(map_wth(1):map_wth(1)+nlayers),        &
                                 m_v(map_wth(1):map_wth(1)+nlayers),          &
                                 exner_wth(map_wth(1):map_wth(1)+nlayers),    &
                                 p_wth,                                       &
                                 height_wth(map_wth(1):map_wth(1)+nlayers),   &
                                 height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
                                 gravity, rd, cp, kappa, one_over_kappa,      &
                                 epsilon, column_cape, column_cin,            &
                                 lnb_height )

    if ( .not. associated( cape, empty_real_data ) ) then
      cape(map_2d(1)) = column_cape
    end if

    if ( .not. associated( cin, empty_real_data ) ) then
      cin(map_2d(1)) = column_cin
    end if

    if ( .not. associated( equilibrium_height, empty_real_data ) ) then
      if ( column_cape > parcel_min_cape ) then
        equilibrium_height(map_2d(1)) = lnb_height
      else
        equilibrium_height(map_2d(1)) = rmdi
      end if
    end if

  end subroutine cape_cin_mixed_layer_code

end module cape_cin_mixed_layer_kernel_mod
