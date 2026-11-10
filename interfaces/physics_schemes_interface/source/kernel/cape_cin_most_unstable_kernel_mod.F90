!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code.
!-----------------------------------------------------------------------------
!> @brief CAPE and CIN of the most unstable parcel and the effective inflow
!>        layer.
module cape_cin_most_unstable_kernel_mod

  use argument_mod,      only: arg_type,                                      &
                               GH_FIELD, GH_SCALAR, GH_REAL,                  &
                               GH_READ, GH_WRITE, CELL_COLUMN,                &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def, l_def
  use empty_data_mod,    only: empty_real_data
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata type.
  type, public, extends(kernel_type) :: cape_cin_most_unstable_kernel_type
    private
    type(arg_type) :: meta_args(20) = (/                                      &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! theta
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! m_v
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                      & ! height_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                          & ! height_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),   & ! t1p5m
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),   & ! q1p5m
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! cape
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! cin
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! start_height
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! equilibrium_height
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! inflow_base
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! inflow_top
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
    procedure, nopass :: cape_cin_most_unstable_code
  end type

  public :: cape_cin_most_unstable_code

contains

  !> @brief CAPE and CIN of the most unstable parcel and the effective inflow
  !>        layer.
  !> @details Parcels are lifted without mixing from screen level and from
  !>          each model level below 500 hPa. The most unstable parcel is the
  !>          one with the largest CAPE; the effective inflow layer contains
  !>          the starting levels whose parcels have enough CAPE and little
  !>          enough CIN. Heights are above sea level and undefined values are
  !>          missing data. Outputs that were not requested are not written.
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     theta              Potential temperature (K)
  !> @param[in]     m_v                Water vapour mixing ratio
  !> @param[in]     exner_wth          Exner pressure on Wtheta
  !> @param[in]     height_wth         Height of Wtheta levels above sea
  !>                                   level (m)
  !> @param[in]     height_w3          Height of W3 levels above sea level (m)
  !> @param[in]     t1p5m              Temperature at 1.5 m (K)
  !> @param[in]     q1p5m              Specific humidity at 1.5 m (kg/kg)
  !> @param[in,out] cape               CAPE of the most unstable parcel (J/kg)
  !> @param[in,out] cin                CIN of the most unstable parcel (J/kg)
  !> @param[in,out] start_height       Starting height of the most unstable
  !>                                   parcel (m)
  !> @param[in,out] equilibrium_height Equilibrium level height of the most
  !>                                   unstable parcel (m)
  !> @param[in,out] inflow_base        Base of the effective inflow layer (m)
  !> @param[in,out] inflow_top         Top of the effective inflow layer (m)
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
  subroutine cape_cin_most_unstable_code( nlayers,                            &
                                          theta, m_v, exner_wth, height_wth,  &
                                          height_w3, t1p5m, q1p5m,            &
                                          cape, cin, start_height,            &
                                          equilibrium_height,                 &
                                          inflow_base, inflow_top,            &
                                          gravity, rd, cp, kappa,             &
                                          one_over_kappa, p_zero, epsilon,    &
                                          ndf_wth, undf_wth, map_wth,         &
                                          ndf_w3, undf_w3, map_w3,            &
                                          ndf_2d, undf_2d, map_2d )

    use parcel_ascent_cape_cin_mod,     only: parcel_ascent_cape_cin,         &
                                              most_unstable_cape_cin
    use science_aviation_constants_mod, only: parcel_surface_layer_depth

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
    real(kind=r_def), intent(in), dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in), dimension(undf_2d)  :: t1p5m, q1p5m

    real(kind=r_def), pointer, intent(inout) :: cape(:)
    real(kind=r_def), pointer, intent(inout) :: cin(:)
    real(kind=r_def), pointer, intent(inout) :: start_height(:)
    real(kind=r_def), pointer, intent(inout) :: equilibrium_height(:)
    real(kind=r_def), pointer, intent(inout) :: inflow_base(:)
    real(kind=r_def), pointer, intent(inout) :: inflow_top(:)

    real(kind=r_def), intent(in) :: gravity, rd, cp, kappa
    real(kind=r_def), intent(in) :: one_over_kappa, p_zero, epsilon

    ! Local variables
    real(kind=r_def) :: p_wth(0:nlayers)
    real(kind=r_def) :: mv_parcel
    real(kind=r_def) :: surface_cape, surface_cin, surface_lnb_height
    real(kind=r_def) :: mu_cape, mu_cin, mu_start_height, mu_lnb_height
    real(kind=r_def) :: column_inflow_base, column_inflow_top

    p_wth = p_zero * exner_wth(map_wth(1):map_wth(1)+nlayers)**one_over_kappa

    ! Ascent from screen level, with the specific humidity as a mixing ratio
    ! and the first layer starting at the top of the surface layer.
    mv_parcel = q1p5m(map_2d(1)) / ( 1.0_r_def - q1p5m(map_2d(1)) )
    call parcel_ascent_cape_cin( nlayers, 1_i_def, .true._l_def,              &
                                 parcel_surface_layer_depth,                  &
                                 parcel_surface_layer_depth,                  &
                                 t1p5m(map_2d(1)),                            &
                                 t1p5m(map_2d(1)) / exner_wth(map_wth(1)),    &
                                 mv_parcel, p_wth(0),                         &
                                 theta(map_wth(1):map_wth(1)+nlayers),        &
                                 m_v(map_wth(1):map_wth(1)+nlayers),          &
                                 exner_wth(map_wth(1):map_wth(1)+nlayers),    &
                                 p_wth,                                       &
                                 height_wth(map_wth(1):map_wth(1)+nlayers),   &
                                 height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
                                 gravity, rd, cp, kappa, one_over_kappa,      &
                                 epsilon, surface_cape, surface_cin,          &
                                 surface_lnb_height )

    ! Ascents from each model level below 500 hPa.
    call most_unstable_cape_cin( nlayers, surface_cape, surface_cin,          &
                                 surface_lnb_height,                          &
                                 theta(map_wth(1):map_wth(1)+nlayers),        &
                                 m_v(map_wth(1):map_wth(1)+nlayers),          &
                                 exner_wth(map_wth(1):map_wth(1)+nlayers),    &
                                 p_wth,                                       &
                                 height_wth(map_wth(1):map_wth(1)+nlayers),   &
                                 height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
                                 gravity, rd, cp, kappa, one_over_kappa,      &
                                 epsilon, mu_cape, mu_cin, mu_start_height,   &
                                 mu_lnb_height, column_inflow_base,           &
                                 column_inflow_top )

    if ( .not. associated( cape, empty_real_data ) ) then
      cape(map_2d(1)) = mu_cape
    end if

    if ( .not. associated( cin, empty_real_data ) ) then
      cin(map_2d(1)) = mu_cin
    end if

    if ( .not. associated( start_height, empty_real_data ) ) then
      start_height(map_2d(1)) = mu_start_height
    end if

    if ( .not. associated( equilibrium_height, empty_real_data ) ) then
      equilibrium_height(map_2d(1)) = mu_lnb_height
    end if

    if ( .not. associated( inflow_base, empty_real_data ) ) then
      inflow_base(map_2d(1)) = column_inflow_base
    end if

    if ( .not. associated( inflow_top, empty_real_data ) ) then
      inflow_top(map_2d(1)) = column_inflow_top
    end if

  end subroutine cape_cin_most_unstable_code

end module cape_cin_most_unstable_kernel_mod
