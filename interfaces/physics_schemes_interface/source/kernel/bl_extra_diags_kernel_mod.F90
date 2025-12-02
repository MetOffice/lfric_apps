!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to various extra bl diagnostics that need to be calculated
!>
module bl_extra_diags_kernel_mod

  use argument_mod,       only : arg_type,                                 &
                                 GH_FIELD, GH_REAL,                        &
                                 GH_READ, GH_WRITE,                        &
                                 CELL_COLUMN,                              &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_def, i_def, i_um, r_um, l_def
  use empty_data_mod,     only : empty_real_data
  use fs_continuity_mod,  only : Wtheta, W3
  use kernel_mod,         only : kernel_type
  use microphysics_config_mod, only : microphysics_casim

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_extra_diags_kernel_type
    private
    type(arg_type) :: meta_args(41) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        & ! rho_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        & ! wetrho_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        & ! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        & ! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! mci
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! mr
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! nr_mphys
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! ns_mphys
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    & ! murk
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! zh
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! t1p5m
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! q1p5m
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! qcl1p5m
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! t1p5m_ssi
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! q1p5m_ssi
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! qcl1p5m_ssi
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! t1p5m_land
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! q1p5m_land
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! qcl1p5m_land
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! wspd10m
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! z0m_eff
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! bl_weight_1dbl
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! ls_rain_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! ls_snow_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! lsca_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! conv_rain_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! conv_snow_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! cca_2d_in
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! ustar_implicit
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! wind_gust
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! scale_dep_wind_gust
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! fog_fraction
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! fog_fraction_ssi
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! fog_fraction_land
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! vis_prob_5km
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! dew_point
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! dew_point_ssi
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! dew_point_land
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! visibility_with_precip
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! visibility_no_precip
                                      /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_extra_diags_code
  end type

  public :: bl_extra_diags_code

contains

  !> @brief Interface to derived boundary layer diagnostics.
  !> @details Calculation of various boundary layer diagnostics.
  !>
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     rho_in_w3              Density field in density space
  !> @param[in]     wetrho_in_w3           Wet density field in w3 space
  !> @param[in]     heat_flux_bl           Vertical heat flux on BL levels
  !> @param[in]     moist_flux_bl          Vertical moisture flux on BL levels
  !> @param[in]     exner_in_wth           Exner
  !> @param[in]     mci                    Cloud ice mixing ratio
  !> @param[in]     mr                     Rain  mixing ratio
  !> @param[in]     nr_mphys               Rain number mixing ratio
  !> @param[in]     ns_mphys               Snow number mixing ratio
  !> @param[in]     zh                     Boundary layer depth
  !> @param[in]     t1p5m                  Diagnostic: 1.5m temperature
  !> @param[in]     q1p5m                  Diagnostic: 1.5m specific humidity
  !> @param[in]     qcl1p5m                Diagnostic: 1.5m specific cloud water content
  !> @param[in]     t1p5m_ssi              Diagnostic: 1.5m temperature over sea and sea-ice
  !> @param[in]     q1p5m_ssi              Diagnostic: 1.5m specific humidity over sea and sea-ice
  !> @param[in]     qcl1p5m_ssi            Diagnostic: 1.5m specific cloud water over sea and sea-ice
  !> @param[in]     t1p5m_land             Diagnostic: 1.5m temperature over land
  !> @param[in]     q1p5m_land             Diagnostic: 1.5m specific humidity over land
  !> @param[in]     qcl1p5m_land           Diagnostic: 1.5m specific cloud water over land
  !> @param[in]     wspd10m                Windspeed at 10m
  !> @param[in]     z0m_eff                Effective roughness length
  !> @param[in]     bl_weight_1dbl         Blending weight to 1D BL scheme in the BL
  !> @param[in]     ls_rain_2d             Surface large-scale  rainfall rate
  !> @param[in]     ls_snow_2d             Surface large-scale snowfall rate
  !> @param[in]     lsca_2d                2D large scale precip fraction
  !> @param[in]     conv_rain_2d           Surface convective rainfall rate
  !> @param[in]     conv_snow_2d           Surface convective snowfall rate
  !> @param[in]     cca_2d_in              2D convective cloud fraction
  !> @param[in]     ustar_implicit         Implicit friction velocity
  !> @param[in,out] wind_gust              Wind gust
  !> @param[in,out] scale_dep_wind_gust    Scale dependent wind gust
  !> @param[in,out] fog_fraction           Fog_fraction
  !> @param[in,out] fog_fraction_ssi       Fog_fraction over sea/sea-ice
  !> @param[in,out] fog_fraction_land      Fog_fraction over land
  !> @param[in,out] vis_prob_5km           vis_prob_5km
  !> @param[in,out] dew_point              Dew point temperature
  !> @param[in,out] dew_point_ssi          Dew point temperature over sea/sea-ice
  !> @param[in,out] dew_point_land         Dew point temperature over land
  !> @param[in,out] visibility_with_precip Visibility with precip included
  !> @param[in,out] visibility_no_precip   Visibility without including precip
  !> @param[in]     ndf_w3                 Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3                Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_wth                Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth               Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth                Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d                 Number of degrees of freedom per cell for 2D fields
  !> @param[in]     undf_2d                Number unique of degrees of freedom  for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields

  subroutine bl_extra_diags_code( nlayers,                  &
                                  rho_in_w3,                &
                                  wetrho_in_w3,             &
                                  heat_flux_bl,             &
                                  moist_flux_bl,            &
                                  exner_in_wth,             &
                                  mci, mr,                  &
                                  nr_mphys, ns_mphys, murk, &
                                  zh,                       &
                                  t1p5m, q1p5m, qcl1p5m,    &
                                  t1p5m_ssi, q1p5m_ssi,     &
                                  qcl1p5m_ssi, t1p5m_land,  &
                                  q1p5m_land, qcl1p5m_land, &
                                  wspd10m,                  &
                                  z0m_eff, bl_weight_1dbl,  &
                                  ls_rain_2d, ls_snow_2d,   &
                                  lsca_2d,                  &
                                  conv_rain_2d,             &
                                  conv_snow_2d, cca_2d_in,  &
                                  ustar_implicit, wind_gust,&
                                  scale_dep_wind_gust,      &
                                  fog_fraction,             &
                                  fog_fraction_ssi,         &
                                  fog_fraction_land,        &
                                  vis_prob_5km, dew_point,  &
                                  dew_point_ssi,            &
                                  dew_point_land,           &
                                  visibility_with_precip,   &
                                  visibility_no_precip,     &
                                  ndf_w3,                   &
                                  undf_w3,                  &
                                  map_w3,                   &
                                  ndf_wth,                  &
                                  undf_wth,                 &
                                  map_wth,                  &
                                  ndf_2d,                   &
                                  undf_2d,                  &
                                  map_2d                  )

    use aerosol_config_mod,   only : murk_visibility
    use beta_precip_mod,      only : beta_precip
    use blayer_config_mod,    only : c_gust
    use cloud_inputs_mod,     only : rhcrit
    use dewpnt_mod,           only : dewpnt
    use fog_fr_mod,           only : fog_fr
    use mphys_constants_mod,  only : mprog_min
    use nlsizes_namelist_mod, only : row_length, rows
    use planet_config_mod,    only : p_zero, kappa, gravity, cp
    use planet_constants_mod, only : vkman, c_virtual
    use vis_precip_mod,       only : vis_precip
    use visbty_constants_mod, only : n_vis_thresh, vis_thresh
    use visbty_mod,           only : visbty
    use variable_precision,   only : wp

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in)     :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in)     :: ndf_2d, undf_2d

    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    real(kind=r_def), intent(in), dimension(undf_w3)    :: rho_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3)    :: wetrho_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3)    :: heat_flux_bl
    real(kind=r_def), intent(in), dimension(undf_w3)    :: moist_flux_bl
    real(kind=r_def), intent(in), dimension(undf_wth)   :: exner_in_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: mci
    real(kind=r_def), intent(in), dimension(undf_wth)   :: mr
    real(kind=r_def), intent(in), dimension(undf_wth)   :: nr_mphys
    real(kind=r_def), intent(in), dimension(undf_wth)   :: ns_mphys
    real(kind=r_def), intent(in), dimension(undf_wth)   :: murk
    real(kind=r_def), intent(in), dimension(undf_2d)    :: zh
    real(kind=r_def), intent(in), dimension(undf_2d)    :: bl_weight_1dbl
    real(kind=r_def), intent(in), dimension(undf_2d)    :: ls_rain_2d
    real(kind=r_def), intent(in), dimension(undf_2d)    :: ls_snow_2d
    real(kind=r_def), intent(in), dimension(undf_2d)    :: lsca_2d
    real(kind=r_def), intent(in), dimension(undf_2d)    :: conv_rain_2d
    real(kind=r_def), intent(in), dimension(undf_2d)    :: conv_snow_2d
    real(kind=r_def), intent(in), dimension(undf_2d)    :: cca_2d_in
    real(kind=r_def), intent(in),    pointer :: t1p5m(:), q1p5m(:), qcl1p5m(:)
    real(kind=r_def), intent(in),    pointer :: t1p5m_ssi(:), q1p5m_ssi(:), qcl1p5m_ssi(:)
    real(kind=r_def), intent(in),    pointer :: t1p5m_land(:), q1p5m_land(:), qcl1p5m_land(:)
    real(kind=r_def), intent(in),    pointer :: wspd10m(:), z0m_eff(:)
    real(kind=r_def), intent(inout), pointer :: ustar_implicit(:)
    real(kind=r_def), intent(inout), pointer :: wind_gust(:), scale_dep_wind_gust(:)
    real(kind=r_def), intent(inout), pointer :: fog_fraction(:), vis_prob_5km(:)
    real(kind=r_def), intent(inout), pointer :: fog_fraction_ssi(:), fog_fraction_land(:)
    real(kind=r_def), intent(inout), pointer :: dew_point(:)
    real(kind=r_def), intent(inout), pointer :: dew_point_ssi(:), dew_point_land(:)
    real(kind=r_def), intent(inout), pointer :: visibility_with_precip(:)
    real(kind=r_def), intent(inout), pointer :: visibility_no_precip(:)

    real(kind=r_def), parameter :: one_third   = 1.0_r_def/3.0_r_def

    ! Tunable parameters used in the calculation of the wind gust
    real(kind=r_def), parameter :: c_ws        = 1.0_r_def/24.0_r_def
    real(kind=r_def), parameter :: gust_const  = 2.29_r_def

    ! Switches needed for visibility calculations
    logical(l_def),      parameter :: pct = .false.  ! Cloud amounts are in %
    logical(l_def),      parameter :: avg = .true.   ! Precip=local*prob
    integer(kind=i_def), parameter :: fog_thres=1
    integer(kind=i_def), parameter :: vis5km_thres=2
    real(kind=r_um),     parameter :: calc_prob_of_vis = 0.5_r_um

    ! single level real fields input
    real(r_um), dimension(row_length,rows) ::                                &
         ls_rain, ls_snow, conv_rain, conv_snow, cca_2d, p_star, rho1, qcf1, &
         qrain1, aerosol1, plsp, t1p5m_loc, q1p5m_loc, qcl1p5m_loc

    ! single level real fields calculated
    real(r_um), dimension(row_length,rows) ::                                &
         beta_ls_rain, beta_ls_snow, beta_c_rain, beta_c_snow,               &
         vis, vis_ls_precip, vis_c_precip, vis_no_precip, dew_pnt

    real(r_um), dimension(row_length,rows,0:1) :: snownumber, rainnumber

    ! fog_fr works for n levels, we want 1
    real(r_um), dimension(row_length,rows,1,n_vis_thresh) :: vis_threshold
    real(r_um), dimension(row_length,rows,n_vis_thresh)   :: pvis

    ! Local scalars
    real(kind=r_def) :: ftl_surf, fqw_surf, &
                        wstar3_imp, std_dev, gust_contribution

    integer(kind=i_def) :: k, icode, i,j

    if ( .not. associated(wind_gust, empty_real_data) .or.                   &
         .not. associated(scale_dep_wind_gust, empty_real_data) ) then
      ftl_surf = heat_flux_bl(map_w3(1)) / cp
      fqw_surf = moist_flux_bl(map_w3(1))
      wstar3_imp = zh(map_2d(1)) * gravity * ( ftl_surf/t1p5m(map_2d(1)) +   &
                                               fqw_surf*c_virtual ) /        &
                                             rho_in_w3(map_w3(1))
      if ( wstar3_imp > 0.0_r_def ) then
        ! Include the stability dependence
        std_dev = gust_const * ( ustar_implicit(map_2d(1))**3.0_r_def +      &
                                 vkman * c_ws * wstar3_imp )**one_third
      else
        std_dev = gust_const * ustar_implicit(map_2d(1))
      end if
      gust_contribution = std_dev * (1.0_r_def/vkman) *                      &
              LOG( (5.0_r_def * EXP(vkman * c_gust) + z0m_eff(map_2d(1)) ) / &
                   (5.0_r_def + z0m_eff(map_2d(1)) ) )

      if ( .not. associated(wind_gust, empty_real_data) ) then
        ! Original scale-independent gust diagnostic
        ! Add the whole gust contribution to the mean wind speed
        wind_gust(map_2d(1)) = wspd10m(map_2d(1)) + gust_contribution
      end if

      if ( .not. associated(scale_dep_wind_gust, empty_real_data) ) then
        ! Scale-dependent gust diagnostic
        ! Note that bl_weight_1dbl is weight_1dbl from the bottom grid level
        scale_dep_wind_gust(map_2d(1)) = wspd10m(map_2d(1)) +                &
                                  bl_weight_1dbl(map_2d(1))*gust_contribution
      end if

    end if

    ! map main input fields
    if (.not. associated(visibility_no_precip, empty_real_data)   .or.       &
        .not. associated(visibility_with_precip, empty_real_data) .or.       &
        .not. associated(fog_fraction, empty_real_data)           .or.       &
        .not. associated(fog_fraction_ssi, empty_real_data)       .or.       &
        .not. associated(fog_fraction_land, empty_real_data)      .or.       &
        .not. associated(vis_prob_5km, empty_real_data)           .or.       &
        .not. associated(dew_point, empty_real_data)              .or.       &
        .not. associated(dew_point_ssi, empty_real_data)          .or.       &
        .not. associated(dew_point_land, empty_real_data) ) then
      ! surface pressure
      p_star(1,1)    = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
      ! level 1 of aerosol (using the standard default of 10 for now)
      if (murk_visibility) aerosol1(1,1)  = murk(map_wth(1)+0)
      ! copy of screen variables
      t1p5m_loc(1,1)   = t1p5m(map_2d(1))
      q1p5m_loc(1,1)   = q1p5m(map_2d(1))
      qcl1p5m_loc(1,1) = qcl1p5m(map_2d(1))
    end if

    ! Visibility
    if ( .not. associated(visibility_no_precip, empty_real_data) .or.        &
         .not. associated(visibility_with_precip, empty_real_data) ) then
      call visbty(                                                           &
                  ! inputs
                  p_star, t1p5m_loc, q1p5m_loc, qcl1p5m_loc, aerosol1,       &
                  calc_prob_of_vis, rhcrit(1), murk_visibility, 1,           &
                  ! output
                  vis_no_precip )
      visibility_no_precip(map_2d(1)) = vis_no_precip(1,1)

      ! Visibility at 1.5 m including precipitation
      if ( .not. associated(visibility_with_precip, empty_real_data) ) then
        ! map additional input fields
        ! level 1 rho
        rho1(1,1)      = wetrho_in_w3(map_w3(1))
        ! level 1 cloud ice mixing ratio
        qcf1(1,1)      = mci(map_wth(1) + 1)
        ! level 1 rain mixing ratio
        qrain1(1,1)    = mr(map_wth(1) + 1)
        ! surface rain and snow rates from large-scale microphysics
        ls_rain(1,1)   = ls_rain_2d(map_2d(1))
        ls_snow(1,1)   = ls_snow_2d(map_2d(1))
        ! surface rain and snow rates from convection
        conv_rain(1,1) = conv_rain_2d(map_2d(1))
        conv_snow(1,1) = conv_snow_2d(map_2d(1))
        ! cca_2d
        cca_2d(1,1)    = cca_2d_in(map_2d(1))
        ! prob of ls precip - just use existing rain area fraction
        plsp(1,1)      = lsca_2d(map_2d(1))

        ! number prognostics used in the visibility calculation
        ! We only copy these if casim is enabled. Otherwise they will
        ! not be used.
        if (microphysics_casim) then
           rainnumber(1,1,1) = nr_mphys(map_wth(1) + 1)
           snownumber(1,1,1) = ns_mphys(map_wth(1) + 1)
        end if

        call beta_precip( ls_rain, ls_snow,                                    &
                          conv_rain, conv_snow, qcf1, qrain1,                  &
                          rho1, t1p5m_loc, p_star, snownumber, rainnumber,     &
                          plsp,cca_2d,pct,avg,                                 &
                          1, 1, 1,                                             &
                          beta_ls_rain, beta_ls_snow,                          &
                          beta_c_rain, beta_c_snow )
        call vis_precip( vis_no_precip,                                        &
                         plsp,cca_2d,pct,                                      &
                         beta_ls_rain, beta_ls_snow,                           &
                         beta_c_rain, beta_c_snow,                             &
                         1, 1, 1,                                              &
                         vis,vis_ls_precip,vis_c_precip,                       &
                         icode )
        visibility_with_precip(map_2d(1)) = vis(1,1)

      end if ! vis with precip
    end if ! any vis

    ! fog fraction
    if ( .not. associated(fog_fraction, empty_real_data) .or.                  &
         .not. associated(vis_prob_5km, empty_real_data) ) then
      do k = 1, n_vis_thresh
        vis_threshold(1,1,1,k)=vis_thresh(k)
      end do
      call fog_fr( p_star, rhcrit, 1, 1,                                       &
                   t1p5m_loc, aerosol1, murk_visibility,                       &
                   q1p5m_loc, qcl1p5m_loc,                                     &
                   vis_threshold, pvis, n_vis_thresh )
      if ( .not. associated(fog_fraction, empty_real_data) )                   &
                fog_fraction(map_2d(1)) = pvis(1,1,fog_thres)
      if ( .not. associated(vis_prob_5km, empty_real_data) )                   &
                vis_prob_5km(map_2d(1)) = pvis(1,1,vis5km_thres)
    end if

    ! dew point
    if ( .not. associated(dew_point, empty_real_data) ) then
      if (q1p5m_loc(1,1) > mprog_min) then
        call dewpnt(q1p5m_loc, p_star, t1p5m_loc, 1, dew_pnt)
      else
        dew_pnt(1,1) = 0.0_r_def  ! no water
      end if
      dew_point(map_2d(1)) = dew_pnt(1,1)
    end if

    ! sea and sea-ice diagnostics
    if (.not. associated(fog_fraction_ssi, empty_real_data)           .or.       &
        .not. associated(dew_point_ssi, empty_real_data) ) then
      ! copy of screen variables
      t1p5m_loc(1,1)   = t1p5m_ssi(map_2d(1))
      q1p5m_loc(1,1)   = q1p5m_ssi(map_2d(1))
      qcl1p5m_loc(1,1) = qcl1p5m_ssi(map_2d(1))
    end if

    if ( .not. associated(fog_fraction_ssi, empty_real_data) ) then
      do k = 1, 1
        vis_threshold(1,1,1,k)=vis_thresh(k)
      end do
      call fog_fr( p_star, rhcrit, 1, 1,                                       &
                   t1p5m_loc, aerosol1, murk_visibility,                       &
                   q1p5m_loc, qcl1p5m_loc,                                     &
                   vis_threshold, pvis, 1 )
      fog_fraction_ssi(map_2d(1)) = pvis(1,1,fog_thres)
    end if

    if ( .not. associated(dew_point_ssi, empty_real_data) ) then
      if (q1p5m_loc(1,1) > mprog_min) then
        call dewpnt(q1p5m_loc, p_star, t1p5m_loc, 1, dew_pnt)
      else
        dew_pnt(1,1) = 0.0_r_def  ! no water
      end if
      dew_point_ssi(map_2d(1)) = dew_pnt(1,1)
    end if

    ! land diagnostics
    if (.not. associated(fog_fraction_land, empty_real_data)           .or.    &
        .not. associated(dew_point_land, empty_real_data) ) then
      ! copy of screen variables
      t1p5m_loc(1,1)   = t1p5m_land(map_2d(1))
      q1p5m_loc(1,1)   = q1p5m_land(map_2d(1))
      qcl1p5m_loc(1,1) = qcl1p5m_land(map_2d(1))
    end if

    if ( .not. associated(fog_fraction_land, empty_real_data) ) then
      do k = 1, 1
        vis_threshold(1,1,1,k)=vis_thresh(k)
      end do
      call fog_fr( p_star, rhcrit, 1, 1,                                       &
                   t1p5m_loc, aerosol1, murk_visibility,                       &
                   q1p5m_loc, qcl1p5m_loc,                                     &
                   vis_threshold, pvis, 1 )
      fog_fraction_land(map_2d(1)) = pvis(1,1,fog_thres)
    end if

    if ( .not. associated(dew_point_land, empty_real_data) ) then
      if (q1p5m_loc(1,1) > mprog_min) then
        call dewpnt(q1p5m_loc, p_star, t1p5m_loc, 1, dew_pnt)
      else
        dew_pnt(1,1) = 0.0_r_def  ! no water
      end if
      dew_point_land(map_2d(1)) = dew_pnt(1,1)
    end if

  end subroutine bl_extra_diags_code

end module bl_extra_diags_kernel_mod
