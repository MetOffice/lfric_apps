!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code.
!-----------------------------------------------------------------------------
!> @brief Column routines for undilute parcel ascents giving CAPE and CIN.
!> @details A parcel is lifted dry adiabatically to its lifting condensation
!>          level and then along a saturated adiabat. Its buoyancy relative
!>          to the environment is summed layer by layer into the convective
!>          available potential energy (CAPE, positive part) and the
!>          convective inhibition (CIN, negative part) until it reaches its
!>          level of neutral buoyancy. A second buoyant layer above a capping
!>          inversion is included. Column arrays on Wtheta levels run from 0
!>          (the surface) to nlayers, and on W3 levels from 0 to nlayers - 1.
module parcel_ascent_cape_cin_mod

  use constants_mod,                  only: r_def, i_def, l_def, rmdi
  use driver_water_constants_mod,     only: latent_heat_h2o_condensation
  use qsat_mod,                       only: qsat_wat_mix
  use science_aviation_constants_mod, only: parcel_mixed_layer_depth,         &
                                            parcel_surface_layer_depth,       &
                                            parcel_min_cape, parcel_min_cin,  &
                                            parcel_qsat_delta_t,              &
                                            parcel_screen_height,             &
                                            parcel_top_start_pressure,        &
                                            inflow_min_cape, inflow_min_cin,  &
                                            bolton_a, bolton_b, bolton_c,     &
                                            bolton_d
  use science_conversions_mod,        only: pascals_to_hectopascals

  implicit none

  private
  public :: mixed_layer_parcel_start, lifting_condensation_level,             &
            parcel_ascent_cape_cin, most_unstable_cape_cin,                   &
            highest_start_level, initialise_most_unstable,                    &
            update_most_unstable

  real(kind=r_def), parameter :: half = 0.5_r_def

contains

  !> @brief Start state of a parcel representing the near-surface mixed layer.
  !> @details Temperature and water vapour are averaged over the mixed layer,
  !>          weighted by wet and dry air mass respectively, using values
  !>          interpolated to W3 levels. The parcel starts at the pressure of
  !>          the middle of the mixed layer, and its ascent starts at the first
  !>          Wtheta level above the mixed-layer top; the depth of the layer
  !>          around that level lying below the mixed-layer top is returned.
  !> @param[in]  nlayers       Number of layers
  !> @param[in]  theta         Potential temperature on Wtheta levels (K)
  !> @param[in]  m_v           Water vapour mixing ratio on Wtheta levels
  !> @param[in]  height_wth    Height of Wtheta levels above sea level (m)
  !> @param[in]  exner_w3      Exner pressure on W3 levels
  !> @param[in]  p_w3          Pressure on W3 levels (Pa)
  !> @param[in]  p_surf        Surface pressure (Pa)
  !> @param[in]  rho_dry       Dry density on W3 levels (kg m-3)
  !> @param[in]  rho_wet       Wet density on W3 levels (kg m-3)
  !> @param[in]  height_w3     Height of W3 levels above sea level (m)
  !> @param[out] t_parcel      Parcel temperature (K)
  !> @param[out] mv_parcel     Parcel water vapour mixing ratio
  !> @param[out] p_parcel      Parcel pressure (Pa)
  !> @param[out] k_start       First Wtheta level above the mixed layer
  !> @param[out] partial_depth Depth of the layer around k_start that lies
  !>                           below the mixed-layer top (m)
  subroutine mixed_layer_parcel_start( nlayers, theta, m_v, height_wth,       &
                                       exner_w3, p_w3, p_surf, rho_dry,       &
                                       rho_wet, height_w3, t_parcel,          &
                                       mv_parcel, p_parcel, k_start,          &
                                       partial_depth )

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers
    real(kind=r_def),    intent(in)  :: theta(0:nlayers)
    real(kind=r_def),    intent(in)  :: m_v(0:nlayers)
    real(kind=r_def),    intent(in)  :: height_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: exner_w3(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: p_w3(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: p_surf
    real(kind=r_def),    intent(in)  :: rho_dry(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: rho_wet(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: height_w3(0:nlayers-1)
    real(kind=r_def),    intent(out) :: t_parcel, mv_parcel, p_parcel
    integer(kind=i_def), intent(out) :: k_start
    real(kind=r_def),    intent(out) :: partial_depth

    integer(kind=i_def) :: k
    real(kind=r_def)    :: z_bottom, z_top, dz, weight, theta_w3, mv_w3
    real(kind=r_def)    :: t_sum, mv_sum, wet_mass, dry_mass
    real(kind=r_def)    :: z_lower, z_upper, p_lower, z_mid, z_layer_base

    z_mid = half * parcel_mixed_layer_depth

    ! Mass-weighted averages over the W3 layers below the mixed-layer top.
    ! Layer k - 1 lies between Wtheta levels k - 1 and k; the lowest layer
    ! takes the values at the first level above the surface.
    t_sum    = 0.0_r_def
    mv_sum   = 0.0_r_def
    wet_mass = 0.0_r_def
    dry_mass = 0.0_r_def
    do k = 1, nlayers
      z_bottom = height_wth(k-1) - height_wth(0)
      if ( z_bottom >= parcel_mixed_layer_depth ) exit
      z_top = height_wth(k) - height_wth(0)
      dz = min( z_top, parcel_mixed_layer_depth ) - z_bottom
      if ( k == 1 ) then
        theta_w3 = theta(1)
        mv_w3    = m_v(1)
      else
        weight = ( height_w3(k-1) - height_wth(k-1) )                         &
               / ( height_wth(k) - height_wth(k-1) )
        theta_w3 = weight * theta(k) + ( 1.0_r_def - weight ) * theta(k-1)
        mv_w3    = weight * m_v(k)   + ( 1.0_r_def - weight ) * m_v(k-1)
      end if
      t_sum    = t_sum    + theta_w3 * exner_w3(k-1) * rho_wet(k-1) * dz
      wet_mass = wet_mass + rho_wet(k-1) * dz
      mv_sum   = mv_sum   + mv_w3 * rho_dry(k-1) * dz
      dry_mass = dry_mass + rho_dry(k-1) * dz
    end do
    t_parcel  = t_sum / wet_mass
    mv_parcel = mv_sum / dry_mass

    ! First Wtheta level above the mixed-layer top. The layer around Wtheta
    ! level k is bounded below by W3 level k - 1, or by the surface for k = 1.
    k_start = nlayers
    do k = 1, nlayers
      if ( height_wth(k) - height_wth(0) > parcel_mixed_layer_depth ) then
        k_start = k
        exit
      end if
    end do
    if ( k_start == 1 ) then
      z_layer_base = 0.0_r_def
    else
      z_layer_base = height_w3(k_start-1) - height_wth(0)
    end if
    partial_depth = parcel_mixed_layer_depth - z_layer_base

    ! Pressure at the middle of the mixed layer, interpolated linearly in
    ! height between the surface and W3 levels.
    z_lower  = 0.0_r_def
    p_lower  = p_surf
    p_parcel = p_surf
    do k = 0, nlayers - 1
      z_upper = height_w3(k) - height_wth(0)
      if ( z_upper > z_mid ) then
        weight   = ( z_upper - z_mid ) / ( z_upper - z_lower )
        p_parcel = weight * p_lower + ( 1.0_r_def - weight ) * p_w3(k)
        exit
      end if
      z_lower = z_upper
      p_lower = p_w3(k)
    end do

  end subroutine mixed_layer_parcel_start

  !> @brief Temperature and pressure of a parcel at its lifting condensation
  !>        level.
  !> @details Uses the empirical fit of Bolton (1980). A parcel with no water
  !>          vapour is given its own temperature and pressure.
  !> @param[in]  t_parcel       Parcel temperature (K)
  !> @param[in]  mv_parcel      Parcel water vapour mixing ratio
  !> @param[in]  p_parcel       Parcel pressure (Pa)
  !> @param[in]  one_over_kappa Ratio of specific heat to gas constant
  !> @param[in]  epsilon        Ratio of gas constants of dry air and vapour
  !> @param[out] t_lcl          Temperature at the condensation level (K)
  !> @param[out] p_lcl          Pressure at the condensation level (Pa)
  subroutine lifting_condensation_level( t_parcel, mv_parcel, p_parcel,       &
                                         one_over_kappa, epsilon,             &
                                         t_lcl, p_lcl )

    implicit none

    real(kind=r_def), intent(in)  :: t_parcel, mv_parcel, p_parcel
    real(kind=r_def), intent(in)  :: one_over_kappa, epsilon
    real(kind=r_def), intent(out) :: t_lcl, p_lcl

    real(kind=r_def) :: vapour_pressure

    vapour_pressure = pascals_to_hectopascals * mv_parcel * p_parcel          &
                    / ( epsilon + mv_parcel )
    if ( vapour_pressure > 0.0_r_def ) then
      t_lcl = bolton_a + bolton_b / ( bolton_c * log( t_parcel )              &
                                      - log( vapour_pressure ) - bolton_d )
      p_lcl = p_parcel * ( t_lcl / t_parcel )**one_over_kappa
    else
      t_lcl = t_parcel
      p_lcl = p_parcel
    end if

  end subroutine lifting_condensation_level

  !> @brief Undilute parcel ascent giving CAPE, CIN and the height of the
  !>        level of neutral buoyancy.
  !> @details The ascent begins at Wtheta level k_start. A parcel that starts
  !>          below that level (a mixed-layer or surface parcel) is first
  !>          lifted to it, and only the part of the layer around k_start
  !>          above partial_depth contributes there. The ascent of a column
  !>          stops once its CIN falls below the minimum or a second level of
  !>          neutral buoyancy is found. The height is zero if no level of
  !>          neutral buoyancy is found.
  !> @param[in]  nlayers           Number of layers
  !> @param[in]  k_start           Wtheta level at which the ascent starts
  !> @param[in]  l_start_below     Parcel starts below level k_start
  !> @param[in]  start_depth       Depth of the layer the parcel represents (m)
  !> @param[in]  partial_depth     Depth of the layer around k_start lying
  !>                               below the parcel's start layer top (m)
  !> @param[in]  t_start           Parcel start temperature (K)
  !> @param[in]  theta_start       Parcel start potential temperature (K)
  !> @param[in]  mv_start          Parcel start water vapour mixing ratio
  !> @param[in]  p_start           Parcel start pressure (Pa)
  !> @param[in]  theta             Potential temperature on Wtheta levels (K)
  !> @param[in]  m_v               Water vapour mixing ratio on Wtheta levels
  !> @param[in]  exner_wth         Exner pressure on Wtheta levels
  !> @param[in]  p_wth             Pressure on Wtheta levels (Pa)
  !> @param[in]  height_wth        Height of Wtheta levels above sea level (m)
  !> @param[in]  height_w3         Height of W3 levels above sea level (m)
  !> @param[in]  gravity           Gravitational acceleration (m s-2)
  !> @param[in]  rd                Gas constant for dry air (J kg-1 K-1)
  !> @param[in]  cp                Specific heat of dry air (J kg-1 K-1)
  !> @param[in]  kappa             Ratio of rd to cp
  !> @param[in]  one_over_kappa    Ratio of cp to rd
  !> @param[in]  epsilon           Ratio of gas constants of dry air and vapour
  !> @param[out] cape              Convective available potential energy (J/kg)
  !> @param[out] cin               Convective inhibition (J/kg)
  !> @param[out] lnb_height        Height of the level of neutral buoyancy
  !>                               above sea level (m)
  subroutine parcel_ascent_cape_cin( nlayers, k_start, l_start_below,         &
                                     start_depth, partial_depth,              &
                                     t_start, theta_start, mv_start, p_start, &
                                     theta, m_v, exner_wth, p_wth,            &
                                     height_wth, height_w3,                   &
                                     gravity, rd, cp, kappa, one_over_kappa,  &
                                     epsilon, cape, cin, lnb_height )

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers, k_start
    logical(kind=l_def), intent(in)  :: l_start_below
    real(kind=r_def),    intent(in)  :: start_depth, partial_depth
    real(kind=r_def),    intent(in)  :: t_start, theta_start, mv_start
    real(kind=r_def),    intent(in)  :: p_start
    real(kind=r_def),    intent(in)  :: theta(0:nlayers)
    real(kind=r_def),    intent(in)  :: m_v(0:nlayers)
    real(kind=r_def),    intent(in)  :: exner_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: p_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: height_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: height_w3(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: gravity, rd, cp, kappa
    real(kind=r_def),    intent(in)  :: one_over_kappa, epsilon
    real(kind=r_def),    intent(out) :: cape, cin, lnb_height

    integer(kind=i_def) :: k
    logical(kind=l_def) :: l_lnb, l_second_lnb
    real(kind=r_def)    :: t_parcel, theta_parcel, mv_parcel
    real(kind=r_def)    :: t_lcl, p_lcl, q_lcl, qsat_parcel
    real(kind=r_def)    :: thv_env, dthv, dthv_buoyant, dcape, dz
    real(kind=r_def)    :: dz_reduction
    real(kind=r_def)    :: cape2, cin2, lnb_height2

    cape  = 0.0_r_def
    cin   = 0.0_r_def
    cape2 = 0.0_r_def
    cin2  = 0.0_r_def
    lnb_height   = 0.0_r_def
    lnb_height2  = 0.0_r_def
    dthv_buoyant = 0.0_r_def
    l_lnb        = .false.
    l_second_lnb = .false.

    ! Part of the first layer that lies below the parcel's start layer top.
    dz_reduction = 0.0_r_def
    if ( l_start_below ) dz_reduction = partial_depth

    t_parcel     = t_start
    theta_parcel = theta_start
    mv_parcel    = mv_start

    call lifting_condensation_level( t_parcel, mv_parcel, p_start,            &
                                     one_over_kappa, epsilon, t_lcl, p_lcl )
    call qsat_wat_mix( q_lcl, t_lcl, p_lcl )

    ! Lift a parcel starting below level k_start up to that level. Any
    ! supersaturation is removed first, its buoyancy relative to the parcel
    ! before condensation acting over the depth of the start layer.
    if ( l_start_below ) then
      call qsat_wat_mix( qsat_parcel, t_parcel, p_start )
      thv_env   = virtual_theta( theta_parcel, mv_parcel, epsilon )
      mv_parcel = min( mv_parcel, qsat_parcel )
      dcape = ( virtual_theta( theta_parcel, mv_parcel, epsilon ) - thv_env ) &
            * gravity * start_depth / thv_env
      if ( dcape > 0.0_r_def ) then
        cape = cape + dcape
      else
        cin = cin + dcape
      end if
      call lift_parcel( t_parcel, theta_parcel, qsat_parcel, p_start,         &
                        p_wth(k_start), exner_wth(k_start), t_lcl, p_lcl,     &
                        q_lcl, rd, cp, kappa )
    end if

    do k = k_start, nlayers - 1

      if ( cin + cin2 < parcel_min_cin .or. l_second_lnb ) exit

      ! Parcel water vapour is limited to its saturation value.
      call qsat_wat_mix( qsat_parcel, t_parcel, p_wth(k) )
      mv_parcel = min( mv_parcel, qsat_parcel )

      ! Buoyancy of the parcel over the layer around level k.
      thv_env = virtual_theta( theta(k), m_v(k), epsilon )
      dz = theta_layer_depth( nlayers, k, height_wth, height_w3 )
      if ( k == k_start ) dz = dz - dz_reduction
      dthv  = virtual_theta( theta_parcel, mv_parcel, epsilon ) - thv_env
      dcape = dthv * gravity * dz / thv_env

      if ( .not. l_lnb ) then
        ! First buoyant layer.
        if ( dcape < 0.0_r_def .and. cape > parcel_min_cape ) then
          lnb_height = neutral_buoyancy_height( dthv, dthv_buoyant,           &
                                                height_wth(k),                &
                                                height_wth(k-1) )
          l_lnb = .true.
          cin2 = cin2 + dcape
        else if ( dcape > 0.0_r_def ) then
          cape = cape + dcape
          dthv_buoyant = dthv
        else if ( cin > parcel_min_cin ) then
          cin = cin + dcape
        end if
      else
        ! Possible second buoyant layer above a capping inversion.
        if ( dcape < 0.0_r_def .and. cape2 > parcel_min_cape ) then
          lnb_height2 = neutral_buoyancy_height( dthv, dthv_buoyant,          &
                                                 height_wth(k),               &
                                                 height_wth(k-1) )
          l_second_lnb = .true.
        else if ( dcape > 0.0_r_def ) then
          cape2 = cape2 + dcape
          dthv_buoyant = dthv
        else if ( cin2 > parcel_min_cin ) then
          cin2 = cin2 + dcape
        end if
      end if

      if ( .not. l_second_lnb ) then
        call lift_parcel( t_parcel, theta_parcel, qsat_parcel, p_wth(k),      &
                          p_wth(k+1), exner_wth(k+1), t_lcl, p_lcl, q_lcl,    &
                          rd, cp, kappa )
      end if

    end do

    ! Include a second buoyant layer if the inhibition below it is not too
    ! large.
    if ( cape2 > 0.0_r_def .and. cin + cin2 > parcel_min_cin ) then
      cin  = cin + cin2
      cape = cape + cape2
      lnb_height = lnb_height2
    end if
    cin = max( cin, parcel_min_cin )

  end subroutine parcel_ascent_cape_cin

  !> @brief CAPE and CIN of the most unstable parcel, and the effective
  !>        inflow layer, from ascents starting at screen level and at each
  !>        model level below a pressure limit.
  !> @details The most unstable parcel is the one with the largest CAPE. The
  !>          effective inflow layer runs from the lowest starting level whose
  !>          parcel has enough CAPE and little enough CIN up to the first
  !>          level above it whose parcel does not. Heights are above sea
  !>          level; undefined values are missing data.
  !> @param[in]  nlayers            Number of layers
  !> @param[in]  surface_cape       CAPE of the screen-level parcel (J/kg)
  !> @param[in]  surface_cin        CIN of the screen-level parcel (J/kg)
  !> @param[in]  surface_lnb_height Height of the screen-level parcel's level
  !>                                of neutral buoyancy (m)
  !> @param[in]  theta              Potential temperature on Wtheta levels (K)
  !> @param[in]  m_v                Water vapour mixing ratio on Wtheta levels
  !> @param[in]  exner_wth          Exner pressure on Wtheta levels
  !> @param[in]  p_wth              Pressure on Wtheta levels (Pa)
  !> @param[in]  height_wth         Height of Wtheta levels above sea level (m)
  !> @param[in]  height_w3          Height of W3 levels above sea level (m)
  !> @param[in]  gravity            Gravitational acceleration (m s-2)
  !> @param[in]  rd                 Gas constant for dry air (J kg-1 K-1)
  !> @param[in]  cp                 Specific heat of dry air (J kg-1 K-1)
  !> @param[in]  kappa              Ratio of rd to cp
  !> @param[in]  one_over_kappa     Ratio of cp to rd
  !> @param[in]  epsilon            Ratio of gas constants of dry air and
  !>                                vapour
  !> @param[out] mu_cape            CAPE of the most unstable parcel (J/kg)
  !> @param[out] mu_cin             CIN of the most unstable parcel (J/kg)
  !> @param[out] mu_start_height    Starting height of the most unstable
  !>                                parcel (m)
  !> @param[out] mu_lnb_height      Height of the most unstable parcel's level
  !>                                of neutral buoyancy (m)
  !> @param[out] inflow_base        Base of the effective inflow layer (m)
  !> @param[out] inflow_top         Top of the effective inflow layer (m)
  subroutine most_unstable_cape_cin( nlayers, surface_cape, surface_cin,      &
                                     surface_lnb_height, theta, m_v,          &
                                     exner_wth, p_wth, height_wth, height_w3, &
                                     gravity, rd, cp, kappa, one_over_kappa,  &
                                     epsilon, mu_cape, mu_cin,                &
                                     mu_start_height, mu_lnb_height,          &
                                     inflow_base, inflow_top )

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers
    real(kind=r_def),    intent(in)  :: surface_cape, surface_cin
    real(kind=r_def),    intent(in)  :: surface_lnb_height
    real(kind=r_def),    intent(in)  :: theta(0:nlayers)
    real(kind=r_def),    intent(in)  :: m_v(0:nlayers)
    real(kind=r_def),    intent(in)  :: exner_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: p_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: height_wth(0:nlayers)
    real(kind=r_def),    intent(in)  :: height_w3(0:nlayers-1)
    real(kind=r_def),    intent(in)  :: gravity, rd, cp, kappa
    real(kind=r_def),    intent(in)  :: one_over_kappa, epsilon
    real(kind=r_def),    intent(out) :: mu_cape, mu_cin
    real(kind=r_def),    intent(out) :: mu_start_height, mu_lnb_height
    real(kind=r_def),    intent(out) :: inflow_base, inflow_top

    integer(kind=i_def) :: k, k_top
    real(kind=r_def)    :: cape, cin, lnb_height

    call initialise_most_unstable( surface_cape, surface_cin,                 &
                                   surface_lnb_height, height_wth(0),         &
                                   mu_cape, mu_cin, mu_start_height,          &
                                   mu_lnb_height, inflow_base, inflow_top )

    k_top = highest_start_level( nlayers, p_wth )

    do k = 1, k_top
      call parcel_ascent_cape_cin( nlayers, k, .false._l_def, 0.0_r_def,      &
                                   0.0_r_def, theta(k) * exner_wth(k),        &
                                   theta(k), m_v(k), p_wth(k), theta, m_v,    &
                                   exner_wth, p_wth, height_wth, height_w3,   &
                                   gravity, rd, cp, kappa, one_over_kappa,    &
                                   epsilon, cape, cin, lnb_height )
      call update_most_unstable( k, k_top, cape, cin, lnb_height,             &
                                 height_wth(k), height_wth(k+1),              &
                                 mu_cape, mu_cin, mu_start_height,            &
                                 mu_lnb_height, inflow_base, inflow_top )
    end do

    ! No model level lies below the pressure limit: an inflow layer found
    ! from the screen-level parcel extends one surface-layer depth above its
    ! base.
    if ( k_top == 0 .and. inflow_base /= rmdi ) then
      inflow_top = inflow_base + parcel_surface_layer_depth
    end if

  end subroutine most_unstable_cape_cin

  !> @brief Highest model level that is below the pressure limit for
  !>        starting a most unstable parcel ascent.
  !> @param[in] nlayers Number of layers
  !> @param[in] p_wth   Pressure on Wtheta levels (Pa)
  !> @return    k_top   Highest starting level, or 0 if there is none
  function highest_start_level( nlayers, p_wth ) result( k_top )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: p_wth(0:nlayers)
    integer(kind=i_def) :: k_top

    integer(kind=i_def) :: k

    k_top = 0_i_def
    do k = 1, nlayers - 1
      if ( p_wth(k) <= parcel_top_start_pressure ) exit
      k_top = k
    end do

  end function highest_start_level

  !> @brief Initialise the most unstable parcel and the effective inflow
  !>        layer from the screen-level parcel.
  !> @param[in]  surface_cape       CAPE of the screen-level parcel (J/kg)
  !> @param[in]  surface_cin        CIN of the screen-level parcel (J/kg)
  !> @param[in]  surface_lnb_height Height of its level of neutral buoyancy (m)
  !> @param[in]  surface_height     Height of the surface above sea level (m)
  !> @param[out] mu_cape            CAPE of the most unstable parcel (J/kg)
  !> @param[out] mu_cin             CIN of the most unstable parcel (J/kg)
  !> @param[out] mu_start_height    Starting height of the most unstable
  !>                                parcel (m)
  !> @param[out] mu_lnb_height      Height of the most unstable parcel's level
  !>                                of neutral buoyancy (m)
  !> @param[out] inflow_base        Base of the effective inflow layer (m)
  !> @param[out] inflow_top         Top of the effective inflow layer (m)
  subroutine initialise_most_unstable( surface_cape, surface_cin,             &
                                       surface_lnb_height, surface_height,    &
                                       mu_cape, mu_cin, mu_start_height,      &
                                       mu_lnb_height, inflow_base,            &
                                       inflow_top )

    implicit none

    real(kind=r_def), intent(in)  :: surface_cape, surface_cin
    real(kind=r_def), intent(in)  :: surface_lnb_height, surface_height
    real(kind=r_def), intent(out) :: mu_cape, mu_cin
    real(kind=r_def), intent(out) :: mu_start_height, mu_lnb_height
    real(kind=r_def), intent(out) :: inflow_base, inflow_top

    mu_cape    = surface_cape
    mu_cin     = surface_cin
    inflow_top = rmdi

    if ( surface_cape > parcel_min_cape ) then
      mu_start_height = surface_height + parcel_screen_height
      mu_lnb_height   = surface_lnb_height
      if ( surface_cape > inflow_min_cape .and.                               &
           surface_cin > inflow_min_cin ) then
        inflow_base = surface_height + parcel_screen_height
      else
        inflow_base = rmdi
      end if
    else
      mu_start_height = rmdi
      mu_lnb_height   = rmdi
      inflow_base     = rmdi
    end if

  end subroutine initialise_most_unstable

  !> @brief Update the most unstable parcel and the effective inflow layer
  !>        with the ascent from one model level.
  !> @param[in]     k               Starting level of the ascent
  !> @param[in]     k_top           Highest starting level
  !> @param[in]     cape            CAPE of the ascent (J/kg)
  !> @param[in]     cin             CIN of the ascent (J/kg)
  !> @param[in]     lnb_height      Height of its level of neutral buoyancy
  !>                                (m)
  !> @param[in]     start_height    Height of level k above sea level (m)
  !> @param[in]     next_height     Height of level k + 1 above sea level (m)
  !> @param[in,out] mu_cape         CAPE of the most unstable parcel (J/kg)
  !> @param[in,out] mu_cin          CIN of the most unstable parcel (J/kg)
  !> @param[in,out] mu_start_height Starting height of the most unstable
  !>                                parcel (m)
  !> @param[in,out] mu_lnb_height   Height of the most unstable parcel's level
  !>                                of neutral buoyancy (m)
  !> @param[in,out] inflow_base     Base of the effective inflow layer (m)
  !> @param[in,out] inflow_top      Top of the effective inflow layer (m)
  subroutine update_most_unstable( k, k_top, cape, cin, lnb_height,           &
                                   start_height, next_height,                 &
                                   mu_cape, mu_cin, mu_start_height,          &
                                   mu_lnb_height, inflow_base, inflow_top )

    implicit none

    integer(kind=i_def), intent(in)    :: k, k_top
    real(kind=r_def),    intent(in)    :: cape, cin, lnb_height
    real(kind=r_def),    intent(in)    :: start_height, next_height
    real(kind=r_def),    intent(inout) :: mu_cape, mu_cin
    real(kind=r_def),    intent(inout) :: mu_start_height, mu_lnb_height
    real(kind=r_def),    intent(inout) :: inflow_base, inflow_top

    if ( cape > parcel_min_cape ) then
      if ( cape > mu_cape ) then
        mu_cape         = cape
        mu_cin          = cin
        mu_lnb_height   = lnb_height
        mu_start_height = start_height
      end if
      if ( inflow_base == rmdi .and. cape > inflow_min_cape .and.             &
           cin > inflow_min_cin ) then
        inflow_base = start_height
      end if
    end if

    ! The inflow layer ends at the first level whose parcel fails the inflow
    ! test, or above the highest starting level.
    if ( inflow_base /= rmdi .and. inflow_top == rmdi ) then
      if ( cape <= inflow_min_cape .or. cin <= inflow_min_cin ) then
        inflow_top = start_height
      else if ( k == k_top ) then
        inflow_top = next_height
      end if
    end if

  end subroutine update_most_unstable

  !> @brief Lift a parcel from one pressure to the next.
  !> @details Below the lifting condensation level the parcel conserves its
  !>          potential temperature. Above it, the temperature change along a
  !>          saturated adiabat is found from a second-order expansion in the
  !>          fractional pressure change, using a finite-difference estimate
  !>          of the rate of change of saturation mixing ratio with
  !>          temperature. A step crossing the condensation level starts the
  !>          saturated ascent from that level.
  !> @param[in,out] t_parcel     Parcel temperature (K)
  !> @param[in,out] theta_parcel Parcel potential temperature (K)
  !> @param[in]     qsat_parcel  Saturation mixing ratio at the parcel's
  !>                             current temperature and pressure
  !> @param[in]     p_now        Current parcel pressure (Pa)
  !> @param[in]     p_next       Pressure to lift the parcel to (Pa)
  !> @param[in]     exner_next   Exner pressure at p_next
  !> @param[in]     t_lcl        Temperature at the condensation level (K)
  !> @param[in]     p_lcl        Pressure at the condensation level (Pa)
  !> @param[in]     q_lcl        Saturation mixing ratio at the condensation
  !>                             level
  !> @param[in]     rd           Gas constant for dry air (J kg-1 K-1)
  !> @param[in]     cp           Specific heat of dry air (J kg-1 K-1)
  !> @param[in]     kappa        Ratio of rd to cp
  subroutine lift_parcel( t_parcel, theta_parcel, qsat_parcel, p_now, p_next, &
                          exner_next, t_lcl, p_lcl, q_lcl, rd, cp, kappa )

    implicit none

    real(kind=r_def), intent(inout) :: t_parcel, theta_parcel
    real(kind=r_def), intent(in)    :: qsat_parcel, p_now, p_next, exner_next
    real(kind=r_def), intent(in)    :: t_lcl, p_lcl, q_lcl, rd, cp, kappa

    real(kind=r_def) :: t_base, p_base, q_base, q_cooler, dp_over_p
    real(kind=r_def) :: dqsat_dt, dt_dry, dt

    if ( p_next <= p_lcl ) then
      if ( p_now > p_lcl ) then
        t_base = t_lcl
        p_base = p_lcl
        q_base = q_lcl
      else
        t_base = t_parcel
        p_base = p_now
        q_base = qsat_parcel
      end if
      dp_over_p = ( p_base - p_next ) / p_base
      call qsat_wat_mix( q_cooler, t_base - parcel_qsat_delta_t, p_base )
      dqsat_dt = ( q_base - q_cooler ) / parcel_qsat_delta_t
      dt_dry = -rd * t_base * dp_over_p                                       &
             * ( 1.0_r_def - ( kappa - 1.0_r_def ) * dp_over_p * half )
      dt = ( dt_dry - latent_heat_h2o_condensation * q_base * dp_over_p )     &
         / ( cp + dqsat_dt * latent_heat_h2o_condensation )
      t_parcel     = t_base + dt
      theta_parcel = t_parcel / exner_next
    else
      t_parcel = theta_parcel * exner_next
    end if

  end subroutine lift_parcel

  !> @brief Depth of the layer around a Wtheta level, between the W3 levels
  !>        either side of it (from the surface for the first level).
  !> @param[in] nlayers    Number of layers
  !> @param[in] k          Wtheta level, 1 to nlayers - 1
  !> @param[in] height_wth Height of Wtheta levels (m)
  !> @param[in] height_w3  Height of W3 levels (m)
  !> @return    dz         Layer depth (m)
  function theta_layer_depth( nlayers, k, height_wth, height_w3 ) result( dz )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers, k
    real(kind=r_def),    intent(in) :: height_wth(0:nlayers)
    real(kind=r_def),    intent(in) :: height_w3(0:nlayers-1)
    real(kind=r_def) :: dz

    if ( k == 1 ) then
      dz = height_w3(1) - height_wth(0)
    else
      dz = height_w3(k) - height_w3(k-1)
    end if

  end function theta_layer_depth

  !> @brief Height at which the parcel becomes neutrally buoyant, by linear
  !>        interpolation of its buoyancy between two levels.
  !> @param[in] dthv         Virtual potential temperature excess at the
  !>                         upper level (K, negative)
  !> @param[in] dthv_buoyant Virtual potential temperature excess at the last
  !>                         buoyant level (K, positive)
  !> @param[in] z_upper      Height of the upper level (m)
  !> @param[in] z_lower      Height of the level below (m)
  !> @return    z_lnb        Height of neutral buoyancy (m)
  function neutral_buoyancy_height( dthv, dthv_buoyant, z_upper, z_lower )    &
                                    result( z_lnb )

    implicit none

    real(kind=r_def), intent(in) :: dthv, dthv_buoyant, z_upper, z_lower
    real(kind=r_def) :: z_lnb

    z_lnb = z_upper - dthv / ( dthv - dthv_buoyant ) * ( z_upper - z_lower )

  end function neutral_buoyancy_height

  !> @brief Virtual potential temperature.
  !> @param[in] theta   Potential temperature (K)
  !> @param[in] m_v     Water vapour mixing ratio
  !> @param[in] epsilon Ratio of gas constants of dry air and vapour
  !> @return    thv     Virtual potential temperature (K)
  function virtual_theta( theta, m_v, epsilon ) result( thv )

    implicit none

    real(kind=r_def), intent(in) :: theta, m_v, epsilon
    real(kind=r_def) :: thv

    thv = theta * ( 1.0_r_def + m_v / epsilon ) / ( 1.0_r_def + m_v )

  end function virtual_theta

end module parcel_ascent_cape_cin_mod
