!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to UKCA Photolysis scheme

module photol_kernel_mod

use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_INTEGER,        &
                             GH_REAL, GH_READWRITE, GH_READ,                   &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4,                        &
                             DOMAIN

use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type

use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_ERROR, &
                             LOG_LEVEL_INFO

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: photol_kernel_type
  private
  type(arg_type) :: meta_args(37) = (/                      &
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1 ), & ! Photol_rates
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_year
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_month
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_day
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_hour
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_minute
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_second
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_daynum       
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_hour
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_minute
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_second       
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! o3     
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! theta_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! exner_in_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! height_in_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! m_cl_n
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! m_cf_n
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! rel_humid_frac       
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! conv_cloud_amount
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! area_cloud_frac       
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! sulp_accum
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! sulp_aitken
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! aod_sulp_accum
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! aod_sulp_aitken
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! photol_rad_jo2
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! photol_rad_jo2b
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! height_in_w3
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! exner_in_w3
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3 ), & ! tile_fraction
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! latitude
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! longitude       
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! sin_stellar_declination_rts
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! stellar_eqn_of_time_rts       
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! conv_cloud_lwp
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! surf_albedo
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! conv_cloud_base
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 )  & ! conv_cloud_top
       /)
  integer :: operates_on = DOMAIN

contains
  procedure, nopass :: photol_code
end type

public photol_code

contains

!> @brief UKCA Photolysis scheme time step.
!> @details Copy or create the environmental driver fields required for the 
!>          current Photolysis configuration, perform a photolysis timestep
!>          and finally copy the rates array back as a LFRic field.

!> @param[in]     nlayers             Number of layers
!> @param[in,out] photol_rates        Photolysis rates (s-1)
!> @param[in]     current_time_year   Current model year
!> @param[in]     current_time_month  Current model month
!> @param[in]     current_time_day    Current model day
!> @param[in]     current_time_hour   Current model hour
!> @param[in]     current_time_minute Current model minute
!> @param[in]     current_time_second Current model second
!> @param[in]     current_time_daynum Current model day number
!> @param[in]     previous_time_hour  Model hour at previous time step
!> @param[in]     previous_time_minute Model minute at previous time step
!> @param[in]     previous_time_second Model second at previous time step
!> @param[in]     o3                  Ozone m.m.r.
!> @param[in]     theta_wth           Potential temperature field (K)
!> @param[in]     exner_in_wth        Exner pressure in theta space
!> @param[in]     height_wth          Height of theta space above surface (m)
!> @param[in]     m_cl_n              Cloud liq mixing ratio at time level n
!> @param[in]     m_cf_n              Cloud frozen mixing ratio at time level n
!> @param[in]     rel_humid_frac      Relative Humidity fraction (0.0-1.0)
!> @param[in]     conv_cloud_amount   Convective clound amount (kg m-2)
!> @param[in]     area_cloud_frac     Area cloud fraction
!> @param[in]     sulp_accum          Sulphate accumulation mode m.m.r (kg/kg)
!> @param[in]     sulp_aiken          Sulphate aikten mode m.m.r (kg/kg)
!> @param[in]     aod_sulp_accum      Optical depth for Sulphate accumulation mode (1)
!> @param[in]     aod_sulp_aiken      Optical depth for Sulphate aikten mode (1)
!> @param[in]     photol_rad_jo2      Rate : O2 => O3P + O3P reaction from Radiation (s-1)
!> @param[in]     photol_rad_jo2b     Rate : O2 => O3P + O1D reaction from Radiation (s-1)
!> @param[in]     height_w3           Height of density space above surface (m)
!> @param[in]     exner_in_w3         Exner pressure in w3 space
!> @param[in]     tile_fraction       Surface tile fractions
!> @param[in]     latitude            Latitude field (radians)
!> @param[in]     longitude           Longitude field (radians)
!> @param[in]     sin_stellar_declination_rts Sin of stellar declination
!> @param[in]     stellar_eqn_of_time Stellar equation of time (radians)
!> @param[in]     conv_cloud_lwp      Convective cloud Liq water path
!> @param[in]     surf_albedo         Surface albedo (0.0-1.0)
!> @param[in]     conv_cloud_base     (Model level for) convective cloud base
!> @param[in]     conv_cloud_top      (Model level for) convective cloud top
!> @param[in]     ndf_nphot           Number of DOFs per cell for photolysis species
!> @param[in]     undf_nphot          Number of unique DOFs for photolysis species
!> @param[in]     map_nphot           Dofmap for cell for photolysis species
!> @param[in]     ndf_wth             Number of DOFs per cell for potential temperature space
!> @param[in]     undf_wth            Number of unique DOFs for potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the column for potential temperature space
!> @param[in]     ndf_w3              Number of DOFs per cell for density space
!> @param[in]     undf_w3             Number of unique DOFs for density space
!> @param[in]     map_w3              Dofmap for the cell at the base of the column for density space
!> @param[in]     ndf_tile            Number of DOFs per cell for surface tiles
!> @param[in]     undf_tile           Number of unique DOFs for surface tiles
!> @param[in]     map_tile            Dofmap for cell for surface tiles
!> @param[in]     ndf_2d              Number of DOFs per cell for 2D fields
!> @param[in]     undf_2d             Number of unique DOFs for 2D fields
!> @param[in]     map_2d              Dofmap for cell for 2D fields


subroutine photol_code( nlayers,                                               &
                        seg_len,                                               &
                        photol_rates,                                          &
                        current_time_year,                                     &
                        current_time_month,                                    &
                        current_time_day,                                      &
                        current_time_hour,                                     &
                        current_time_minute,                                   &
                        current_time_second,                                   &
                        current_time_daynum,                                   &
                        previous_time_hour,                                    &
                        previous_time_minute,                                  &
                        previous_time_second,                                  &
                        o3,                                                    &
                        theta_wth,                                             &
                        exner_in_wth,                                          &
                        height_wth,                                            &
                        m_cl_n,                                                &
                        m_cf_n,                                                &
                        rel_humid_frac,                                        &
                        conv_cloud_amount,                                     &
                        area_cloud_frac,                                       &
                        sulp_accum,                                            &
                        sulp_aitken,                                           &
                        aod_sulp_accum,                                        &
                        aod_sulp_aitken,                                       &
                        photol_rad_jo2,                                        &
                        photol_rad_jo2b,                                       &
                        height_w3,                                             &
                        exner_in_w3,                                           &
                        tile_fraction,                                         &
                        latitude,                                              &
                        longitude,                                             &
                        sin_stellar_declination_rts,                           &
                        stellar_eqn_of_time_rts,                               &
                        conv_cloud_lwp,                                        &
                        surf_albedo,                                           &
                        conv_cloud_base,                                       &
                        conv_cloud_top,                                        &
                        ndf_nphot, undf_nphot, map_nphot,                      &
                        ndf_wth, undf_wth, map_wth,                            &
                        ndf_w3, undf_w3, map_w3,                               &
                        ndf_tile, undf_tile, map_tile,                         &
                        ndf_2d, undf_2d, map_2d)

  use constants_mod,    only: r_def, i_def, r_um, i_um, i_timestep,            &
                              radians_to_degrees
  use conversions_mod,  only: rsec_per_day, rsec_per_hour
  use jules_control_init_mod, only: n_surf_tile, n_land_tile
  use um_ukca_init_mod, only: n_phot_spc, n_phot_flds_req,                     &
                              ratj_data, ratj_varnames,                        &
                              photol_fldnames_scalar_real,                     &
                              photol_fldnames_flat_integer,                    &
                              photol_fldnames_flat_real,                       &
                              photol_fldnames_fullht_real,                     &
                              photol_fldnames_fullht0_real,                    &
                              photol_fldnames_fullhtphot_real
  
  use chemistry_config_mod, only: chem_scheme, chem_scheme_strattrop,          &
                              photol_scheme, photol_scheme_fastjx

  ! Photolysis API module
  USE photol_api_mod,     ONLY: photol_fastjx, photol_fieldname_len,           &
                              photol_fldname_aod_sulph_aitk,                   &
                              photol_fldname_aod_sulph_accum,                  &
                              photol_fldname_area_cloud_fraction,              &
                              photol_fldname_conv_cloud_amount,                &
                              photol_fldname_conv_cloud_base,                  &
                              photol_fldname_conv_cloud_lwp,                   &
                              photol_fldname_conv_cloud_top,                   &
                              photol_fldname_cos_latitude,                     &
                              photol_fldname_equation_of_time,                 &
                              photol_fldname_land_fraction,                    &
                              photol_fldname_longitude,                        &
                              photol_fldname_ozone_mmr,                        &
                              photol_fldname_p_layer_boundaries,               &
                              photol_fldname_p_theta_levels,                   &
                              photol_fldname_qcf, photol_fldname_qcl,          &
                              photol_fldname_rad_ctl_jo2,                      &
                              photol_fldname_rad_ctl_jo2b,                     &
                              photol_fldname_r_rho_levels,                     &
                              photol_fldname_r_theta_levels,                   &
                              photol_fldname_sec_since_midnight,               &
                              photol_fldname_sin_declination,                  &
                              photol_fldname_sin_latitude,                     &
                              photol_fldname_so4_aitken,                       &
                              photol_fldname_so4_accum,                        &
                              photol_fldname_surf_albedo,                      &
                              photol_fldname_t_theta_levels,                   &
                              photol_fldname_tan_latitude,                     &
                              photol_fldname_z_top_of_model,                   &
                              photol_step_control
  ! UM modules  
  use planet_constants_mod, only: p_zero, kappa, planet_radius

  ! UKCA API module
  use ukca_error_mod,         only: ukca_maxlen_message => maxlen_message,     &
                                    ukca_maxlen_procname => maxlen_procname

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: seg_len
  integer(kind=i_def), intent(in) :: ndf_nphot
  integer(kind=i_def), intent(in) :: undf_nphot
  integer(kind=i_def), dimension(ndf_nphot, seg_len), intent(in) :: map_nphot
  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth, seg_len), intent(in) :: map_wth
  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), dimension(ndf_w3, seg_len), intent(in) :: map_w3
  integer(kind=i_def), intent(in) :: ndf_tile
  integer(kind=i_def), intent(in) :: undf_tile
  integer(kind=i_def), dimension(ndf_tile, seg_len), intent(in) :: map_tile  
  integer(kind=i_def), intent(in) :: ndf_2d
  integer(kind=i_def), intent(in) :: undf_2d
  integer(kind=i_def), dimension(ndf_2d, seg_len), intent(in) :: map_2d  

  real(kind=r_def), intent(in out), dimension(undf_nphot) :: photol_rates

  integer(kind=i_def), intent(in) :: current_time_year
  integer(kind=i_def), intent(in) :: current_time_month
  integer(kind=i_def), intent(in) :: current_time_day
  integer(kind=i_def), intent(in) :: current_time_hour
  integer(kind=i_def), intent(in) :: current_time_minute
  integer(kind=i_def), intent(in) :: current_time_second
  integer(kind=i_def), intent(in) :: current_time_daynum
  integer(kind=i_def), intent(in) :: previous_time_hour
  integer(kind=i_def), intent(in) :: previous_time_minute
  integer(kind=i_def), intent(in) :: previous_time_second

  real(kind=r_def), intent(in), dimension(undf_wth) :: o3
  real(kind=r_def), intent(in), dimension(undf_wth) :: theta_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: exner_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: m_cl_n
  real(kind=r_def), intent(in), dimension(undf_wth) :: m_cf_n
  real(kind=r_def), intent(in), dimension(undf_wth) :: rel_humid_frac
  real(kind=r_def), intent(in), dimension(undf_wth) :: conv_cloud_amount
  real(kind=r_def), intent(in), dimension(undf_wth) :: area_cloud_frac
  real(kind=r_def), intent(in), dimension(undf_wth) :: sulp_accum
  real(kind=r_def), intent(in), dimension(undf_wth) :: sulp_aitken
  real(kind=r_def), intent(in), dimension(undf_wth) :: aod_sulp_accum
  real(kind=r_def), intent(in), dimension(undf_wth) :: aod_sulp_aitken
  real(kind=r_def), intent(in), dimension(undf_wth) :: photol_rad_jo2
  real(kind=r_def), intent(in), dimension(undf_wth) :: photol_rad_jo2b
  real(kind=r_def), intent(in), dimension(undf_w3) :: height_w3
  real(kind=r_def), intent(in), dimension(undf_w3) :: exner_in_w3

  real(kind=r_def), intent(in), dimension(undf_tile) :: tile_fraction  
  real(kind=r_def), intent(in), dimension(undf_2d) :: latitude
  real(kind=r_def), intent(in), dimension(undf_2d) :: longitude  
  real(kind=r_def), intent(in), dimension(undf_2d) ::                          &
    sin_stellar_declination_rts
  real(kind=r_def), intent(in), dimension(undf_2d) :: stellar_eqn_of_time_rts  
  real(kind=r_def), intent(in), dimension(undf_2d) :: conv_cloud_lwp
  real(kind=r_def), intent(in), dimension(undf_2d) :: surf_albedo
  real(kind=r_def), intent(in), dimension(undf_2d) :: conv_cloud_base
  real(kind=r_def), intent(in), dimension(undf_2d) :: conv_cloud_top
  
  ! Local variables for the kernel

  ! Current model time (year, month, day, hour, minute, second, day number)
  integer(i_um) :: current_time(7)
  ! Seconds since midnight
  real(r_um) :: sec_since_midnight

  ! Groups of environment driving fields 
  ! Scalars, dimension (m_fields)
  real(r_um), allocatable :: fldgroup_scalar_real(:)

  ! Flat; dimensions (seg_len, 1, m_fields)
  integer(i_um), allocatable :: fldgroup_flat_integer(:,:,:)
  real(r_um), allocatable :: fldgroup_flat_real(:,:,:)

  ! Dimensions: (seg_len, 1, nlayers, m_fields)
  real(r_um), allocatable :: fldgroup_fullht_real(:,:,:,:)
  ! Dimensions: (seg_len, 1, 0:nlayers, m_fields)
  real(r_um), allocatable :: fldgroup_fullht0_real(:,:,:,:)
  ! Dimensions: (seg_len, 1, nlayers, n_phot, M=1)
  real(r_um), allocatable :: fldgroup_fullhtphot_real(:,:,:,:,:)  

  ! Environmental driver fields to be calculated
  real(r_um) :: r_theta_levels(seg_len,1,0:nlayers)
  real(r_um) :: r_rho_levels(seg_len,1,nlayers)
  real(r_um) :: t_theta_levels(seg_len,1,nlayers)
  real(r_um) :: p_theta_levels(seg_len,1,nlayers)
  real(r_um) :: p_layer_boundaries(seg_len,1,0:nlayers)

  real(r_um) :: frac_land(seg_len,1)   ! Land fraction in cell
  
  ! Photolysis rates in UM type, to be passed to control routine
  real(r_um), allocatable :: photol_rates_loc(:,:,:,:)

  real(r_um) :: z_top_of_model ! Height (distance from surface) at
                               ! model top
  ! Working variables
  integer(i_um) :: m_fields   ! Number of fields in a group
  integer(i_um) :: n_land_pts ! Number of land points
  integer(i_um) :: i          ! Model horizontal loop counter  
  integer(i_um) :: k          ! Model level loop counter
  integer(i_um) :: m          ! Fields Loop counter
  integer(i_um) :: n          ! Fields counter
  integer(i_um) :: jp1        ! Loop counter for photolytic species
  integer(i_um) :: error_code

  ! Array to hold names of fields actually being passed to photolysis
  character(len=photol_fieldname_len), allocatable :: env_fldnames_avail(:)

  ! UKCA error reporting variables
  character(len=ukca_maxlen_message)  :: photol_errmsg  ! Error return message
  character(len=ukca_maxlen_procname) :: photol_errproc ! Routine in which error
                                                      ! was trapped
  !-----------------------------------------------------------------------
  ! Derive required Environmental Drivers and add to field groups
  !-----------------------------------------------------------------------
  ! list of environment drivers being provided.
  ! n_phot_flds_req is count of fields actually required for this configuration.
  if ( .not. allocated(env_fldnames_avail)) then
    allocate(env_fldnames_avail(n_phot_flds_req))
    env_fldnames_avail(:) = ''
  end if

  ! Collate time variables
  current_time(1) = int( current_time_year, i_um )
  current_time(2) = int( current_time_month, i_um )
  current_time(3) = int( current_time_day, i_um )
  current_time(4) = int( current_time_hour, i_um )
  current_time(5) = int( current_time_minute, i_um )
  current_time(6) = int( current_time_second, i_um )
  current_time(7) = int( current_time_daynum, i_um )
  
  sec_since_midnight = real( (previous_time_hour * 3600_i_def +                &
                              previous_time_minute * 60_i_def +                &
                              previous_time_second), r_um)

  ! Derive fields required for more than one photolysis scheme
  ! Temperature on theta levels
  do i = 1, seg_len
    do k = 1, nlayers
      t_theta_levels(i,1,k) = real( exner_in_wth( map_wth(1,i) + k ), r_um )   &
                              * real( theta_wth( map_wth(1,i) + k ), r_um )
    end do
  end do
  ! Pressure on theta levels
  do i = 1, seg_len
    do k = 1, nlayers
      p_theta_levels(i,1,k) = p_zero *                                         &
          real( exner_in_wth( map_wth(1,i) + k ), r_um )**( 1.0_r_um / kappa )
    end do
  end do
  ! Pressure at level boundaries
  do i = 1, seg_len
    ! At surface == pstar
    p_layer_boundaries(i,1,0) = p_zero *                                       &
          real( exner_in_wth(map_wth(1,i) + 0), r_um ) ** ( 1.0_r_um / kappa )

    ! At levels < top == p on rho levels
    do k = 1, nlayers-1
      p_layer_boundaries(i,1,k) = p_zero *                                     &
          real( exner_in_w3( map_w3(1,i) + k-1 ), r_um )**( 1.0_r_um / kappa )
    end do

    ! At topmost level == p_theta_level
    p_layer_boundaries(i,1,nlayers) = p_theta_levels(i,1,nlayers)
  end do

  z_top_of_model = real( height_wth( map_wth(1,1) + nlayers ), r_um )

  ! processing required for fast-jx scheme inputs
  if ( photol_scheme == photol_scheme_fastjx ) then
    ! height of theta levels from centre of planet
    do i = 1, seg_len
      do k = 0, nlayers
        r_theta_levels( i, 1, k ) =                                            &
           real( height_wth( map_wth(1,i) + k ), r_um ) + planet_radius
      end do
    end do    

    ! height of rho levels from centre of planet
    do i = 1, seg_len
      do k = 1, nlayers
        r_rho_levels( i, 1, k ) =                                              &
           real( height_w3( map_w3(1,i) + k - 1 ), r_um ) + planet_radius
      end do
    end do

    ! land fraction (by summing over tiles)
    frac_land(:,:) = 0.0_r_um
    do i = 1, seg_len
      do m = 1, n_land_tile
        frac_land( i, 1 ) = frac_land( i, 1 ) +                                &
          real( tile_fraction( map_tile(1,i) + m - 1 ), r_um )
      end do
    end do

  else
    r_theta_levels(:,:,:) = 0.0_r_um
    r_rho_levels(:,:,:) = 0.0_r_um
    frac_land(:,:) = 0.0_r_um
  end if  ! if fast-jx scheme
 
  ! Add the driving fields into field groups, based on list of names received
  ! via the photolysis api (in um_ukca_init). Also, append the name of the
  ! added field to list of env_fieldnames passed separately.
  n = 0   ! Overall field counter

  ! Populate FLDGROUP_SCALAR_REAL
  ! Note : sin_stellar_declination and eqn_of_time are scalars
  ! but have been indexed with map_2d(1,1)

  m_fields = size(photol_fldnames_scalar_real)
  if (m_fields > 0 ) then
    allocate(fldgroup_scalar_real(m_fields))
    fldgroup_scalar_real(:) = 0.0_r_um
    do m = 1, m_fields      
      select case(photol_fldnames_scalar_real(m))
      case(photol_fldname_equation_of_time)
        ! Eqn of time
        fldgroup_scalar_real(m) = real(stellar_eqn_of_time_rts(map_2d(1,1)), r_um)        
        call append_fieldname(n, photol_fldname_equation_of_time,              &
                              env_fldnames_avail)
      case(photol_fldname_sec_since_midnight)
        ! Seconds since midnight
        fldgroup_scalar_real(m) = sec_since_midnight
        call append_fieldname(n,photol_fldname_sec_since_midnight,             &
                              env_fldnames_avail)
      case(photol_fldname_sin_declination)
        ! Sin declination angle
        fldgroup_scalar_real(m) =                                              &
            real( sin_stellar_declination_rts( map_2d(1,1) ), r_um )        
        call append_fieldname(n,photol_fldname_sin_declination,                &
                              env_fldnames_avail)
      case(photol_fldname_z_top_of_model)
        ! Top of model (height from surface)
        fldgroup_scalar_real(m) = z_top_of_model        
        call append_fieldname(n,photol_fldname_z_top_of_model,                 &
                              env_fldnames_avail)
      case default
        write(log_scratch_space, '(a,a)')                                      &
          'Unknown environment field in FLDGROUP_SCALAR_REAL',                 &
          photol_fldnames_scalar_real(m)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end select

    end do ! m_fields
  else
    ! Allocate to minimal size; last dimension = 0
    allocate(fldgroup_scalar_real(0))
  end if

  ! Populate FLDGROUP_FLAT_INTEGER
  m_fields = size(photol_fldnames_flat_integer)
  if (m_fields > 0 ) then
    allocate(fldgroup_flat_integer(seg_len,1,m_fields))
    fldgroup_flat_integer(:,:,:) = 0_i_um
    do m = 1, m_fields      
      select case(photol_fldnames_flat_integer(m))
      case(photol_fldname_conv_cloud_base)
        ! Base (model level) for convective cloud
        do i = 1, seg_len
          fldgroup_flat_integer(i,1,m) = int(conv_cloud_base(map_2d(1,i)), i_um)
        end do
        call append_fieldname(n,photol_fldname_conv_cloud_base,              &
                              env_fldnames_avail)        
      case(photol_fldname_conv_cloud_top)
        ! Top (model level) for convective cloud
        do i = 1, seg_len
          fldgroup_flat_integer(i,1,m) = int(conv_cloud_top(map_2d(1,i)), i_um)
        end do
        call append_fieldname(n,photol_fldname_conv_cloud_top,              &
                              env_fldnames_avail)
      case default
        write(log_scratch_space, '(a,a)')     &
          'Unknown environment field in FLDGROUP_FLAT_INTEGER',  &
          photol_fldnames_flat_integer(m)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end select

    end do ! m_fields
  else
    allocate(fldgroup_flat_integer(1,1,0))  
  end if

  ! Populate FLDGROUP_FLAT_REAL
  m_fields = size(photol_fldnames_flat_real)
  if (m_fields > 0 ) then
    allocate(fldgroup_flat_real(seg_len,1,m_fields))
    fldgroup_flat_real(:,:,:) = 0.0_r_um
    do m = 1, m_fields      
      select case(photol_fldnames_flat_real(m))
      case(photol_fldname_conv_cloud_lwp)
        ! Convective cloud liquid water path
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = real(conv_cloud_lwp(map_2d(1,i)), r_um)
        end do
        call append_fieldname(n,photol_fldname_conv_cloud_lwp,              &
                              env_fldnames_avail)
      case(photol_fldname_cos_latitude)
        ! COSine of latitude
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = real( cos( latitude(map_2d(1,i)) ), r_um )
        end do
        call append_fieldname(n,photol_fldname_cos_latitude,              &
                              env_fldnames_avail)
      case(photol_fldname_land_fraction)
        ! Fraction of land in grid box
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = frac_land(i, 1)          
        end do
        call append_fieldname(n,photol_fldname_land_fraction,              &
                              env_fldnames_avail)
      case(photol_fldname_longitude)
        ! Longitude (degrees)
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) =                                          &
               real( radians_to_degrees * longitude(map_2d(1,i)), r_um )
          if (fldgroup_flat_real(i,1,m) < 0.0_r_um) then
            fldgroup_flat_real(i,1,m) = fldgroup_flat_real(i,1,m) + 360.0_r_um
          end if
        end do
        call append_fieldname(n,photol_fldname_longitude,              &
                              env_fldnames_avail)
      case(photol_fldname_sin_latitude)
        ! Sine of latitude
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = real( sin( latitude(map_2d(1,i)) ), r_um )
        end do
        call append_fieldname(n,photol_fldname_sin_latitude,              &
                              env_fldnames_avail)
      case(photol_fldname_surf_albedo)
        ! Surface albedo
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = real(surf_albedo(map_2d(1,i)), r_um)
        end do
        call append_fieldname(n,photol_fldname_surf_albedo,              &
                              env_fldnames_avail)
      case(photol_fldname_tan_latitude)
        ! TAN of latitude
        do i = 1, seg_len
          fldgroup_flat_real(i,1,m) = real( tan( latitude(map_2d(1,i)) ), r_um )
        end do
        call append_fieldname(n,photol_fldname_tan_latitude,              &
                              env_fldnames_avail)
      case default
        write(log_scratch_space, '(a,a)')     &
          'Unknown environment field in FLDGROUP_FLAT_REAL',  &
          photol_fldnames_flat_real(m)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end select

    end do ! m_fields
  else
    allocate(fldgroup_flat_real(1,1,0))  
  end if

  ! Populate FLDGROUP_FULLHT_REAL
  m_fields = size(photol_fldnames_fullht_real)
  if (m_fields > 0 ) then
    allocate(fldgroup_fullht_real(seg_len,1,nlayers,m_fields))
    fldgroup_fullht_real(:,:,:,:) = 0.0_r_um
    do m = 1, m_fields      
      select case(photol_fldnames_fullht_real(m))
      case(photol_fldname_aod_sulph_aitk)
        ! Optical depth from Sulphate in aitken mode
        do i = 1, seg_len
           do k = 1, nlayers
          fldgroup_fullht_real(i,1,k,m) =                                      &
                              real(aod_sulp_aitken( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_aod_sulph_aitk,              &
                              env_fldnames_avail)
      case(photol_fldname_aod_sulph_accum)
        ! Optical depth from Sulphate in accumulation mode
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                    &
                             real(aod_sulp_accum( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_aod_sulph_accum,              &
                              env_fldnames_avail)
      case(photol_fldname_area_cloud_fraction)
        ! Area Cloud fraction
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                   &
                              real(area_cloud_frac( map_wth(1,i) + k ), r_um)          
          end do
        end do
        call append_fieldname(n,photol_fldname_area_cloud_fraction,            &
                              env_fldnames_avail)
      case(photol_fldname_conv_cloud_amount)
        ! Convective cloud amount
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                      &
                             real(conv_cloud_amount( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_conv_cloud_amount,              &
                              env_fldnames_avail)
      case(photol_fldname_ozone_mmr)
        ! Ozone mass mixing ratio
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = real(o3( map_wth(1,i) + k ), r_um)            
          end do
        end do
        call append_fieldname(n,photol_fldname_ozone_mmr,              &
                              env_fldnames_avail)
      case(photol_fldname_p_theta_levels)
        ! Pressure on theta levels
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = p_theta_levels(i, 1, k )            
          end do
        end do
        call append_fieldname(n,photol_fldname_p_theta_levels,              &
                              env_fldnames_avail)
      case(photol_fldname_qcf)
        ! Cloud frozen (ice) mixing ratio
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = real(m_cf_n( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_qcf,              &
                              env_fldnames_avail)
      case(photol_fldname_qcl)
        ! Cloud liquid mixing ratio
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = real(m_cl_n( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_qcl,env_fldnames_avail)
      case(photol_fldname_rad_ctl_jo2)
        ! Rate of O2 => O3P + O3P Photolysis reaction from Radiation scheme
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                &
                          real(photol_rad_jo2( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_rad_ctl_jo2,              &
                              env_fldnames_avail)
      case(photol_fldname_rad_ctl_jo2b)
        ! Rate of O2 => O1D + O3P Photolysis reaction from Radiation scheme
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                 &
                       real(photol_rad_jo2b( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_rad_ctl_jo2b,              &
                              env_fldnames_avail)
      case(photol_fldname_r_rho_levels)
        ! Distance of middle of level from planet centre
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = r_rho_levels( i, 1, k )
          end do
        end do
        call append_fieldname(n,photol_fldname_r_rho_levels,              &
                              env_fldnames_avail)      
      case(photol_fldname_so4_accum)
        ! Mass mix ratio of sulphate aeroosol in accumulation mode
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                 &
                                real(sulp_accum( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_so4_accum,                 &
                              env_fldnames_avail)
      case(photol_fldname_so4_aitken)
        ! Mass mix ratio of sulphate aeroosol in aitken mode
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) =                                 &
                              real(sulp_aitken( map_wth(1,i) + k ), r_um)
          end do
        end do
        call append_fieldname(n,photol_fldname_so4_aitken,              &
                              env_fldnames_avail)
      case(photol_fldname_t_theta_levels)
        ! Temperature on theta levels
        do i = 1, seg_len
          do k = 1, nlayers
            fldgroup_fullht_real(i,1,k,m) = t_theta_levels( i, 1, k )           
          end do
        end do
        call append_fieldname(n,photol_fldname_t_theta_levels,              &
                              env_fldnames_avail)
      case default
        write(log_scratch_space, '(a,a)')     &
          'Unknown environment field in FLDGROUP_FULLHT_REAL',  &
          photol_fldnames_fullht_real(m)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end select

    end do ! m_fields
  else
    allocate(fldgroup_fullht_real(1,1,1,0))  
  end if

  ! Populate FLDGROUP_FULLHT0_REAL
  m_fields = size(photol_fldnames_fullht0_real)
  if (m_fields > 0 ) then
    allocate(fldgroup_fullht0_real(seg_len,1,0:nlayers,m_fields))
    fldgroup_fullht0_real(:,:,:,:) = 0.0_r_um
    do m = 1, m_fields      
      select case(photol_fldnames_fullht0_real(m))
      case(photol_fldname_p_layer_boundaries)
        ! Pressure on level boundaries
        do i = 1, seg_len
          do k = 0, nlayers
            fldgroup_fullht0_real(i,1,k,m) = p_layer_boundaries( i, 1, k )            
          end do
        end do
        call append_fieldname(n,photol_fldname_p_layer_boundaries,              &
                              env_fldnames_avail)
      case(photol_fldname_r_theta_levels)
        ! Distance of level from planet centre
        do i = 1, seg_len
          do k = 0, nlayers
            fldgroup_fullht0_real(i,1,k,m) = r_theta_levels( i, 1, k )          
          end do
        end do
        call append_fieldname(n,photol_fldname_r_theta_levels,              &
                              env_fldnames_avail)
      case default
        write(log_scratch_space, '(a,a)')     &
          'Unknown environment field in FLDGROUP_FULLHT0_REAL',  &
          photol_fldnames_fullht0_real(m)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end select

    end do ! m_fields
  else
    allocate(fldgroup_fullht0_real(1,1,1,0))
  end if

  ! call wrapper routine for photolysis and pass required field groups
  ! Initialise UM copy of rates array first
  if (.not. allocated(photol_rates_loc)) then
    allocate(photol_rates_loc(seg_len,1,nlayers,n_phot_spc))
  end if
  photol_rates_loc(:,:,:,:) = 0.0_r_um
 
  call photol_step_control(current_time, seg_len, 1, nlayers, n_phot_spc,      &
                         ratj_data, ratj_varnames, error_code, photol_rates_loc,&
                         ! names of environment fields provided
                         envfield_names_in=env_fldnames_avail,                 &
                         ! environment field groups
                         envgroup_scalar_real=fldgroup_scalar_real,            &
                         envgroup_flat_integer=fldgroup_flat_integer,          &
                         envgroup_flat_real=fldgroup_flat_real,                &
                         envgroup_fullht_real=fldgroup_fullht_real,            &
                         envgroup_fullht0_real=fldgroup_fullht0_real,          &
                         ! error return info
                    error_message=photol_errmsg, error_routine=photol_errproc)

  if (error_code > 0) then
    write( log_scratch_space, '(A)' )                                          &
     trim(photol_errmsg) // ' in Photolyis routine ' // trim(photol_errproc)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Transfer rates from UM type to LFRic array
  do i = 1, seg_len
    do k = 1, nlayers
      do jp1 = 1, n_phot_spc
        photol_rates( map_nphot(1, i) + ( (jp1-1)*(nlayers+1) ) + k ) =        &
            photol_rates_loc(i, 1, k, jp1)
      end do
    end do
  end do        
  ! deallocate field groups and local arrays
  if (allocated(photol_rates_loc)) deallocate(photol_rates_loc)
  if (allocated(fldgroup_fullht0_real)) deallocate(fldgroup_fullht0_real)
  if (allocated(fldgroup_fullht_real)) deallocate(fldgroup_fullht_real)
  if (allocated(fldgroup_flat_real)) deallocate(fldgroup_flat_real)
  if (allocated(fldgroup_flat_integer)) deallocate(fldgroup_flat_integer)
  if (allocated(fldgroup_scalar_real)) deallocate(fldgroup_scalar_real)

  if ( allocated(env_fldnames_avail)) deallocate(env_fldnames_avail)

return

end subroutine photol_code

! ----------------------------------------------------------------------
! Helper routine to append name of added driving field to names array
! ----------------------------------------------------------------------
subroutine append_fieldname(x,fldname_in, fldnames_avail)

implicit none

integer, intent(in out) :: x
character(len=*), intent(in) :: fldname_in
character(len=*), intent(out) :: fldnames_avail(:)

x = x + 1
if ( x > size(fldnames_avail) ) then
  write(log_scratch_space, '(a,a)')                                            &
    'ERROR from photol_kernel_mod: append_filename: ',                         &
    ' Index exceeds env_fldnames_avail arraysize '
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
ELSE  
  ! Append to list of fieldnames being provided.
  fldnames_avail(x) = fldname_in
END IF

return
end subroutine append_fieldname

end module photol_kernel_mod
