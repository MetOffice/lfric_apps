!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM spectral gravity wave drag scheme

module spectral_gwd_kernel_mod

  use argument_mod,             only: arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_READ, GH_WRITE, &
                                      DOMAIN,            &
                                      ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,            only: r_def, i_def, r_um, i_um
  use empty_data_mod,           only: empty_real_data
  use fs_continuity_mod,        only: W3, Wtheta
  use kernel_mod,               only: kernel_type
  use spectral_gwd_config_mod,  only: ussp_heating

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spectral_gwd_kernel_type
    private
    type(arg_type) :: meta_args(16) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! du_spectral_gwd, u wind increment
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! dv_spectral_gwd, v wind increment
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! dtemp_spectral_gwd, temperature increment
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! u_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! v_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! wetrho_in_w3, wet density in w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta, theta in wtheta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! height_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wtheta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! totalppn
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! latitude
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! tau_east_spectral_gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! tau_south_spectral_gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! tau_west_spectral_gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta)                     & ! tau_north_spectral_gwd
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: spectral_gwd_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: spectral_gwd_code

contains

  !> @brief Call the UM spectral gravity wave drag scheme
  !> @details This code calls the UM USSP spectral gravity wave drag scheme, which
  !>          calculates the zonal and meridional winds and temperature increments
  !>          from parametrized non-orographic gravity wave drag.
  !> @param[in]     nlayers               Integer the number of layers
  !> @param[in,out] du_spectral_gwd       'Zonal' wind increment from spectral gravity wave drag
  !> @param[in,out] dv_spectral_gwd       'Meridional' wind increment from spectral gravity wave drag
  !> @param[in,out] dtemp_spectral_gwd    Theta increment from spectral gravity wave drag
  !> @param[in]     u_in_w3               'Zonal' wind in density space
  !> @param[in]     v_in_w3               'Meridional' wind in density space
  !> @param[in]     wetrho_in_w3           Wet density in density space
  !> @param[in]     exner                  Exner pressure in density space
  !> @param[in]     theta_in_wth           Theta in density space
  !> @param[in]     height_w3              Height of theta space levels above surface
  !> @param[in]     height_wtheta          Height of density space levels above surface
  !> @param[in]     totalppn_2d            Total precipitation from twod fields
  !> @param[in]     latitude               Latitude field
  !> @param[in,out] tau_east_spectral_gwd  Eastward flux
  !> @param[in,out] tau_south_spectral_gwd Southward flux
  !> @param[in,out] tau_west_spectral_gwd  Westward flux
  !> @param[in,out] tau_north_spectral_gwd Northward flux
  !> @param[in]     ndf_w3                 Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3                Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_wth                Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth               Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth                Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d                 Number of degrees of freedom per cell for 2D fields
  !> @param[in]     undf_2d                Number unique of degrees of freedom  for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields
  !>
  subroutine spectral_gwd_code(nlayers, seg_len,                             &
                               du_spectral_gwd, dv_spectral_gwd,             &
                               dtemp_spectral_gwd, u_in_w3, v_in_w3,         &
                               wetrho_in_w3, exner_in_wth, theta_in_wth,     &
                               height_w3, height_wth, totalppn_2d, latitude, &
                               tau_east_spectral_gwd, tau_south_spectral_gwd,&
                               tau_west_spectral_gwd, tau_north_spectral_gwd,&
                               ndf_w3, undf_w3, map_w3, ndf_wth, undf_wth,   &
                               map_wth, ndf_2d, undf_2d, map_2d)

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use planet_constants_mod,       only: p_zero, kappa, planet_radius
    use timestep_mod,               only: timestep
    use gw_ussp_mod,                only: gw_ussp

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_w3, ndf_wth, ndf_2d
    integer(kind=i_def), intent(in) :: undf_w3, undf_wth, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_w3, seg_len)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_wth, seg_len) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d, seg_len)  :: map_2d

    real(kind=r_def), intent(inout), dimension(undf_w3)  :: du_spectral_gwd
    real(kind=r_def), intent(inout), dimension(undf_w3)  :: dv_spectral_gwd
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtemp_spectral_gwd
    real(kind=r_def), intent(in), dimension(undf_w3)   :: u_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3)   :: v_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3)   :: wetrho_in_w3
    real(kind=r_def), intent(in), dimension(undf_wth)  :: exner_in_wth
    real(kind=r_def), intent(in), dimension(undf_wth)  :: theta_in_wth
    real(kind=r_def), intent(in), dimension(undf_w3)   :: height_w3
    real(kind=r_def), intent(in), dimension(undf_wth)  :: height_wth
    real(kind=r_def), intent(in), dimension(undf_2d)   :: totalppn_2d
    real(kind=r_def), intent(in), dimension(undf_2d)   :: latitude

    ! Diagnostics
    real(kind=r_def), pointer, intent(inout) :: tau_east_spectral_gwd(:)
    real(kind=r_def), pointer, intent(inout) :: tau_south_spectral_gwd(:)
    real(kind=r_def), pointer, intent(inout) :: tau_west_spectral_gwd(:)
    real(kind=r_def), pointer, intent(inout) :: tau_north_spectral_gwd(:)

    ! Local variables for the kernel
    ! Pressure at theta levels (but 0 at top)
    real(r_um), dimension(seg_len,1,0:nlayers) :: p_theta_levels, &
                                                  r_theta_levels

    real(r_um), dimension(seg_len,1) :: totalppn, & ! total percipitation
                                        sin_theta_latitude

    ! Local variables for the kernel
    real(r_um), dimension(seg_len,1,nlayers) ::   u_on_p, v_on_p, & ! u and v at p points
                                                  theta,          &
                                                  dtemp_on_t,     & ! temperature increment
                                                  wetrho_rsq,     & ! wetrho*r*r
                                                  r_rho_levels,   &
                                                  du_on_p, dv_on_p  ! wind increments
    ! Local variables that are diagnostics
    real(r_um), dimension(:,:,:), allocatable :: gwspec_eflux, &
                                                 gwspec_sflux, &
                                                 gwspec_wflux, &
                                                 gwspec_nflux
    ! Unused diagnostics
    real(r_um), dimension(1,1,1) :: gwspec_ewacc, gwspec_nsacc

    integer(i_um) :: k, i

    integer(i_um), parameter :: global_row_length = 1

    ! These are flags for diagnostics that are used in LFRic
    logical ::                   &
               gwspec_eflux_on,  &
               gwspec_sflux_on,  &
               gwspec_wflux_on,  &
               gwspec_nflux_on

    ! These are flags for diagnostics that are not used in LFRic
    logical, parameter ::                   &
               gwspec_eflux_p_on = .false., &
               gwspec_wflux_p_on = .false., &
               gwspec_ewacc_on = .false.,   &
               gwspec_ewacc_p_on = .false., &
               gwspec_nsacc_on = .false.

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------
    ! Initialise some variables (and diagnostics) to zero

    do k = 1, nlayers
      do i = 1, seg_len
        dtemp_on_t(i,1,k) = 0.0_r_um
        du_on_p(i,1,k) = 0.0_r_um
        dv_on_p(i,1,k) = 0.0_r_um
      end do
    end do

    ! This assumes that map_wth(1) points to level 0
    ! and map_w3(1) points to level 1

    do i = 1, seg_len
      do k = 0, nlayers-1
        ! Pressure on layer boundaries (note, top layer is set to zero below)
        p_theta_levels(i,1,k) = real(p_zero*(exner_in_wth(map_wth(1,i) + k)) &
                                      **(1.0_r_um/kappa), r_um)
      end do   ! k
      p_theta_levels(i,1,nlayers) = 0.0_r_um
    end do ! i

    do i = 1, seg_len
      r_theta_levels(i,1,0) = real(height_wth(map_wth(1,i) + 0) &
                                 + planet_radius, r_um)
      do k = 1, nlayers

        r_theta_levels(i,1,k) = real(height_wth(map_wth(1,i) + k) + planet_radius, r_um)

        u_on_p(i,1,k) = real(u_in_w3(map_w3(1,i) + k-1), r_um)
        v_on_p(i,1,k) = real(v_in_w3(map_w3(1,i) + k-1), r_um)

        theta(i,1,k) = real(theta_in_wth(map_wth(1,i) + k), r_um)

        r_rho_levels(i,1,k)   = real(height_w3(map_w3(1,i) + k-1) + planet_radius, r_um)
        wetrho_rsq(i,1,k) = real(wetrho_in_w3(map_w3(1,i) + k-1) * &
                               (r_rho_levels(i,1,k)**2), r_um)

      end do   ! k
    end do ! i

    do i = 1, seg_len
      totalppn(i,1) = real(totalppn_2d(map_2d(1,i)), r_um)
      sin_theta_latitude(i,1) = real(sin( latitude(map_2d(1,i)) ) ,r_um)
    end do

    ! Set stash flags and arrays
    if (.not. associated(tau_east_spectral_gwd, empty_real_data) ) then
      gwspec_eflux_on = .true.
      allocate(gwspec_eflux(seg_len,1,nlayers))
    else
      gwspec_eflux_on = .false.
      allocate(gwspec_eflux(1,1,1))
    end if
    if (.not. associated(tau_south_spectral_gwd, empty_real_data) ) then
      gwspec_sflux_on = .true.
      allocate(gwspec_sflux(seg_len,1,nlayers))
    else
      gwspec_sflux_on = .false.
      allocate(gwspec_sflux(1,1,1))
    end if
    if (.not. associated(tau_west_spectral_gwd, empty_real_data) ) then
      gwspec_wflux_on = .true.
      allocate(gwspec_wflux(seg_len,1,nlayers))
    else
      gwspec_wflux_on = .false.
      allocate(gwspec_wflux(1,1,1))
    end if
    if (.not. associated(tau_north_spectral_gwd, empty_real_data) ) then
      gwspec_nflux_on = .true.
      allocate(gwspec_nflux(seg_len,1,nlayers))
    else
      gwspec_nflux_on = .false.
      allocate(gwspec_nflux(1,1,1))
    end if

    ! call USSP code from UM
    call gw_ussp(nlayers, 1, seg_len,                                        &
                 global_row_length,                                          &
                 r_rho_levels, r_theta_levels, p_theta_levels,               &
                 sin_theta_latitude, theta, wetrho_rsq, u_on_p, v_on_p,      &
                 totalppn, timestep, du_on_p, dv_on_p, dtemp_on_t,           &
                 ussp_heating, gwspec_eflux,gwspec_sflux,gwspec_wflux,       &
                 gwspec_nflux, gwspec_ewacc,gwspec_nsacc,                    &
                 gwspec_eflux_on, gwspec_eflux_p_on, gwspec_sflux_on,        &
                 gwspec_wflux_on, gwspec_wflux_p_on, gwspec_nflux_on,        &
                 gwspec_ewacc_on, gwspec_ewacc_p_on, gwspec_nsacc_on)

    do k = 1, nlayers
      do i = 1, seg_len
        du_spectral_gwd(map_w3(1,i) + k-1) = real(du_on_p(i,1,k), r_def)
        dv_spectral_gwd(map_w3(1,i) + k-1) = real(dv_on_p(i,1,k), r_def)
        dtemp_spectral_gwd(map_wth(1,i) + k) = real(dtemp_on_t(i,1,k), r_def)
      end do ! i
    end do   ! k

    ! Set level 0 increment such that theta increment will equal level 1
    do i = 1, seg_len
      dtemp_spectral_gwd(map_wth(1,i) + 0) = real(dtemp_on_t(i,1,1), r_def) &
                                         * exner_in_wth(map_wth(1,i) + 0)   &
                                         / exner_in_wth(map_wth(1,i) + 1)
    end do

    ! Map diagnostics back
    if (.not. associated(tau_east_spectral_gwd, empty_real_data) ) then
      do k = 1, nlayers
        do i = 1, seg_len
          tau_east_spectral_gwd(map_wth(1,i) + k) = real(gwspec_eflux(i,1,k), r_def)
        end do
      end do
      do i = 1, seg_len
        tau_east_spectral_gwd(map_wth(1,i) + 0) = real(gwspec_eflux(i,1,1), r_def)
      end do
    end if
    if (.not. associated(tau_south_spectral_gwd, empty_real_data) ) then
      do k = 1, nlayers
        do i = 1, seg_len
          tau_south_spectral_gwd(map_wth(1,i) + k) = real(gwspec_sflux(i,1,k), r_def)
        end do
      end do
      do i = 1, seg_len
        tau_south_spectral_gwd(map_wth(1,i) + 0) = real(gwspec_sflux(i,1,1), r_def)
      end do
    end if
    if (.not. associated(tau_west_spectral_gwd, empty_real_data) ) then
      do k = 1, nlayers
        do i = 1, seg_len
          tau_west_spectral_gwd(map_wth(1,i) + k) = real(gwspec_wflux(i,1,k), r_def)
        end do
      end do
      do i = 1, seg_len
        tau_west_spectral_gwd(map_wth(1,i) + 0) = real(gwspec_wflux(i,1,1), r_def)
      end do
    end if
    if (.not. associated(tau_north_spectral_gwd, empty_real_data) ) then
      do k = 1, nlayers
        do i = 1, seg_len
          tau_north_spectral_gwd(map_wth(1,i) + k) = real(gwspec_nflux(i,1,k), r_def)
        end do
      end do
      do i = 1, seg_len
        tau_north_spectral_gwd(map_wth(1,i) + 0) = real(gwspec_nflux(i,1,1), r_def)
      end do
    end if

    deallocate(gwspec_eflux)
    deallocate(gwspec_sflux)
    deallocate(gwspec_wflux)
    deallocate(gwspec_nflux)

  end subroutine spectral_gwd_code


end module spectral_gwd_kernel_mod
