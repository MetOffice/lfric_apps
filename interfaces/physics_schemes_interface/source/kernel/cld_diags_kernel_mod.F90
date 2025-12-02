!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud diagnostics.
!>
module cld_diags_kernel_mod

  use argument_mod,       only : arg_type,                                 &
                                 GH_FIELD, GH_REAL,                        &
                                 GH_READ, GH_READWRITE,                    &
                                 GH_WRITE, CELL_COLUMN,                    &
                                 GH_SCALAR, GH_LOGICAL,                    &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um, l_def
  use empty_data_mod,     only : empty_real_data
  use fs_continuity_mod,  only : Wtheta, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: cld_diags_kernel_type
    private
    type(arg_type) :: meta_args(20) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                         & ! combined_cld_amount_wth
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! cld_amount_max
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! cld_amount_rnd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! cld_amount_maxrnd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! ceil_cld_amount_maxrnd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! cld_base_altitude
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! low_cld_base_altitude
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! very_low_cld_amount
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! low_cld_amount
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! medium_cld_amount
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! high_cld_amount
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! very_high_cld_amount
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! cloud_fraction_below_1000feet_asl
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                         & ! mi_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                         & ! cf_frozen
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                         & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                         & ! rho_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                             & ! height_w3
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                                & ! opt_depth_thresh
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              & ! filter_optical_depth
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cld_diags_code
  end type

  public :: cld_diags_code

contains

  !> @brief Interface to the cloud diagnostics.
  !> @details Calculation of various cloud diagnostics.
  !>
  !> @param[in]     nlayers                    Number of layers
  !> @param[in]     combined_cld_amount_wth    Combined cloud amount
  !> @param[in,out] cld_amount_max             Cloud amount maximum overlap
  !> @param[in,out] cld_amount_rnd             Cloud amount random overlap
  !> @param[in,out] cld_amount_maxrnd          Cloud amount maximum-random overlap
  !> @param[in,out] ceil_cld_amount_maxrnd     Ceilometer filtered cloud amount maximum-random overlap
  !> @param[in,out] cld_base_altitude          Cloud base altitude wrt sea level
  !> @param[in,out] low_cld_base_altitude      Cloud base altitude wrt sea level for very low amount of cloud
  !> @param[in,out] very_low_cld_amount        Maximum cloud amount below 111m
  !> @param[in,out] low_cld_amount             Maximum cloud amount between 111 and 1949m above sea level
  !> @param[in,out] medium_cld_amount          Maximum cloud amount between 1949 and 5574m above sea level
  !> @param[in,out] high_cld_amount            Maximum cloud amount between 5574 and 13608m above sea level
  !> @param[in,out] very_high_cld_amount       Maximum cloud amount above 13608m above sea level
  !> @param[in,out] cloud_fraction_below_1000feet_asl Cloud fraction below 1000 feet above sea level
  !> @param[in]     mi_wth                     Ice water content
  !> @param[in]     cf_frozen                  Frozen fraction
  !> @param[in]     height_wth                 Height above sea level in wtheta
  !> @param[in]     rho_wth                    Dry air density in wtheta
  !> @param[in]     height_w3                  Height above sea level in w3
  !> @param[in]     opt_depth_thresh           Optical depth filter threshold
  !> @param[in]     filter_optical_depth       Apply optical depth filter for cloud diagnostics
  !> @param[in]     ndf_2d                     Number of degrees of freedom per cell for 2D fields
  !> @param[in]     undf_2d                    Number unique of degrees of freedom  for 2D fields
  !> @param[in]     map_2d                     Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_wth                    Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth                   Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth                    Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3                     Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3                    Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3                     Dofmap for the cell at the base of the column for density space


 ! Note that combined_cld_amount_wth needs to be the first argument below, since
 ! this field is used to determine the value of nlayers in psyclone and hence
 ! needs to be a 3D field.
  subroutine cld_diags_code( nlayers,                 &
                             combined_cld_amount_wth, &
                             cld_amount_max,          &
                             cld_amount_rnd,          &
                             cld_amount_maxrnd,       &
                             ceil_cld_amount_maxrnd,  &
                             cld_base_altitude,       &
                             low_cld_base_altitude,   &
                             very_low_cld_amount,     &
                             low_cld_amount,          &
                             medium_cld_amount,       &
                             high_cld_amount,         &
                             very_high_cld_amount,    &
                             cloud_fraction_below_1000feet_asl, &
                             mi_wth,                  &
                             cf_frozen,               &
                             height_wth,              &
                             rho_wth,                 &
                             height_w3,               &
                             opt_depth_thresh,        &
                             filter_optical_depth,    &
                             ndf_wth,                 &
                             undf_wth,                &
                             map_wth,                 &
                             ndf_2d,                  &
                             undf_2d,                 &
                             map_2d,                  &
                             ndf_w3,                  &
                             undf_w3,                 &
                             map_w3                   )

    use science_conversions_mod, only: feet_to_metres

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in)     :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in)     :: ndf_w3, undf_w3

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3

    real(kind=r_def), pointer, intent(inout)            :: cld_amount_max(:)
    real(kind=r_def), pointer, intent(inout)            :: cld_amount_rnd(:)
    real(kind=r_def), pointer, intent(inout)            :: cld_amount_maxrnd(:)
    real(kind=r_def), pointer, intent(inout)            :: ceil_cld_amount_maxrnd(:)
    real(kind=r_def), pointer, intent(inout)            :: cld_base_altitude(:)
    real(kind=r_def), pointer, intent(inout)            :: low_cld_base_altitude(:)
    real(kind=r_def), pointer, intent(inout)            :: very_low_cld_amount(:)
    real(kind=r_def), pointer, intent(inout)            :: low_cld_amount(:)
    real(kind=r_def), pointer, intent(inout)            :: medium_cld_amount(:)
    real(kind=r_def), pointer, intent(inout)            :: high_cld_amount(:)
    real(kind=r_def), pointer, intent(inout)            :: very_high_cld_amount(:)
    real(kind=r_def), pointer, intent(inout)            :: cloud_fraction_below_1000feet_asl(:)
    real(kind=r_def), intent(in), dimension(undf_wth)   :: combined_cld_amount_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: mi_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: cf_frozen
    real(kind=r_def), intent(in), dimension(undf_wth)   :: height_wth
    real(kind=r_def), intent(in), dimension(undf_wth)   :: rho_wth
    real(kind=r_def), intent(in), dimension(undf_w3)    :: height_w3

    real(kind=r_def), intent(in)                        :: opt_depth_thresh
    logical(kind=l_def), intent(in)                     :: filter_optical_depth

    integer(kind=i_def) :: k, k_top, k_bot

    real(kind=r_def), dimension(nlayers) :: area_cld_fraction,            &
         combined_cld_amount, ceil_combined_cld_amount, z_agl_centre_of_levels,&
         z_asl_base_of_levels, z_asl_centre_of_levels

    ! Certain cloud amount diagnostics are values between certain height ranges.
    ! Height of the transition between two ranges are defined as
    ! e.g. low_to_medium which is transition from low cloud amount to medium.
    ! Height are defined based on temporally-fixed ICAO (International Civil
    ! Aviation Organization) standard atmosphere pressure levels for:
    ! ... 1000 hPa:
    real(kind=r_def), parameter :: very_low_to_low   =   111.0_r_def ! metres
    ! ...  800 hPa:
    real(kind=r_def), parameter :: low_to_medium     =  1949.0_r_def  ! metres
    ! ...  500 hPa:
    real(kind=r_def), parameter :: medium_to_high    =  5574.0_r_def  ! metres
    ! ...  150 hPa:
    real(kind=r_def), parameter :: high_to_very_high = 13608.0_r_def  ! metres
    ! When looking for cloud base: how much cloud cover defines cloud base.
    real(kind=r_def), parameter :: cld_cover_for_cld_base = 2.5_r_def/8.0_r_def
    ! When looking for cloud base sometimes want much smaller amount of cloud.
    real(kind=r_def), parameter :: low_cld_cover_for_cld_base = 0.05_r_def
    ! For tracking whether cloud has been found.
    logical(kind=l_def) :: found
    ! Cloud-base height diagnostic is often used by the aviation community,
    ! so the units need to be converted from metres to kilofeet
    real(kind=r_def), parameter :: m_to_kfeet = 0.001_r_def / feet_to_metres
    ! Max range of ceilometers used to ignore any cloud beyond range.
    real(kind=r_def), parameter :: ceilometer_range = 6.0e3_r_def
    ! If no model levels are found in a certain height range category, set the
    ! cloud fraction to missing data but don't use rmdi=-huge as that will
    ! make it hard to do quickviews as the dynamic range of the data will get
    ! compressed. Instead use unphysical value of the same order of magnitude.
    real(kind=r_def), parameter :: cld_mdi = -0.999_r_def

    ! Small number: Tolerance of cloud fraction for resetting of cloud
    !               or minimum total thickness of layers to divide by.
    real(kind=r_def), parameter :: cld_tol = 0.001_r_def

    ! Scalar field for working out calculations
    real(kind=r_def) :: work_scalar, work_counter, dz, top

    ! Scalars for optical depth filtering
    real(kind=r_def) :: sigma, tau, filtered_cloud, ice_fraction

    do k = 1, nlayers
      combined_cld_amount(k) = combined_cld_amount_wth(map_wth(1) + k)
    end do

    if (filter_optical_depth) then
      ! Filter the combined cloud amount by removing sub-visual
      ! (optically thin) ice cloud before calculating other diagnostics
      ! such as total cloud amount.
      do k = 1, nlayers
        if (combined_cld_amount(k) > cld_tol) then
          ! Extinction (m-1) taken from Heymsfield et al 2003.
          sigma = ((( mi_wth(map_wth(1) + k) *                                 &
                      rho_wth(map_wth(1) + k) * 1000.0_r_def ) /               &
                      combined_cld_amount(k))**0.9_r_def ) * 0.02_r_def
          ! Optical depth=integral of sigma over depth of model level.
          tau = sigma * (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + k-1))
          if (tau < opt_depth_thresh) then
            filtered_cloud = 0.0_r_def
          else
            filtered_cloud = combined_cld_amount(k)
          end if

          ice_fraction = cf_frozen(map_wth(1) + k) / combined_cld_amount(k)
          combined_cld_amount(k) = ( ice_fraction * filtered_cloud ) +         &
             ((1.0_r_def - ice_fraction) * combined_cld_amount(k))

        else
          ! If there is practically no cloud
          ! set filtered combined cloud to zero.
          combined_cld_amount(k) = 0.0_r_def
        end if
      end do
    end if

    ! cld_amount_max
    if (.not. associated(cld_amount_max, empty_real_data) ) then
      cld_amount_max(map_2d(1)) = maxval(combined_cld_amount(:))
    end if

    ! cld_amount_rnd
    if (.not. associated(cld_amount_rnd, empty_real_data) ) then
      work_scalar = 1.0_r_def
      do k = 1, nlayers
        work_scalar = work_scalar * ( 1.0_r_def - combined_cld_amount(k) )
      end do
      cld_amount_rnd(map_2d(1)) = 1.0_r_def - work_scalar
    end if

    ! cld_amount_maxrnd
    if (.not. associated(cld_amount_maxrnd, empty_real_data) ) then
      ! from Raisanen (1998), MWR, vol 126, eqn 2
      work_scalar = 1.0_r_def - combined_cld_amount(nlayers)
      do k = 1, nlayers-1
        if ( combined_cld_amount(k+1) < 1.0_r_def  ) then
          work_scalar = work_scalar * &
            ( 1.0_r_def - max(combined_cld_amount(k+1),combined_cld_amount(k)) )&
            / ( 1.0_r_def - combined_cld_amount(k+1) )
        else
          ! Overcast, overlapped cloud cover must be too. But need to avoid
          ! divide by zero. Set this to 0.0, so next line gives 1.0
          work_scalar = 0.0_r_def
        end if
      end do
      cld_amount_maxrnd(map_2d(1)) = 1.0_r_def - work_scalar
    end if

    ! ceil_cld_amount_maxrnd
    if (.not. associated(ceil_cld_amount_maxrnd, empty_real_data) ) then
      do k = 1, nlayers
        ! NB: We are considering the range of a ceilometer, so the height
        ! needs to be above ground level (agl) rather than above sea level (asl).
        z_agl_centre_of_levels(k)=height_wth(map_wth(1)+k)-height_wth(map_wth(1)+0)
        if ( z_agl_centre_of_levels(k) <= ceilometer_range ) then
          ceil_combined_cld_amount(k) = combined_cld_amount(k)
        else
          ! Simple filtering where any model cloud above a certain height is
          ! assumed not detectable by ceilometer and hence removed for fairer comparison.
          ceil_combined_cld_amount(k) = 0.0_r_def
        end if
      end do
      ! from Raisanen (1998), MWR, vol 126, eqn 2
      work_scalar = 1.0_r_def - ceil_combined_cld_amount(nlayers)
      do k = 1, nlayers-1
        if ( ceil_combined_cld_amount(k+1) < 1.0_r_def  ) then
          work_scalar = work_scalar * &
            ( 1.0_r_def - max( ceil_combined_cld_amount(k+1), ceil_combined_cld_amount(k) )) &
            / ( 1.0_r_def - ceil_combined_cld_amount(k+1) )
        else
          ! Overcast, overlapped cloud cover must be too. But need to avoid
          ! divide by zero. Set this to 0.0, so next line gives 1.0
          work_scalar = 0.0_r_def
        end if
      end do
      ceil_cld_amount_maxrnd(map_2d(1)) = 1.0_r_def - work_scalar
    end if

    ! Find heights above sea level (asl) if required.
    if (.not. associated(cld_base_altitude, empty_real_data)     .or. &
        .not. associated(low_cld_base_altitude, empty_real_data) .or. &
        .not. associated(very_low_cld_amount, empty_real_data)   .or. &
        .not. associated(low_cld_amount, empty_real_data)        .or. &
        .not. associated(medium_cld_amount, empty_real_data)     .or. &
        .not. associated(high_cld_amount, empty_real_data)       .or. &
        .not. associated(very_high_cld_amount, empty_real_data)  .or. &
        .not. associated(cloud_fraction_below_1000feet_asl, empty_real_data) ) then
      z_asl_base_of_levels(1)   = height_wth(map_wth(1) + 0 )
      z_asl_centre_of_levels(1) = height_wth(map_wth(1) + 1 )

      do k = 2, nlayers
        ! NB: these height are above sea level (asl) not above ground level.
        z_asl_base_of_levels(k)   = height_w3(map_w3(1)   + k-1 )
        z_asl_centre_of_levels(k) = height_wth(map_wth(1) + k   )
      end do
    end if

    ! cld_base_altitude (in kilofeet)
    if (.not. associated(cld_base_altitude, empty_real_data) ) then
      ! As a default, set cloud-base to beyond top of model
      cld_base_altitude(map_2d(1)) = 1.1_r_def *                               &
                                     z_asl_centre_of_levels(nlayers) *         &
                                     m_to_kfeet

      do k = 1, nlayers
        if ( combined_cld_amount(k) >= cld_cover_for_cld_base ) then
          cld_base_altitude(map_2d(1)) = z_asl_base_of_levels(k) * m_to_kfeet
          exit
        end if
      end do
    end if

    ! low_cld_base_altitude (in feet) for low threshold cloud amounts
    if (.not. associated(low_cld_base_altitude, empty_real_data) ) then
      ! As a default, set cloud-base to beyond top of model
      low_cld_base_altitude(map_2d(1)) = 1.1_r_def * &
                            z_asl_centre_of_levels(nlayers) / feet_to_metres
      do k = 1, nlayers
        if ( combined_cld_amount(k) >= low_cld_cover_for_cld_base ) then
          low_cld_base_altitude(map_2d(1)) = z_asl_base_of_levels(k) / feet_to_metres
          exit
        end if
      end do
    end if

    ! very_low_cld_amount
    if (.not. associated(very_low_cld_amount, empty_real_data) ) then
      found = .false.
      do k = nlayers, 1, -1
        if ( z_asl_centre_of_levels(k) < very_low_to_low .and. &
            .not. found ) then
          k_top = k
          found = .true.
        end if
      end do
      if ( found ) then
        very_low_cld_amount(map_2d(1)) = maxval(combined_cld_amount(1:k_top))
      else
        ! Orography at this location is high enough that no levels are in
        ! desired height range. So set to missing data but don't use -huge
        ! as that will make quick view of the data hard due to
        ! compressing of dynamic range
        very_low_cld_amount(map_2d(1)) = cld_mdi
      end if
    end if

    ! low_cld_amount
    if (.not. associated(low_cld_amount, empty_real_data) ) then
      k_bot = 1
      found = .false.
      do k = 1, nlayers
        if ( z_asl_centre_of_levels(k) >= very_low_to_low .and. &
            .not. found ) then
          k_bot = k
          found = .true.
        end if
      end do
      found = .false.
      do k = nlayers, 1, -1
        if ( z_asl_centre_of_levels(k) < low_to_medium .and. &
            .not. found ) then
          k_top = k
          found = .true.
        end if
      end do
      if ( found ) then
        ! A layer has been found whose top is within required height range.
        ! However, note k_bot could be 1 if orog is above 111m above sea level.
        low_cld_amount(map_2d(1)) = maxval(combined_cld_amount(k_bot:k_top))
      else
        ! Orography at this location is high enough that no levels are in
        ! desired height range. So set to missing data
        low_cld_amount(map_2d(1)) = cld_mdi
      end if
    end if

    ! medium_cld_amount
    if (.not. associated(medium_cld_amount, empty_real_data) ) then
      k_bot = 1
      found = .false.
      do k = 1, nlayers
        if ( z_asl_centre_of_levels(k) >= low_to_medium .and. &
            .not. found ) then
          k_bot = k
          found = .true.
        end if
      end do
      found = .false.
      do k = nlayers, 1, -1
        if ( z_asl_centre_of_levels(k) < medium_to_high .and. &
            .not. found ) then
          k_top = k
          found = .true.
        end if
      end do
      if ( found ) then
        ! A layer has been found whose top is within required height range.
        ! However, note k_bot could be 1 if orog is above 1949m above sea level.
        medium_cld_amount(map_2d(1)) = maxval(combined_cld_amount(k_bot:k_top))
      else
        ! Orography at this location is high enough that no levels are in
        ! desired height range. So set to missing data
        medium_cld_amount(map_2d(1)) = cld_mdi
      end if
    end if

    ! high_cld_amount
    if (.not. associated(high_cld_amount, empty_real_data) ) then
      k_bot = 1
      found = .false.
      do k = 1, nlayers
        if ( z_asl_centre_of_levels(k) >= medium_to_high .and. &
            .not. found ) then
          k_bot = k
          found = .true.
        end if
      end do
      found = .false.
      do k = nlayers, 1, -1
        if ( z_asl_centre_of_levels(k) < high_to_very_high .and. &
            .not. found ) then
          k_top = k
          found = .true.
        end if
      end do
      if ( found ) then
        ! A layer has been found whose top is within required height range.
        ! However, note k_bot could be 1 if orog is above 5574m above sea level.
        high_cld_amount(map_2d(1)) = maxval(combined_cld_amount(k_bot:k_top))
      else
        ! Orography at this location is high enough that no levels are in
        ! desired height range. So set to missing data
        high_cld_amount(map_2d(1)) = cld_mdi
      end if
    end if

    ! very_high_cld_amount
    if (.not. associated(very_high_cld_amount, empty_real_data) ) then
      k_bot = 1
      found = .false.
      do k = 1, nlayers
        if ( z_asl_centre_of_levels(k) >= high_to_very_high .and. &
            .not. found ) then
          k_bot = k
          found = .true.
        end if
      end do
      very_high_cld_amount(map_2d(1)) = maxval(combined_cld_amount(k_bot:nlayers))
    end if

    !  cloud_fraction_below_1000feet_asl
    if (.not. associated(cloud_fraction_below_1000feet_asl, empty_real_data) ) then

      if ( z_asl_base_of_levels(1) * m_to_kfeet <= 1.0_r_def ) then
        ! Bottom of lowest layer is below 1000 feet.
        work_scalar  = 0.0_r_def
        work_counter = 0.0_r_def

        do k = 1, nlayers-1
          ! Find top of this layer, or the height threshold, whichever is lower.
          top = min(z_asl_base_of_levels(k+1) * m_to_kfeet, 1.0_r_def)

          ! Hence find thickness of the region being added to running sum.
          dz = top - (z_asl_base_of_levels(k) * m_to_kfeet)
          ! If layer is entirely above height of interest, dz will be negative.
          if (dz <= 0.0_r_def) exit

          ! Calculate a running sum of the cloud fraction and the layers depths.
          work_scalar  = work_scalar  + ( dz * combined_cld_amount(k) )
          work_counter = work_counter +   dz

        end do

        ! Find a mean cloud fraction up to that height.
        cloud_fraction_below_1000feet_asl(map_2d(1)) = work_scalar / work_counter

      else
        ! Orography at this location is high enough that no levels are below
        ! 1000 feet ASL. So set to missing data
        cloud_fraction_below_1000feet_asl(map_2d(1)) = cld_mdi
      end if

    end if

  end subroutine cld_diags_code

end module cld_diags_kernel_mod
