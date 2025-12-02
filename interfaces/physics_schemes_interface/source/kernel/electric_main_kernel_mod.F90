!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the microphysics-based electric (lightning) scheme.

module electric_main_kernel_mod

use argument_mod, only: arg_type,                                              &
                        GH_FIELD, GH_SCALAR, GH_REAL,                          &
                        GH_READ, GH_READWRITE,                                 &
                        ANY_DISCONTINUOUS_SPACE_1,                             &
                        CELL_COLUMN

use empty_data_mod,    only: empty_real_data
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: electric_main_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! mi_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! ms_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! mg_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! temp_in_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! w_in_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ, WTHETA),       & ! rhodz_in_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ),               & ! dt
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: electric_main_code
end type

public :: electric_main_code

contains

!> @brief Interface to the main UM electric (or lightning) scheme
!> @param[in]     nlayers               Number of layers
!> @param[in]     mi_wth                Ice cloud mass mixing ratio
!> @param[in]     ms_wth                Snow mass mixing ratio
!> @param[in]     mg_wth                Graupel mass mixing ratio
!> @param[in]     temp_in_wth           Temperature field (K)
!> @param[in]     w_in_wth              'Vertical' wind in theta space
!> @param[in]     rhodz_in_wth          'Air density * layer thickness
!> @param[in]     dt                    The model timestep length (s)
!> @param[in,out] flash_potential       Flash potential, accumulated over
!!                                       model timesteps
!> @param[inout]  num_flashes_2d        Number of lightning flashes discharged
!!                                       on this timestep due to total lightning
!!                                       activity
!> @param[inout]  total_flash_rate_2d   Total flash rate (s-1), for all
!!                                       lightning types
!> @param[inout]  storm_field_2d        Flag to indicate whether a grid column
!!                                       contains a thunderstorm
!> @param[inout]  fr1_mc_2d             Lightning flash rate (s-1) due to the
!!                                       upward flux of graupel at the -15
!!                                       Celsius level in the atmosphere.
!> @param[inout]  fr2_mc_2d             Lightning flash rate (s-1) due to the
!!                                       total ice water path in the model
!!                                       column
!> @param[inout]  gwp_2d                Graupel water path in the model column
!!                                       (kg m-2)
!> @param[inout]  tiwp_2d               Total ice water path in the model
!!                                       column (kg m-2)
!> @param[in]     ndf_wth               Number of degrees of freedom per cell
!!                                       for potential temperature space
!> @param[in]     undf_wth              Number unique of degrees of freedom in
!!                                       segment for potential temperature space
!> @param[in]     map_wth               Dofmap for a segment of any potential
!!                                       temperature space field
!> @param[in]     ndf_2d                Number of degrees of freedom per cell for
!!                                       2D fields
!> @param[in]     undf_2d               Number unique of degrees of freedom in
!!                                       segment for 2D fields
!> @param[in]     map_2d                Dofmap for a segment of any 2D field

subroutine electric_main_code(nlayers, mi_wth, ms_wth, mg_wth, temp_in_wth,    &
                              w_in_wth, rhodz_in_wth,                          &
                              dt, flash_potential, num_flashes_2d,             &
                              total_flash_rate_2d, storm_field_2d,             &
                              fr1_mc_2d, fr2_mc_2d, gwp_2d, tiwp_2d,           &
                              ndf_wth, undf_wth, map_wth,                      &
                              ndf_2d, undf_2d, map_2d )

    use constants_mod, only: r_def, r_um, i_um, i_def

    !------------------
    ! UM modules
    !------------------
    use mphys_inputs_mod,     only: l_mcr_qcf2

    ! Modules from UM 'electric' directory for subroutine calls.
    use calc_wp_mod,      only: calc_wp
    use define_storm_mod, only: define_storm
    use flash_rate_mod,   only: flash_rate

    implicit none

    !-----------------------
    ! Subroutine arguments
    !-----------------------
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wth, undf_2d, ndf_wth, ndf_2d
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    real(kind=r_def), intent(in), dimension(undf_wth) :: mi_wth ! ice mmr
    real(kind=r_def), intent(in), dimension(undf_wth) :: ms_wth ! snow mmr
    real(kind=r_def), intent(in), dimension(undf_wth) :: mg_wth ! graupel mmr
    real(kind=r_def), intent(in), dimension(undf_wth) :: temp_in_wth
                                                         ! Temperature [K]
    real(kind=r_def), intent(in), dimension(undf_wth) :: w_in_wth
                                                         ! Vertical velocity
    real(kind=r_def), intent(in), dimension(undf_wth) :: rhodz_in_wth
    real(kind=r_def), intent(inout), dimension(undf_2d) :: flash_potential, &
                                                           num_flashes_2d

    ! Diagnostic arguments
    real(kind=r_def), pointer, intent(inout) :: total_flash_rate_2d(:)
    real(kind=r_def), pointer, intent(inout) :: storm_field_2d(:)
    real(kind=r_def), pointer, intent(inout) :: fr1_mc_2d(:)
    real(kind=r_def), pointer, intent(inout) :: fr2_mc_2d(:)
    real(kind=r_def), pointer, intent(inout) :: gwp_2d(:)
    real(kind=r_def), pointer, intent(inout) :: tiwp_2d(:)

    ! The model timestep length
    real(kind=r_def), intent(in) :: dt

    !-------------------------------------
    ! Local variables for the kernel
    !-------------------------------------
    logical, dimension(1,1) :: storm_field_work

    real(r_um), dimension(1,1,nlayers) :: rhodz_work
    real(r_um), dimension(1,1,nlayers) :: qcf_tot_work
    real(r_um), dimension(1,1,nlayers) :: qgraup_work
    real(r_um), dimension(1,1,nlayers) :: t_work
    real(r_um), dimension(1,1,0:nlayers) :: w_work

    real(r_um), dimension(1,1) :: fr1_mc_work, fr2_mc_work
    real(r_um), dimension(1,1) :: gwp_work, tiwp_work
    real(r_um), dimension(1,1) :: total_fr_work

    ! Minimum flash potential carried
    real(r_um), parameter :: flashmin = 1.0e-10_r_def
    real(r_um) :: flash_pot

    integer(i_um) :: k, nspts

    !-------------------------------------------------------------
    ! Convert inputs from the algorithm layer into UM-dimensions
    !-------------------------------------------------------------
    do k = 1, nlayers
      t_work(1,1,k) = temp_in_wth(map_wth(1) + k)
      w_work(1,1,k) = w_in_wth(map_wth(1) + k)
      qgraup_work(1,1,k) = mg_wth(map_wth(1) + k)
      rhodz_work(1,1,k) = rhodz_in_wth(map_wth(1) + k)
    end do

    !-------------------------------------------------------------
    ! Initialise key variables:
    ! In the UM this is done via a call to electric_init. However
    ! this largely sets variables to zero and contains a lot of
    ! specific quantity code not relevant to LFRic. So let's do
    ! the initialisation over just a few lines here.
    !-------------------------------------------------------------
    storm_field_work(1,1) = .false.
    total_fr_work(1,1) = 0.0_r_um
    fr1_mc_work(1,1) = 0.0_r_um
    fr2_mc_work(1,1) = 0.0_r_um

    ! Add up the total ice mixing ratio in the column
    if (l_mcr_qcf2) then
      ! Two ice prognostics in use
      do k = 1, nlayers
        qcf_tot_work(1,1,k) = ms_wth(map_wth(1) + k) + mi_wth(map_wth(1) + k)
      end do
    else
      ! One ice prognostic
      do k = 1, nlayers
        qcf_tot_work(1,1,k) = ms_wth(map_wth(1) + k)
      end do
    end if ! l_mcr_qcf2

    ! Calculate graupel water path (gwp_work)
    call calc_wp(qgraup_work, rhodz_work, gwp_work)

    ! Calculate total ice water path (tiwp_work)
    call calc_wp(qcf_tot_work, rhodz_work, tiwp_work)

    ! Locate where storms are present within the model fields
    call define_storm(gwp_work, tiwp_work, storm_field_work, nspts)

    ! Calculate the flash rate of storms
    if (nspts > 0) then
      call flash_rate( storm_field_work, qgraup_work, t_work, w_work,          &
                       rhodz_work, gwp_work, tiwp_work, total_fr_work,         &
                       fr1_mc_work, fr2_mc_work)
    end if

    !---------------------------------------------------------------------------
    ! Need to distribute the flash into a 2-D flash potential and flashes.
    ! In the UM this is done in the distribute_flash subroutine, but for LFRic
    ! the process can be made much simpler and is done in the loops below.
    !---------------------------------------------------------------------------
    flash_pot = flash_potential(map_2d(1)) + ( total_fr_work(1,1) * dt )

    if (storm_field_work(1,1)) then

      if (flash_pot >= 1.0_r_def) then
        ! Flash rate is sufficient to discharge lightning flashes
        num_flashes_2d(map_2d(1)) = real(floor(flash_pot))
        flash_pot = flash_pot - num_flashes_2d(map_2d(1))
      else
        ! Flash rate is zero
        num_flashes_2d(map_2d(1)) = 0.0_r_def
      end if

    else

      ! Flash rate is zero
      num_flashes_2d(map_2d(1)) = 0.0_r_def

      ! If this point is not defined as a storm, then by definition
      ! no charging can take place. Therefore, set the flash potential
      ! to zero. This prevents any unused flash from older storms which
      ! passed over this grid box many timesteps ago contributing to the
      ! present lightning flash.
      flash_pot = 0.0_r_def

    end if

    ! Update prognostic flash potential
    if (flash_pot >= flashmin) then
      flash_potential(map_2d(1)) = flash_pot
    else
      flash_potential(map_2d(1)) = 0.0_r_def
    end if

    !--------------------------------------------------
    ! Calculate diagnostics to pass back to the kernel
    !--------------------------------------------------
    ! Lightning flash rates
    if (.not. associated(total_flash_rate_2d, empty_real_data) ) then
      total_flash_rate_2d(map_2d(1)) = total_fr_work(1,1)
    end if

    if (.not. associated(fr1_mc_2d, empty_real_data) ) then
      fr1_mc_2d(map_2d(1)) = fr1_mc_work(1,1)
    end if

    if (.not. associated(fr2_mc_2d, empty_real_data) ) then
      fr2_mc_2d(map_2d(1)) = fr2_mc_work(1,1)
    end if

    if (.not. associated(storm_field_2d, empty_real_data) ) then
      ! Convert storm field from logical to real type
      if (storm_field_work(1,1)) then
        storm_field_2d(map_2d(1)) = 1.0_r_def
      else
        storm_field_2d(map_2d(1)) = 0.0_r_def
      end if
    end if

    ! Water paths for diagnostics
    if (.not. associated(gwp_2d, empty_real_data) ) then
      gwp_2d(map_2d(1))  = gwp_work(1,1)
    end if

    if (.not. associated(tiwp_2d, empty_real_data) ) then
      tiwp_2d(map_2d(1)) = tiwp_work(1,1)
    end if

end subroutine electric_main_code
end module electric_main_kernel_mod




