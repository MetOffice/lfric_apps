!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Source time axis dimensions from ancil files
!>
module time_dimensions_mod

  use constants_mod,             only: i_def, l_def, str_def, cmdi
  use log_mod,                   only: log_event,                             &
                                       log_scratch_space,                     &
                                       log_level_error
#ifdef UM_PHYSICS
  ! This import split to support fparser which gets confused by FPP directives
  ! in the middle of a syntactic unit.
  !
  use initialization_config_mod, only: ancil_option,                          &
                                       ancil_option_fixed,                    &
                                       ancil_option_updating
#endif
  use initialization_config_mod, only: lbc_option,                            &
                                       lbc_option_gungho_file,                &
                                       lbc_option_um2lfric_file
#ifdef UM_PHYSICS
  ! This import split to support fparser which gets confused by FPP directives
  ! in the middle of a syntactic unit.
  !
  use files_config_mod,          only: ancil_dir => ancil_directory,          &
                                       sst_ancil_path,                        &
                                       emiss_bc_biofuel_ancil_path,           &
                                       emiss_bc_fossil_ancil_path,            &
                                       emiss_bc_biomass_ancil_path,           &
                                       emiss_bc_biomass_hi_ancil_path,        &
                                       emiss_bc_biomass_lo_ancil_path,        &
                                       emiss_om_biofuel_ancil_path,           &
                                       emiss_om_fossil_ancil_path,            &
                                       emiss_om_biomass_ancil_path,           &
                                       emiss_om_biomass_hi_ancil_path,        &
                                       emiss_om_biomass_lo_ancil_path,        &
                                       emiss_so2_low_ancil_path,              &
                                       emiss_so2_high_ancil_path
#endif
  use files_config_mod,          only: lbc_dir => lbc_directory,              &
                                       lbc_filename
#ifdef UM_PHYSICS
  use aerosol_config_mod,        only: glomap_mode,                           &
                                       glomap_mode_ukca
  use chemistry_config_mod,      only: chem_scheme, chem_scheme_strattrop,    &
                                       chem_scheme_strat_test
#endif
  implicit none

  private
  public :: sync_time_dimensions

  contains

  !> @brief Source time dimension from netcdf file
  !> @param[in] path Fully qualified name of netcdf file
  !> @result    Time dimension
  function get_netcdf_time_dim(path) result(time_dim)
    use netcdf, only: nf90_open, nf90_close, nf90_nowrite,                    &
                      nf90_inq_dimid, nf90_inquire_dimension
    implicit none

    character(*), intent(in) :: path
    integer(i_def) :: ncid
    integer(i_def) :: dimid
    integer(i_def) :: ierr
    character(str_def) :: time_name
    integer(i_def) :: time_dim

    ierr = nf90_open(path, NF90_NOWRITE, ncid)
    if (ierr /= 0) then
      write(log_scratch_space,'(A, A, A, I3)')                                &
        'error opening file ', trim(path),                                    &
        ' for input, error code: ', ierr
      call log_event(log_scratch_space, log_level_error)
    end if

    ierr = nf90_inq_dimid(ncid, 'time', dimid)
    if (ierr /= 0) then
      write(log_scratch_space, '(A, A, A, I3)')                               &
        'error inquiring time dimension id in file ',                         &
        trim(path),                                                           &
        ', error code: ', ierr
      call log_event(log_scratch_space, log_level_error)
    end if

    ierr = nf90_inquire_dimension(ncid, dimid, time_name, time_dim)
    if (ierr /= 0) then
      write(log_scratch_space, '(A, A, I3)')                                  &
        'error inquiring time dimension, ',                                   &
        'error code: ', ierr
      call log_event(log_scratch_space, log_level_error)
    end if

    if (time_name /= 'time') then
      write(log_scratch_space, '(A, A)')                                      &
        'unexpected time dimension name, ',                                   &
        time_name
      call log_event(log_scratch_space, log_level_error)
    end if

    ierr = nf90_close(ncid)
    if (ierr /= 0) then
      write(log_scratch_space, '(A, A, A, I3)')                               &
        'error closing file, ', trim(path),                                   &
        ', error code: ', ierr
      call log_event(log_scratch_space, log_level_error)
    end if

  end function get_netcdf_time_dim

  !> @brief Get time dimension of ancil file
  !> @param[in]  dir  Ancil file directory
  !> @param[in]  file Ancil file name (may include subdirectories)
  !> @param[out] tdim Time dimension of ancil file
  !> @result     True if and only if the dimension could be obtained
  function get_ancil_dim(dir, file, tdim) result(status)
    implicit none
    character(*), intent(in) :: dir
    character(*), intent(in) :: file
    integer(i_def), intent(out) :: tdim

    logical(l_def) :: status

    if (file == cmdi) then
      status = .false.
    else if (file(1:1) == '/') then
      tdim = get_netcdf_time_dim(trim(file) // '.nc')
      status = .true.
    else
      tdim = get_netcdf_time_dim(trim(dir) // '/' // trim(file) // '.nc')
      status = .true.
    end if
  end function get_ancil_dim

  !> @brief Synchonise XIOS time axis dimensions with the files being read.
  !>
  subroutine sync_time_dimensions()
    use lfric_xios_diag_mod, only: set_axis_dimension
    implicit none

    logical(l_def), parameter :: tolerate_missing_axes = .true.
    integer(i_def) :: tdim

    ! determine time axis dimensions from ancil file being read
    tdim = get_lbc_dim()
    if (tdim /= 0) &
      call set_axis_dimension('lbc_axis', tdim, tolerate_missing_axes)
    tdim = get_reynolds_dim()
    if (tdim /= 0) &
      call set_axis_dimension('reynolds_timeseries', tdim, &
        tolerate_missing_axes)
    tdim = get_emiss_dim()
    if (tdim /= 0) &
      call set_axis_dimension('emiss_axis', tdim, tolerate_missing_axes)
  end subroutine sync_time_dimensions

  !> @brief Source the dimension of the lbc_axis from the lbc file.
  !>
  !> @result The dimension of the lbc_axis
  !>         or zero if there is no lbc file.
  !>
  function get_lbc_dim() result(tdim)
    implicit none

    integer(i_def) :: tdim
    tdim = 0
    if (lbc_option == lbc_option_gungho_file .or.                             &
        lbc_option == lbc_option_um2lfric_file) then

      if (.not. get_ancil_dim(lbc_dir, lbc_filename, tdim)) return

    end if
  end function get_lbc_dim

  !> @brief Source the dimension of the reynolds_timeseries axis
  !> from the ancil file being read.
  !>
  !> @result The dimension of the reynolds_timeseries axis
  !>         or zero if the relevant ancil file is not enabled.
  !>
  function get_reynolds_dim() result(tdim)
    implicit none

    integer(i_def) :: tdim
    tdim = 0
#ifdef UM_PHYSICS
    if(ancil_option == ancil_option_fixed .or.                                &
       ancil_option == ancil_option_updating ) then

        if (.not. get_ancil_dim(ancil_dir, sst_ancil_path, tdim)) return

    end if
#endif
  end function get_reynolds_dim

  !> @brief Source the dimension of the emiss_axis from the
  !> ancil file being read.
  !>
  !> @result The dimension of the emiss_axis
  !>         or zero if there are no enabled emiss files.
  !>
  function get_emiss_dim() result(tdim)
    implicit none
    integer(i_def) :: tdim
    tdim = 0
#ifdef UM_PHYSICS
    ! conditions and list of emiss files from gungho_setup_io_mod;
    ! to be safe, we try them all in sequence
    if (glomap_mode == glomap_mode_ukca   .and.                               &
        ancil_option == ancil_option_updating) then

      if (get_ancil_dim(ancil_dir, emiss_bc_biofuel_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_bc_fossil_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_bc_biomass_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_bc_biomass_hi_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_bc_biomass_lo_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_om_biofuel_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_om_fossil_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_om_biomass_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_so2_low_ancil_path, tdim)) return
      if (get_ancil_dim(ancil_dir, emiss_so2_high_ancil_path, tdim)) return

    end if
#endif
  end function get_emiss_dim

end module time_dimensions_mod
