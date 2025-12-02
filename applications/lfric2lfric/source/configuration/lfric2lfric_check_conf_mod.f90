!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief     Checks lfric2lfric configuration namelist.

!> @details   Imports the lfric2lfric configuration namelist and checks that
!!            no configuration options conflict. Specifically, checks that the
!!            regridding method and the interpolation direction are supported.
module lfric2lfric_check_conf_mod

  use constants_mod,          only : i_def
  use log_mod,                only : log_event, log_scratch_space, &
                                     LOG_LEVEL_ERROR
  use namelist_mod,           only : namelist_type

  ! Configuration modules
  use lfric2lfric_config_mod, only : ORIGIN_DOMAIN_LAM,         &
                                     TARGET_DOMAIN_LAM,         &
                                     regrid_method_lfric2lfric, &
                                     key_from_target_domain


  implicit none

  private
  public lfric2lfric_check_configuration

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief     Checks namelist configuration unsupported combinations of valid
  !!            options.
  !> @details   Currently uses a pointer to the lfric2lfric namelist to check
  !!            the lfric2lfric configuration options for any invalid
  !!            combinations of config options.
  !> @param [in]   lfric2lfricnml   Pointer to the lfric2lfric namelist
  !>                                configuration object.
  subroutine lfric2lfric_check_configuration( lfric2lfric_nml )

    implicit none

    type(namelist_type), pointer, intent(in) :: lfric2lfric_nml

    ! Namelist options are enumerated so become randomized integers
    integer(kind=i_def) :: origin_domain
    integer(kind=i_def) :: target_domain
    integer(kind=i_def) :: regrid_method

    ! Extract target and origin domains from config namelist
    call lfric2lfric_nml%get_value( 'origin_domain', origin_domain )
    call lfric2lfric_nml%get_value( 'target_domain', target_domain )
    call lfric2lfric_nml%get_value( 'regrid_method', regrid_method )

    ! Check our origin and target domains are compatible, currently the
    ! interpolations NOT allowed are lam => lbc and lam => global
    if( origin_domain == ORIGIN_DOMAIN_LAM .and. &
        target_domain /= TARGET_DOMAIN_LAM ) then
        write( log_scratch_space, '(A)' ) &
               'Invalid target domain "'                                // &
               trim(key_from_target_domain(target_domain))              // &
               '" detected for lfric2lfric origin_domain "lam", valid ' // &
               'options are: [ '                                        // &
               trim(key_from_target_domain(TARGET_DOMAIN_LAM))          // &
               ' ]'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Check that the regridding method is supported by lfric2lfric, currently
    ! we do not support weights-lfric2lfric regridding
    if( regrid_method == regrid_method_lfric2lfric ) then
      write( log_scratch_space, '(A)' ) &
             'Weights-lfric2lfric regridding is not currently supported in lfric2lfric'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine lfric2lfric_check_configuration

end module lfric2lfric_check_conf_mod
