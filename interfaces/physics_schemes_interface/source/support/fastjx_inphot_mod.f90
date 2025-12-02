!----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Call reads for fastjx input files
!
!  Description:
!   Fast-jx is an updated routine for calculating online photolysis rates
!   This calls routines to initialise various quantities. Based on the routine
!   fastj_inphot
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
module fastjx_inphot_mod

    use constants_mod,        only : r_def, r_um, i_um
    use log_mod,              only : log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR,   &
                                     LOG_LEVEL_INFO

    implicit none

    private
    public :: fastjx_inphot

  contains

    !>@brief    Calls read subroutines for fastjx input files
    !>@details  Calls read subroutines for fastjx input files

    subroutine fastjx_inphot(                                                  &
                 ! (Max) data dimensions in files
                  max_miesets, n_solcyc_av, sw_band_aer, sw_phases, max_wvl,   &
                  wvl_intervals, max_crossec, num_tvals, n_phot_spc,           &
                 ! Variables used for/ set from cross-section data
                  njval, nw1, nw2, fl, q1d, qo2, qo3, qqq, qrayl, tqq, wl,     &
                  titlej, jlabel, jfacta, jind,                                &
                 ! Variables used for/ set from scatterrer data
                  jtaumx, naa, atau, atau0, daa, paa, qaa, raa, saa, waa,      &
                 ! Variables used for/ set from solar cycle data
                  n_solcyc_ts, solcyc_av, solcyc_quanta, solcyc_ts,        &
                  solcyc_spec )
      
      use fastjx_specs_mod,        only: fastjx_rd_xxx,                        &
                                     fastjx_rd_sol, fastjx_rd_mie
      use chemistry_config_mod,    only: fastjx_dir, fjx_spec_file,            &
                                     fjx_scat_file, fjx_solar_file,            &
                                     fastjx_numwavel, fjx_solcyc_type   

      use errormessagelength_mod,  only: errormessagelength
      use io_utility_mod,          only: claim_io_unit, release_io_unit
      
      implicit none

      ! Arguments - dimensions
      integer(kind=i_um), intent(in) :: max_miesets, n_solcyc_av
      integer(kind=i_um), intent(in) :: sw_band_aer, sw_phases
      integer(kind=i_um), intent(in) :: max_wvl, wvl_intervals
      integer(kind=i_um), intent(in) :: max_crossec, num_tvals
      integer(kind=i_um), intent(in) :: n_phot_spc

      ! Variables set in/ related to x-section data (fjx_rd_xxx)
      integer(kind=i_um), intent(out) :: njval, nw1, nw2
      real(kind=r_um), intent(out) :: fl(max_wvl)
      real(kind=r_um), intent(out) :: q1d(max_wvl, num_tvals)
      real(kind=r_um), intent(out) :: qo2(max_wvl, num_tvals)
      real(kind=r_um), intent(out) :: qo3(max_wvl, num_tvals)
      real(kind=r_um), intent(out) :: qqq(max_wvl, 2, max_crossec)
      real(kind=r_um), intent(out) :: qrayl(max_wvl+1)
      real(kind=r_um), intent(out) :: tqq(num_tvals, max_crossec)
      real(kind=r_um), intent(out) :: wl(max_wvl)

      character(len=*), INTENT(out) :: titlej(max_crossec)

      character(len=*), INTENT(in) :: jlabel(n_phot_spc)
      real(kind=r_um), intent(in) :: jfacta(n_phot_spc)
      integer(kind=i_um), intent(out) :: jind(n_phot_spc)

      ! Variables used for/ set from scatterrer data
      integer(kind=i_um), intent(out) :: jtaumx, naa
      real(kind=r_um), intent(out) :: atau, atau0
      real(kind=r_um), intent(out) :: daa(max_miesets)
      real(kind=r_um), intent(out) :: paa(sw_phases, sw_band_aer, max_miesets)
      real(kind=r_um), intent(out) :: qaa(sw_band_aer, max_miesets)
      real(kind=r_um), intent(out) :: raa(max_miesets)
      real(kind=r_um), intent(out) :: saa(sw_band_aer, max_miesets)
      real(kind=r_um), intent(out) :: waa(sw_band_aer, max_miesets)

      ! Variables used for/ set from solar cycle data
      integer(kind=i_um), intent(in) :: n_solcyc_ts
      real(kind=r_um), intent(out) :: solcyc_av(n_solcyc_av)
      real(kind=r_um), intent(out) :: solcyc_quanta(wvl_intervals)
      real(kind=r_um), intent(out) :: solcyc_ts(n_solcyc_ts)
      real(kind=r_um), intent(out) :: solcyc_spec(max_wvl)

      character(len=255)       :: jv_fullpath
      logical                  :: l_exist
      integer                  :: errorstatus
      integer                  :: ukcafjxx_unit
      integer                  :: ukcafjsc_unit
      integer                  :: ukcafjsol_unit

      character(len=errormessagelength) :: cmessage   ! Error message
      
      write( log_scratch_space, '(A,I6)' ) 'fastjx_numwl=', fastjx_numwavel
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      !--------------------------------------------------
      ! Read in Fast-JX photolysis cross sections
      l_exist = .false.
      jv_fullpath = trim( fastjx_dir )//'/'//trim( fjx_spec_file )
      inquire( file = trim( jv_fullpath ), exist = l_exist )
      if (.not. l_exist) then
        errorstatus = 1
        cmessage = trim( jv_fullpath )//                                      &
                   ': Fast-JX spectral file does not exist'
        write( log_scratch_space, '(A)' ) cmessage
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Read in Fast-JX photolysis cross sections
      ukcafjxx_unit = claim_io_unit()
      call fastjx_rd_xxx( ukcafjxx_unit, jv_fullpath, n_phot_spc, njval,      &
                          nw1, nw2, fl, q1d, qo2, qo3, qqq, qrayl,tqq, wl,    &
                          jlabel, jfacta, jind, titlej )
      call release_io_unit( ukcafjxx_unit )

      !--------------------------------------------------
      ! Read in Fast-JX photolysis solar cycle
      if (fjx_solcyc_type > 0) then
        l_exist = .false.
        jv_fullpath = trim( fastjx_dir )//'/'//trim( fjx_solar_file )
        inquire( file = trim( jv_fullpath ), exist = l_exist )
        if (.not. l_exist) then
          errorstatus = 1
          cmessage = trim( jv_fullpath )//                                    &
                     ': Fast-JX solar cycle file does not exist'
          write( log_scratch_space, '(A)' ) cmessage
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if

        ! Read in Fast-JX photolysis solar cycle
        ukcafjsol_unit = claim_io_unit()
        call fastjx_rd_sol( ukcafjsol_unit, jv_fullpath, n_solcyc_ts,         &
                      solcyc_av, solcyc_quanta, solcyc_ts, solcyc_spec )
        call release_io_unit( ukcafjsol_unit )
      end if

      !----------------------------------------------
      ! Read in Fast-JX scattering cross sections
      l_exist = .false.
      jv_fullpath = trim( fastjx_dir )//'/'//trim( fjx_scat_file )
      inquire( file = trim( jv_fullpath ), exist = l_exist )
      if (.not. l_exist) then
        errorstatus = 1
        cmessage = trim( jv_fullpath )//                                      &
                   ': Fast-JX scattering file does not exist'
        write( log_scratch_space, '(A)' ) cmessage
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Read in Fast-JX photolysis cross sections
      ukcafjsc_unit = claim_io_unit()
      call fastjx_rd_mie( ukcafjsc_unit, jv_fullpath, jtaumx, naa, atau,      &
                         atau0, daa, paa, qaa, raa, saa, waa )
      call release_io_unit( ukcafjsc_unit )
   
      return
    end subroutine fastjx_inphot
  end module fastjx_inphot_mod