!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief UKCA initialisation subroutine for UM science configuration
!
!  Description:
!   Fast-jx is an updated routine for calculating online photolysis rates
!   This module contains routines that read in the phase factors etc from file
!   Based upon the fast_specs routine, though with large differences caused by
!   differences between fast-j and fast-jx
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
module fastjx_specs_mod

    use photol_api_mod, only: fastjx_rd_mie_file => photol_fjx_rd_mie_file,   &
          fastjx_rd_xxx_file => photol_fjx_rd_xxx_file,                       &
          fastjx_rd_sol_file => photol_fjx_rd_sol_file, photol_jlabel_len,    &
          !  Max sizes for spectral file data
          a_ => photol_max_miesets, n_solcyc_av => photol_n_solcyc_av,        &
          sw_band_aer => photol_sw_band_aer, sw_phases => photol_sw_phases,   &
          wx_ => photol_max_wvl, x_ => photol_max_crossec,                    &
          jpwav => photol_wvl_intervals
    
    use chemistry_config_mod, only: fastjx_numwavel
    
    use log_mod, only : log_event,                                            &
                        log_scratch_space,                                    &
                        LOG_LEVEL_ERROR, LOG_LEVEL_INFO,                      &
                        LOG_LEVEL_DEBUG, LOG_LEVEL_WARNING

    use constants_mod, only: r_def, r_um, i_um

    use errormessagelength_mod, only: errormessagelength

    use lfric_mpi_mod,       only: global_mpi

    implicit none

    ! *********************************************

    public fastjx_rd_mie, fastjx_rd_xxx, fastjx_rd_sol

    contains
    
    ! ######################################################################
    !-----------------------------------------------------------------------
    !-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
    !  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
    !-----------------------------------------------------------------------
    subroutine fastjx_rd_mie( nj1,namfil, jtaumx, naa, atau, atau0,           &
                             daa, paa, qaa, raa, saa, waa )

    implicit none

    integer, intent(in) :: nj1  ! Channel number for reading data file
    character(len=*), intent(in) :: namfil ! Name of scattering data file
                                           ! (e.g., FJX_scat.dat)
    integer(kind=i_um), intent(out) :: jtaumx ! Max num cloud sub-layers
    integer(kind=i_um), intent(out) :: naa    ! Num aerosol/ cloud data types
    real(kind=r_um), intent(out) :: atau      ! Cloud sub-layer factor
    real(kind=r_um), intent(out) :: atau0     ! min atau
                                              ! Properties of scattering types:
    real(kind=r_um), intent(out) :: daa(a_)   ! Density
    real(kind=r_um), intent(out) :: paa(sw_phases, sw_band_aer, a_)  ! Phases
    real(kind=r_um), intent(out) :: qaa(sw_band_aer, a_)  ! Q
    real(kind=r_um), intent(out) :: raa(a_)   ! Effective radius
    real(kind=r_um), intent(out) :: saa(sw_band_aer, a_)  ! Single Scat Albedo
    real(kind=r_um), intent(out) :: waa(sw_band_aer, a_)  ! Wavelengths

    integer :: j, k                      ! Loop variables

    ! String containing cloud/aerosol scattering
    character(len=photol_jlabel_len)  ::  aerosol_cloud_title(a_)

    ! ***********************************
    ! End of Header

    if (global_mpi%get_comm_rank() == 0) then
      call fastjx_rd_mie_file( nj1,namfil, jtaumx, naa, atau, atau0, daa,     &
                             paa, qaa, raa, saa, waa, aerosol_cloud_title )
    end if

    call global_mpi%broadcast( naa, 0 )
    call global_mpi%broadcast( jtaumx, 0 )
    call global_mpi%broadcast( atau, 0 )
    call global_mpi%broadcast( atau0, 0 )

    call global_mpi%broadcast( aerosol_cloud_title, photol_jlabel_len*a_, 0 )

    call global_mpi%broadcast( raa, a_, 0 )
    call global_mpi%broadcast( daa, a_, 0 )
    call global_mpi%broadcast( waa, sw_band_aer*a_, 0 )
    call global_mpi%broadcast( qaa, sw_band_aer*a_, 0 )
    call global_mpi%broadcast( saa, sw_band_aer*a_, 0 )
    call global_mpi%broadcast( paa, sw_phases*sw_band_aer*a_, 0 )

    ! Output some information
    write( log_scratch_space, '(A,5F8.1)' )                                   &
            ' Aerosol optical: r-eff/rho/Q(@wavel):', (waa(k,1),k=1,5)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Loop over aerosol types writing radius, density and q
    do j=1,naa
      write( log_scratch_space, '(i3,1x,a8,7f8.3)' ) j,                       &
              aerosol_cloud_title(j), raa(j), daa(j), (qaa(k,j),k=1,5)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    return
    end subroutine fastjx_rd_mie


    !########################################################################
    !-----------------------------------------------------------------------
    !  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
    !
    !>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
    !           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
    !           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
    !           if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
    !           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
    !
    !-----------------------------------------------------------------------
    subroutine fastjx_rd_xxx( nj1, namfil, jppj, njval, nw1, nw2, fl,         &
                         q1d, qo2, qo3, qqq, qrayl, tqq, wl, jlabel,          &
                         jfacta, jind, titlej )

    implicit none

    integer, intent(in)       :: nj1 ! Channel number for reading data file
    character(len=*), intent(in)  :: namfil    ! Name of spectral data file
                                               ! (JX_spec.dat)
    integer(kind=i_um), intent(in)  :: jppj    ! Num of photolytic species
    integer(kind=i_um), intent(out) :: njval   ! Num of species with x-section
    integer(kind=i_um), intent(out) :: nw1     ! Min wavelength bin
    integer(kind=i_um), intent(out) :: nw2     ! Max wavelength bin
    real(kind=r_um), intent(out) :: fl(wx_)    ! TOA solar flux
                                               ! Photolysis rates for:
    real(kind=r_um), intent(out) :: q1d(wx_,3)  ! O1D (3 temperatures)
    real(kind=r_um), intent(out) :: qo2(wx_,3)  ! O2
    real(kind=r_um), intent(out) :: qo3(wx_,3)  ! O3
    real(kind=r_um), intent(out) :: qqq(wx_,2,x_) ! All other species

    real(kind=r_um), intent(out) :: qrayl(wx_+1)  ! Rayleigh parameters
    real(kind=r_um), intent(out) :: tqq(3,x_)     ! temp corresponding to rates
    real(kind=r_um), intent(out) :: wl(wx_)       ! effective wavelengths

    character(len=photol_jlabel_len), intent(in) :: jlabel(jppj)
    real(kind=r_um), intent(in) :: jfacta(jppj)   ! Quantum yields  
    integer(kind=i_um), intent(out) :: jind(jppj)  ! Index of species from file

    character (len=errormessagelength)        :: cmessage
                                               ! String for error handling

    integer ::  j, k                ! Loop variables

    ! List of species being photolysed, as read from file
     character(len=photol_jlabel_len), intent(out) :: titlej(x_)

    ! *****************************

    ! Initialise the temperatures to 0
    do j = 1,x_
      do k = 1,3
        tqq(k,j) = 0.0e0_r_def
      end do
    end do

    !----------spectral data----set for new format data J-ver8.3------------------
    !         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
    !         for 2005a data, NJVAL = 62 (including a spare XXXX) and
    !              NQQQ = 64 so that 4 wavelength datasets read in for acetone
    !         note NQQQ is not used outside this subroutine!
    ! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

    if (global_mpi%get_comm_rank() == 0) then
      call fastjx_rd_xxx_file( nj1, namfil, fastjx_numwavel, njval, nw1, nw2, &
                              fl, q1d, qo2, qo3, qqq, qrayl, tqq, wl, titlej )
    end if

    call global_mpi%broadcast( njval, 0 )
    call global_mpi%broadcast( nw1, 0 )
    call global_mpi%broadcast( nw2, 0 )

    call global_mpi%broadcast( wl, wx_, 0 )
    call global_mpi%broadcast( fl, wx_, 0 )
    call global_mpi%broadcast( qrayl, wx_+1, 0 )

    call global_mpi%broadcast( titlej, x_*photol_jlabel_len, 0 )
    call global_mpi%broadcast( tqq, 3*x_, 0 )
    call global_mpi%broadcast( qo2, wx_*3, 0 )
    call global_mpi%broadcast( qo3, wx_*3, 0 )
    call global_mpi%broadcast( q1d, wx_*3, 0 )
    call global_mpi%broadcast( qqq, wx_*2*x_, 0 )


    ! *************************************************************
    ! Map local indices to UKCA ones

    do j = 1,njval
      do k = 1,jppj
        if ( k == 1 ) then
          write( log_scratch_space, '(A,A,A)' ) 'FASTJX Compare titles ',     &
                  titlej(j), jlabel(k)
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
        end if
        if (jlabel(k) == titlej(j)) jind(k) = j
      end do
    end do

    do k = 1,jppj

      write( log_scratch_space, '(A,A,A)' )                                   &
              'Comparing J-rate for photolysis reaction ', jlabel(k), ' ?'
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write( log_scratch_space, '(A,I0)' ) 'Using index ', jind(k)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if (jfacta(k) < 1.0e-20_r_def) then
        write( log_scratch_space, '(A,A)' ) 'Not using photolysis reaction ',  &
                jlabel(k)
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
      end if

      if (jind(k) == 0) then

        if (jfacta(k) < 1.0e-20_r_def) then
          jind(k) = 1
        else
          cmessage =                                                          &
           ' Which J-rate for photolysis reaction '//jlabel(k)//' ?'
          write( log_scratch_space, '(A)' ) cmessage
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

    return
    end subroutine fastjx_rd_xxx

    !########################################################################
    !-----------------------------------------------------------------------
    !  Read in solar cycle data
    !-----------------------------------------------------------------------
    subroutine fastjx_rd_sol( nj1, namfil, n_solcyc_ts, solcyc_av,            &
                             solcyc_quanta, solcyc_ts, solcyc_spec )   

    implicit none

    integer, intent(in)       :: nj1 ! Channel number for reading data file
    character(len=*), intent(in)  :: namfil    ! Name of solar cycle data file

    integer(kind=i_um), intent(in) :: n_solcyc_ts  ! Total months of solcyc val
    real(kind=r_um), intent(out) :: solcyc_av(n_solcyc_av)  ! Avg solar cycle
    real(kind=r_um), intent(out) :: solcyc_quanta(jpwav)    ! Quantum component
    real(kind=r_um), intent(out) :: solcyc_ts(n_solcyc_ts)  ! Observed tseries
    real(kind=r_um), intent(out) :: solcyc_spec(wx_)       ! Spectral component

    integer(kind=i_um)           :: n_solcyc_ts_file    

    ! Pointer to pass to read routine
    real(kind=r_um), pointer :: solcyc_ts_ptr(:)

    if (global_mpi%get_comm_rank() == 0) then
      call fastjx_rd_sol_file( nj1, namfil, n_solcyc_ts_file, solcyc_av,      &
                             solcyc_quanta, solcyc_ts_ptr, solcyc_spec )
      if (n_solcyc_ts_file /= n_solcyc_ts) then
        write( log_scratch_space, '(2(A,I0),A)' )                              &
          'Solar Cycle file records do not match specified value: expected ',  &
          n_solcyc_ts, ', found: ', n_solcyc_ts_file,                          &
          '. Check chemistry namelist value fjx_solcyc_months '
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      ! Copy data to array, to enable broadcast      
      solcyc_ts(:) = solcyc_ts_ptr(:)
    end if

    if (associated( solcyc_ts_ptr )) nullify( solcyc_ts_ptr )

    call global_mpi%broadcast( n_solcyc_ts_file, 0 )
    call global_mpi%broadcast( solcyc_av, n_solcyc_av, 0 )
    call global_mpi%broadcast( solcyc_ts, n_solcyc_ts, 0 )
    call global_mpi%broadcast( solcyc_spec, wx_, 0 )
    call global_mpi%broadcast( solcyc_quanta, jpwav, 0 )

    end subroutine fastjx_rd_sol


    !#######################################################################

    end module fastjx_specs_mod
