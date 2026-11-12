! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Reads in the pre-computed variables for use in the
!  interaction between UKCA aerosols and the radiation
!  code.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
module ukca_radaer_read_precalc_mod

use ukca_radaer_precalc, only:                                                 &
    npd_integ_pts,                                                             &
    npd_ukca_band,                                                             &
    npd_ukca_spectrum,                                                         &
    npd_ukca_maxcomptype,                                                      &
    npd_ukca_aod_wavel,                                                        &
    precalc

implicit none

type :: my_namelist
  sequence
  integer :: n_integ_pts
  integer :: n_ukca_aod_wavel
  real :: wavelength(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)
  real :: irrad(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)
  real :: weight(npd_integ_pts)
  real :: flux(npd_ukca_band, npd_ukca_spectrum)
  real :: realrefr(npd_ukca_maxcomptype, npd_integ_pts, npd_ukca_band,         &
                   npd_ukca_spectrum)
  real :: imagrefr(npd_ukca_maxcomptype, npd_integ_pts, npd_ukca_band,         &
                   npd_ukca_spectrum)
  real :: aod_wavel(npd_ukca_aod_wavel)
  real :: aod_realrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)
  real :: aod_imagrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)
end type my_namelist

type (my_namelist) :: my_nml


character(len=*), parameter, private ::                                        &
  ModuleName = 'UKCA_RADAER_READ_PRECALC_MOD'

contains



subroutine ukca_radaer_read_precalc( filename,                                 &
                                     sw_wavelength_short,                      &
                                     sw_wavelength_long,                       &
                                     lw_wavelength_short,                      &
                                     lw_wavelength_long,                       &
                                     nbr_band_sw,                              &
                                     nbr_band_lw )

use ereport_mod,            only: ereport
use errormessagelength_mod, only: errormessagelength
use file_manager,           only: assign_file_unit, release_file_unit
use filenamelength_mod,     only: filenamelength
use parkind1,               only: jpim, jprb
use spcrg3a_mod,            only: ip_solar, ip_infra_red
use umprintmgr,             only: newline
use yomhook,                only: lhook, dr_hook
use ukca_mode_setup,        only: cp_bc
use ukca_option_mod,        only:                                              &
    ! Logical and real to apply scaling to BC refractive indices
    l_ukca_scale_ppe, bc_refrac_im_scaling

implicit none

!
! Arguments
!

! Name and path of precalc file
character(len=filenamelength), intent(in) :: filename

! Number of wavebands in shortwave spectrum decomposition
integer, intent(in) :: nbr_band_sw

! Number of wavebands in longwave spectrum decomposition
integer, intent(in) :: nbr_band_lw

! Lower boundary of shortwave waveband
real, intent(in) :: sw_wavelength_short( npd_ukca_band )

! Lower boundary of longwave waveband
real, intent(in) :: lw_wavelength_short( npd_ukca_band )

! Upper boundary of shortwave waveband
real, intent(in) :: sw_wavelength_long(  npd_ukca_band )

! Upper boundary of longwave waveband
real, intent(in) :: lw_wavelength_long(  npd_ukca_band )

!
! Local variables
!

!
! Error indicator
!
integer :: ierr
!
! Error message
!
character (len=errormessagelength) :: cmessage
character (len=errormessagelength) :: iomessage
character (len=*), parameter :: RoutineName='UKCA_RADAER_READ_PRECALC'
!
! Local variables
!
integer :: ios


character (len=9) :: spectrum_name (npd_ukca_spectrum)

integer :: n_band

integer :: aerosol_index

real :: wl_short
real :: wl_long
!
! Number of bands and AOD wavelengths in the spectral files.
! These must match the values read from the precalc file.
!
integer :: specf_n_band(npd_ukca_spectrum)
real :: wavelength_short(npd_ukca_spectrum, npd_ukca_band)
real :: wavelength_long(npd_ukca_spectrum, npd_ukca_band)

integer :: funit

character (len=80) :: line ! buffer for neglected lines (headers,...)

!
! Loop indices
!
integer :: i, j, k, n

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

! New real for applying PPE scaling to the imaginary part of the BC refractive
! index
real :: bc_refrac_im_tmp

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

wavelength_short(ip_infra_red,1:nbr_band_lw)=lw_wavelength_short(1:nbr_band_lw)
wavelength_long( ip_infra_red,1:nbr_band_lw)=lw_wavelength_long( 1:nbr_band_lw)
wavelength_short(ip_solar,    1:nbr_band_sw)=sw_wavelength_short(1:nbr_band_sw)
wavelength_long( ip_solar,    1:nbr_band_sw)=sw_wavelength_long( 1:nbr_band_sw)
specf_n_band(ip_solar)           = nbr_band_sw
specf_n_band(ip_infra_red)       = nbr_band_lw

! Check precalc_filename is not blank
if (len_trim(filename) == 0) then
  ierr = 1
  cmessage = 'UKCA pre-computed filename is blank'
  call ereport(routinename, ierr, cmessage)
end if

call assign_file_unit(filename, funit, handler="fortran")
open(unit=funit, file=filename, action='read', iostat=ios, iomsg=iomessage)

if (ios /= 0) then
  ierr = 1
  cmessage = 'Error opening UKCA pre-computed file:' //newline//               &
              trim(filename) //newline//                                       &
              'IoMsg: '//trim(iomessage)
  call ereport(routinename,ierr,cmessage)
end if

! Read the formatted file.

! Number of integration points and weight of each point in the integration.
read(unit=funit, fmt='(29x,i3)') precalc%n_integ_pts
if (precalc%n_integ_pts > npd_integ_pts) then
  ierr = 1
  cmessage = 'Increase the maximum number of integration points in ' //        &
             'UKCA_RADAER pcalc file ' // filename
  close(funit)
  call release_file_unit(funit, handler="fortran")
  call ereport(routinename,ierr,cmessage)
end if

read(unit=funit, fmt='(a)') line ! column header

do n = 1, precalc%n_integ_pts
  read(unit=funit, fmt='(4x,E12.5)') precalc%weight(n)
end do ! n

! For each spectrum
do i = 1, npd_ukca_spectrum

  read(unit=funit, fmt='(a)') line ! spectrum header

  ! Number of wavebands

  spectrum_name = [ 'shortwave', 'longwave ' ]

  read(unit=funit, fmt='(20x,i3)') n_band
  if (n_band /= specf_n_band(i)) then
    ierr = 1
    cmessage = 'UKCA RADAER pcalc file ' // filename // newline //             &
               ' is not compatible with the ' //                               &
               spectrum_name(i) // ' spectral file.' // newline //             &
               'i = '  //              trim(char( i )) // newline //           &
               'n_band = ' //          trim(char( n_band )) // newline //      &
               'specf_n_band(i) = ' // trim(char( specf_n_band(i) ))

    close(funit)
    call release_file_unit(funit, handler="fortran")
    call ereport(routinename,ierr,cmessage)
  end if

  if (n_band > npd_ukca_band) then
    ierr = 1
    cmessage = 'Too many wavebands. Increase the maximum number of wavebands.'
    close(funit)
    call release_file_unit(funit, handler="fortran")
    call ereport(routinename,ierr,cmessage)
  end if

  ! Loop on all wavebands
  do j = 1, n_band

    read(unit=funit, fmt='(14x,E12.5,3x,E12.5,2x)') wl_short, wl_long
    if ( (abs(wl_short-wavelength_short(i,j)) / wl_short  >=  1.0e-5) .or.     &
         (abs(wl_long -wavelength_long(i,j) ) / wl_long   >=  1.0e-5)) then
      ierr = 1
      cmessage =                                                               &
      'Wavelength limits in RADAER PCALC.ukca file do not match those in ' //  &
                   spectrum_name(i) // ' spectral file.'
      close(funit)
      call release_file_unit(funit, handler="fortran")
      call ereport(routinename,ierr,cmessage)
    end if

    read(unit=funit, fmt='(a)') line ! column header

    ! Wavelength at integration point and monochromatic solar
    ! irradiance or Planckian at this wavelength
    do n = 1, precalc%n_integ_pts

      read(unit=funit, fmt='(6x,2(E12.5,1x))')                                 &
          precalc%wavelength(n, j, i), precalc%irrad(n, j, i)

    end do ! n (integration points)

    ! Waveband-integrated irradiance or Planckian
    read(unit=funit, fmt='(27x,E12.5)') precalc%flux(j, i)

    ! Complex refractive index of aerosols at integration points.
    read(unit=funit, fmt='(a)') line ! header

    do k = 1, npd_ukca_maxcomptype

      read(unit=funit, fmt='(25x,i2)') aerosol_index

      read(unit=funit, fmt='(a)') line ! column header

      do n = 1, precalc%n_integ_pts

        ! If requested, scale the imaginary part of the refractive index
        ! for black carbon
        if (l_ukca_scale_ppe .and.                                             &
            aerosol_index == cp_bc .and.                                       &
            spectrum_name(i) == 'shortwave') then

          read(unit=funit, fmt='(8x,2(E12.5,1x))')                             &
            precalc%realrefr(aerosol_index, n, j, i),                          &
            bc_refrac_im_tmp

          precalc%imagrefr(aerosol_index, n, j, i) = bc_refrac_im_tmp *        &
                                                     bc_refrac_im_scaling

          ! Otherwise just read normally
        else

          read(unit=funit, fmt='(8x,2(E12.5,1x))')                             &
            precalc%realrefr(aerosol_index, n, j, i),                          &
            precalc%imagrefr(aerosol_index, n, j, i)

        end if

      end do ! n

    end do ! k

  end do ! j (wavebands)

end do ! i (spectrum)

! Number and values of wavelengths for the aerosol optical depth
read(unit=funit, fmt='(a)') line ! header

read(unit=funit, fmt='(22x,i3)') precalc%n_ukca_aod_wavel
if (precalc%n_ukca_aod_wavel > npd_ukca_aod_wavel) then
  ierr = 1
  cmessage = 'Increase npd_ukca_aod_wavel in UKCAPCALC.'
  close(funit)
  call release_file_unit(funit, handler="fortran")
  call ereport(routinename,ierr,cmessage)
end if

read(unit=funit, fmt='(a)') line ! column header

do i = 1, precalc%n_ukca_aod_wavel
  read(unit=funit, fmt='(4x,E12.5)') precalc%aod_wavel(i)
end do ! i

! Complex refractive index of aerosols at AOD wavelengths.
read(unit=funit, fmt='(a)') line ! header

do k = 1, npd_ukca_maxcomptype

  read(unit=funit, fmt='(25x,i2)') aerosol_index

  read(unit=funit, fmt='(a)') line ! column header

  do n = 1, precalc%n_ukca_aod_wavel

    ! If requested, scale the imaginary part of the refractive index
    ! for black carbon
    if (l_ukca_scale_ppe .and.                                                 &
        aerosol_index == cp_bc) then

      read(unit=funit, fmt='(8x,2(E12.5,1x))')                                 &
        precalc%aod_realrefr(aerosol_index, n),                                &
        bc_refrac_im_tmp

      precalc%aod_imagrefr(aerosol_index, n) = bc_refrac_im_tmp *              &
                                               bc_refrac_im_scaling
      ! Otherwise just read normally
    else

      read(unit=funit, fmt='(8x,2(E12.5,1x))')                                 &
        precalc%aod_realrefr(aerosol_index, n),                                &
        precalc%aod_imagrefr(aerosol_index, n)

    end if

  end do ! n

end do ! k

close(funit)
call release_file_unit(funit, handler="fortran")

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

return
end subroutine ukca_radaer_read_precalc

end module ukca_radaer_read_precalc_mod
