! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module dms_flux_mod_4A

use um_types, only: real_umphys

implicit none

! Indicators for sea-air flux calculation methods
integer, parameter :: i_liss_merlivat = 1
integer, parameter :: i_wanninkhof = 2
integer, parameter :: i_nightingale = 3

character(len=*), parameter, private :: ModuleName='DMS_FLUX_MOD_4A'

contains

subroutine dms_flux_4A(                                                        &
  ! Arguments in
  row_length, rows,                                                            &
  wind_10m, tstar, land_fract, dms_conc,                                       &
  i_dms_flux,                                                                  &
  ! Arguments out
  f_dms )
!---------------------------------------------------------------------
! Purpose: To calculate the flux of DMS (as kg m-2 s-1 of sulphur)
!          from the ocean surface as a function of its concentration
!          in seawater and of windspeed. The sea-air exchange can
!          be determined according to one of three commonly-used
!          parametrization schemes, those of Liss & Merlivat (1986),
!          Wanninkhof (1992) or Nightingale et al. (2000). The routine
!          is called by Aero_Ctl.

! Method:  The Schmidt number
!          for DMS is calculated as in Saltzman et al. (1993), and
!          used with the windspeed to determine the mass transfer (or
!          "piston") velocity according to the desired parametrization.
!          This is then used to determine the sea-air mass flux of DMS
!          as a function of sea-water DMS concentration. High surface
!          temperatures (caused by the land portion of a gridbox when
!          coastal tiling is not active) cause negative Sc values which
!          would give a floating-point error in the k_DMS calculation,
!          so the Tstar values are capped. This shouldn't be a problem
!          when coastal tiling is on as then the Tstar values passed in
!          are those for sea only.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP 003 programming standards

!---------------------------------------------------------------------

use conversions_mod, only: zerodegc
use water_constants_mod, only: tfs

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use umPrintMgr, only: umPrint, umMessage
use ereport_mod, only: EReport

implicit none

! Arguments with intent in:
integer :: row_length
integer :: rows

real(kind=real_umphys)    :: wind_10m(row_length, rows)   ! 10m windspeed (ms-1)
real(kind=real_umphys)    :: tstar(row_length, rows)
                                        ! Surface temperature (K)
real(kind=real_umphys)    :: land_fract(row_length, rows)
                                        ! Fraction of land in gridbox
real(kind=real_umphys)    :: dms_conc(row_length, rows)
                                  ! Concentration of DMS in seawater (nmol l-1)

! Switch to determine which scheme to use to calculate mass transfer velocity
integer :: i_dms_flux

! Arguments with intent out:

real(kind=real_umphys)    :: f_dms(row_length, rows)
                                         ! Sea-air flux of DMS (kg[S] m-2 s-1)

! Local variables:
integer  :: i,j                          ! Loop counters
real(kind=real_umphys)     :: sc(row_length, rows)         ! Schmidt number
real(kind=real_umphys)     :: k_dms(row_length, rows)
                                         ! Piston velocity of DMS (cm h-1)
real(kind=real_umphys)     :: t_c        ! Surface temperature degrees Celsius
! Piston velocities for gases with Schmidt numbers of 600 & 660 resp. (cm h-1)
real(kind=real_umphys)     :: k_600
real(kind=real_umphys)     :: k_660
real(kind=real_umphys)     :: n          ! Schmidt number exponent
real(kind=real_umphys), parameter :: t_max=47.0
                                         ! Max T to avoid breaking Sc fit (C)

integer  :: errcode

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DMS_FLUX_4A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the Schmidt number (Sc):

do j = 1, rows
  do i = 1, row_length
    t_c = min((max(tstar(i, j), tfs) - zerodegc), t_max)
    sc(i, j) = 2674.0 - (147.12*t_c) + (3.726*t_c**2)                          &
      - (0.038*t_c**3)
  end do
end do

! Determine the mass transfer (or "piston") velocity (k_DMS) over sea
! according to the specified parametrization scheme:

select case (i_dms_flux)

case (i_liss_merlivat)
  do j = 1, rows
    do i = 1, row_length
      if (wind_10m(i, j)  <=  3.6) then
        k_600 = 0.17 * wind_10m(i, j)
        n = -2.0/3.0
      end if
      if (wind_10m(i, j)  >   3.6 .and.                                        &
        wind_10m(i, j)  <=  13.0) then
        k_600 = (2.85*wind_10m(i, j)) - 9.65
        n = -0.5
      end if
      if (wind_10m(i, j)  >   13.0) then
        k_600 = (5.9*wind_10m(i, j)) - 49.3
        n = -0.5
      end if
      if (land_fract(i, j)  <   1.0) then
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      else
        k_dms(i, j) = 0.0
      end if
    end do
  end do

case (i_wanninkhof)
  do j = 1, rows
    do i = 1, row_length
      k_660 = 0.31 * wind_10m(i, j)**2
      n = -0.5
      if (land_fract(i, j)  <   1.0) then
        k_dms(i, j) = k_660 * (sc(i, j)/660.0)**n
      else
        k_dms(i, j) = 0.0
      end if
    end do
  end do

case (i_nightingale)
  do j = 1, rows
    do i = 1, row_length
      k_600 = (0.222*wind_10m(i, j)**2) + (0.333*wind_10m(i, j))
      n = -0.5
      if (land_fract(i, j)  <   1.0) then
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      else
        k_dms(i, j) = 0.0
      end if
    end do
  end do

case DEFAULT
  write(umMessage,'(A,I4,A)') 'DMS_FLUX_4A, i_dms_flux=',i_dms_flux,           &
                      ' is not a valid choice.'
  call umPrint(umMessage,src='dms_flux_4a')
  errcode = 1
  call EReport('dms_flux_4a', errcode, 'i_dms_flux not valid')

end select

! Finally, calculate the sea-air flux of DMS as a function of k_DMS
! and dissolved DMS concentration. The former requires a conversion
! from cm hour-1 to ms-1, and the latter from nanomoles per litre to
! kg[S] m-3, to return the flux in kg[S] m-2 sec-1.

do j = 1, rows
  do i = 1, row_length
    f_dms(i, j) = (k_dms(i, j) / 3.6e5)                                        &
      * (dms_conc(i, j) * 32.0e-9)
  end do
end do

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine dms_flux_4A

end module
