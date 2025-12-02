! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   Module to treat the growth of CLASSIC aerosol particles with humidity.
!
! Method:
!   The subroutine changes the logarithms of the geometric mean radius
!   (or median radius) and geometric standard deviation of a log-normal
!   distribution of particle sizes, to account for growth due to humidity.
!   The routine has been adapted from an analogous one present in the
!   Edwards-Slingo Offline Radiation Code and in SOCRATES (Suite Of
!   Community RAdiative Transfer codes based on Edwards and Slingo).
!   Note also that the subroutine hygro_fact in this code section (i.e.
!   CLASSIC aerosols) does a similar job but only for ammonium sulphate.
!
!   In the case of ammonium sulphate, sea-salt and ammonium nitrate:
!   * For humidities above the deliquescence point, the scheme applied
!     is the one due to J. W. Fitzgerald, J. Appl. Met., 14, 1044-1049,
!     1975. It states that: wr = alpha * r ** beta
!     where:
!       wr          = wet radius (microns)
!       r           = dry radius (microns)
!       alpha, beta = parameters calculated as a function of humidity
!   * There is no growth for humidities below the efflorescence point
!     (i.e. both alpha and beta are 1).
!   * As an approximation, the hysteresis effect (i.e. the dependence
!     on the past state) is modelled by a region where growth is linear.
!     For that it is assumed that Fitzgerald's formulation is valid for
!     humidities between the efflorescence and deliquescence points, with
!     alpha and beta parameters being interpolated between the values at
!     these two points. The justification for this is that near the
!     deliquescence point, the probablility that a particle has previously
!     been exposed to a humidity in excess of that needed for deliquescence
!     is high, so that on average the upper branch of the growth curve
!     should receive a higher weighting than the lower branch. At low
!     humidities the reverse applies.
!
!   Note that Fitzgerald's formulation has been derived for radii in
!   microns. However, the natural logarithm of the dry radius (input
!   argument) is given for a radius in metres. As a consequence the
!   parameter alpha needs to account for some unit conversions as
!   indicated below.
!
!   Fitzgerald's formulation has been derived for ammonium sulphate.
!   The same formulation is valid for sea-salt aerosol and ammonium
!   nitrate, with the following modifications:
!   * The efflorescence / deliquescence points are different.
!   * Alpha for sea-salt (ammonium nitrate) is 1.35 (1.06) times its
!     value for ammonium sulphate, while beta has the same value
!     for all substances.
!   (For sea-salt see Tang et al., J. Aerosol Sci., 8,149-159, 1977)
!
!   Extra comments regarding  biomass burning, fossil-fuel organic carbon
!   (OCFF) and biogenic aerosols: Hygroscopic growth is prescribed using
!   arrays of growth factors as a function of relative humidity. OCFF
!   aerosols use the same factors as biomass burning aerosols.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: aerosols
!
! Code description:
!   Language: FORTRAN 2003.
!   This code is written to UMDP3 programming standards.
!----------------------------------------------------------------------

module grow_particles_mod

use rad_pcf,     only: ip_accum_sulphate, ip_aitken_sulphate,                  &
  ip_fresh_soot,  ip_aged_soot,  ip_biomass_1,     ip_biomass_2,               &
  ip_ocff_fresh,  ip_ocff_aged,  ip_seasalt_film,  ip_seasalt_jet,             &
  ip_biogenic,    ip_nitrate
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim
use ereport_mod,            only: ereport
use errormessagelength_mod, only: errormessagelength

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='GROW_PARTICLES_MOD'

contains

subroutine grow_particles (aerosol_type, humidity, ln_r0, ln_sigma,            &
                           ln_wr0,  ln_wsigma)

use umPrintMgr,             only: umPrint, umMessage

implicit none

! Subroutine arguments
integer, intent(in) :: aerosol_type     ! indicator for aerosol type

real(kind=real_umphys),    intent(in) :: humidity
                                        ! relative humidity (0-1)

! Natural logarithm of geometric mean radius of log-normal distribution   (dry)
! Note that radius is expressed in metres.
real(kind=real_umphys), intent(in) :: ln_r0

! Natural logarithm of geometric standard deviation of log-normal distrib (dry)
real(kind=real_umphys), intent(in) :: ln_sigma

! Natural logarithm of geometric mean radius        (with hygroscopic growth)
! Note that radius is expressed in metres.
real(kind=real_umphys), intent(out) :: ln_wr0

! Natural logarithm of geometric standard deviation (with hygroscopic growth)
real(kind=real_umphys), intent(out) :: ln_wsigma

! Local variables
real(kind=real_umphys) :: growth_factor ! growth factor for dilution
real(kind=real_umphys) :: alpha         ! Fitzgerald's alpha parameter
real(kind=real_umphys) :: beta          ! Fitzgerald's beta parameter
real(kind=real_umphys) :: alpha_deliq   ! value of alpha at RH (deliquescence)
real(kind=real_umphys) :: beta_deliq    ! value of beta  at RH (deliquescence)
real(kind=real_umphys) :: power         ! exponent for units conversion factor

! RH capped to 0.995 (i.e. 99.5%) if Fitzgerald used,
! or to 1.0 (i.e. 100%) if lookup table.
real(kind=real_umphys) :: humidity_eff


! Relative humidity of the deliquescence and efflorescence points:
! deliq & effl:
!
! for ammonium sulphate
real(kind=real_umphys), parameter :: deliq_so4 = 0.81
real(kind=real_umphys), parameter :: effl_so4  = 0.30
!
! for sea-salt
real(kind=real_umphys), parameter :: deliq_ss  = 0.75
real(kind=real_umphys), parameter :: effl_ss   = 0.42
!
! for ammonium nitrate
real(kind=real_umphys), parameter :: deliq_no3 = 0.61
real(kind=real_umphys), parameter :: effl_no3  = 0.30
!
! Generic values of RH(deliquescence) and RH(efflorescence)
! for any aerosol type:
real(kind=real_umphys):: deliq
real(kind=real_umphys):: effl

integer                            :: ierr      ! error flag
character (len=errormessagelength) :: cmessage  ! error message

! Growth factors for biomass-burning aerosols as a function of relative
! humidity. They were computed to match water uptake properties observed
! during the SAFARI campaign:
! Magi, B. I., and P. V. Hobbs, Effects of humidity on aerosols in
! southern Africa during the biomass burning season, J. Geophys. Res.,
! 108(D13), 8495, doi:10.1029/2002JD002144, 2003.
! They are given for 21 values of relative humidity, ranging from 0 to
! 100% with a step of 5% (but note that humidity here is expressed as
! 0-1 instead of 0-100).
!
integer, parameter :: nhum = 21 ! total no. of RH values in the table
integer            :: i_hum     ! no. of element currently being used
!
real(kind=real_umphys),    parameter :: humidity_step = 0.05
!
! Fresh biomass
real(kind=real_umphys), parameter ::                                           &
                   biom_gf1(nhum) = [ 1.0000, 1.0000, 1.0000,                  &
                                       1.0000, 1.0000, 1.0025,                 &
                                       1.0045, 1.0070, 1.0085,                 &
                                       1.0100, 1.0180, 1.0280,                 &
                                       1.0430, 1.0600, 1.0850,                 &
                                       1.1100, 1.1450, 1.1870,                 &
                                       1.2350, 1.2850, 1.3430 ]
!
! Aged biomass
real(kind=real_umphys), parameter :: biom_gf2(nhum) = [ 1.000, 1.000, 1.000,   &
                                       1.000, 1.000, 1.005,                    &
                                       1.015, 1.025, 1.040,                    &
                                       1.060, 1.080, 1.115,                    &
                                       1.150, 1.185, 1.230,                    &
                                       1.280, 1.333, 1.390,                    &
                                       1.447, 1.507, 1.573 ]
!
! Growth factors for biogenic aerosols as a function of relative humidity.
! They are taken from Varutbangkul et al.: Hygroscopicity of secondary
! organic aerosols formed by oxidation of cycloalkenes, monoterpenes,
! sesquiterpenes, and related compounds, Atmos. Chem. Phys., 6, 2367-2388,
! doi:10.5194/acp-6-2367-2006, 2006. They are given for 21 values of
! relative humidity, ranging from 0 to 100% with a step of 5%.
real(kind=real_umphys) :: biogenic_gf(nhum) =   [ 1.00000, 1.00055, 1.00113,   &
                                 1.00186, 1.00281, 1.00401,                    &
                                 1.00553, 1.00744, 1.00983,                    &
                                 1.01279, 1.01646, 1.02100,                    &
                                 1.02663, 1.03366, 1.04257,                    &
                                 1.05409, 1.06959, 1.09195,                    &
                                 1.12896, 1.21524, 1.69764 ]

integer (kind=jpim), parameter :: zhook_in  = 0
integer (kind=jpim), parameter :: zhook_out = 1
real    (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='GROW_PARTICLES'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise error flag
ierr = 0

! Select type of aerosol
!
! -----------------------------------------------------------------
!
! Use Fitzgerald formulation for ammonium sulphate, sea-salt and
! ammonium nitrate.
if ( ( aerosol_type  ==  ip_accum_sulphate )  .or.                             &
     ( aerosol_type  ==  ip_aitken_sulphate ) .or.                             &
     ( aerosol_type  ==  ip_seasalt_film )    .or.                             &
     ( aerosol_type  ==  ip_seasalt_jet )     .or.                             &
     ( aerosol_type  ==  ip_nitrate ) ) then

  ! Setting RH of deliquescescence and efflorescence points
  ! for ammonium sulphate.
  if ( ( aerosol_type  ==  ip_accum_sulphate ) .or.                            &
       ( aerosol_type  ==  ip_aitken_sulphate ) ) then
    deliq = deliq_so4
    effl  = effl_so4

    ! Setting RH of deliquescescence and efflorescence points
    ! for sea-salt.
  else if ( ( aerosol_type  ==  ip_seasalt_film ) .or.                         &
            ( aerosol_type  ==  ip_seasalt_jet ) ) then
    deliq = deliq_ss
    effl  = effl_ss

    ! Setting RH of deliquescescence and efflorescence points
    ! for ammonium nitrate.
  else if ( aerosol_type  ==  ip_nitrate ) then
    deliq = deliq_no3
    effl  = effl_no3

    ! Trap incorrect aerosol types here. Even if this looks redundant it may
    ! be safe to keep it in case any new aerosol type is added above.
  else
    ierr     = 1
    cmessage = 'Unavailable aerosol composition type for Fitzgerald formulation'
    call ereport ('GROW_PARTICLES_MOD', ierr, cmessage)
  end if

  ! Cap humidity to the max value allowed by Fitzgerald: 99.5%.
  ! Note that this ceiling is needed because of the highly non-
  ! linear growth as saturation is approached; in any case, at
  ! such high humidities the impact of aerosol will be swamped
  ! by the effect of cloud.
  humidity_eff = min (humidity, 0.995)

  ! Initialise to default values of alpha and beta (1.0 means that
  ! there is no growth, which will be the case for humidities
  ! below the efflorescence point)
  alpha = 1.0
  beta  = 1.0

  ! Specify Fitzgerald's alpha and beta parameters.
  !
  ! Note that in Eq. 8 of Fitzgerald:
  !    wr = alpha * r ** beta
  ! where:
  !   wr = wet radius (microns)
  !   r  = dry radius (microns)
  !
  ! If we take natural logarithms:
  !   ln (wr) = ln (alpha) + beta * ln (r)
  !
  ! In Fitzgerald's fomulation alpha was given in a form derived for particle
  ! radii in microns. However ln_r0 (i.e. natural logarithm of the radius) has
  ! been passed to this routine for radius in metres (1 metre = 1.0e6 microns).
  ! Some conversions are then needed. These result in keeping beta unchanged
  ! but increasing alpha by a factor which turns out to be:
  !    10.0**( 6.0*(beta-1) )
  ! See demonstration here:
  !
  ! Convert wet/dry radii from metres to microns and introduce them in Eq. 8:
  !   wr * 1.0e6 = alpha * (r * 1.0e6) ** beta
  !
  ! Rearranging:
  !   wr =  10.0 ** [6.0 * (beta -1)] * alpha * r ** beta
  !
  ! Finally, wr can be expressed as:
  !   wr    = 10.0 ** power * alpha * r ** beta    (Eq. 8.1)
  !   power = 6 * (beta -1)
  !
  ! Thus, as mentioned above, beta remains unchanged and alpha needs to
  ! be multiplied by the following factor in all the calculations below:
  !    10.0 ** [6 * (beta -1)]
  !
  ! Or if we use logarithms in Eq. 8.1:
  !   ln (wr) = ln (alpha) + beta * ln (r) + power * ln (10.0)

  ! Note that no growth takes place below the efflorescence point.
  if ( humidity_eff >= effl ) then

    ! To get the values of alpha and beta for RH < RH (deliq) we
    ! need to interpolate between the known values of alpha/beta
    ! at the efflorescence point (1, no growth) and those at the
    ! deliquescence point (alpha_deliq, beta_deliq).
    if ( humidity_eff < deliq ) then

      ! First, call the functions bfunc and afunc to get the values of
      ! beta and alpha at RH(deliq). Also take into account the unit
      ! conversion factor for alpha mentioned above.
      beta_deliq  = bfunc (deliq)
      power       = 6.0 * (beta_deliq - 1.0)
      alpha_deliq = afunc (aerosol_type, deliq) * (10.0**power)

      ! Now interpolate both beta and alpha
      beta  = 1.0 + (humidity_eff - effl) * (beta_deliq  - 1.0) / (deliq - effl)
      alpha = 1.0 + (humidity_eff - effl) * (alpha_deliq - 1.0) / (deliq - effl)

      ! Above deliquescence use Fitzgerald's formulation, which is valid
      ! up to RH = 0.995.
    else if ( humidity_eff  <=  0.995 ) then
      beta  = bfunc (humidity_eff)

      ! With the following 2 lines one could do all conversions affecting
      ! alpha when humidity is at or above RH(deliq):
      !   power = 6.0 * (beta - 1.0)
      !   alpha = afunc (aerosol_type, humidity_eff) * (10.0**power)
      !
      ! However that didn't preserve model results compared with
      ! the first implementation that followed the Edwards-Slingo
      ! Offline Radiation Code. As a consequence, the corrections
      ! are not applied here but later on after calculating the
      ! logarigthmic term.
      alpha = afunc (aerosol_type, humidity_eff)

    else
      write (umMessage,'(A8,1X,3E12.4)') 'Humidity', humidity_eff
      call umPrint(umMessage,src='grow_particles')

      ierr     = 2
      cmessage = 'humidity_eff used in Fitzgerald formulation ' //             &
                 'should be lower than 0.995'
      call ereport ('GROW_PARTICLES', ierr, cmessage)
    end if
  end if

  ! Increase natural logarithms of radius and standard deviation.
  ln_wr0        = beta * ln_r0 + log (alpha)
  ln_wsigma     = beta * ln_sigma

  ! Now do the unit conversion which affects alpha when humidity
  ! is at or above RH(deliq).  We do not make this correction
  ! below RH(deliq), because values of alpha already corrected
  ! were used in the interpolation.
  if ( humidity_eff >= deliq) then
    power         = 6.0 * (beta - 1.0)
    ln_wr0        = ln_wr0 + power * log (10.0)
  end if

  ! ----------------------------------------------------------------------------
  ! Biomass-burning aerosol: Use a pre-defined array of growth factors.
  ! We suppose that the standard deviation of the log-normal distrib does
  ! not vary with humidity (i.e. all particle sizes grow at the same rate).
  !
  ! Fossil-fuel organic carbon aerosols also use the same growth factors
  ! as biomass burning aerosols.
  !
else if ( ( aerosol_type == ip_biomass_1  ) .or.  &  ! biomass fresh
          ( aerosol_type == ip_biomass_2  ) .or.  &  ! biomass aged
          ( aerosol_type == ip_ocff_fresh ) .or.                               &
          ( aerosol_type == ip_ocff_aged  ) ) then

  ! Not using Fitzgerald: cap humidity to max value of 100%
  humidity_eff = min (humidity, 1.0)

  ! Verify that 0.0 <= humidity_eff <= 1.0.
  if ( humidity_eff < 0.0 .or. humidity_eff > 1.0) then
    write (umMessage,'(A8,1X,3E12.4)') 'Humidity', humidity_eff
    call umPrint(umMessage,src='grow_particles')

    ierr     = 3
    cmessage = 'Humidity has an unexpected value'

    ! Get the index in the growth factor array and calculate
    ! growth of the aerosol radius.
  else
    i_hum = nint (humidity_eff / humidity_step) + 1

    ! biomass fresh or OCFF fresh
    if (aerosol_type == ip_biomass_1 .or.                                      &
        aerosol_type == ip_ocff_fresh ) then
      growth_factor = biom_gf1 (i_hum)
      ! biomass aged or OCFF aged
    else
      growth_factor = biom_gf2 (i_hum)
    end if

    ! The line below is the same as: ln(wr0) = ln (r0 * growth_factor)
    ln_wr0    = log ( exp (ln_r0) * growth_factor )
    ln_wsigma = ln_sigma
  end if


  ! Biogenic aerosol. Pre-defined growth factors are used.
else if ( aerosol_type == ip_biogenic ) then

  ! Not using Fitzgerald: cap humidity to max value of 100%
  humidity_eff = min (humidity, 1.0)

  ! Verify that 0.0 <= humidity_eff <= 1.0.
  if ( humidity_eff < 0.0 .or. humidity_eff > 1.0) then
    write (umMessage,'(A8,1X,3E12.4)') 'Humidity', humidity_eff
    call umPrint(umMessage,src='grow_particles')

    ierr     = 3
    cmessage = 'Humidity has an unexpected value'

    ! Get the index in the growth factor array and calculate
    ! growth of the aerosol radius.
  else
    i_hum         = nint (humidity_eff / humidity_step) + 1
    growth_factor = biogenic_gf (i_hum)
    ln_wr0        = log ( exp (ln_r0) * growth_factor )
    ln_wsigma     = ln_sigma
  end if

else
  ! Trap incorrect aerosol types here.
  ierr     = 4
  cmessage = 'Unavailable aerosol composition type'
end if

if (ierr > 0) then
  call ereport ('GROW_PARTICLES', ierr, cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

return

end subroutine grow_particles


!----------------------------------------------------------------------
! Function to calculate the alpha parameter of Fitzgerald (1975).
! Only valid for RH at the deliquescence point and above up to 0.995.
!----------------------------------------------------------------------
real(kind=real_umphys) function afunc (aerosol_type, h)

implicit none

! Function arguments
integer, intent(in) :: aerosol_type

real(kind=real_umphys),    intent(in) :: h
                                     ! relative humidity (0-0.995)

! Local variables
real(kind=real_umphys) :: phi           ! Fitzgerald's phi parameter

integer                            :: ierr      ! error flag
character (len=errormessagelength) :: cmessage  ! error message

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='AFUNC'

!
if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! According to Eq. 10 of Fitzgerald the calculations below are
! valid for ammonium sulphate.
!
! First calculate the phi parameter depending on RH fraction
! as in Eq. 10 of Fitzgerald.
if ( h  <=  0.97 ) then
  phi = 1.058
else if ( h  <=  0.995 ) then
  phi = 1.058 - 0.0155 * (h - 0.97) / (1.02 - h**1.4)
end if
!
! Now calculate alpha, also following Eq. 10 of Fitzgerald.
afunc = 1.2 * exp ( 0.066 * h / (phi - h) )

! For sea-salt (ammonium nitrate), alpha is 1.35 (1.06) times
! its value for ammonium sulphate.
if ( ( aerosol_type  ==  ip_seasalt_film )  .or.                               &
     ( aerosol_type  ==  ip_seasalt_jet ) ) then
  afunc = afunc * 1.35

else if ( aerosol_type  ==  ip_nitrate ) then
  afunc = afunc * 1.06

  ! Trap incorrect aerosol types here.
else if ( ( aerosol_type  /=  ip_accum_sulphate )  .and.                       &
          ( aerosol_type  /=  ip_aitken_sulphate ) ) then
  ierr     = 1
  cmessage = 'Unavailable aerosol composition type for Fitzgerald formulation'
  call ereport ('GROW_PARTICLES:AFUNC', ierr, cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

return

end function afunc

!----------------------------------------------------------------------
! Function to calculate the beta parameter of Fitzgerald (1975).
! Only valid for RH at the deliquescence point and above up to 0.995.
!----------------------------------------------------------------------

real(kind=real_umphys) function bfunc (h)

implicit none

! Function arguments
real(kind=real_umphys), intent(in) :: h  ! relative humidity (0-0.995)

! Local variables
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BFUNC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Calculation of beta following Eq. 9 of Fitzgerald.
bfunc = exp ( 7.7e-4 * h / (1.009 - h) )

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

return

end function bfunc


end module grow_particles_mod
