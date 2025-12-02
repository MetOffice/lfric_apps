! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates constants used in large-scale precipitation scheme.
module lspcon_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LSPCON_MOD'

contains

subroutine lspcon()


  ! General modules
use stochastic_physics_run_mod, only: l_rp2, i_rp_scheme, i_rp2b,              &
                                      x1r_rp, m_ci_rp, rp_idx
use conversions_mod,            only: pi
use planet_constants_mod,       only: g, repsilon, r
use water_constants_mod,        only: rho_water, lc, lf
use missing_data_mod,           only: rmdi

!   Microphysics modules
use lsp_dif_mod,         only: air_density0, air_viscosity0,                   &
                               air_conductivity0, air_diffusivity0,            &
                               air_pressure0, apb1, apb2, apb3, apb4,          &
                               apb5, apb6, sc, vent_ice1, vent_ice2,           &
                               vent_rain1, vent_rain2

use mphys_psd_mod,       only: l_calcfall, cr, dr, c1r, d1r, h1r,              &
                               c2r, d2r, h2r, x2i, x3i, x4i, ri, si,           &
                               ci0, di0, x2ic, x3ic, x4ic, ric,                &
                               sic, cic0, dic0, x1g, x2g, x4g,                 &
                               ag, bg, cg, dg, x1gl, x2gl, x4gl,               &
                               x1gf, x2gf, x4gf

use mphys_constants_mod, only: cx, constp, m_ci_sav, rho_q_veloc,              &
                               l_calc_mp_const,                                &
                               x4r, x1i, x1ic, lsp_ei, lsp_fi, lsp_eic,        &
                               lsp_fic, n0_murk, m0_murk, lam_evap_enh,        &
                               qclmin_rime, area_ratio_prefac,                 &
                               area_ratio_expn, timestep_mp

use mphys_inputs_mod,    only: x1r, x2r, l_psd, ai, bi, ar, arc,               &
                               l_diff_icevt, cic_input, dic_input,             &
                               ci_input, di_input, l_mcr_qrain,                &
                               l_mcr_qgraup,                                   &
                               i_mcr_iter, i_mcr_iter_none,                    &
                               i_mcr_iter_niters, i_mcr_iter_tstep,            &
                               niters_mp, timestep_mp_in,                      &
                               l_shape_rime, qclrime, a_ratio_fac,             &
                               a_ratio_exp, z_surf, l_droplet_tpr,             &
                               l_autoconv_murk,                                &
                               l_mcr_arcl, l_use_sulphate_autoconv,            &
                               graupel_option, gr_field_psd,                   &
                               l_ice_shape_parameter

use timestep_mod,         only: timestep

use mphys_bypass_mod,     only: l_crystals, mphys_mod_top

use lsp_autoc_consts_mod, only: n0_clark, m0_clark,                            &
                                z_low_nd, eta_peak,                            &
                                eta_low_nd, level_peak, level_surf,            &
                                vala_fac1, vala_fac2, eta_before_taper
use level_heights_mod,    only: eta_theta_levels
use atm_fields_bounds_mod,only: tdims
use murk_inputs_mod,      only: l_murk, m0_murk_scale
use ukca_option_mod,      only: l_ukca_aie2
use glomap_clim_option_mod, only: l_glomap_clim_aie2
use easyaerosol_option_mod,only: l_easyaerosol_autoconv
use mphys_radar_mod,      only: rho_g, rho_i, rho_i2

!   Dr hook modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

use gammaf_mod,           only: gammaf

!   Printing and error reporting
use umPrintMgr, only: umprint, ummessage, newline, prstatus_oper
use um_parcore,           only: mype
use ereport_mod,          only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none

!  Description:
!     Calculates constants used within the LSP_ICE routine.

!  Method:
!     Calculate powers, gamma functions and constants from combinations
!     of physical parameters.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

!  Code description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! Declarations:

!  Local scalars:
integer :: k
! Counter to print out the contents of CX and CONSTP
real(kind=real_umphys) :: temp,                                                &
! Forms input to the GAMMAF routine which calculates gamma functions.
       g1, g2, g3,                                                             &
       gb1, gb2, gb3,                                                          &
       gbc1, gbd1, gbdc1,                                                      &
       gc1, gc2, gc3,                                                          &
       gd3, gdc3, gd52, gdc52,                                                 &
       gdr3, gdr4, gdr52,                                                      &
       gd1r3, gd2r3, gd1r4, gd2r4,                                             &
       gr2, gr4, gr5, gr6,                                                     &
       gg1, ggb1, ggbd1,                                                       &
       gg2, gdg52, gdg3, gg3, g7,                                              &
       gref1x4g, gref1x4ic, gref1x4i,                                          &
! Represents the gamma function of BI+DI+1 etc.
! Fall speed of ice particles parameters
       ci,di,cic,dic,aic,bic

! Height (m) below which to taper to surface droplet number
real(kind=real_umphys), parameter :: z_peak_nd = 150.0
                                      ! Removed from mphys_input_mod

! Eta value at which to taper to surface droplet number
real(kind=real_umphys) :: eta_surf

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSPCON'

integer :: errorstatus
character(len=errormessagelength) :: cmessage

!- End of header
!---------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First initialise cx and constp to rmdi where terms are
! not being used to prevent non-initialisation
cx(:)     = rmdi
constp(:) = rmdi

! Mass diameter relationship for ice crystals:  m(D) = AI D^BI
! Recommended values aic,bis for the generic particle size distribution
! from Brown and Francis are identical to ai,bi.  Copy here
aic = ai
bic = bi

!----------------------------------------------------------------------
! If the generic ice psd is being used and l_diff_icevt = .true.
! then different vt-diameter relations can be used for crystals
! and aggregates.
! Otherwise, a single fallspeed curve is used for all ice particles.
!----------------------------------------------------------------------
if ( .not. l_diff_icevt .or. .not. l_psd ) then
  !----------------------------------------------------------------------
  ! Use a single ice fallspeed curve determined from either a
  ! Re-X relation or a directly specified powerlaw
  !----------------------------------------------------------------------
  ! Do we need to calculate fall speeds?
  if (l_calcfall) then
    ! Define fall speeds
    ci=lsp_ei*air_viscosity0**(1.0-2.0*lsp_fi)                                 &
       *air_density0**(lsp_fi-1.0)                                             &
       *(2.0*g)**lsp_fi*(ai/ri)**lsp_fi
    di=lsp_fi*(bi+2.0-si)-1.0

    if (l_crystals) then
      cic=lsp_eic*air_viscosity0**(1.0-2.0*lsp_fic)                            &
         *air_density0**(lsp_fic-1.0)                                          &
         *(2.0*g)**lsp_fic*(aic/ric)**lsp_fic

      dic=lsp_fic*(bic+2.0-sic)-1.0

    else
      cic = 1.0
      dic = 1.0
    end if ! l_crystals

    ! Modify fallspeeds for random parameters 2
    if (l_rp2) then
      ci  = ci  * m_ci_rp(rp_idx)
      cic = cic * m_ci_rp(rp_idx)
    end if

  else ! l_calcfall
    ! Use preset parameters
    ci  = ci0
    di  = di0
    cic = cic0
    dic = dic0
  end if  ! Calculation of fall speeds (l_calcfall)

else ! not difficevt  / l_psd
  !----------------------------------------------------------------------
  ! Generic psd is used and crystals and aggregates can
  ! have different fallspeed relations
  !----------------------------------------------------------------------
  ! Fallspeed parameters for crystals
  cic = cic_input
  dic = dic_input
  ! Fallspeed parameters for aggregates
  ci = ci_input
  di = di_input

  ! Modify fallspeeds for random parameters 2
  if (l_rp2) then
    ci  = ci  * m_ci_rp(rp_idx)
    cic = cic * m_ci_rp(rp_idx)
  end if

end if   ! Use different fallspeeds with generic psd

if (l_shape_rime) then

   ! Set parameters to use shape-dependent
   ! riming rate with a threshold LWC
   !
  qclmin_rime = qclrime
  area_ratio_prefac = a_ratio_fac
  area_ratio_expn = a_ratio_exp

else

  ! Set parameters so that riming is
  ! not shape-dependent and all LWC
  ! is available to rime
  !
  qclmin_rime = 0.0
  area_ratio_prefac = 1.0
  area_ratio_expn = 0.0

end if

if ( graupel_option == gr_field_psd ) then
  ! Switch to new graupel PSD parameters
  x1g = x1gf
  x2g = x2gf
  x4g = x4gf
else
  ! Use original graupel PSD parameters
  x1g = x1gl
  x2g = x2gl
  x4g = x4gl
end if

! CX values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-62 are for rain. 63-80 are for graupel.
! 81-99 are for the generic ice size distribution.

! Crystals

if (l_crystals) then

  cx(1)=(bic+1.0+x4ic-x2ic)/bic
  cx(2)=-(x4ic+1.0-x2ic)/bic
  cx(3)=dic/(bic+1.0+x4ic-x2ic)
  cx(4)=(2.0+x4ic-x2ic)/(bic+1.0+x4ic-x2ic)
  cx(5)=(5.0+dic+2.0*x4ic-2.0*x2ic)*0.5/(bic+1.0+x4ic-x2ic)
  cx(6)=(3.0+dic+x4ic-x2ic)/(bic+1.0+x4ic-x2ic)
  cx(7)=1.0/(x2ic-x4ic-1.0-bic)
  cx(8)=1.0+x4ic
  cx(9)=2.0+x4ic
  cx(10)=3.0+x4ic
  cx(11)=x2ic
  cx(12)=x3ic
  cx(13)=1.0+x4ic+bic
  cx(14)=bic

end if

! Aggregates
cx(23)=di/(bi+1.0+x4i-x2i)
cx(24)=(2.0+x4i-x2i)/(bi+1.0+x4i-x2i)
cx(25)=(5.0+di+2.0*x4i-2.0*x2i)*0.5/(bi+1.0+x4i-x2i)
cx(26)=(3.0+di+x4i-x2i)/(bi+1.0+x4i-x2i)
cx(27)=1.0/(x2i-x4i-1.0-bi)
cx(28)=1.0+x4i
cx(29)=2.0+x4i
cx(30)=3.0+x4i
cx(31)=x2i
cx(32)=x3i
cx(33)=3.0+x4i+bi
cx(34)=2.0+x4i+bi
cx(35)=1.0+x4i+bi
! Rain
cx(41)=dr/(4.0+dr-x2r+x4r)
cx(42)=1.0/(4.0+dr-x2r+x4r)
cx(43)=6.0+x4r
cx(44)=5.0+x4r
cx(45)=4.0+x4r
cx(46)=x2r
cx(47)=2.0+x4r-x2r
cx(48)=1.0/(4.0+x4r+dr-x2r)
cx(49)=(dr+5.0)*0.5-x2r+x4r
cx(50)=(3.0+dr-x2r+x4r)/(4.0+dr-x2r+x4r)
! Rain mixing ratio
cx(51)=dr/(4.0-x2r+x4r)
cx(52)=1.0/(4.0-x2r+x4r)
cx(53)=(3.0+dr-x2r+x4r)/(4.0-x2r+x4r)
! Abel & Shipway Terms
cx(56)=h1r
cx(57)=h2r
cx(59)=4.0+d1r+x4r
cx(60)=4.0+d2r+x4r
cx(61)=3.0+d1r+x4r
cx(62)=3.0+d2r+x4r

!Note there is no space between rain and graupel.

! Graupel

if (l_mcr_qgraup) then
  cx(63)=dg/(bg+1.0+x4g-x2g)
  cx(64)=(2.0+x4g-x2g)/(bg+1.0+x4g-x2g)
  cx(65)=(5.0+dg+2.0*x4g-2.0*x2g)*0.5/(bg+1.0+x4g-x2g)
  cx(66)=(3.0+dg+x4g-x2g)/(bg+1.0+x4g-x2g)
  cx(67)=1.0/(x2g-x4g-1.0-bg)
  cx(68)=1.0+x4g
  cx(69)=2.0+x4g
  cx(70)=3.0+x4g
  cx(71)=x2g

  cx(73)=3.0+x4g+bg
  cx(74)=2.0+x4g+bg
  cx(75)=1.0+x4g+bg
end if


! Generic ice particle size distribution
cx(81)=2.0+di
cx(82)=di+bi
cx(83)=1.0+bi
cx(84)=1.0
cx(85)=1.0+0.5*(di+1.0)

! Radar reflectivity
cx(90) = 0.224 * ((6.0 * ( (ai/pi)/900.0) )**2)
cx(91) = 6.0/(pi * rho_water)
cx(92) = x2r - (1.0 + (2.0*3.0) + x4r )

if (l_mcr_qgraup) then
  cx(93) = 1.0 / (bg  + 1.0 + x4g  - x2g)
  cx(96) = -1.0 * (1.0 + x4g  + (2.0*bg)  - x2g)
end if

if (l_crystals) then
  cx(94) = 1.0 / (bic + 1.0 + x4ic - x2ic)
  cx(97) = -1.0 * (1.0 + x4ic + (2.0*bic) - x2ic)
end if

cx(95) = 1.0 / (bi  + 1.0 + x4i  - x2i)
cx(98) = -1.0 * (1.0 + x4i  + (2.0*bi)  - x2i)

!-----------------------------------------------------------------------
! Additional constants for generic psd crystals.
! Only used if crystals and aggregates can have different vt-D relations
!-----------------------------------------------------------------------
cx(181)=2.0+dic
cx(182)=dic+bi
cx(185)=1.0+0.5*(dic+1.0)

! Define gamma values
! Crystals
if (l_crystals) then
  temp=1.0+x4ic
  call gammaf(temp,gc1)
  temp=bic+1.0+x4ic
  call gammaf(temp,gbc1)
  temp=bic+dic+1.0+x4ic
  call gammaf(temp,gbdc1)
  temp=2.0+x4ic
  call gammaf(temp,gc2)
  temp=(dic+5.0+2.0*x4ic)*0.5
  call gammaf(temp,gdc52)
  temp=dic+3.0+x4ic
  call gammaf(temp,gdc3)
  temp=3.0+x4ic
  call gammaf(temp,gc3)
end if ! l_crystals

! Aggregates
temp=1.0+x4i
call gammaf(temp,g1)
temp=bi+1.0+x4i
call gammaf(temp,gb1)
temp=bi+2.0+x4i
call gammaf(temp,gb2)
temp=bi+3.0+x4i
call gammaf(temp,gb3)
temp=bi+di+1.0+x4i
call gammaf(temp,gbd1)
temp=2.0+x4i
call gammaf(temp,g2)
temp=(di+5.0+2.0*x4i)*0.5
call gammaf(temp,gd52)
temp=di+3.0+x4i
call gammaf(temp,gd3)
temp=3.0+x4i
call gammaf(temp,g3)
! Rain
temp=dr+4.0+x4r
call gammaf(temp,gdr4)
temp=4.0+x4r
call gammaf(temp,gr4)
temp=6.0+x4r
call gammaf(temp,gr6)
temp=5.0+x4r
call gammaf(temp,gr5)
temp=2.0+x4r
call gammaf(temp,gr2)
temp=(dr+5.0+2.0*x4r)*0.5
call gammaf(temp,gdr52)
temp=dr+3.0+x4r
call gammaf(temp,gdr3)
!Abel & Shipway Rain Terms
temp= 3.0+d1r+x4r
call gammaf(temp,gd1r3)
temp= 3.0+d2r+x4r
call gammaf(temp,gd2r3)
temp= 4.0+d1r+x4r
call gammaf(temp,gd1r4)
temp= 4.0+d2r+x4r
call gammaf(temp,gd2r4)

! Graupel
if (l_mcr_qgraup) then
  temp=1.0+x4g
  call gammaf(temp,gg1)
  temp=bg+1.0+x4g
  call gammaf(temp,ggb1)
  temp=bg+dg+1.0+x4g
  call gammaf(temp,ggbd1)
  temp=2.0+x4g
  call gammaf(temp,gg2)
  temp=(dg+5.0+2.0*x4g)*0.5
  call gammaf(temp,gdg52)
  temp=dg+3.0+x4g
  call gammaf(temp,gdg3)
  temp=3.0+x4g
  call gammaf(temp,gg3)
end if ! l_mcr_qgraup

! Radar reflectivity
temp = 1.0 + (2.0 * 3.0) +x4r
call gammaf(temp, g7)
! called g7 as it will be gamma(7) unless x4r not equal to 0

if (l_mcr_qgraup) then
  temp = 1.0 + x4g + (2.0 * bg)
  call gammaf(temp, gref1x4g)
end if

if (l_crystals) then
  temp = 1.0 + x4ic + (2.0 * bic)
  call gammaf(temp, gref1x4ic)
end if

temp = 1.0 + x4i  + (2.0 * bi)
call gammaf(temp, gref1x4i)

! CONSTP values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-60 are for rain. 61-80 are for graupel.
! 81-99 are for the generic ice size distribution

! First update value of x1r if random parameters 2b
! scheme is switched on
if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
  x1r = x1r_rp(rp_idx)
end if

! Crystals
if (l_crystals) then
  constp(1)=x1ic
  constp(2)=1.0/gc1
  constp(3)=1.0/(aic*gbc1)
  constp(4)=cic*gbdc1/gbc1
  constp(5)=1.0/(aic*x1ic*gbc1)
  constp(6)=2.0*pi*x1ic
  constp(7)=vent_ice1*gc2
  constp(8)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                        &
            *sqrt(cic)*gdc52
  constp(9)=pi*0.25*x1ic*cic*gdc3
  constp(10)=5.0*gc1
  constp(11)=2.0*gc2
  constp(12)=0.25*gc3
  constp(13)=pi**2*rho_water*x1ic*x1r
  constp(14)=2.0*pi*air_conductivity0/lf*x1ic
  ! Capacitance relative to spheres of same maximum dimension
  ! Formula depends on value of axial ratio
  constp(15)=arc
  if (arc  >   1.0) then
    ! Prolate
    constp(15)=(1.0-(1.0/constp(15))**2)**0.5 /                                &
               log( constp(15) + (constp(15)**2-1.0)**0.5 )
  else if (arc  ==  1.0) then
    ! Spherical
    constp(15)=1.0
  else
    ! Oblate
    constp(15)=(1.0-constp(15)**2)**0.5                                        &
               /asin((1.0-constp(15)**2)**0.5)
  end if
  ! Now adjust diffusional growth constants for capacitance
  constp(6)=constp(6)*constp(15)
  constp(14)=constp(14)*constp(15)
  !            constp(20) is reserved for lsp_collection

end if ! l_crystals

! Values 16 to 23 are unused
! Aggregates

constp(24)=ci*gbd1/gb1
constp(25)=1.0/(ai*x1i*gb1)
constp(26)=2.0*pi*x1i
constp(27)=vent_ice1*g2
constp(28)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                         &
           *sqrt(ci)*gd52
constp(29)=pi*0.25*x1i*ci*gd3
constp(30)=5.0*g1
constp(31)=2.0*g2
constp(32)=0.25*g3
constp(33)=pi**2*rho_water*x1i*x1r
constp(34)=2.0*pi*air_conductivity0/lf*x1i
! Capacitance relative to spheres of same maximum dimension
! Formula depends on value of axial ratio
if (l_ice_shape_parameter) then
  constp(35)=0.5 !ar won't be used
else
  constp(35)=ar
  if (ar  >   1.0) then
      ! Prolate
    constp(35)=(1.0-(1.0/constp(35))**2)**0.5 /                                &
        log( constp(35) + (constp(35)**2-1.0)**0.5 )
  else if (ar  ==  1.0) then
      ! Spherical
    constp(35)=1.0
  else
      ! Oblate
    constp(35)=(1.0-constp(35)**2)**0.5                                        &
        / asin((1.0-constp(35)**2)**0.5)
  end if
end if
! Now adjust diffusional growth constants for capacitance
constp(26)=constp(26)*constp(35)
constp(34)=constp(34)*constp(35)

if ( l_crystals ) then
  constp(36) = gc1*gb3
  constp(37) = 2.0*gc2*gb2
  constp(38) = gc3*gb1
  constp(39) = ai*x1ic*x1i*pi/4.0
end if

!            constp(40) is reserved for lsp_collection
! Rain
constp(41)=6.0*cr*gdr4/gr4
constp(42)=pi*rho_water/6.0*x1r*gdr4*cr
constp(43)=1.0/120.0*gr6
constp(44)=1.0/24.0*gr5
constp(45)=1.0/6.0*gr4
constp(46)=2.0*pi*x1r
constp(47)=vent_rain1*gr2
constp(48)=vent_rain2*sc**(1.0/3.0)/air_viscosity0**0.5                        &
           *gdr52*sqrt(cr)
constp(49)=pi*0.25*x1r*cr*gdr3
constp(50)=1.0/(pi*rho_water*x1r*gr4/6.0)

! Abel & Shipway Rain Section
!-----------------------------------------------------------
!First overwrite lam_evap_enh if maximum rain rate specified
!-----------------------------------------------------------

constp(51)=c1r*x1r*pi*0.25*gd1r3
constp(52)=c2r*x1r*pi*0.25*gd2r3
constp(53)=gr4
constp(54)=c1r*gd1r4
constp(55)=c2r*gd2r4
constp(56)=(pi/6.0)*rho_water*x1r*(lam_evap_enh**x2r)
constp(57)=(pi/6.0)*rho_water*x1r
constp(58)=20.0*pi*pi*rho_water*x1r

! Values 59 to 60 are unused

! Graupel

if (l_mcr_qgraup) then
  constp(64) = cg*ggbd1/ggb1
  constp(65) = 1.0/(ag*x1g*ggb1)
  constp(67) = vent_rain1*gg2
  constp(68) = vent_rain2*sc**(1.0/3.0)/air_viscosity0**0.5                    &
               *sqrt(cg)*gdg52
  constp(69) = pi*0.25*x1g*cg*gdg3
  constp(74) = 2.0*pi*air_conductivity0*x1g/lf
  constp(76) = gg1*gb3
  constp(77) = 2.0*gg2*gb2
  constp(78) = gg3*gb1
  constp(79) = ai*x1g*x1i*pi/4.0
  constp(80) = pi*pi/24.0*x1g*ag

end if ! l_mcr_qgraup

! Generic ice particle size distribution
constp(81)=pi*0.25*ci
constp(82)=ci*ai
constp(83)=2.0*pi*constp(35)
constp(84)=vent_ice1
constp(85)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                         &
           *sqrt(ci)
constp(86)=pi*pi/24.0*rho_water*x1r
constp(87)=gr6
constp(88)=2.0*gr5
constp(89)=gr4
constp(90)=2.0*pi*constp(35)*air_conductivity0/lf
constp(91)=gb3
constp(92)=2.0*gb2
constp(93)=gb1

! Radar reflectivity
if (l_mcr_qgraup) constp(100) = ag  * x1g  * ggb1
if (l_crystals)   constp(101) = aic * x1ic * gbc1

constp(102) = ai  * x1i  * gb1

if (l_mcr_qgraup) then
  constp(103) = x1g  * (ag**2)  * gref1x4g  * ( 6.0 /(pi * rho_g))**2
end if

if (l_crystals) then
  constp(104) = x1ic * (aic**2) * gref1x4ic * ( 6.0 /(pi * rho_i2))**2
end if

constp(105) = x1i  * (ai**2)  * gref1x4i  * ( 6.0 /(pi * rho_i))**2

constp(106) = 1.0 / constp(50)
constp(107) = g7 * x1r


!-----------------------------------------------------------------------
! Additional constants for generic psd crystals.
! Only used if crystals and aggregates can have different vt-D relations
!-----------------------------------------------------------------------
constp(181)=pi*0.25*cic
constp(182)=cic*ai
constp(185)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                        &
           *sqrt(cic)

!-----------------------------------------------------------------------
! Determine Abel & Shipway Prognostic Fall Velocity at the point where
! it diverges from the traditional UM parametrization
! Note: already have lambda from module
!-----------------------------------------------------------------------
! Determine rhoqv = droplet velocity multiplied by any given rho and q.
! rho_q_veloc is a constant and
! droplet velocity = rho_q_veloc * air density correction / (rho q)

if ( .not. l_mcr_qrain ) then

  rho_q_veloc = constp(56) *  (                                                &
  ( constp(54) / ((lam_evap_enh+cx(56))**cx(59) ) ) +                          &
  ( constp(55) / ((lam_evap_enh+cx(57))**cx(60) ) )    )

end if

! Use Clark et al (2008, QJRMS) constants for Murk coupling.
! Haywood et al (2008, QJRMS) scheme has now been retired.
n0_murk = n0_clark
m0_murk = m0_murk_scale * m0_clark

!-----------------------------------------------------------------
! Values of diffusional growth parameters
!-----------------------------------------------------------------
! Terms in deposition and sublimation
apb1=(lc+lf)**2 * repsilon /(r*air_conductivity0)
apb2=(lc+lf) / air_conductivity0
apb3=r/(repsilon*air_pressure0*air_diffusivity0)
! Terms in evap of melting snow and rain
apb4=lc**2*repsilon/(r*air_conductivity0)
apb5=lc /air_conductivity0
apb6=r/(repsilon*air_pressure0*air_diffusivity0)

!-----------------------------------------------------------------
! Set up microphysics iterations
! These is only calculated at the beginning of the run (or CRUN)
!-----------------------------------------------------------------
if (l_calc_mp_const) then

  select case (i_mcr_iter)

  case ( i_mcr_iter_none )
    !     If multiple iterations are not selected, then ensure the
    !     number of all iterations is always exactly 1.
    niters_mp = 1
    timestep_mp = timestep

  case ( i_mcr_iter_niters )
    !     niters_mp is set in the namelist so no need to calculate
    timestep_mp = timestep/niters_mp
    write(umMessage,fmt='(A,I0,A,F14.2,A)')                                    &
         'Microphysics is using ', niters_mp, ' iterations:'                   &
         // newline // ' timestep_mp = ', timestep_mp, ' seconds'
    call umPrint(umMessage,pe=0,level=prstatus_oper)

  case ( i_mcr_iter_tstep )
    !     Calculate niters_mp from user requested microsphysics timestep
    !     Note we use the nearest integer so that this can round up or down
    if (timestep_mp_in < timestep) then
      niters_mp = nint(timestep/timestep_mp_in)
      timestep_mp = timestep/niters_mp
      write(umMessage,fmt='(A,I0,A,I0,A,F14.2,A)')                             &
           'User has requested a ', timestep_mp_in,                            &
           ' second microphysics iterative timestep ' // newline //            &
           'Using ', niters_mp, ' iterations, such that'   //                  &
           ' timestep_mp = ', timestep_mp, ' seconds'
      call umPrint(umMessage,pe=0,level=prstatus_oper)
    else
      !       If requested microsphysics timestep is longer than the model
      !       timestep then run with a single iteration and warn the user
      !       Also warn the user that this is what we have done
      niters_mp = 1
      timestep_mp = timestep
      if (mype == 0) then
        errorstatus = -100
        write(cmessage,fmt='(A,I0,A,F14.2,A,F14.2,A)')                         &
           'User has requested a ', timestep_mp_in,                            &
           ' second microphysics iterative '                                   &
           // newline // ' timestep ' // newline //                            &
           'This is longer than the model timestep (',                         &
           timestep, ' seconds) ' // newline //                                &
           'Running with a single microphysics iteration such that'            &
           // newline // ' timestep_mp = ', timestep_mp, ' seconds'
        call ereport('lspcon', errorstatus, cmessage)
      end if ! mype == 0
    end if ! timestep_mp_in < timestep

  end select

  ! If droplet tapering is on, pre-calculate the height of the taper
  ! and some constants related to it
  if (l_droplet_tpr) then

    eta_peak    = z_peak_nd / mphys_mod_top
    eta_low_nd  = z_low_nd  / mphys_mod_top
    eta_surf    = z_surf    / mphys_mod_top

    ! First find index of level of eta_peak
    do k = 1, tdims%k_end
      if (eta_theta_levels(k) >= eta_peak) then
        level_peak = k
        eta_before_taper = eta_theta_levels(k)
        exit
      end if
    end do
    ! Find index of level of eta_surf
    ! start loop from 2 so that level 1 is always found as a minimum
    do k = 2, level_peak
      if (eta_theta_levels(k) >= eta_surf) then
        level_surf = k-1
        exit
      end if
    end do

    ! Calculate constants for drop taper curve
    ! One of these is probably a bug, but it is not entirely clear which
    ! one. Since both are used in operational models, we maintain the
    ! behaviour of both. The difference is very small.
    if ((l_murk .and. l_autoconv_murk) .or.                                    &
         l_mcr_arcl .or.                                                       &
         (l_use_sulphate_autoconv .or. l_ukca_aie2 .or.                        &
          l_glomap_clim_aie2 .or. l_easyaerosol_autoconv) ) then
      vala_fac1 = log( eta_before_taper / eta_theta_levels( 1 ) )
      vala_fac1 = 1.0 / vala_fac1
      vala_fac2 = log( eta_before_taper / eta_theta_levels(level_surf) )
      vala_fac2 = 1.0 / vala_fac2
    else
      vala_fac1 = log( eta_peak / eta_theta_levels( 1 ) )
      vala_fac1 = 1.0 / vala_fac1
      vala_fac2 = log( eta_peak / eta_theta_levels(level_surf) )
      ! need to do vala_fac2 = 1.0 / vala_fac2 in lsp_taper_ndrop
    end if

  end if ! l_droplet_tpr

end if ! l_calc_mp_const

!-----------------------------------------------------------------
! Constants are now calculated, no need to recalculate
l_calc_mp_const = .false.

!-----------------------------------------------------------------
! Update m_ci_sav to the new value defined by m_ci_rp

m_ci_sav = m_ci_rp(rp_idx)
!-----------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
! End the subroutine
end subroutine lspcon

end module lspcon_mod
