! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change phase for points where no downdraught occurring.
!
module chg_phse_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'CHG_PHSE_MOD'
contains

subroutine chg_phse (npnts,rain,snow,dthbydt_km1,                              &
                     exk,exkm1,delpkm1,the_k,the_km1,timestep,                 &
                     cca, wtrac_conv_old)


use cv_run_mod, only: t_melt_snow, pr_melt_frz_opt
use cv_param_mod, only:                                                        &
     precip_area_fac, t_frez_rain, pr_melt_tdep, pr_melt_frz_tdep

use planet_constants_mod, only: g, cp
use water_constants_mod,  only: lf, tm
use wtrac_conv_mod,       only: l_wtrac_conv
use wtrac_conv_store_mod, only: conv_old_wtrac_type

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ------------------------------------------------------------------------------
! Description:
! Change phase for points where no downdraught occurring. (The code of this
! routine is almost identical to crs_frzl. At some future date it would be
! better to have one common routine if possible.)
! Updates potential temperature of layer k as precipitation changes phase
! in situ.
! Add latent heating where precipitation crosses a melting or freezing level.
!
! See UM Documentation paper 27.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                  ! Number of points

real(kind=real_umphys), intent(in) ::                                          &
  timestep               ! Timestep

real(kind=real_umphys), intent(in) ::                                          &
  exk(npnts)        & ! Exner ratio for layer k
 ,exkm1(npnts)      & ! Exner ratio for layer k-1
 ,delpkm1(npnts)    & ! Pressure difference across layer k-1 (Pa)
 ,the_k(npnts)      & ! potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! potential temperature of environment in layer k-1 (K)
 ,cca(npnts)          ! Convective cloud amount


!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
real(kind=real_umphys), intent(in out) ::                                      &
  rain(npnts)          & ! in Amount of falling rain (kg/m**2/s)
                         ! out Updated amount of falling rain (kg/m**2/s)
 ,snow(npnts)          & ! in Amount of falling Snowfall (kg/m**2/s)
                         ! out Updated amount of falling snow (kg/m**2/s)
 ,dthbydt_km1(npnts)     ! in Increment to model potential temperature in
                         !    layer k-1
                         ! out Updated increment to model potential
                         !     temperature in layer k-1 due to change of phase

type(conv_old_wtrac_type), optional, intent(in out) :: wtrac_conv_old
                         ! Store of amount of precip freezing/melting
                         ! optional as only used when called from 6A routines

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i                 ! loop counter

logical ::                                                                     &
  bppnwt_k       & ! Mask where precipitation is liquid in layer k
 ,bppnwt_km1       ! Mask where precipitation is liquid in layer k-1

real(kind=real_umphys) ::                                                      &
  factor          &   ! Used in the calculation of change of phase of falling
                      ! precipitation.
 ,precip_limit    &   ! Limit on amount of precip melting/freezing to avoid
                      ! a phase change
 ,precip_area     &   ! Area fraction of precip = cca * precip_area_fac
 ,the_km1_new         ! the_km1 updated with all increments prior to this
                      ! subroutine
real(kind=real_umphys) ::                                                      &
  precip_frez(npnts)     &   ! Amount of precipitation freezing
                             ! (This field is stored for water tracer use)
 ,precip_melt(npnts)         ! Amount of precipitation melting
                             ! (This field is stored for water tracer use)
real(kind=real_umphys) :: melt_tscale
                      ! e-folding temperature thickness for snow melt
real(kind=real_umphys) :: r_melt_tscale2
                      ! reciprocal squared of e-folding temperature thickness
                      ! for snow melt
real(kind=real_umphys) :: frez_tscale
                      ! e-folding temperature thickness for freezing rain
real(kind=real_umphys) :: r_frez_tscale2
                      ! reciprocal squared of e-folding temperature thickness
                      ! for freezing rain

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CHG_PHSE'
!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! ----------------------------------------------------------------------
!   Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Section (11), equation (42)
! ----------------------------------------------------------------------
if (pr_melt_frz_opt == pr_melt_tdep .or.                                       &
    pr_melt_frz_opt == pr_melt_frz_tdep) then
  ! Melt snow at at rate proportional to the temperature
  ! above the melting temperature (usually 0C) and (if switched on) freeze
  ! rain below the melting temperature

  !Melting scale
  melt_tscale   = t_melt_snow - tm
  r_melt_tscale2= 1.0/melt_tscale**2

  !Freezing scale
  frez_tscale   = t_frez_rain - tm
  r_frez_tscale2= 1.0/frez_tscale**2

  do i=1,npnts

    precip_melt(i) = 0.0
    precip_frez(i) = 0.0

    the_km1_new = the_km1(i)  + timestep*dthbydt_km1(i)
    bppnwt_km1  = the_km1_new >   tm/exkm1(i)
    factor      = lf*g/(exkm1(i)*cp*delpkm1(i))
    precip_area = precip_area_fac*cca(i)

    if ( bppnwt_km1 .and. (snow(i) > 0.0) ) then

      if ( (3.0*melt_tscale) > (the_km1_new*exkm1(i)-tm) ) then
        ! If temperature still below emergency melting level
        !then calculate melting rate
        precip_melt(i)   = snow(i)*(the_km1_new*exkm1(i)-tm)                   &
                      *(the_km1_new*exkm1(i) - the_k(i)*exk(i))*r_melt_tscale2
      else
        ! If above emergency melting level melt all snow
        precip_melt(i)   = snow(i)
      end if

      ! Check that precip_melt is not negative.
      ! to prevent refreezing of rain
      precip_melt(i)     = max(0.0,precip_melt(i))

      ! Limit the amount of melting so that is cannot
      ! lower the temperature to below tm.
      precip_limit    = precip_area*(the_km1_new - tm/exkm1(i)) /              &
                                                        (timestep * factor)
      precip_melt(i)     = min(precip_limit,precip_melt(i))

      ! Limit the melting to the amount of snow
      precip_melt(i)     = min(snow(i),precip_melt(i))

      !Calculate temperature increment and update rain and snow
      dthbydt_km1(i)  = dthbydt_km1(i) - precip_melt(i) * factor
      rain(i)         = rain(i)+precip_melt(i)
      snow(i)         = snow(i)-precip_melt(i)

    end if  !Test on Melting

    if ( pr_melt_frz_opt == pr_melt_frz_tdep .and.                             &
         (.not. bppnwt_km1 .and. (rain(i) > 0.0))) then
      ! Calculate the freezing rate. Note no emergency freeze.
      precip_frez(i)     = rain(i)*(tm-the_km1_new*exkm1(i))                   &
                      *(the_km1_new*exkm1(i) - the_k(i)*exk(i))*r_frez_tscale2

      ! Check that precip_frez is not negative.
      ! to prevent melting of snow
      precip_frez(i)     = max(0.0,precip_frez(i))

      ! Limit the amount of melting so that it cannot
      ! raise the temperature to above tm.
      precip_limit    = precip_area*(tm/exkm1(i) - the_km1_new) /              &
                                                        (timestep * factor)
      precip_frez(i)     = min(precip_limit,precip_frez(i))

      ! Limit the freezing to the amount of rain
      precip_frez(i)     = min(rain(i),precip_frez(i))

      !Calculate temperature increment and update rain and snow
      dthbydt_km1(i)  = dthbydt_km1(i) + precip_frez(i) * factor
      rain(i)         = rain(i)-precip_frez(i)
      snow(i)         = snow(i)+precip_frez(i)

    end if  !Test on freezing

  end do

else
  ! original code melts all snow at freezing level and immediately
  ! freezes any supercooled rain.
  do i=1,npnts

    precip_melt(i) = 0.0
    precip_frez(i) = 0.0

    the_km1_new = the_km1(i) + timestep*dthbydt_km1(i)
    bppnwt_k    = the_k(i) >  tm/exk(i)
    bppnwt_km1  = the_km1_new  >   tm/exkm1(i)
    factor      = lf*g/(exkm1(i)*cp*delpkm1(i))
    precip_area = precip_area_fac*cca(i)

    if (.not. bppnwt_km1 .and. (bppnwt_k .or. rain(i) >  0.0)) then
      ! Freeze
      precip_frez(i) = min( rain(i), precip_area*(tm/exkm1(i) - the_km1_new)/  &
                                                         (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) + precip_frez(i) * factor
      snow(i) = snow(i) + precip_frez(i)
      rain(i) = rain(i) - precip_frez(i)
    end if


    if (bppnwt_km1 .and. (.not. bppnwt_k .or. snow(i) >  0.0)) then
      ! Melt
      precip_melt(i) = min( snow(i),                                           &
               precip_area * (the_km1_new - tm/exkm1(i)) /                     &
                  (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) - precip_melt(i) * factor
      rain(i) = rain(i) + precip_melt(i)
      snow(i) = snow(i) - precip_melt(i)

    end if ! test on melting /freezing
  end do

end if

! If water tracers, save amount of freezing or melting
! (Note, freeze/melt fields are combined for efficiency only - there is either
!  freeze or melt at each point, not both.)
if (l_wtrac_conv) then
  do i=1,npnts
    wtrac_conv_old%precip_frez(i) = precip_frez(i) - precip_melt(i)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine chg_phse

end module chg_phse_mod
