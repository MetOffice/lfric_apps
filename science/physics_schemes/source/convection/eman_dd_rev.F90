! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module eman_dd_rev_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='EMAN_DD_REV_MOD'

contains

!  Emanuel Downdraught code - revised version

subroutine eman_dd_rev(n_dp,kmax_term,nlev,trlev,ntra,kterm,l_tracer           &
,                      exner_layer_centres                                     &
,                      p, ph, scale_factor                                     &
,                      th, q, qse, tracer, precip                              &
,                      dthbydt, dqbydt, dtrabydt                               &
,                      rain, snow ,down_flux,dt_dd,dqbydt_dd          )

use cv_run_mod, only:                                                          &
  t_melt_snow

use cv_derived_constants_mod, only:                                            &
  ls, lsrcp, lcrcp, lfrcp

use planet_constants_mod, only: r, cp, c_virtual, g
use water_constants_mod,  only: lc, tm

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Description ;
!   See Convection Documentation paper 27 which provides a reference to the
!   paper by Emanuel describing his convection scheme.
!   Calculation of unsaturated Down draughts using Emanuel code.
!   Note at present the convectin code deals in specific humidity
!   etc but the Emanuel code says it requires mixing ratios.
!   William Ingram appears to have taken the 4.3 version of the Emanuel
!   code and imported this into the old dynamics.
!   William Ingram chose to set various constants in a way consistent
!   with the UM. These choices mean that Emanuel's original choice
!   of a varying latent heat lv=lc+(cpv-cl)*(T-273.15)is reduced
!   to lv=lc.
!  Revised version to treat rain and snow separately and allow a mix
!  over a range of environmental temperatures.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards 8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  n_dp                 & ! No. of deep convection points
 ,nlev                 & ! No. of model layers
 ,trlev                & ! No. of tracer levels
 ,ntra                 & ! No. of tracer fields
 ,kmax_term            & ! highest level reached by convection
 ,kterm(n_dp)            ! level reached by convection

logical ,intent(in) ::                                                         &
  l_tracer               ! true in tracers present

real(kind=real_umphys) ,intent(in)    ::                                       &
  exner_layer_centres(n_dp,0:nlev)    & ! Exner
 ,p(n_dp,0:nlev)                      & ! Pressure at theta levels (Pa)
 ,ph(n_dp,0:nlev)                     & ! Pressure at half level (Pa)
 ,scale_factor(n_dp)                  & ! factor to scale sigd by. Comes from
                                        ! updraught closure.
 ,qse(n_dp,nlev)                      & ! Saturation specific humidity of cloud
                                        ! environment (kg/kg)
 ,q(n_dp,nlev)                        & ! specific humidity of cloud
                                        ! environment(kg/kg)
 ,th(n_dp,nlev)                       & ! theta of cloud environment(K)
 ,precip(n_dp,nlev)                   & ! Precip from updraught G-R scheme
                                        ! (kg/m2/s). Not the same units as
                                        !  wdtrain.
 ,tracer(n_dp,trlev,ntra)               !  Tracer on model levels  (kg/kg)


real(kind=real_umphys), intent(in out)    ::                                   &
  rain(n_dp)              & ! rainfall at surface (kg/m**2/s)
 ,snow(n_dp)                ! snowfall at surface (kg/m**2/s)

! increments
real(kind=real_umphys), intent(in out)    ::                                   &
  dqbydt(n_dp,nlev)        & ! increments to q (kg/kg/s)
 ,dthbydt(n_dp,nlev)       & ! increments to potential temperature (K/s)
 ,dtrabydt(n_dp,nlev,ntra)   ! increments to model tracers (kg/kg/s)

real(kind=real_umphys), intent(out)    ::                                      &
 down_flux(n_dp,nlev)        ! Downdraught mass flux (Pa/s) diagnostic

real(kind=real_umphys), intent(in out)   ::                                    &
  dt_dd(n_dp,nlev)         & ! dT/dt from convective downdraughts (K/s)
 ,dqbydt_dd(n_dp,nlev)       ! dq/dt from convective downdraughts (kg/kg/s)
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer ::                                                                     &
 i,j,k             & ! loop counters
, jtt(n_dp)          ! level k at which the pressure is > pcrit*p(1)

logical ::                                                                     &
  l_dd_present(n_dp)  ! DD present as precip from updraught

real(kind=real_umphys) ::                                                      &
  wdtrain(n_dp)             & ! detrained precip in level k (kg/m/s3)
 ,dpinv(n_dp)               & ! 1/dp
 ,sigd(n_dp)                & ! DD fractional area of a grid box
 ,sigs(n_dp)                & ! fraction of precip falling through unsat DD
 ,qsum(n_dp)                & ! Used for conservation check - column rate of
                              ! change of q
 ,qsum2(n_dp)               & ! Second array for conservation checking
 ,qloss(n_dp)               & ! rate of loss of q from the atmosphere
 ,surface_prec(n_dp)        & ! surface prec rate from this call to Emanuel
 ,cfactor(n_dp)               ! correction factor

real(kind=real_umphys) ::                                                      &
  evap(n_dp,nlev)           & ! evaporation from Downdraught (kg/kg/s)
 ,water(n_dp,nlev)          & ! rainwater in Downdraught (lp in paper)(kg/kg)
 ,wt(n_dp,nlev)             & ! fall speed of precipitate (Pa/s)
 ,t(n_dp,nlev)              & ! temperature on model levels (K)
 ,h(n_dp,nlev)              & ! static energy (J ?)
 ,gz(n_dp,nlev)             & ! gz a form of height measured from level 1
 ,dqbydt_dd1(n_dp,nlev)     & ! increments to q as a rate (kg/kg/s)
 ,dqbydt_dd2(n_dp,nlev)     & ! increments to q as a rate (kg/kg/s)
 ,md(n_dp,nlev)             & ! Downdraught mass flux for Emanuel scheme
                              ! (units kg/m2/s not Pa/s)
 ,qp_down (n_dp,nlev)       & ! specific humidity of downdraught (kg/kg)
 ,fraction_rain(n_dp,nlev)  & ! fraction of water(lp) in DD which is rain
 ,trap_down(n_dp,nlev,ntra)   ! tracer values of downdraught (kg/kg)

real(kind=real_umphys) ::                                                      &
 qsm          & ! initial estimate of q in downdraught (kg/kg)
,qstm         & ! Estimate of qsat in DD  (kg/kg)
,afac         & ! Evap = afac*revap (units unclear)
,b6,c6        & ! Coefficients for quadratic equation in (sqrt(lp))
,revap        & ! Proportional to reevaporation (units unclear)
,dhdp         & ! dz/dp
,rat          & ! Holds ratio of mass fluxes
,ginv         & ! 1/g
,tvx,tvy      & ! Virtual temperature (K)
,dp           & ! Layer thickness in pressure (Pa)
,coeff        & ! Temporary coeficient
,fac          & ! Temporary factor
,inc_t        & ! increment to temperature (K/s)
,lv_temp      & ! value between lc and ls depending on proportion of rain/snow
,rlv          & ! 1/lc or 1/ls
,denom                                                                         &
,dfraction_rain  &  ! change in fraction of rain water
,inc_rain                                                                      &
,inc_snow

! parameters for Emanuel scheme

real(kind=real_umphys), parameter ::                                           &
  coeffr = 1.0     & ! coefficient govering rate of evaporation of rain
, coeffs = 0.8     & ! coefficient govering rate of evaporation of snow
, omtrain = 50.0   & ! the assumed fall speed of rain (kg/m/s3 = Pa/s)
, omtsnow =  5.5   & ! the assumed fall speed of snow (kg/m/s3 = Pa/s)
, pcrit = 0.900    & ! point at which DD mass flux forced to decrease
                     ! with height (Note original Emanuel scheme had a value
                     ! of pcrit=0.949)
, sigd0 = 0.1      & ! fractional area covered by unsaturated DD (Note orig
                     ! value in the Emanuel code was 0.05).
, sigs0 = 0.12       ! fraction of precipitation falling outside cloud
                     ! i.e. through unsaturated DD

real(kind=real_umphys), parameter :: minloss = tiny(qloss)
                                          ! minimum allowable moisture loss

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EMAN_DD_REV'

!-----------------------------------------------------------------------
! 1.0 Initialisation of arrays
!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ginv = 1.0/g

do i= 1,n_dp
  qsum(i)  = 0.0
  qsum2(i) = 0.0
  sigs(i)  = sigs0
  jtt(i)   = 2
  l_dd_present(i) = .false.
  surface_prec(i) = 0.0

  ! Scale DD area using scaling from updraught closure

  sigd(i) = sigd0*scale_factor(i)
  ! Want to restrict sigd, cannot be > 1. but probably a problem if a
  ! very large fraction of area so choosing to limit to 0.7
  if (sigd(i) > 0.7) then
    sigd(i) = 0.7
  end if
end do

do k=1,nlev
  do i= 1,n_dp
    md(i,k) = 0.0
    wt(i,k) = 0.0
    down_flux(i,k) = 0.0
    evap(i,k)  = 0.0
    water(i,k) = 0.0
    fraction_rain(i,k) = 0.0

    dqbydt_dd1(i,k)  = 0.0
    dqbydt_dd2(i,k)  = 0.0

    ! convert from theta to T
    t(i,k)  = th(i,k) *exner_layer_centres(i,k)
  end do
end do

! Want the level k at which the pressure is > pcrit*p(k=1)
! Default value is 2
! Cannot use kterm here as jtt could be above top of convection for very
! shallow cases.

do k=nlev,2,-1
  do i= 1,n_dp
    if (p(i,k) >  (pcrit*p(i,1)) .and. jtt(i) == 2) then
      jtt(i) = k
    end if

    ! Is there precip from the updraught if so there is a downdraught & evap
    if (precip(i,k) > 0.0 .and. kterm(i) /= 0 .and. k <= kterm(i)+1) then
      l_dd_present(i) = .true.
    end if

  end do
end do

! Appears to be deriving some form of height using
!
!      p=rhoRT and dp/dz=-rhog    =>   gdz = -dp RT/p
!
!  equations here are using Tv  not T
!  Taking lowest model T level as zero height.
! Also initialising q in downdraught - note will only be used for highest
! level of downdraught. Lower levels will be calculated using value of qp_down
! in level above.

k = 1
do i= 1,n_dp
  gz(i,k) = 0.0
  h(i,k)  = cp*t(i,k)+gz(i,k)
  qp_down(i,k) = q(i,k)
end do
do k=2,nlev
  do i= 1,n_dp
    tvx=t(i,k)  *(1.0+c_virtual*q(i,k))
    tvy=t(i,k-1)*(1.0+c_virtual*q(i,k-1))
    gz(i,k)=gz(i,k-1)+0.5*r*(tvx+tvy)*(p(i,k-1)-p(i,k))/ph(i,k)

    ! dry static energy but gz term taking account of moisture
    h(i,k)  = cp*t(i,k)+gz(i,k)
    qp_down(i,k) = q(i,k-1)
  end do
end do

! initialise tracer in downdraught

if (l_tracer) then
  do j=1,ntra
    k=1
    do i=1,n_dp
      trap_down(i,k,j) = tracer(i,k,j)
    end do
    do k=2,nlev
      do i=1,n_dp
        trap_down(i,k,j) = tracer(i,k-1,j)
      end do
    end do
  end do
end if

!-----------------------------------------------------------------------
! 2.0 Main Down draught calculation
!-----------------------------------------------------------------------
! Assuming the array precip on model levels from G-R scheme corresponds
! to the Emanuel variable wdtrain - detrained precipitation.
!  Note In Emanuel scheme
!  wdtrain = g*ep*M*clw
!  ep  - precip efficiency
!  M   - up draught mass flux
!  clw - condensed cloud water
!-----------------------------------------------------------------------
!  In Emanuel's paper rainwater lp
!
!   -lp(k)+b6(k)sqrt(lp(k))+c6(k) = 0
!  The quadratic equation is solved to find the value of sqrt(lp)=revap(k)
!
!  Note both evaporation and reevaporation are proportional to sqrt(lp).
!-----------------------------------------------------------------------

do k=kmax_term+1,1,-1          ! level loop working Downwards

  do i=1,n_dp           ! grid point loop

    ! start at level above termination
    if (k <= kterm(i)+1 .and. l_dd_present(i)) then   ! test for Downdraughts

      ! From the updraught detrained precipitation (rain+snow) for level k
      wdtrain(i) = precip(i,k)*g

      if (t(i,k) <=  tm) then  ! All snow
        fraction_rain(i,k) = 0.0
      else if (t(i,k) > tm .and. t(i,k) < t_melt_snow) then
        fraction_rain(i,k) = (t(i,k) -tm)/(t_melt_snow - tm)
      else                     ! All rain
        fraction_rain(i,k) = 1.0
      end if

      wt(i,k) = fraction_rain(i,k)*omtrain +(1.0-fraction_rain(i,k))*omtsnow
      coeff   = fraction_rain(i,k)*coeffr  +(1.0-fraction_rain(i,k))*coeffs

      ! qsm - first guess for water vapour of parcel at level k
      ! For the top level of the downdraught this will just be q(i,k)

      qsm=0.5*(q(i,k)+qp_down(i,k+1))

      ! calculation for rain/snow water
      ! afac looks like expression for evaporation without sqrt(lp(k))

      afac=coeff*ph(i,k)*(qse(i,k)-qsm)                                        &
                       /(1.0e6+2.0e3*ph(i,k)*qse(i,k))
      afac=max(afac,0.0)
      b6=(ph(i,k)-ph(i,k+1))*sigs(i)*afac/wt(i,k)

      c6=( water(i,k+1)*wt(i,k+1) + wdtrain(i)/sigd(i) )/wt(i,k)

      if (c6 <= 0.0) then
        revap =0.0     ! set because of numeric problems
      else
        revap=0.5*(-b6+sqrt(b6*b6+4.0*c6))
      end if
      if (revap <  0.0) then
        revap =0.0        ! Cannot have negative value
      end if

      evap(i,k)=afac*revap

      ! water - rain water (lp in Emanuel paper)

      water(i,k)=revap*revap

      ! Calculate precipitating Downdraught mass flux under hydrostatic approx

      if (k >  1) then
        dhdp=(h(i,k)-h(i,k-1))/(p(i,k-1)-p(i,k))
        dhdp=max(dhdp,0.1)      ! correct units

        lv_temp = lc*fraction_rain(i,k)+ls*(1.0-fraction_rain(i,k))
        md(i,k)=ginv*sigd(i)*(lv_temp*sigs0*evap(i,k))/dhdp
        md(i,k)=max(md(i,k),0.0)

        ! Add a small amount of inertia to Downdraught
        ! (This is not referred to in Emanuel paper but is in the code with
        ! no explantion of why he had to do this.)
        ! Using layer thickness of layer k

        fac=2000.0/(ph(i,k-1)-ph(i,k))

        md(i,k)=(fac*md(i,k+1)+md(i,k))/(1.0+fac)


        ! Force Downdraught mass flux to decrease linearly to zero between about
        ! 950mb and the surface.
        ! Actually coded as decreasing linearly to zero between
        ! pcrit*lowest level pressure and model T level and the surface.
        ! (Is this what we really want it to do? Did the original Emanuel code
        ! assume level 1 was the surface?)
        ! Note Gregory-Rowntree DD are forced to start to decrease when within
        ! 100hPa of the surface. i.e. this imples a pcrit~0.9

        if (p(i,k) >  (pcrit*p(i,1))) then
          md(i,k)=md(i,jtt(i))*(p(i,1)-p(i,k))/(p(i,1)-p(i,jtt(i)))
        end if

      end if   ! k>1

      ! Find mixing ratio of precipitation Downdraught
      ! qstm - Estimated qsat value for downdraught used to apply a limit
      !        on qp_down

      if (k == 1) then
        qstm=qse(i,k)
      else
        qstm=qse(i,k-1)
      end if

      ! See Emanuel's paper for why the calculation of qp_down varies according
      ! to the mass flux.
      if (md(i,k) >  md(i,k+1)) then
        rat=md(i,k+1)/md(i,k)
        qp_down(i,k)=qp_down(i,k+1)*rat+q(i,k)*(1.0-rat)+ginv*                 &
                  sigd(i)*(ph(i,k-1)-ph(i,k))*(sigs0*evap(i,k)/md(i,k))

      else
        if (md(i,k+1) >  0.0) then
          if (t(i,k) > tm) then
            rlv = 1.0/lc
          else
            rlv = 1.0/ls
          end if
          ! This appears to be important in ensuring qp_down has a sensible
          ! value.
          qp_down(i,k)=(gz(i,k+1)-gz(i,k)+cp*(t(i,k+1)-t(i,k)))*rlv            &
                                       + qp_down(i,k+1)
        end if
      end if
      ! Not allowed to be greater that qsat for the level below
      ! qp down cannot be negative
      qp_down(i,k)=min(qp_down(i,k),qstm)
      qp_down(i,k)=max(qp_down(i,k),0.0)

    end if   ! test on DD present

  end do   ! loop over points
end do     ! loop over levels


! calculate tracers in downdraught plume

if (l_tracer) then
  do j=1,ntra
    do k=kmax_term+1,1,-1   ! level loop working Downwards
      do i=1,n_dp           ! grid point loop

        if (l_dd_present(i)) then
          if (k <= kterm(i)+1) then

            if (md(i,k) >  md(i,k+1)) then
              rat=md(i,k+1)/md(i,k)
              trap_down(i,k,j)=trap_down(i,k+1,j)*rat+tracer(i,k,j)*(1.0-rat)
            else
              if (md(i,k+1) >  0.0) then
                trap_down(i,k,j)=trap_down(i,k+1,j)
              end if
            end if

          end if
        end if

      end do   ! loop over points
    end do     ! loop over levels
  end do     ! loop over tracers
end if      ! test on tracers

! water corresponds to lp in paper
! Emanuel scheme assumes all surface precip rain. Note wt depends on
! whether rain or snow should we use this information?
!  now in mm/s = kg/m2/s  (wt in Pa/s)

! Add to any existing surface prec i.e. in the case of mid called after deep

do i=1,n_dp
  if (l_dd_present(i)) then
    inc_rain =     fraction_rain(i,1) *omtrain*sigd(i)*water(i,1)*ginv
    inc_snow = (1.0-fraction_rain(i,1))*omtsnow*sigd(i)*water(i,1)*ginv
    rain(i) = rain(i) + inc_rain
    snow(i) = snow(i) + inc_snow
    ! Surface precipitation from this call
    surface_prec(i) = inc_rain + inc_snow
  end if
end do

! ----------------------------------------------------------------------
! Calculation of increments
! ----------------------------------------------------------------------

!  level 1
k=1
do i=1,n_dp
  if (l_dd_present(i)) then
    dpinv(i) = 1.0/(ph(i,0)-ph(i,k+1))

    ! increments to temperature    dtheta = dt/exner
    ! evaporation of precipitation
    inc_t = (lcrcp*fraction_rain(i,k)+lsrcp*(1.0-fraction_rain(i,k)))          &
                           *sigd(i)*sigs(i)*evap(i,k)

    dthbydt(i,k) = dthbydt(i,k) - inc_t/exner_layer_centres(i,k)
    dt_dd(i,k)   = dt_dd(i,k)   - inc_t

    ! Cooling due to melting of falling snow.
    ! Calculation based on the increase in fraction of rain water

    dfraction_rain = fraction_rain(i,k) - fraction_rain(i,k+1)
    if (dfraction_rain > 0.0) then  ! increase in rain
      inc_t = sigd(i)*lfrcp*dfraction_rain*water(i,k)*omtsnow*dpinv(i)
      dthbydt(i,k) = dthbydt(i,k) - inc_t/exner_layer_centres(i,k)
      dt_dd(i,k)   = dt_dd(i,k)   - inc_t
    end if

    ! increments to q as rates
    ! change in q due to transport by Downdraught

    dqbydt_dd2(i,k)= g*md(i,k+1)*(qp_down(i,k+1)-q(i,k))*dpinv(i)

    ! increase in q due to evaporation of precip

    dqbydt_dd1(i,k) = sigd(i)*sigs(i)*evap(i,k)

    ! total increment to q

    dqbydt_dd(i,k) = dqbydt_dd1(i,k)+dqbydt_dd2(i,k)

  end if
end do      ! loop over gridpoints


!  Levels above k=1

do k=2,kmax_term+1        ! level loop
  do i=1,n_dp            ! grid point loop

    if (k <= kterm(i)+1 .and. l_dd_present(i) ) then

      dpinv(i)=1.0/(ph(i,k)-ph(i,k+1))

      ! increments to temperature    dtheta = dt/exner
      ! evaporation of precipitation
      inc_t = (lcrcp*fraction_rain(i,k)+lsrcp*(1.0-fraction_rain(i,k)))*       &
                       sigd(i)*sigs(i)*evap(i,k)

      dthbydt(i,k) = dthbydt(i,k) - inc_t/exner_layer_centres(i,k)
      dt_dd(i,k)   = dt_dd(i,k)   - inc_t

      ! Cooling due to melting of falling snow.
      ! Calculation based on the increase in fraction of rain water

      dfraction_rain = fraction_rain(i,k) - fraction_rain(i,k+1)
      if (dfraction_rain > 0.0) then  ! increase in rain
        inc_t = sigd(i)*lfrcp*dfraction_rain*water(i,k)*omtsnow*dpinv(i)
        dthbydt(i,k) = dthbydt(i,k) - inc_t/exner_layer_centres(i,k)
        dt_dd(i,k)   = dt_dd(i,k)   - inc_t
      end if

      ! increments to q as rates
      ! change in q due to transport by Downdraught

      dqbydt_dd2(i,k) = g*dpinv(i)*(md(i,k+1)*(qp_down(i,k+1)-q(i,k))          &
                                   -md(i,k)*(qp_down(i,k)-q(i,k-1)))

      ! increase in q due to evaporation of precip

      dqbydt_dd1(i,k) = sigd(i)*sigs(i)*evap(i,k)

      ! Total q increment

      dqbydt_dd(i,k) = dqbydt_dd1(i,k) +  dqbydt_dd2(i,k)

    end if           ! test on whether level in convecting column
  end do             ! gridpoint loop
end do               ! level loop

if (l_tracer) then
  do j = 1, ntra
    k = 1
    do i=1,n_dp
      if (l_dd_present(i)) then
        dpinv(i) = 1.0/(ph(i,1)-ph(i,2))
        dtrabydt(i,k,j) = dtrabydt(i,k,j)                                      &
                    +  g*dpinv(i)*md(i,k+1)*(trap_down(i,k+1,j)-tracer(i,k,j))
      end if
    end do
  end do

  do j = 1, ntra
    do k= 2, kmax_term+1
      do i=1,n_dp
        if (k <= kterm(i) .and. l_dd_present(i)) then

          dpinv(i)=1.0/(ph(i,k)-ph(i,k+1))
          dtrabydt(i,k,j) = dtrabydt(i,k,j)                                    &
                   +g*dpinv(i)*(md(i,k+1)*(trap_down(i,k+1,j)-tracer(i,k,j))-  &
                                md(i,k)*(trap_down(i,k,j)-tracer(i,k-1,j)))

        end if       ! test on whether level in convecting column
      end do         ! gridpoint loop
    end do           ! level loop
  end do             ! tracer loop

end if

! --------------------------------------------------------------------------
! Conservation checks on moisture - check assumes conserving for dp integral
! --------------------------------------------------------------------------
do k= 1, kmax_term+1
  do i=1,n_dp
    if (k <= kterm(i)+1 .and. l_dd_present(i)) then
      if (k == 1) then
        dp=(ph(i,0)-ph(i,k+1))       ! correct for surface layer
      else
        dp=(ph(i,k)-ph(i,k+1))       ! correct for layer k
      end if
      qsum(i) = qsum(i)  + dqbydt_dd(i,k)*dp*ginv
      qsum2(i)= qsum2(i) + precip(i,k)
    end if
  end do           ! gridpoint loop
end do             ! level


! Need to correct dq increments as moisture not conserved and non
! conservation is systematic leading to loss of q from atmosphere.

do i = 1, n_dp
  if (l_dd_present(i)) then
    qloss(i) = qsum2(i)-surface_prec(i)-qsum(i)
    ! Try to correct if qloss non zero and divisor not small
    denom = qsum(i)
    if (abs(qloss(i)) > minloss .and. abs(denom) > 10.0*minloss) then
      cfactor(i) = (qsum2(i)-surface_prec(i))/denom
      ! Check not a large change
      if (cfactor(i) < 0.5 .or. cfactor(i) > 1.5) then
        cfactor(i) = 1.0  ! don't try to correct
      end if
    else
      cfactor(i) = 1.0  ! don't try to correct as error very small
    end if
  end if
end do

! Scale q increments and add to output increments

do k= 1, kmax_term+1
  do i=1,n_dp
    if (k <= kterm(i)+1 .and. l_dd_present(i)) then
      dqbydt_dd(i,k) = dqbydt_dd(i,k)*cfactor(i)

      dqbydt(i,k) = dqbydt(i,k) + dqbydt_dd(i,k)
    end if
  end do           ! gridpoint loop
end do             ! level

! --------------------------------------------------------------------------
! change units of Emanuel Downdraught mass flux for diagnostic output
! --------------------------------------------------------------------------

do k= 1, kmax_term+1
  do i=1,n_dp
    if (k <= kterm(i)+1 .and. l_dd_present(i)) then
      down_flux(i,k) = g*md(i,k)
    end if
  end do
end do

! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine eman_dd_rev

end module eman_dd_rev_mod
