! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate wet bulb temperatute and potential temperature
!
! Subroutine Interface:
module thetaw_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='THETAW_MOD'

contains

subroutine Thetaw(                                                             &
! In:
       n,                                                                      &
       t,q,                                                                    &
       pressure, l_potential,                                                  &
! Out:
       tw)

use c_rmol,              only: rmol
use conversions_mod,     only: zerodegc
use qsat_mod,            only: qsat
use water_constants_mod, only: lc

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! Description:
!   Thetaw generates the wet bulb temperature or wet bulb potential
!   temperature for a field of points which on a surface defined
!   by array of pressures p, with input temperature and humidity
!   already interpolated on the pressure surface.
!
! Method:
! ---------------------------------------------------------------------
! 1. Calculate wet bulb temperature at requested pressure level
! 1.1 Determine g function for first guess of wet bulb T
! 1.2 Solve for wet bulb T using Newton iteration
! [Note that solution is to within a convergence criterion for each
! point individually, using flags so that results don't change if there
! is a different set of accompanying points.]
! ---------------------------------------------------------------------
! 2. Descend wet bulb line to 1000mb to determine wet bulb theta
! 2.1 Loop over delta(pressure) intervals
! 2.2 Integrate dT/dp from requested pressure level to 1000mb using
!     Runge-Kutta method:
!     a1= Dp dT/dp(p0     ,T0     ) a2= Dp dT/dp(p0+Dp/2,T0+a1/2)
!     a3= Dp dT/dp(p0+Dp/2,T0+a2/2) a4= Dp DT/dp(p0+Dp  ,T0+a3  )
!     T1 = T0 + a1/6 + a2/3 + a3/3 + a4/6
! ---------------------------------------------------------------------

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: physics_diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!   Documentation: UMDP 80

! Subroutine arguments
!   Scalar arguments with intent(in):
integer, intent(in) ::                                                         &
 n                   ! number of points
!   Array  arguments with intent(in):
real(kind=real_umphys), intent(in) ::                                          &
 t(n)                                                                          &
                     ! temperature on defined surface
,q(n)                                                                          &
                     ! humidity on defined surface
,pressure(n)         ! pressure on defined surface
logical ::                                                                     &
 l_potential         ! T=Output wet bulb potential temperature
                     ! F=Output wet bulb temperature
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):
real(kind=real_umphys), intent(out) ::                                         &
 tw(n)               ! wet bulb potential temperature

! Local parameters:
real(kind=real_umphys), parameter ::                                           &
 L_coeff = 2.34e3                                                              &
                     ! latent heat temperature dependence
,Cp_v    = 1850.0                                                              &
                     ! specific heat of water vapour
,Cp_d    = 1005.0                                                              &
                     ! specific heat of dry air
,m_v     = 0.01801                                                             &
                     ! Mol wt of water vapour KG/MOL
,m_d     = 0.02896                                                             &
                     ! Mol wt of dry air      KG/MOL
,Rstar   = rmol                                                                &
                     ! Universal gas constant R = rmol
,R_d     = Rstar/m_d                                                           &
                     ! R for dry air
,R_v     = Rstar/m_v                                                           &
                     ! R for water vapour
,delta_TW_tol = 0.005                                                          &
                     ! tolerance for TW convergence (degrees)
,p_1000mb = 1.0e5    ! pressure at 1000mb
integer, parameter ::                                                          &
 loop_max = 10                                                                 &
                     ! max loops for TW iterations
,loop_RK  = 5        ! 5 Iterations of Runge-Kutta are
                     ! sufficient for convergence upto 250Mb.
                     ! Significantly more are required for
                     ! higher levels, but are expensive.

! Local scalars:
integer ::                                                                     &
 i,loop
logical ::                                                                     &
 all_converged       ! test for TW convergence of Newton method
real(kind=real_umphys) ::                                                      &
 dgbydt                                                                        &
                     ! derivative of g with respect to T
,delta_TW                                                                      &
                     ! increment wet bulb T
,g                                                                             &
                     ! function of wet bulb T
,dTbydP                                                                        &
                     ! derivative of T with respect to p
,lh                                                                            &
                     ! latent heat scalar
,a1,a2,a3,a4         ! Runge-Kutta coefficients

! Local dynamic arrays:
logical ::                                                                     &
 converged(n)        ! test each TW convergence of Newton method
real(kind=real_umphys) ::                                                      &
 p(n)                                                                          &
                     ! pressure
,delta_p(n)                                                                    &
                          ! Pressure Increment
,qs(n)                                                                         &
                     ! saturated humidity
,g_TW(n)                                                                       &
                     ! g(TW)
,l(n)                                                                          &
                     ! latent heat
,Cp_moist(n)                                                                   &
                     ! specific heat of moist air
,tw_0(n),tw_1(n)     ! TW on pressure surface temporary arrays

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='THETAW'

!- End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ---------------------------------------------------------------------
! 1. Calculate wet bulb temperature at requested pressure level
! ---------------------------------------------------------------------

! Populate pressure array for this (and subsequent) qsat calls
do i=1,n
  p(i) = pressure(i)
end do

! Get saturated humidity with respect to water/ice
call qsat(qs,t,p,n)

! ---------------------------------------------------------------------
! 1.1 Determine g function for first guess of wet bulb T
!     g(T) = qs(T) L + C T
! ---------------------------------------------------------------------

do i=1,n
  tw(i) = t(i)                                         ! first guess
  l(i) = lc - L_coeff*(t(i) - zerodegc)     ! latent heat
  Cp_moist(i) = (1.0-q(i))*Cp_d + q(i)*Cp_v            ! specific heat
  g_TW(i) = l(i)*q(i)  + Cp_moist(i) * t(i)
  converged(i) = .false.
end do

! ---------------------------------------------------------------------
! 1.2 Solve for wet bulb T using Newton iteration,
!     where g(TW) = Cp T0 + L q = constant
! ---------------------------------------------------------------------

do loop=1,loop_max
  all_converged=.true.

  do i=1,n
    if (.not. converged(i)) then
      g      = l(i)*qs(i) + Cp_moist(i) * tw(i)
      dgbydt = l(i)*l(i)*qs(i) /                                               &
               max(0.0+epsilon(0.0),R_v*tw(i)*tw(i)) + Cp_moist(i)
      delta_TW = (g_TW(i) - g) / dgbydt
      tw(i) = tw(i) + delta_TW
      if (abs(delta_TW) >  delta_TW_tol) then
        all_converged=.false.
      else
        converged(i)=.true.
      end if
    end if ! not converged

  end do ! i

  if (all_converged) exit

  call qsat(qs,tw,p,n)     ! revise qs values

end do ! loop

if (l_potential) then
  ! ---------------------------------------------------------------------
  ! 2. Descend wet bulb line to 1000mb to determine wet bulb theta
  ! ---------------------------------------------------------------------

  do i=1,n
    delta_p(i)= (p_1000mb - pressure(i))/loop_RK
  end do

  ! ---------------------------------------------------------------------
  ! 2.1 Loop over delta(pressure) intervals
  ! ---------------------------------------------------------------------

  do loop=1,loop_RK

    ! Determine first guess for dT/dp at wet bulb temperature
    do i=1,n
      tw_0(i)=tw(i)
    end do

    ! ---------------------------------------------------------------------
    ! 2.2 Integrate dT/dp from requested pressure level to 1000mb using
    !     Runge-Kutta method:
    !     a1= Dp dT/dp(p0     ,T0     ) a2= Dp dT/dp(p0+Dp/2,T0+a1/2)
    !     a3= Dp dT/dp(p0+Dp/2,T0+a2/2) a4= Dp DT/dp(p0+Dp  ,T0+a3  )
    !     T1 = T0 + a1/6 + a2/3 + a3/3 + a4/6
    ! ---------------------------------------------------------------------

    ! A1 = dp * DTbyDP (p0      ,T0  )
    call qsat(qs,tw,p,n)
    do i=1,n
      lh = lc - L_coeff*(tw(i) - zerodegc)
      dTbydp = (lh*qs(i) + R_d*tw(i)) /                                        &
      (p(i) * (Cp_moist(i) + lh*lh*qs(i)/(R_v*tw(i)*tw(i)) ))
      p(i) = p(i) + delta_p(i)*0.5
      a1 = delta_p(i)*dTbydp
      tw(i) = tw_0(i) + a1*0.5
      tw_1(i) = tw_0(i) + a1/6.0 ! build up Runge-Kutta estimate
    end do

    ! A2 = dp * DTbyDP (p0 + dp/2,T0 + a1/2  )
    call qsat(qs,tw,p,n)
    do i=1,n
      lh = lc - L_coeff*(tw(i) - zerodegc)
      dTbydp = (lh*qs(i) + R_d*tw(i)) /                                        &
      (p(i) * (Cp_moist(i) + lh*lh*qs(i)/(R_v*tw(i)*tw(i)) ))
      a2 = delta_p(i)*dTbydp
      tw(i) = tw_0(i) + a2*0.5
      tw_1(i) = tw_1(i) + a2/3.0 ! build up Runge-Kutta estimate
    end do

    ! A3 = dp * DTbyDP (p0 + dp/2,T0 + a2/2  )
    call qsat(qs,tw,p,n)
    do i=1,n
      lh = lc - L_coeff*(tw(i) - zerodegc)
      dTbydp = (lh*qs(i) + R_d*tw(i)) /                                        &
      (p(i) * (Cp_moist(i) + lh*lh*qs(i)/(R_v*tw(i)*tw(i)) ))
      p(i) = p(i) + delta_p(i)*0.5
      a3 = delta_p(i)*dTbydp
      tw(i) = tw_0(i) + a3
      tw_1(i) = tw_1(i) + a3/3.0 ! build up Runge-Kutta estimate
    end do

    ! A4 = dp * DTbyDP (p0 + dp  ,T0 + a3  )
    call qsat(qs,tw,p,n)
    do i=1,n
      lh = lc - L_coeff*(tw(i) - zerodegc)
      dTbydp = (lh*qs(i) + R_d*tw(i)) /                                        &
      (p(i) * (Cp_moist(i) + lh*lh*qs(i)/(R_v*tw(i)*tw(i)) ))
      a4 = delta_p(i)*dTbydp
      tw(i) = tw_1(i) + a4/6.0 ! final Runge-Kutta estimate
    end do

  end do ! Loop over Runge-Kutta

end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine Thetaw
end module thetaw_mod
