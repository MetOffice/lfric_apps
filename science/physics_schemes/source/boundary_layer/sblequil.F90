! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: Estimate full and partial diffusivities
!           in an equilibrium Stable Boundary Layer

!  Programming standard: UMDP 3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
! Abbreviations used in variable descriptions:
! --------------------------------------------
! z:      height
! MO_loc: local Monin-Obukhov length
! ll:     turbulence mixing length
! zet:    ll/MO_loc (greek symbol zeta)
! kappa:  von Karman constant
! us_loc: local u_star  estimate (|uw|^0.5)
! WT_loc: local <w*THv> estimate
! KM/KH:  diffusivities for momentum and heat
!       [For X=U,THv: <w*X>=-KX*(dX/dz)]
! KUU/KTT/KUT/KTU: partial diffusivities.
!       [For X=U,THv: <w*X>=-KXU*(dU/dz)-KXT*(dTHv/dz)]
! PHI_KX: Normalized diffusivities
!----------------------------------------------------------------------
module sblequil_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'SBLEQUIL_MOD'
contains

subroutine sblequil (                                                          &
  Zhat,                                                                        &
  zetkz,pkmzet,pkhzet,pkuu,pktt,pkut,pktu,                                     &
  PHIe,phiww,ri,ce,rpow,rcb,cn,ierr                                            &
  )

use planet_constants_mod, only: vkman
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

!  Arguments passed in / out
real(kind=real_umphys), intent(in) ::                                          &
 Zhat
               ! in  z/MO_loc (>0.0)
real(kind=real_umphys), intent(out) ::                                         &
 zetkz,                                                                        &
               ! out l_hat=zeta/(kappa*Zhat)=ll/(kappa*z)
 pkmzet,                                                                       &
               ! out PHI_KM/zet=KM/(ll*us_loc)
 pkhzet,                                                                       &
               ! out PHI_KH/zet=KH/(ll*us_loc)
 pkuu,                                                                         &
               ! out KUU/(ll*us_loc)
 pktt,                                                                         &
               ! out KTT/(ll*us_loc)
 pkut,                                                                         &
               ! out KUT/(ll^2*g/Thv)
 pktu,                                                                         &
               ! out KTU/(ll^2*(g/Thv)*-1.0*(WT_loc/us_loc^2)^2)
 PHIe,                                                                         &
               ! out TKE/stress ratio
 phiww,                                                                        &
               ! out w-variance/stress ratio
 ri,                                                                           &
               ! out Equilib. SBL Richardson no.
 ce,                                                                           &
               ! out Turbulent dissipation timescl con
 rpow,                                                                         &
               ! out Copy of parameter npow (see below)
 rcb,                                                                          &
               ! out Copy of parameter CB (see below)
 cn        ! out Constant in length scale formula

integer, intent(out) ::                                                        &
 ierr      ! out Error status. Five-digit integer "MLKJI",
               !     where each digit has the value 0 (no-error)
               !     or 1 (error). Errors detected are as follows:
               !       I: Input Zhat < 0 (set to absolute value)
               !       J: Outer (zetkz,A1) loop exceeded max iteration
               !       K: Inner (PHIww) loop exceeded max iterations
               !          at least once
               !       L: Error in Newton's method for calculating
               !          PHIww (under/over-flow or wrong root found).
               !          Set PHIww to neutral value and continue.
               !       M: Either or both of PKHzet and PKMzet are too
               !          big/small, or negative. Return neutral soln.

! local constants
!  Tunable parameters
real(kind=real_umphys) :: npow,gmma,c1,c3,cb
parameter (                                                                    &
 npow=1.0,                                                                     &
                       !1/power in mixing-length formula
 gmma=4.0,                                                                     &
                       !Constant in 1st guess mix-len form
 c1=-0.002,                                                                    &
                       !Constant in param. of press-stress cov
 c3= -0.5,                                                                     &
                      !Constant in param. of press-stress cov
 cb= 1.05                                                                      &
                       !Constant in mixing length formula
)

!  Matching conditions for homogeneous shear (HS) layers
real(kind=real_umphys) :: pehsn,pehsw,puhsn,puhsw,pvhsn,pvhsw,pwhsn,pwhsw
real(kind=real_umphys) :: prhsn,prhsw,ptuhsw,ptthsw
parameter (                                                                    &
 pehsn= 1.0/0.170,                                                             &
                        !PHI_e (no wall)
 pehsw=  1.0/0.0862,                                                           &
                         !PHI_e (wall)
 puhsn= 0.480*pehsn,                                                           &
                        !PHI_u (no wall)
 puhsw=  0.460*pehsw,                                                          &
                         !PHI_u (wall)
 pvhsn= 0.260*pehsn,                                                           &
                        !PHI_v (no wall)
 pvhsw=  0.340*pehsw,                                                          &
                         !PHI_v (wall)
 pwhsn= 0.260*pehsn,                                                           &
                        !PHI_w (no wall)
 pwhsw=  0.200*pehsw,                                                          &
                         !PHI_w (wall)
 prhsn= 0.67,                                                                  &
                        !Prandtl no. (no wall)
 prhsw=  0.85,                                                                 &
                         !Prandtl no. (wall)
 ptuhsw= -2.1,                                                                 &
                         !PHI_ut (wall)
 ptthsw=  7.7                                                                  &
                         !PHI_tt (wall)

)

!  Parameters for iteration loops
real(kind=real_umphys) :: epszet,epsa1,epspw,vsml,vbig,pwlim
integer :: maxit
parameter (                                                                    &
 epszet =1.0e-3,                                                               &
                        !Tolerance for convergence of zeta
 epsa1  =1.0e-4,                                                               &
                        !Tolerance for convergence of A1
 epspw  =1.0e-2,                                                               &
                        !Tolerance for convergence of PHIww
 maxit  =100,                                                                  &
                        !Maximum permissible number of iterations
 vsml   =1.0e-20,                                                              &
                        !Limit for small numbers
 vbig   =1.0e+20,                                                              &
                        !Limit for big numbers
 pwlim  =1.0                                                                   &
                        !Lower limit in PHIw iteration
                        !(check for convergence towards wrong root)
)

!  Local storage

!  Dependent constants
real(kind=real_umphys) ::                                                      &
 c2,                                                                           &
                 !Constant in pressure-stress cov
 cc,                                                                           &
                 !Constant in pressure-stress cov
 ccn,                                                                          &
                 !no-wall value for CC
 d1,                                                                           &
                 !Constant in wall-effect terms
 d2,                                                                           &
                 !Constant in wall-effect terms
 d1t,                                                                          &
                 !Constant in wall-effect terms
 dd,                                                                           &
                 !Constant in param. of press-temp cov
 a1,                                                                           &
                 !Constant in param. of press-temp cov
 a2,                                                                           &
                 !Constant in param. of press-temp cov
 ct,                                                                           &
                 !Temperature variance  timescl con
 h2,                                                                           &
                 !Constant in PHI_KM expression
 ews,                                                                          &
                  !Intermediate constants
 ewb,                                                                          &
                 !Intermediate constants
 ewsn,                                                                         &
                 !no-wall value for EWS
 f1,f2,f3,f4,                                                                  &
                 !Intermediate constants
 f1n,f3n,f4n,                                                                  &
                 !no-wall values for F1, F3, F4
 h1,                                                                           &
                 !Constant in PHI_KH expression
 m1,m2,                                                                        &
                       !Constants in PHIww equation
 n0,n2,n3,n5,n6,n9 !Constants in PHIww equation

! Other local quantities
real(kind=real_umphys) ::                                                      &
 phitt,phiut,                                                                  &
                   !Norm T-var & horiz heat flx
 kc,                                                                           &
                   !kappa*CE
 kz,zeta,                                                                      &
                   !kappa*Zhat, zeta
 fwt,                                                                          &
                   !Wall-effect weighting fn, f=ll/(kappa*z)
 n_zetkz,n_pkhzet,n_pkmzet,n_phie,n_phiww !Neutral quantities

! Temporary variables
real(kind=real_umphys) ::                                                      &
 rcespw,                                                                       &
                       !1.0/(CE*sqrt(PHIww))
 tmp,                                                                          &
                       !temporary computational variables
 zetprev,zetchange,                                                            &
                       !temp values for zetkz
 a1prev,a1change,                                                              &
                       !temp values for A1
 pwprev,pwchange,                                                              &
                       !temp values for PHIww
 yy,py,dpy         !temp variables in PHIww-iteration

integer ::                                                                     &
 iloop,jloop,                                                                  &
                       !loop counters
 newterr,kxerr,                                                                &
                       !error flags
 ie1,ie2,ie3       !temporary error indicators

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SBLEQUIL'

!-----------------------------------------------------------------------
! Start of code
!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

rpow=npow
rcb =cb
ierr=0
ie1 =0
ie2 =0
ie3 =0
if (Zhat <  0) ierr=1

!-----------------------------------------------------------------------
! 0. Calculate dependent constants (C2,CC,D1,D2,DD,D1T,A2,CE,CT,H2)
!    and other useful local quantities
!-----------------------------------------------------------------------
tmp=2.0*(pwhsn-2.0/(puhsn-pvhsn))/pehsn
c2=1.0-c1/tmp
ccn=2.0*pwhsn*(1.0-c2)/(puhsn-pvhsn)
cc =2.0*pwhsw*(1.0-c2)/(puhsw-pvhsw)
tmp=pwhsw/(3.0*pwhsw*pwhsw+2.0)
d1=tmp*(2.0*(1.0-c2)*(2.0*pwhsw*puhsw+pwhsw*pvhsw                              &
   -3.0*pwhsw*pwhsw-4.0)/(3.0*(puhsw-pvhsw))-2.0*c1*pehsw/3.0)
tmp=c1-(2.0*(1.0-c2)*(pwhsw-2.0/(puhsw-pvhsw))/pehsw)                          &
   +3.0*d1/(pehsw*pwhsw)
d2=tmp*pehsw/(3.0*c2*pwhsw)
ewsn=ccn/(ccn+2.0*(1.0-c2)/3.0)
f1n=ccn
f3n=1.0-c2
f4n=sqrt((f3n-(3.0*c1/(2.0*ewsn)))/f1n)
ews=cc/(cc+2.0*(1.0-c2)/3.0+2.0*(2.0*c2*d2/3.0+d1))
f1=cc+1.5*d1
f3=1.0-c2+1.5*d2*c2
f4=sqrt((f3-(3.0*c1/(2.0*ews)))/f1)
dd=prhsn/(f4n*f4n)
d1t=prhsw/(f4*f4)-dd
f2=dd+d1t
a2=-1.0-((f2+ptuhsw*pehsw*dd*ews/3.0)*(f3-1.5*c1/ews)/f1)
ce=(pwhsw)**(-1.5)
ct=3.0*f2/(ptthsw*pehsw*ews)
h2=(1.0+c3)/dd
cn=(1.0-((cb*f4/ce)**(-1.0/npow)))**(-npow)

kc=vkman*ce
kz=vkman*Zhat

!-----------------------------------------------------------------------
! Neutral solution
!-----------------------------------------------------------------------
n_zetkz = 1.0
n_pkhzet= 1.0/(f2*ce*sqrt(f4))
n_pkmzet= 1.0
n_phiww = 1.0/f4
n_phie  = 3.0*n_phiww/ews

!-----------------------------------------------------------------------
! 1. Loop to get zetkz (and A1)
!-----------------------------------------------------------------------
zetprev=-vbig
a1prev =-vbig
    !first guesses for zetkz and A1
zetkz=1.0/((1.0+((kz*gmma)**(1.0/npow)))**(npow))
a1=0.5
zetchange=abs(zetprev-zetkz)
a1change =abs(a1prev -a1)
iloop=0
jloop=0

do while (                                                                     &
          ((zetchange >  epszet) .or. (a1change >  epsa1))                     &
          .and. (iloop <  maxit)                                               &
        )
  iloop=iloop+1
  !--------------------------------------------------
  ! 1.1 Parametric functions requiring zetkz (and A1)
  !--------------------------------------------------
  h1=(1.0-a1)/ct
  fwt=zetkz
  tmp=cc+2.0*(1.0-c2)/3.0+fwt*2.0*(2.0*c2*d2/3.0+d1)
  ews=cc/tmp
  ewb=(c2*(1.0-fwt*2.0*d2)-3.0-2.0*c3)/tmp
  f1=cc+fwt*1.5*d1
  f2=dd+fwt*d1t
  f3=1.0-c2+fwt*1.5*d2*c2
  f4=sqrt((f3-(3.0*c1/(2.0*ews)))/f1)
  zeta=zetkz*kz

  !--------------------------------------------------
  ! 1.2 Calculate PHIw, given zetkz
  !--------------------------------------------------
        !Set temporary constants
  m1=f1*f4*f4
  m2=c1*ewb/ews-h2*(1.0+a2)
  n0=zeta*zeta*zeta*(-h1*m2/(vkman*kc))
  n2=zeta*(ce*(h1*f1-h2*f2))
  n3=zeta*zeta*(((1.0-h1)*m2-h1*m1)/vkman)
  n5=-kc*ce*f1
  n6=zeta*(ce*(m2+(1.0-h1)*m1))
  n9=kc*ce*m1

      !Iteration starting with neutral PHIw (Newton's method)
  phiww=n_phiww+0.5
  pwprev=-vbig
  pwchange=abs(pwprev-phiww)
  jloop=0
  newterr=0
  do while ((pwchange >  epspw) .and. (jloop <  maxit))
    jloop=jloop+1
    pwprev=phiww
    yy=sqrt(phiww)
    py= n9*(yy**(9.0))+n6*(yy**(6.0))+n5*(yy**(5.0))                           &
       +n3*(yy**(3.0))+n2*(yy**(2.0))+n0
    dpy= 9.0*n9*(yy**(8.0))+6.0*n6*(yy**(5.0))                                 &
        +5.0*n5*(yy**(4.0))+3.0*n3*(yy**(2.0))+2.0*n2*yy
    if (abs(dpy) <  vsml) then
      phiww=pwprev
      newterr=1
    else
      yy=yy-py/dpy
      phiww=yy*yy
      if ((phiww <  pwlim) .and. (iloop >  1)) then
        phiww=pwprev
        newterr=1
      end if
    end if
    pwchange=abs(pwprev-phiww)
  end do !while loop to get PHIe
  if (jloop >= maxit) ie1=1
  if (newterr == 1) then
    phiww=n_phiww
    ie2=1
  end if

  !--------------------------------------------------
  ! 1.3 Compute other dimensionless functions
  !--------------------------------------------------
  kxerr=0
  PHIe=3.0*(phiww-zeta*2.0*ewb/(3.0*kc*sqrt(phiww)))/ews
  pkhzet=(sqrt(phiww)-zeta*h1/(kc*phiww))/(f2*ce)
  if (abs(pkhzet) <  vsml) pkhzet=vsml
  if (abs(pkhzet) >  vbig) pkhzet=vbig
  if ((pkhzet <  vsml) .and. (iloop == 1)) pkhzet=1.0

  if (pkhzet <  vsml) then
    kxerr=1
  else
    tmp=kc*f1*f4*f4*(phiww**(1.5))                                             &
       +zeta*(c1*ewb/ews-h2*(1.0+a2))
    pkmzet=tmp/(kc*ce*f1*phiww+zeta*h2/pkhzet)
    if (abs(pkmzet) <  vsml) pkmzet=vsml
    if (abs(pkmzet) >  vbig) pkmzet=vbig
    if ((pkmzet <  vsml) .and. (iloop == 1)) pkmzet=1.0
    if (pkmzet <  vsml) kxerr=1
  end if

  !--------------------------------------------------
  ! 1.4 Recalculate A1 & zetkz; return to top of loop
  !--------------------------------------------------
  zetprev=zetkz
  a1prev =a1
  if (kxerr == 0) then
    ri=0.0
    ri=zeta*pkmzet*pkmzet/(vkman*pkhzet)
    a1=0.5+1.5*ri*ri-ri*ri*ri
    if (ri <= 0.0)  a1=0.5
    if (ri >  1.0) then
      a1=1.0
      zetkz=1.0/((                                                             &
            ((cn                             )**(-1.0/npow))                   &
           +((Zhat/(zetkz*pkhzet*cb*cb*phiww))**( 0.5/npow))                   &
                  )**(npow))
    else
      zetkz=1.0/((                                                             &
            ((cn                         )**(-1.0/npow))                       &
           +((pkmzet*cb*sqrt(phiww)*zetkz)**(-1.0/npow))                       &
                  )**(npow))
    end if
  else
    zetkz=zetprev
    a1   =a1prev
  end if
  zetchange=abs(zetprev-zetkz)
  a1change =abs(a1prev -a1)
end do !while main (zetkz,A1) loop

if (kxerr == 1) then !return neutral soln
  zetkz =n_zetkz
  pkhzet=n_pkhzet
  pkmzet=n_pkmzet
  PHIe  =n_phie
  phiww =n_phiww
  ie3=1
end if

if (iloop >= maxit) ierr=ierr+10
if (ie1 /= 0)       ierr=ierr+100
if (ie2 /= 0)       ierr=ierr+1000
if (ie3 /= 0)       ierr=ierr+10000

!-----------------------------------------------------------------------
! 2. Calculate residual output quantities
!-----------------------------------------------------------------------
rcespw=1.0/(ce*sqrt(phiww))
phitt =rcespw/(ct*pkhzet)
tmp   =1.0/pkhzet+(1.0+a2)/pkmzet
phiut =-rcespw*tmp/dd

pkuu  = (f3*phiww-c1*PHIe/2.0)*rcespw/f1
pkut  = (1.0+c3)*phiut*pkhzet *rcespw/f1
pktt  = phiww                 *rcespw/f2
pktu  = (1.0-a1)*phitt*pkmzet *rcespw/f2

!-----------------------------------------------------------------------
! Finish up
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine sblequil
end module sblequil_mod
