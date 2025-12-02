! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the mass flux profile for deep CMT.
!
module cmt_mass_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'CMT_MASS_MOD'
contains

subroutine cmt_mass(np_field, nconv, nlevs, nterm, cu_term,                    &
                    kterm, cu_tend, nlcl, ntop,                                &
                    mb, p_0degc, plcl, ptop, phalf,                            &
                  ! Output arguments
                    mass_up,mass_dwn,visc)

use cv_run_mod, only:                                                          &
    deep_cmt_opt

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use ereport_mod, only: ereport

use errormessagelength_mod, only: errormessagelength
implicit none

!-----------------------------------------------------------------------
! Description :
!  Calculates the mass flux profile for deep convection to be used in CMT
!  calculations. Uses the cloud-base mass flux from the plume scheme, but
!  profile is not the same as used for the thermodynamic part of the
!  convection scheme.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  np_field             & ! Full field length
 ,nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_term(nterm)       & ! Indices for terminating points
 ,kterm(np_field)      & ! Terminating levels
 ,cu_tend(nterm)       & !
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)            ! Top level of convection

real(kind=real_umphys), intent(in)    ::                                       &
  mb(nconv)            & ! Cloud base mass flux
 ,p_0degc(nconv)       & ! Pressure of melting level (hPa)
 ,plcl(nconv)          & ! Pressure of LCL (hPa)
 ,ptop(nconv)          & ! Pressure at top of convection (hPa)
 ,phalf(nlevs,nconv)     ! Pressure on model half levels (hPa)

real(kind=real_umphys), intent(out) ::                                         &
  mass_up(nlevs,nconv)  & ! Updraught mass flux profile (Pa/s)
 ,mass_dwn(nlevs,nconv) & ! Downdarught mass flux
 ,visc(nlevs,nconv)       ! Viscosity


! Local variables

integer ::                                                                     &
  i,j,k,n           ! loop counters
integer ::                                                                     &
  icode             ! error code

real(kind=real_umphys) ::                                                      &
  beta_l                                                                       &
 ,beta_u                                                                       &
 ,a_m(nterm)                                                                   &
 ,alpha_up                                                                     &
 ,alpha_dwn                                                                    &
 ,delta_p                                                                      &
 ,shape_factor

character (len=*), parameter ::  RoutineName = 'CMT_MASS'
character (len=errormessagelength) :: cmessage        ! error message

! Parameters (in future consider putting in a module).

real(kind=real_umphys), parameter ::                                           &
  a_max=1.5         & ! Ratio mass flux at melting level to cloud base flux
 ,a_top=0.2         & ! Ratio mass flux at top to cloud base flux
 ,a_dwn=0.0         & !
 ,p_ref=60000.0     & ! Used to scale ratios as p_0degC increase
 ,p_min=10000.0     & !
 ,top_press=15000.0 & !
 ,alpha_visc=0.30     !

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
!------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Determine shape of mass flux profile. If the depth of convection below
! the zero degree level is large enough, and the depth above is as well,
! the mass flux profile is assumed to have a maximum at the zero degree
! level.

do i=1,nterm
  j=cu_term(i)
  if ((p_0degc(j)-ptop(j)) >= p_min .and. (plcl(j)-p_0degc(j)) >= p_min) then

    ! Mass flux profile has an elevated maximum, the value of the maximum
    ! mass flux increase with decreasing pressure of the zero degree level
    ! while it is below p_ref. Above p_ref the maximum value is fixed.

    if (p_0degc(j) >= p_ref) then
      a_m(i)=1.0+(a_max-1.0)*(p_0degc(j)-(plcl(j)-p_min))/                     &
                             (p_ref     -(plcl(j)-p_min))
    else
      a_m(i)=a_max
    end if
  else
    a_m(i)=1.0
  end if
end do

do i=1,nterm
  j=cu_term(i)
  n=cu_tend(i)
  mass_up(nlcl(j),j)=mb(j)
  mass_dwn(nlcl(j),j)=a_dwn*mb(j)
  beta_l=log(a_m(i))
  beta_u=log(a_m(i)/a_top)
  if ((p_0degc(j)-ptop(j)) >= p_min .and. (plcl(j)-p_0degc(j)) >= p_min) then

    do k=nlcl(j)+1,ntop(j)+1
      if (phalf(k,j) >= p_0degc(j)) then
        mass_up(k,j)=a_m(i)*mb(j)*                                             &
                      exp(-beta_l*((phalf(k,j)-p_0degc(j))/                    &
                                   (plcl(j)-p_0degc(j)))**2)
        mass_dwn(k,j)=a_dwn*mb(j)

      else
        mass_up(k,j)=a_m(i)*mb(j)*                                             &
                      exp(-beta_u*((phalf(k,j)-p_0degc(j))/                    &
                                   (ptop(j)-p_0degc(j)))**2)
        mass_dwn(k,j)=0.0
      end if
    end do
  else
    do k=nlcl(j)+1,ntop(j)+1

      ! Choose either the diagnosed top of convection or the level at which
      ! scheme detrains, whichever is higher (this prevents odd failures
      ! in the CMT scheme at mid-latitudes)

      if (ptop(j) >= phalf(kterm(n)+1,j)) then
        mass_up(k,j)=a_m(i)*mb(j)*                                             &
                      exp(-beta_u*((phalf(k,j)-plcl(j))/                       &
                                   (phalf(kterm(n)+1,j)-plcl(j)))**2)
      else
        mass_up(k,j)=a_m(i)*mb(j)*                                             &
                     exp(-beta_u*((phalf(k,j)-plcl(j))/                        &
                                  (ptop(j)-plcl(j)))**2)
      end if
      mass_dwn(k,j)=0.0
    end do
  end if
end do

! Calculate eddy viscosity (Note viscosity=0 at NLCL and  KTERM+2)

do i=1,nterm
  j=cu_term(i)
  visc(nlcl(j),j)  =0.0
  visc(ntop(j)+2,j)=0.0
  delta_p=-(ptop(j)-plcl(j))

  ! alpha_up is a very crude representation of the non-dimensional profile
  ! of updraught vertical velocity. Tuning of the viscosity was done in the
  ! SCM, comparing momentum fluxes with those derived from CRM simulation
  ! of periods during TOGA-COARE.
  ! Allows for possibility of a downdraught component, but set to zero at
  ! present.

  if (deep_cmt_opt == 0) then

    do k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      if (p_0degc(j) < plcl(j) .and. phalf(k,j) > p_0degc(j)) then
        alpha_dwn =2.0*sqrt((phalf(k,j)-p_0degc(j))/(plcl(j)-p_0degc(j)))
      else
        alpha_dwn =0.0
      end if
      visc(k,j)=(alpha_up*mass_up(k,j)+alpha_dwn*mass_dwn(k,j))*               &
                  alpha_visc*delta_p
    end do

  else if (deep_cmt_opt == 1) then

    ! Added a shape factor to control gradient term
    ! Note as a_dwn=0.0 removed downward terms from calculation as waste of CPU

    do k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      shape_factor = 1.0-0.75*(phalf(k,j)-plcl(j))/(ptop(j)-plcl(j))
      visc(k,j)=alpha_up*mass_up(k,j)*alpha_visc*delta_p*shape_factor
    end do

  else if (deep_cmt_opt == 5) then

    ! Added a quadratic shape factor

    do k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      shape_factor = (1.0-0.75*(phalf(k,j)-plcl(j))/(ptop(j)-plcl(j)))**2
      visc(k,j)=alpha_up*mass_up(k,j)*alpha_visc*delta_p*shape_factor
    end do

  else       ! unallowed option

    icode = 1
    write(cmessage,'(a35,i3)') 'Unacceptable value of deep_cmt_opt ',          &
                deep_cmt_opt
    call ereport(routinename,icode,cmessage)

  end if

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine cmt_mass
end module cmt_mass_mod
