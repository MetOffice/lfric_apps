! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the evaporation of precipitation
!
module evp_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'EVP_MOD'
contains

subroutine evp(npnts, iphase, bevap, area_fac                                  &
               , precip, tevp, cca, rho, delq, delpkm1, pkm1                   &
               , evap, full_evap)

use planet_constants_mod, only: g

use cv_param_mod, only:                                                        &
   p_lq1, p_lq2, p_ice1, p_ice2, rho_lqp2, rho_lqa, rho_lqb,                   &
   rho_icp2, rho_icea, rho_iceb,                                               &
   lq_a, lq_b, lq_c, ice_a, ice_b, ice_c, ice_d

use science_fixes_mod, only: l_fix_conv_precip_evap

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! Description: Calculates the evaporation of precipitation
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
integer, intent(in) ::                                                         &
  npnts                 & ! Vector length
 ,iphase                  ! Indication for rain (1), or snow (2)

logical, intent(in) ::                                                         &
  bevap(npnts)            ! Mask for points where evaporation takes place

real(kind=real_umphys), intent(in) ::                                          &
  area_fac                ! Fraction of convective cloud amount to give
                          ! local cloud area
real(kind=real_umphys), intent(in) ::                                          &
  precip(npnts)         & ! Amount of precipitation (kg/m**2/s)

 ,tevp(npnts)           & ! Temperature of layer k (K)

 ,cca(npnts)            & ! Convective cloud amount (fraction)

 ,rho(npnts)            & ! Density of air

 ,delq(npnts)           & ! Change in humidity mixing ratio across layer k
                          ! (kg/kg)
 ,delpkm1(npnts)        & ! Change in pressure across layer k-1 (Pa)

 ,pkm1(npnts)             ! Pressure at level k-1 (Pa)

real(kind=real_umphys), intent(out) ::                                         &
  evap(npnts)             ! Evaporation

logical, intent(out) ::                                                        &
  full_evap(npnts)        ! True if all precip evaporated

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i               ! Loop counter


real(kind=real_umphys) ::                                                      &
  econ         & ! Quadratic term
 ,c1           & ! Constant
 ,c2           & ! Constant
 ,lrate        & ! Local rate of precipitation
 ,area         & ! Fractional area occupied by precipitation / downdraught
 , tl1,ti1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EVP'

!-----------------------------------------------------------------------
! Start of routine
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

full_evap(:) = .false.

tl1=0.5*p_lq1
ti1=0.5*p_ice1

if (iphase == 1) then        ! Rain

  do i=1,npnts
    if (bevap(i)) then
      if (precip(i)  >   0.0) then

        ! The following code implements the formulae in section 2.8.4.1b
        ! of UM Documentation Paper 027: Convection Schemes:
        econ = ((lq_a*tevp(i)+lq_b)*tevp(i)+lq_c)* (100000.0/pkm1(i))
        area = area_fac*cca(i)
        lrate = precip(i)/area

        if ( l_fix_conv_precip_evap ) then

          ! Correct formula doesn't include scaling by the area
          ! fraction here, as this is done afterwards in the routines that
          ! call evp (pevp_bcb and devap)
          c1 = rho_lqa * (lrate * lrate * rho(i))**tl1
          c2 = rho_lqb * (lrate**p_lq2) * (rho(i)**rho_lqp2)
          evap(i) = econ * (c1 + c2) * delq(i) * delpkm1(i) / g

          ! Limit evaporation to not exceed the total precipitation available
          if ( evap(i) >= lrate ) then
            evap(i) = lrate
            full_evap(i) = .true.
          end if

        else

          ! Original formula with erroneous duplicated factor of area
          c1 = rho_lqa * area * (lrate * lrate * rho(i))**tl1
          c2 = rho_lqb * area * (lrate**p_lq2) * (rho(i)**rho_lqp2)
          evap(i) = min( econ*(c1+c2)*delq(i)*delpkm1(i)/g, lrate )
          ! Possibly numerically dodgy way of assessing whether full evap
          if (evap(i) == lrate) full_evap(i) = .true.

        end if

      else
        evap(i) = 0.0
      end if
    end if
  end do

else if (iphase == 2) then        ! Snow

  do i=1,npnts
    if (bevap(i)) then
      if (precip(i)  >   0.0) then

        ! The following code implements the formulae in section 2.8.4.1b
        ! of UM Documentation Paper 027: Convection Schemes:
        if (tevp(i) <= 243.58) then
          econ = ice_d*(100000.0/pkm1(i))
        else
          econ = ((ice_a*tevp(i)+ice_b)*tevp(i)+ice_c)*(100000.0/pkm1(i))
        end if
        area = area_fac*cca(i)
        lrate = precip(i)/area

        if ( l_fix_conv_precip_evap ) then

          ! Correct formula doesn't include scaling by the area
          ! fraction here, as this is done afterwards in the routines that
          ! call evp (pevp_bcb and devap)
          c1 = rho_icea * (lrate * lrate * rho(i))**ti1
          c2 = rho_iceb * (lrate**p_ice2) * (rho(i)**rho_icp2)
          evap(i) = econ * (c1 + c2) * delq(i) * delpkm1(i) / g

          ! Limit evaporation to not exceed the total precipitation available
          if ( evap(i) >= lrate ) then
            evap(i) = lrate
            full_evap(i) = .true.
            ! Also don't allow negative evaporation (vapour subliming onto falling
            ! convective snow when supersaturated w.r.t. ice).
          else if ( evap(i) < 0.0 ) then
            evap(i) = 0.0
          end if

        else

          ! Original formula with erroneous duplicated factor of area
          c1 = rho_icea * area * (lrate * lrate * rho(i))**ti1
          c2 = rho_iceb * area * (lrate**p_ice2) * (rho(i)**rho_icp2)
          evap(i)=max( 0.0, min( econ*(c1+c2)*delq(i)*delpkm1(i)/g, lrate ) )
          ! Possibly numerically dodgy way of assessing whether full evap
          if (evap(i) == lrate) full_evap(i) = .true.

        end if

      else
        evap(i) = 0.0
      end if
    end if
  end do

end if     ! test on iphase

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine evp

end module evp_mod
