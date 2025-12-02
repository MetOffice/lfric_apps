! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates layer dependent constants

module layer_cn_6a_mod

use um_types, only: real_umphys

implicit none

! mutually exclusive convection indicator types.
integer, parameter ::  deep = 1          ! indicator all points are deep
integer, parameter ::  congestus = 2     ! indicator all points are congestus
integer, parameter ::  shallow = 3       ! indicator all points are shallow
integer, parameter ::  midlevel = 4      ! indicator all points are mid
!
! Description:
!   Calculates the following layer dependent constants:
!   * pressure (k, K+1/2 and k+1)
!   * layer thickness (k, k+1/2, k+1)
!   * entrainment coefficients (k+1/4, k+3/4)
!   * detrainment coefficients (k)
!
! Method:
!      See Unified Model documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='LAYER_CN_6A_MOD'

contains

subroutine layer_cn_6a(k, npnts, nlev,                                         &
                       ntml, ntpar, start_lev,                                 &
                       exner_layer_centres,                                    &
                       p_layer_boundaries, p_layer_centres,                    &
                       z_rho,                                                  &
                       conv_prog_precip,                                       &
                       recip_pstar, entrain_coef, rhum,                        &
                       ccp_strength,                                           &
                       zk, zkp12, zkp1,                                        &
                       wsc_o_mb, qsat_lcl, w_max,                              &
                       conv_indicator,                                         &
                       bconv,                                                  &
                       ! Out
                       pk, pkp1, exk, exkp1,                                   &
                       delpk, delpkp12, delpkp1,                               &
                       delp_uv_k, delp_uv_kp1,                                 &
                       ekp14, ekp34, amdetk                                    &
                       )

use cv_run_mod, only:                                                          &
    ent_fac_sh, ent_fac_dp, ent_fac_md, orig_mdet_fac, amdet_fac, ent_opt_dp,  &
    ent_opt_md, ent_dp_power, ent_md_power, mdet_opt_dp, mdet_opt_md,          &
    icvdiag, cldbase_opt_sh, c_mass_sh, w_cape_limit,                          &
    prog_ent_grad, prog_ent_int, prog_ent_max, prog_ent_min,                   &
    entrain_max, entrain_min
use cv_param_mod, only: ae2, refdepth_dp, refqsat, sh_wstar_closure,           &
    sh_grey_closure, entcoef, ccp_exp
use cv_dependent_switch_mod, only: l_new_det,                                  &
    l_const_ent

use ereport_mod,            only: ereport
use errormessagelength_mod, only: errormessagelength
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------
! Vector lengths and loop counters

integer,intent(in) :: k             ! present model layer
integer,intent(in) :: npnts         ! Number of points
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: ntml(npnts)   ! Number of levels in the surface-based
                                    ! turbulently mixed layer
integer,intent(in) :: ntpar(npnts)  ! Top of initial parcel ascent
integer,intent(in) :: start_lev(npnts)  ! Initiation level

! Field on model levels
real(kind=real_umphys),intent(in) :: exner_layer_centres(npnts,0:nlev)
                                                        ! Exner function
                                                        ! at layer centre
real(kind=real_umphys),intent(in) :: p_layer_centres(npnts,0:nlev)
                                                        ! Pressure
                                                        ! at layer centre (Pa)
real(kind=real_umphys),intent(in) :: p_layer_boundaries(npnts,0:nlev)
                                                        ! Pressure
                                                        ! at layer boundary (Pa)
real(kind=real_umphys),intent(in) :: z_rho(npnts,nlev)  ! height of rho levels
                                                        ! (m)
real(kind=real_umphys),intent(in) :: conv_prog_precip(npnts,nlev)
                                                ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
! Fields on a single level
real(kind=real_umphys),intent(in) :: recip_pstar(npnts)
                                      ! Reciprocal of pstar array (1/Pa)
real(kind=real_umphys),intent(in) :: entrain_coef(npnts)
                                      ! entrainment coefficients
real(kind=real_umphys),intent(in) :: rhum(npnts)
                                      ! Relative humidity at level K
real(kind=real_umphys),intent(in) :: zk(npnts)          ! height on k
real(kind=real_umphys),intent(in) :: zkp12(npnts)       ! height on k+1/2
real(kind=real_umphys),intent(in) :: zkp1(npnts)        ! height on k+1
real(kind=real_umphys),intent(in) :: wsc_o_mb(npnts)
                                      ! Convective velocity scale divided
                                      ! by cloud base mass flux mb
real(kind=real_umphys),intent(in) :: qsat_lcl(npnts)    ! qsat at the LCL
real(kind=real_umphys),intent(in) :: w_max(npnts)
                                      ! maximum large-scale w in column
real(kind=real_umphys),intent(in) :: ccp_strength(npnts)! Cold-pool strength

integer,intent(in) :: conv_indicator  ! type of convection for layer calls.
logical,intent(in) :: bconv(npnts)    ! Mask for points at which
                                      ! convection is occurring

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
real(kind=real_umphys),intent(out) :: pk(npnts)
                                      ! pressure at mid-point of layer k (Pa)
real(kind=real_umphys),intent(out) :: pkp1(npnts)
                                      ! pressure at mid-point of layer k+1 (Pa)
real(kind=real_umphys),intent(out) :: exk(npnts)
                                      ! Exner ratio at mid-point of layer k
real(kind=real_umphys),intent(out) :: exkp1(npnts)
                                      ! Exner ratio at mid-point of layer k+1
real(kind=real_umphys),intent(out) :: delpk(npnts)
                                      ! pressure difference across layer k (Pa)
real(kind=real_umphys),intent(out) :: delpkp12(npnts)
                                      ! pressure diff. across layer k+1/2 (Pa)
real(kind=real_umphys),intent(out) :: delpkp1(npnts)
                                      ! pressure diff. across layer k+1 (Pa)
real(kind=real_umphys),intent(out) :: delp_uv_k(npnts)
                                      ! pressure difference across UV
                                      ! layer k (Pa)
real(kind=real_umphys),intent(out) :: delp_uv_kp1(npnts)
                                      ! pressure difference across UV
                                      ! layer k+1 (Pa)

real(kind=real_umphys),intent(out) :: ekp14(npnts)
                                    ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(out) :: ekp34(npnts)
                                    ! Entrainment coefficient at level k+3/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(out) :: amdetk(npnts)
                                    ! Mixing detrainment coefficient at level k
                                    ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------

integer :: i              ! loop counter

real(kind=real_umphys) :: aekp14
                          ! Used in calculation of entrainment rate
real(kind=real_umphys) :: aekp34
                          ! Used in calculation of entrainment rate
real(kind=real_umphys) :: pkp12(npnts)
                          ! Pressure at upper boundary of layer k (Pa)
real(kind=real_umphys) :: pntml(npnts)
                          ! Pressure at upper boundary of layer ntml (Pa)
real(kind=real_umphys) :: delp_cld(npnts)   ! thickness of cloud layer (Pa)
real(kind=real_umphys) :: delz_cld(npnts)   ! thickness of cloud layer (m)
real(kind=real_umphys) :: delpkp14(npnts)   ! thickness of layer for k+1/4 (Pa)
real(kind=real_umphys) :: delpkp34(npnts)   ! thickness of layer for k+3/4 (Pa)
real(kind=real_umphys) :: tmp_entrain_coef(npnts)! entrainment coefficients
real(kind=real_umphys) :: factor
real(kind=real_umphys) :: factor1
real(kind=real_umphys) :: factor2

real(kind=real_umphys) :: entrain_ccp
                    ! entrainment dependence on cold-pool strength

character (len=errormessagelength) :: cmessage      ! used for ereport
integer                            :: icode         ! used for ereport
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LAYER_CN_6A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Initialise entrainment and detrainment rates to zero and hence they
! will default to 0.0 if not explicitly set.
!----------------------------------------------------------------------
do i=1,npnts
  ekp14(i)  = 0.0
  ekp34(i)  = 0.0
  amdetk(i) = 0.0
end do

!----------------------------------------------------------------------
! Set constant ae used in calculation of entrainment and detrainment
! rates depending upon level.
!----------------------------------------------------------------------
if (conv_indicator == deep) then
  !Deep only
  aekp14 = ent_fac_dp*ae2
  aekp34 = ent_fac_dp*ae2
else
  !Used only for mid-level
  aekp14 = ent_fac_md*ae2
  aekp34 = ent_fac_md*ae2
end if

!---------------------------------------------------------------------
! Calculate pressurea and pressure thicknesses
! NB pressure differences are calculated at k+1 - k to give positive
! thicknesses.
!---------------------------------------------------------------------
do i=1,npnts
  pk(i)           = p_layer_centres(i,k)
  pkp12(i)        = p_layer_boundaries(i,k)
  pkp1(i)         = p_layer_centres(i,k+1)
  exk(i)          = exner_layer_centres(i,k)
  exkp1(i)        = exner_layer_centres(i,k+1)
  delpk(i)        = p_layer_boundaries(i,k-1)   - p_layer_boundaries(i,k)
  delpkp1(i)      = p_layer_boundaries(i,k)     - p_layer_boundaries(i,k+1)
  delpkp12(i)     = pk(i)                       - pkp1(i)
  delpkp14(i)     = pk(i)                       - p_layer_boundaries(i,k)
  delpkp34(i)     = p_layer_boundaries(i,k)     - pkp1(i)
  delp_uv_k(i)    = p_layer_centres(i,k-1)      - p_layer_centres(i,k)
  delp_uv_kp1(i)  = p_layer_centres(i,k)        - p_layer_centres(i,k+1)
end do

if (conv_indicator == shallow) then
  ! Cloud thickness  - only used for shallow convection
  do i=1,npnts
    delp_cld(i)  = p_layer_boundaries(i,ntml(i)) -                             &
                   p_layer_boundaries(i,ntpar(i))
    pntml(i)     = p_layer_boundaries(i,ntml(i))
  end do
end if

if ((conv_indicator == deep) .and.                                             &
        (ent_opt_dp == 4 .or. ent_opt_dp == 5)) then
  ! Cloud thickness  - only used for deep convection with particular options.
  do i=1,npnts
    delz_cld(i)  = z_rho(i,ntpar(i)+1)-z_rho(i,ntml(i)+1)
  end do
end if

! ---------------------------------------------------------------------
! Calculate entrainment coefficients multiplied by approppriate
! layer thickness.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Shallow convection entrainment coefficients
! ---------------------------------------------------------------------
! Exponential decrease of entrainment rate with height above NTPAR
! is stopped to ensure massflux continues to decrease
! significantly with height (characteristic of shallow)

if (conv_indicator == shallow) then
  if (cldbase_opt_sh == sh_wstar_closure) then
    ! Original entrainment rates
    do i=1,npnts
      if (k >= ntml(i)) then
        ! Note: the factor of c_mass just cancels out the factor 1/c_mass
        ! held in wsc_o_mb passed in from shallow_conv.
        ekp14(i) = ent_fac_sh * wsc_o_mb(i) * delpkp14(i) * c_mass_sh          &
                 * exp( -1.0*min( 1.0, (pntml(i)-pk(i)   )/delp_cld(i) ) )     &
                 / delp_cld(i)
        ekp34(i) = ent_fac_sh * wsc_o_mb(i) * delpkp34(i) * c_mass_sh          &
                 * exp( -1.0*min( 1.0, (pntml(i)-pkp12(i))/delp_cld(i) ) )     &
                 / delp_cld(i)
      end if    ! level
      if (k == ntml(i)) ekp14(i) = 0.0
    end do
  else if (cldbase_opt_sh == sh_grey_closure) then
    ! Increased entrainment for grey zone (though decreases faster with height)
    do i=1,npnts
      if (k >= ntml(i)) then
        ekp14(i) = 2.0 * ent_fac_sh * wsc_o_mb(i) * delpkp14(i) * c_mass_sh    &
                 * exp( -1.5*min( 1.0, (pntml(i)-pk(i)   )/delp_cld(i) ) )     &
                 / delp_cld(i)
        ekp34(i) = 2.0 * ent_fac_sh * wsc_o_mb(i) * delpkp34(i) * c_mass_sh    &
                 * exp( -1.5*min( 1.0, (pntml(i)-pkp12(i))/delp_cld(i) ) )     &
                 / delp_cld(i)
      end if    ! level
      if (k == ntml(i)) ekp14(i) = 0.0
    end do
  end if
  ! ---------------------------------------------------------------------
  ! Congestus convection entrainment coefficients
  ! ---------------------------------------------------------------------
  ! 0.5/z rates
else if (conv_indicator == congestus) then
  do i=1,npnts
    if (k > ntml(i)) then
      ! Is there any reason why k=ntml isn't set?
      ekp34(i) = 0.5*2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
      ekp14(i) = 0.5*2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
    end if    ! type of convection and level
  end do

  ! ---------------------------------------------------------------------
  ! Deep convection entrainment coefficients
  ! ---------------------------------------------------------------------
else if (conv_indicator == deep) then

  select case (ent_opt_dp)
  case (0)     ! original Ap/(p*)^2  style entrainment
    do i=1,npnts
      ekp14(i) = entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do
  case (1)     !  n/z style entrainment  where n = ent_fac
    do i=1,npnts
      ekp14(i) = ent_fac_dp *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
      ekp34(i) = ent_fac_dp *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
    end do

  case (2)     ! Higher entrainment near surface (various versions have
              ! been used during GA3.0 plus testing) current code has
              ! version known as New 3.
              ! factor * original Ap/(p*)^2
    do i=1,npnts
      if (pk(i) > 50000.0) then
        factor1 = 1.0 + 1.25*(1.0-(100000.0 - pk(i))/50000.0)
      else
        factor1 = 1.0
      end if
      if (pkp12(i) > 50000.0) then
        factor2 = 1.0 + 1.25*(1.0-(100000.0 - pkp12(i))/50000.0)
      else
        factor2 = 1.0
      end if
      ekp14(i) = factor1 * entcoef * aekp14 * pk(i)    * delpkp14(i) *         &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = factor2 * entcoef * aekp34 * pkp12(i) * delpkp34(i) *         &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment
    do i=1,npnts
      ekp14(i) = entcoef * aekp14 * delpkp14(i) * recip_pstar(i) *             &
                (pk(i)    * recip_pstar(i))**ent_dp_power
      ekp34(i) = entcoef * aekp34 * delpkp34(i) * recip_pstar(i) *             &
                (pkp12(i) * recip_pstar(i))**ent_dp_power
    end do

  case (4)     ! variable n/z style entrainment
    do i=1,npnts
      tmp_entrain_coef(i) = refdepth_dp/delz_cld(i)
    end do
    do i=1,npnts
      ekp34(i) = tmp_entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))                 &
                               / (zkp12(i)+zkp1(i))
      ekp14(i) = tmp_entrain_coef(i) *2.0 * (zkp12(i) - zk(i))                 &
                               / (zkp12(i) + zk(i))
    end do

  case (5)     ! variable p/(p*)^2  style entrainment
    do i=1,npnts
      tmp_entrain_coef(i) = refdepth_dp/delz_cld(i)
    end do
    do i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (6)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the initiation level.
    do i=1,npnts
      if (w_max(i) < w_cape_limit) then
        if (conv_prog_precip(i,start_lev(i)) > 1.0e-10) then
          tmp_entrain_coef(i) = prog_ent_grad *                                &
                                log10(conv_prog_precip(i,start_lev(i)) *       &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = min(                                           &
                                max(tmp_entrain_coef(i), prog_ent_min),        &
                                prog_ent_max)
        else
          tmp_entrain_coef(i) = prog_ent_max
        end if
      else
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      end if
    end do
    do i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (7)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the current level.
    do i=1,npnts
      if (w_max(i) < w_cape_limit) then
        if (conv_prog_precip(i,k) > 1.0e-10) then
          tmp_entrain_coef(i) = prog_ent_grad *                                &
                                log10(conv_prog_precip(i,k) *                  &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = min(                                           &
                                max(tmp_entrain_coef(i), prog_ent_min),        &
                                prog_ent_max)
        else
          tmp_entrain_coef(i) = prog_ent_max
        end if
      else
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      end if
    end do
    do i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (8)   ! Like case 3 but use cold-pool strength to modify entrainment
    do i=1,npnts

      if (ccp_strength(i) > 0.0) then
        entrain_ccp = entrain_min                                              &
                      + (1.0 -  ccp_strength(i)**ccp_exp)                      &
                      * (entrain_max - entrain_min)
      else
        entrain_ccp = entrain_max
      end if

      ekp14(i) = entrain_ccp * delpkp14(i) * recip_pstar(i) *                  &
                (pk(i)    * recip_pstar(i))**ent_dp_power
      ekp34(i) = entrain_ccp * delpkp34(i) * recip_pstar(i) *                  &
                (pkp12(i) * recip_pstar(i))**ent_dp_power
    end do

  case DEFAULT

    icode = 1
    write(cmessage,fmt='(A)')                                                  &
     "Unexpected value of ent_opt_dp."
    call ereport(RoutineName, icode, cmessage)

  end select  ! test on ent_opt_dp

  ! ---------------------------------------------------------------------
  ! Mid-Level convection entrainment coefficients
  ! ---------------------------------------------------------------------
else

  select case (ent_opt_md)
  case (0)     ! orignal Ap/(p*)^2  style entrainment

    do i=1,npnts
      ekp14(i) = entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (1)     !  n/z style entrainment  where n = ent_fac

    do i=1,npnts
      ekp34(i) = ent_fac_md *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
      ekp14(i) = ent_fac_md *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
    end do

  case (2)     ! New 3 profile - higher near surface (option not used
              ! during GA3.0 plus testing).
    do i=1,npnts
      if (pk(i) > 50000.0) then
        factor1 = 1.0 + 1.25*(1.0-(100000.0 - pk(i))/50000.0)
      else
        factor1 = 1.0
      end if
      if (pkp12(i) > 50000.0) then
        factor2 = 1.0 + 1.25*(1.0-(100000.0 - pkp12(i))/50000.0)
      else
        factor2 = 1.0
      end if
      ekp14(i) = factor1 * entcoef * aekp14 * pk(i)    * delpkp14(i) *         &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = factor2 * entcoef * aekp34 * pkp12(i) * delpkp34(i) *         &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment
    do i=1,npnts
      ekp14(i) = entcoef * aekp14 * delpkp14(i) * recip_pstar(i) *             &
                (pk(i)    * recip_pstar(i))**ent_md_power
      ekp34(i) = entcoef * aekp34 * delpkp34(i) * recip_pstar(i) *             &
                (pkp12(i) * recip_pstar(i))**ent_md_power
    end do

  case (6)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the initiation level.
    do i=1,npnts
      if (w_max(i) < w_cape_limit  .and. qsat_lcl(i) > 0.0) then
        if (conv_prog_precip(i,start_lev(i)) > 1.0e-10) then
          tmp_entrain_coef(i) = prog_ent_grad *                                &
                                log10(conv_prog_precip(i,start_lev(i)) *       &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = min(                                           &
                                max(tmp_entrain_coef(i), prog_ent_min),        &
                                prog_ent_max)
        else
          tmp_entrain_coef(i) = prog_ent_max
        end if
      else
        ! If w_cape_limit is exceeded or qstat_lcl is negative
        ! then turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      end if
    end do
    do i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (7)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the current level.
    do i=1,npnts
      if (w_max(i) < w_cape_limit .and. qsat_lcl(i) > 0.0) then
        if (conv_prog_precip(i,k) > 1.0e-10) then
          tmp_entrain_coef(i) = prog_ent_grad *                                &
                                log10(conv_prog_precip(i,k) *                  &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = min(                                           &
                                max(tmp_entrain_coef(i), prog_ent_min),        &
                                prog_ent_max)
        else
          tmp_entrain_coef(i) = prog_ent_max
        end if
      else
        ! If w_cape_limit is exceeded or qstat_lcl is negative
        ! then turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      end if
    end do
    do i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                         &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *                   &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case (8)   ! Like case 0 but use cold-pool strength to modify entrainment

    do i=1,npnts

      if (ccp_strength(i) > 0.0) then
        entrain_ccp = entrain_min                                              &
                      + (1.0 -  ccp_strength(i)**ccp_exp)                      &
                      * (entrain_max - entrain_min)
      else
        entrain_ccp = entrain_max
      end if

      ekp14(i) = entrain_ccp * pk(i)    * delpkp14(i) *                        &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = entrain_ccp * pkp12(i) * delpkp34(i) *                        &
                 recip_pstar(i) * recip_pstar(i)
    end do

  case DEFAULT

    icode = 2
    write(cmessage,fmt='(A)')                                                  &
     "Unexpected value of ent_opt_md."
    call ereport(RoutineName, icode, cmessage)

  end select  ! test on ent_opt_md

end if        ! type of convection


! ---------------------------------------------------------------------
! If variable entrainment then overwrite previously calculated
! entrainment values when there is a valid (positive) entrain_coef
! entrain_coef is calculated in the diagnosis (conv_diag_comp).
! ---------------------------------------------------------------------
if ( (icvdiag == 4) .or. (icvdiag == 5) ) then
  if ( (conv_indicator == shallow) .or.                                        &
       (conv_indicator == congestus) .or.                                      &
       (conv_indicator == deep) ) then
    if (l_const_ent) then         ! no height dependence
      do i=1,npnts
        if (entrain_coef(i)  > 0.0) then
          ekp14(i) = entrain_coef(i) * (zkp12(i) - zk(i))
          ekp34(i) = entrain_coef(i) * (zkp1(i)-zkp12(i))
        end if
        if ((conv_indicator == shallow) .or.                                   &
            (conv_indicator == congestus )) then
          if (k == ntml(i)) ekp14(i) = 0.0
        end if
      end do
    else      ! rates vary with height
      do i=1,npnts
        if (entrain_coef(i)  > 0.0) then

          factor = 1.0

          ekp34(i) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))                 &
                                  *factor / (zkp12(i)+zkp1(i))
          ekp14(i) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))                 &
                                  *factor / (zkp12(i) + zk(i))
        end if    ! test on entrain_coef
        if ((conv_indicator == shallow) .or.                                   &
            (conv_indicator == congestus )) then
          if (k == ntml(i)) ekp14(i) = 0.0
        end if
      end do
    end if    ! test on l_const_ent
  end if      ! conv type
end if


! ---------------------------------------------------------------------
! Calculate mixing detrainment coefficient multiplied by appropriate
! layer thickness.
! ---------------------------------------------------------------------
! NB the mixing detrainment is always zero for the first convecting
! level.

! ---------------------------------------------------------------------
! Shallow convection mixing detrainment coefficients
! ---------------------------------------------------------------------
if (conv_indicator == shallow) then

  if (l_new_det) then    ! Use new relationship for detrainment
    do i=1,npnts
      if ( k > ntml(i) ) then
        if (entrain_coef(i) > 0.0) then    ! alter detrainment
          ! Trying 2* entrainment rates for shallow
          amdetk(i) = (1.0 + 1.0) * (ekp14(i) + ekp34(i))
          !    But set to 1 if RH is greater than 100%
          if (rhum(i)  >   1.0) then
            amdetk(i) = 1.0 * (ekp14(i) + ekp34(i))
          end if
        end if
      end if
    end do    ! npnts
  else if (cldbase_opt_sh == sh_grey_closure) then
    ! Increased entrainment for grey zone
    do i=1,npnts
      if (k  > ntml(i) ) then
        amdetk(i) = 1.5*(ekp14(i) + ekp34(i))
      end if
    end do    ! npnts
  else  !Original detrainment rates
    do i=1,npnts
      if ( k > ntml(i) ) then
        if (rhum(i)  <=  0.85) then
          amdetk(i) = (1.0 + 0.3) * (ekp14(i) + ekp34(i))
        else if (rhum(i)  >   1.0) then
          amdetk(i) =  1.0        * (ekp14(i) + ekp34(i))
        else  ! 0.85 <  rhum <= 1.0
          amdetk(i) = (1.0 + (0.3/0.15) * (1.0-rhum(i)))                       &
                                  * (ekp14(i) + ekp34(i))
        end if
      end if
    end do    ! npnts
  end if


  ! ---------------------------------------------------------------------
  ! Congestus convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
  ! adaptive mixing detrainment used in all cases
else if (conv_indicator == congestus) then
  do i=1,npnts
    if ( k > ntml(i) .and. rhum(i) <= 1.0 ) then
      amdetk(i) = (1.0-rhum(i)) * (ekp14(i) + ekp34(i))
    end if
  end do !npnts

  ! ---------------------------------------------------------------------
  ! Deep convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
else if (conv_indicator == deep) then
  select case (mdet_opt_dp)
  case (0)     ! Original orig_mdet_fac*(1.-1./(ent_fac_dp*ae2))
    do i=1,npnts
      if ( k > ntml(i)) then
        amdetk(i) = orig_mdet_fac*(ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
      end if
    end do !npnts
  case (1)     ! Adaptive mixing detrainment
    do i=1,npnts
      if ( k > ntml(i) .and. rhum(i) <= 1.0 ) then
        amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i))*(1-rhum(i))
      end if
    end do !npnts
  end select  ! mdet_dpt_dp

  ! ---------------------------------------------------------------------
  ! Mid-level convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
else
  select case (mdet_opt_md)
  case (0)     ! Original orig_mdet_fac*(1.-1./(ent_fac_md*ae2))
    do i=1,npnts
      if (bconv(i)) then
        amdetk(i) = orig_mdet_fac*(ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
      end if
    end do !npnts
  case (1)     ! Adaptive mixing detrainment
    do i=1,npnts
      if ( bconv(i) .and. rhum(i) <= 1.0 ) then
        amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i))*(1-rhum(i))
      end if
    end do !npnts
  end select  ! mdet_opt_md
end if  ! type of convection scheme

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return

end subroutine layer_cn_6a

end module layer_cn_6a_mod
