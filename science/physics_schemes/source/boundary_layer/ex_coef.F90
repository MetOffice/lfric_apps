! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate exchange coefficients for boundary layer
!           subroutine KMKH.

!  Programming standard: UMDP 3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module ex_coef_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'EX_COEF_MOD'
contains

subroutine ex_coef (                                                           &
! in levels/logicals
 bl_levels, k_log_layr, BL_diag,                                               &
! in fields
 sigma_h,flandg,dvdzm,ri,rho_wet_tq,z_uv,z_tq,z0m,zhpar,ntpar,                 &
 ntml_nl,ntdsc,nbdsc,l_shallow_cth,rmlmax2,rneutml_sq, delta_smag,             &
! in/out fields
 cumulus,weight_1dbl,                                                          &
! out fields
 lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm,rhokh,                   &
 fm_3d,fh_3d,tke_loc                                                           &
 )

use atm_fields_bounds_mod, only: pdims, tdims, pdims_s
use bl_diags_mod, only: strnewbldiag
use bl_option_mod, only:  WeightLouisToLong, Variable_RiC, cbl_op,             &
   sg_orog_mixing, ricrit_sharp, pr_max,                                       &
   local_fa,Prandtl,ishear_bl,L_SBLco,Muw_SBL,Mwt_SBL,sbl_op,                  &
   LockMailhot2004, lem_stability, lem_std, lem_conven,                        &
   lem_adjust, cbl_mix_fac_nml,                                                &
   off, on, sharpest, sharp_sea_long_land, sharp_sea_mes_land,                 &
   louis_tails, sharp_sea_louis_land, long_tails, mes_tails, ritrans,          &
   neut_cbl, lambda_min_nml, lambda_max_nml,                                   &
   lambda_fac, beta_bl, beta_fa, rlinfac, linear0,                             &
   to_sharp_across_1km, ntml_level_corrn, free_trop_layers, two_thirds,        &
   blending_option, blend_except_cu, blend_gridindep_fa, blend_cth_shcu_only,  &
   extended_tail, zero, one, one_half
use conversions_mod, only: pi => pi_bl
use gen_phys_inputs_mod, only: l_mr_physics

use planet_constants_mod, only: vkman => vkman_bl


use stochastic_physics_run_mod, only:                                          &
   l_rp2, rp_idx, par_mezcla_rp, g0_rp, ricrit_rp, lambda_min_rp,              &
   i_rp_scheme, i_rp2b, cbl_mix_fac_rp
use turb_diff_mod, only: l_subfilter_vert, l_subfilter_horiz

use model_domain_mod, only: model_type, mt_single_column

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use sblequil_mod, only: sblequil
implicit none

integer, intent(in) ::                                                         &
 bl_levels,                                                                    &
                 ! in maximum number of boundary layer levels
 k_log_layr
                 ! in num of levs requiring log-profile correction

integer, intent(in) ::                                                         &
 ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                 ! in Number of model layers in the turbulently
                 !    mixed layer as determined from the non-local
                 !    scheme.
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                 ! in Top level of any decoupled Sc
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                 ! in Bottom level of any decoupled Sc layer.
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! in Top level of parcel ascent

real(kind=r_bl), intent(in) ::                                                 &
 sigma_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                 ! in Standard deviation of subgrid
                 !    orography (m)
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                 ! in density on theta levels;
                 !    used in RHOKM so wet density
 rmlmax2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),       &
                 ! in Square of asymptotic mixing length for Smagorinsky scheme
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
                 ! in Z_UV(K) is height of u level k
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                 ! in Z_TQ(K) is height of T,Q level k
                 !    NOTE: RI(K) is held at Z_TQ(K-1)
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                 ! in Height of top of initial parcel ascent
 z0m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                 ! in Roughness length for momentum (m).
 dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       2:bl_levels),                                                           &
                 ! in Modulus of wind shear.
 ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),          &
                 ! in Local Richardson number.
 flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),          &
                 ! in Land fraction on all tiles.
 rneutml_sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                 ! in Square of the neutral mixing length for Smagorinsky
 delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! in delta_x used by Smagorinsky

logical, intent(in) ::                                                         &
  l_shallow_cth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! in Flag to indicate shallow convection based on cl-top

! Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

! INOUT arguments
logical, intent(in out) ::                                                     &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! INOUT Flag for boundary layer cumulus.
                 !       Can only be changed if ISHEAR_BL=1

real(kind=r_bl), intent(in out) ::                                             &
 weight_1dbl(pdims%i_start:pdims%i_end,                                        &
             pdims%j_start:pdims%j_end ,bl_levels)
                 ! INOUT Weighting applied to 1D BL scheme
                 !       to blend with Smagorinsky scheme,
                 !       index k held on theta level (k-1)
! out arguments
real(kind=r_bl), intent(out) ::                                                &
 rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,            &
       bl_levels),                                                             &
                 ! out Layer K-1 - to - layer K exchange coefficient
                 !       for momentum, on UV-grid with first and last
                 !       levels set to "missing data"
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       bl_levels),                                                             &
                 ! out Layer K-1 - to - layer K exchange coefficient
                 !       for scalars (but currently on th-levels)
                 ! On out: still to be multiplied by rho(if l_mr_physics)
                 !         and, for Ri-based scheme, interpolated to
                 !         rho levels in BDY_EXPL2
 zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! out Mixing layer height (m).

integer, intent(out) ::                                                        &
 ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! out Number of model layers in the turbulently
                 !     mixed layer as determined from the local
                 !     Richardson number profile.

real(kind=r_bl), intent(out) ::                                                &
 lambda_min,                                                                   &
                 ! out Min value of length scale LAMBDA.
 fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                 ! out stability function for momentum transport.
                 !     level 1 value is dummy for use in diagnostics
 fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                 ! out stability function for heat and moisture.
                 !     level 1 value is dummy for use in diagnostics
 tke_loc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
            2:bl_levels),                                                      &
                 ! out Ri-based scheme diagnosed TKE
 elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),         &
                 ! out Mixing length for momentum
 elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),         &
                 ! out Mixing length for scalars on theta levels
 elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels)
                 ! out Mixing length for scalars on rho levels

!-----------------------------------------------------------------------
!    Local and other symbolic constants :-

character(len=*), parameter ::  RoutineName = 'EX_COEF'

real(kind=r_bl) :: eh,em,g0,dh,dm,r_c_tke
real(kind=r_bl) :: subbmin,subbmax,subcmin,subcmax
real(kind=r_bl) :: a_ri,b_ri

parameter (                                                                    &
 eh=25.0_r_bl,                                                                 &
                 ! Used in calc of stability function FH.
 em=4.0_r_bl,                                                                  &
                 ! Used in calc of stability function FM.
 r_c_tke=one/0.41_r_bl                                                         &
                 ! used in calc of TKE (1/stress-energy ratio, see UMDP25)
           )

! Parameters used to set the stability function when the RP2b scheme
! is in use.  The min values correspond to the "conventional" subgrid
! model and the max values correspond to the "standard" subgrid model.
! See UMDP(024) and UMDP(081) for more details.
parameter (                                                                    &
 subbmin=1.43_r_bl,                                                            &
                 ! Used in calc of unstability function FH.
 subbmax=40.0_r_bl,                                                            &
                 ! Used in calc of unstability function FH.
 subcmin=1.43_r_bl,                                                            &
                 ! Used in calc of unstability function FM.
 subcmax=16.0_r_bl                                                             &
                 ! Used in calc of unstability function FM.
           )

!  Define local storage.
real(kind=r_bl) ::                                                             &
 ricrit(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                 ! Critical Richardson number
 func(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                 ! 2D variable for SBL stabiliy function options
 sharp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                 ! 2D variable for SHARP stabiliy function
 prandtl_number(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                 ! = KM/KH (currently only calculated for stable)
 BL_weight(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           bl_levels),                                                         &
                 ! Fractional weight applied to
                 ! BL function, vs free atmos
 turb_length(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
             2:bl_levels),                                                     &
                 ! Turbulent length scale on theta levels,
                 ! indexed as Ri (m)
 weight_bltop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! weight_1dbl at the top of the PBL

logical ::                                                                     &
 topbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! Flag for having reached the
                 ! top of the turbulently mixed layer.

! Variables for stability function tails
real(kind=r_bl) ::                                                             &
 fm_louis,                                                                     &
                 ! FM calculated using Louis
 fm_sharpest,                                                                  &
                 ! FM calculated using SHARPEST
 z_scale,                                                                      &
                 ! Scale height for interpolation
 g0_orog,                                                                      &
                 ! Orog dependent version of G0
 zpr,                                                                          &
                 ! z/sigma_h
 rpr
                 ! !/Pr

real(kind=r_bl) ::                                                             &
 pr_n,                                                                         &
                  ! neutral Prandtl number
 r_pr_n
                  ! 1 / neutral Prandtl number

real(kind=r_bl) ::                                                             &
 subb, subc, subg, ric, ricinv, rifac
                  !Constants for LEM stability functions

real(kind=r_bl) ::                                                             &
 f_log,                                                                        &
                 ! Temporary in calculation of logarithmic correction
 fh,                                                                           &
                 ! (Value of) stability function for heat & moisture.
 fm,                                                                           &
                 ! (Value of) stability function for momentum transport.
 rtmri,                                                                        &
                 ! Temporary in stability function calculation.
 vkz,                                                                          &
                 ! Temporary in calculation of ELH.
 lambdam,                                                                      &
                 ! Asymptotic mixing length for turbulent transport
                 !   of momentum.
 lambdah,                                                                      &
                 ! Asymptotic mixing length for turbulent transport
                 ! of heat/moisture.
 lambdah_rho,                                                                  &
                 ! Asymptotic mixing length for turbulent transport
                 ! of heat/moisture on rho levels
 rlambda_fac,                                                                  &
                 ! reciprocal of lambda_fac
 turb_length_layer,                                                            &
                 ! lengthscale in current layer
 beta,                                                                         &
                 ! empirical factor multiplying zh/delta
 zht,                                                                          &
                 ! top of boundary layer mixing
 zfa,                                                                          &
                 ! height to use beta_fa in blendin
 zz

integer ::                                                                     &
 i,j,                                                                          &
               ! Loop counter (horizontal field index).
 k, kl,                                                                        &
               ! Loop counters (vertical level index).
 kb, kt
               ! Base and top level of unstable Ri layers

logical ::                                                                     &
 subcrit       ! flag for being in a subcritical ri layer

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!   if stochastic physics random parameters is used set the parameter
!   used to vary the stability function to a perturbed value, if not
!   use the standard setting.
!-----------------------------------------------------------------------
if (l_rp2) then
  g0=g0_rp(rp_idx)
  ritrans=one/g0
  a_ri = zero       ! a_ri and b_ri need to be rewritten like this
  b_ri = 2.0_r_bl * g0  ! for bit-reproducibility
else
  g0=10.0_r_bl
  a_ri = (one-g0*ritrans)/(one-g0*ritrans/2.0_r_bl)**2.0_r_bl
  b_ri =    (g0/2.0_r_bl)/(one-g0*ritrans/2.0_r_bl)**2.0_r_bl
end if

dh=g0/eh                 ! Used in calc of stability function FH.
dm=g0/em                 ! Used in calc of stability function FM.
!---------------------------------------------------------------
! Set neutral and default Prandtl number (Pr=KM/KH)
!---------------------------------------------------------------
pr_n = one
if (Prandtl == LockMailhot2004) pr_n = 0.7_r_bl
! Use pr_n=0.7 if any LEM stability function selected
if (sbl_op == lem_stability .or. cbl_op == lem_std                             &
                            .or. cbl_op == lem_conven                          &
                            .or. cbl_op == lem_adjust) pr_n = 0.7_r_bl
r_pr_n = one / pr_n
!-----------------------------------------------------------------------
!  Set LAMBDA_MIN
!-----------------------------------------------------------------------
if (l_rp2) then
  lambda_min = lambda_min_rp(rp_idx)
else
  lambda_min = lambda_min_nml
end if

! Settings for LEM stability functions
if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
  ! Introduce stochastic physics in stability functions through the
  ! random parameter cbl_mix_fac_rp.  This allows the subb and
  ! subc settings for the stability function to vary between the
  ! "conventional" (subb = subbmin, subc = subcmin) and the
  ! "standard" (subb = subbmax, subc = subcmax) subgrid model
  ! (Brown 1999).

  subb = subbmin + (subbmax - subbmin) * cbl_mix_fac_rp(rp_idx)
  subc = subcmin + (subcmax - subcmin) * cbl_mix_fac_rp(rp_idx)

else

  if (cbl_op == lem_conven) then
    ! the "conventional" subgrid model, Brown (1999)
    subb = 1.43_r_bl
    subc = 1.43_r_bl
  else if (cbl_op == lem_adjust) then
    ! Using Parameterised LEM model, where 'cbl_mix_fac_nml' used to set
    ! subb and subc settings for the between the 'conventional' and 'standard'
    ! model. NB: This differs from RP2 parameterisation, as cbl_mix_fac_nml
    ! doesn't vary in time.
    subb = subbmin + (subbmax - subbmin) * cbl_mix_fac_nml
    subc = subcmin + (subcmax - subcmin) * cbl_mix_fac_nml
  else
    ! the "standard" LEM subgrid model, Brown (1999)
    subb = 40.0_r_bl
    subc = 16.0_r_bl
  end if

end if

subg = 1.2_r_bl
ric = 0.25_r_bl
ricinv = one/ric
rlambda_fac=one/lambda_fac

!$OMP PARALLEL DEFAULT(SHARED) private ( i, j )
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    !-----------------------------------------------------------------------
    ! 0. Initialise flag for having reached top of turbulently mixed layer
    !-----------------------------------------------------------------------
    topbl(i,j)     = .false.
    prandtl_number(i,j) = pr_n
    ! initialise blending weight at top of BL to one
    weight_bltop(i,j)   = one
  end do
end do
!$OMP end do NOWAIT
!-----------------------------------------------------------------------
! Set critical Richardson number
!-----------------------------------------------------------------------
if (l_rp2) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ricrit(i,j) = ricrit_rp(rp_idx)
    end do
  end do
!$OMP end do NOWAIT
else
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! Default critical Ri for Long_tails and Louis
      ricrit(i,j) = one
    end do
  end do
!$OMP end do NOWAIT
end if
!$OMP end PARALLEL

if (Variable_RiC == on) then

  select case (sbl_op)

    !--------------------------------------------
    ! SHARP TAILS
    !--------------------------------------------
  case (sharpest)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ricrit(i,j) = ricrit_sharp
      end do
    end do

    !--------------------------------------------
    ! LEM TAILS
    !--------------------------------------------
  case (lem_stability)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ricrit(i,j) = ric
      end do
    end do

    !--------------------------------------------
    ! SHARP over sea; longer tails over land
    !--------------------------------------------
  case (sharp_sea_long_land, sharp_sea_mes_land,                               &
       sharp_sea_louis_land)

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, flandg, ricrit, ricrit_rp, l_rp2 )                        &
!$OMP private( i, j )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (flandg(i,j) < one_half) then
              ! SHARPEST over sea
          ricrit(i,j) = ricrit_sharp
        else
              ! Longer tails over land
          if (l_rp2) then
            ricrit(i,j) = ricrit_rp(rp_idx)
          else
            ricrit(i,j) = one
          end if
        end if
      end do
    end do
!$OMP end PARALLEL do

  end select ! SBL_OP

end if
!-----------------------------------------------------------------------
! Initialise 3D stability functions
!-----------------------------------------------------------------------
if (l_subfilter_vert .or. l_subfilter_horiz) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( bl_levels, pdims, fm_3d, fh_3d )                                 &
!$OMP private( i, j, k )
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        fm_3d(i,j,k) = zero
        fh_3d(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
end if

!-----------------------------------------------------------------------
!  1.1 Use Richardson number profile to calculate BL depth, zh
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none) private(k,j,i)                                    &
!$OMP SHARED(bl_levels,pdims,topbl,ri,ricrit,local_fa,ntml_local,zh_local,z_uv)
do k = 2, bl_levels
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      !------------------------------------------------------------------
      ! If either a stable layer (Ri>RiCrit) or the maximum BL
      ! height has been reached, set boundary layer height (ZH_LOCAL) to
      ! the height of the lower boundary of the current layer
      !------------------------------------------------------------------
      if ( .not. topbl(i,j) .and.                                              &
           (ri(i,j,k) >  ricrit(i,j) .or. k == bl_levels) ) then
        topbl(i,j) = .true.
        if (local_fa >= ntml_level_corrn) then
          ! Ri(k)>RiC => theta-level(k-1) is supercrit => NTML=k-2
          ntml_local(i,j) = max( 1, k-2 )
        else
          ntml_local(i,j) = k-1
        end if
        zh_local(i,j) = z_uv(i,j,ntml_local(i,j)+1)
      end if
    end do  ! Loop over points
  end do  ! Loop over points
!$OMP end do NOWAIT
end do  ! Loop over levels
!$OMP end PARALLEL

! Save original diagnosis
if (BL_diag%l_zhlocal) then
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      BL_diag%zhlocal(i,j)=zh_local(i,j)
    end do  ! Loop over points
  end do  ! Loop over levels
end if

!-----------------------------------------------------------------------
! 1.2 In CUMULUS layers the local scheme is capped at the LCL (given in
!     this case by NTML_NL).
!     If NTML_LOCAL is greater than the top of the parcel ascent (NTPAR)
!     for a cumulus-capped layer, shear driven mixing is allowed to
!     dominate (if ISHEAR_BL=1 selected)
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( pdims, ishear_bl, ntml_local, ntpar, cumulus,                    &
!$OMP         bl_levels, lambda_min, rlambda_fac,                              &
!$OMP         turb_length, blending_option, rmlmax2)                           &
!$OMP private( i, j, k )
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( ishear_bl == 1 .and. ntml_local(i,j) > ntpar(i,j) ) then
      cumulus(i,j) = .false.
    end if
  end do
end do
!$OMP end do
!-----------------------------------------------------------------------
! 1.3 Search for sub-critical layers above the PBL and set the
!      mixing length to scale with these layer depths
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      turb_length(i,j,k) = lambda_min*rlambda_fac
    end do
  end do
end do
!$OMP end do NOWAIT
if (blending_option == blend_cth_shcu_only) then
  ! use Smag mixing length as background length scale if smaller
  ! than lambda_min (ie ignore lambda_min for high res simulations)
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        turb_length(i,j,k) = min( turb_length(i,j,k), sqrt(rmlmax2(i,j,k)) )
      end do
    end do
  end do
!$OMP end do NOWAIT
end if
!$OMP end PARALLEL

if (local_fa == free_trop_layers) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, bl_levels, ntml_local, ri, ricrit, z_uv,                  &
!$OMP         turb_length, rlambda_fac, lambda_min )                           &
!$OMP private( i, j, k, subcrit, kb, kt, kl, turb_length_layer )
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      subcrit = .false.
      do k = 3, bl_levels

        if ( k > ntml_local(i,j)+1  .and.                                      &
             ! we know Ri(ntml_local(i,j)+2) > RiCrit
             ri(i,j,k) < ricrit(i,j) .and. .not. subcrit ) then
          kb      = k   ! first level of subcritical Ri in layer
          subcrit = .true.
        end if
        if (ri(i,j,k) >= ricrit(i,j) .and. subcrit ) then
          kt      = k-1 ! last level of subcritical ri
          subcrit = .false.
          !---------------------------------------------------------
          ! turb_length(k) is held, with Ri(k), on th-level(k-1)
          !---------------------------------------------------------
          turb_length_layer   = z_uv(i,j,kt) - z_uv(i,j,kb-1)
          do kl = kb, kt
            turb_length(i,j,kl) = max( turb_length(i,j,kl),                    &
                        min(turb_length_layer,lambda_max_nml*rlambda_fac)   )
          end do
        end if

      end do
    end do
  end do
!$OMP end PARALLEL do
end if
!-----------------------------------------------------------------------
! When using turb_length, calculate within the BL
! and use the DSC layer depth as the length scale within a DSC layer
! Remember turb_length(k) is held, with Ri(k), on th-level(k-1)
!-----------------------------------------------------------------------
if (blending_option /= off) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( bl_levels, pdims, ntml_nl, ntml_local, turb_length, z_uv,        &
!$OMP         zh_local, nbdsc, ntdsc )                                         &
!$OMP private( i, j, k )
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( k-1 <= max(ntml_nl(i,j),ntml_local(i,j)) ) then
          turb_length(i,j,k) =  max( turb_length(i,j,k),                       &
              max( z_uv(i,j,ntml_nl(i,j)+1), zh_local(i,j) ) )
        end if
        if ( k-1 >= nbdsc(i,j) .and. k-1 <= ntdsc(i,j) ) then
          turb_length(i,j,k) = max( turb_length(i,j,k),                        &
                  ( z_uv(i,j,ntdsc(i,j)+1)-z_uv(i,j,nbdsc(i,j)) ) )
        end if
      end do
    end do
  end do
!$OMP end PARALLEL do
end if
!-----------------------------------------------------------------------
! 2.0 Loop over levels; calculate the mixing lengths
!-----------------------------------------------------------------------
do k = 2, bl_levels
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  PRIVATE(z_scale,j,i,lambdam,lambdah,                                    &
!$OMP  lambdah_rho,vkz,f_log,zz,zht,zfa,beta)                                  &
!$OMP  SHARED(k,pdims,ri,ricrit,flandg,ntml_local,ntml_nl,z_tq,                &
!$OMP  l_rp2,lambda_min,par_mezcla_rp,zh_local,turb_length,k_log_layr,         &
!$OMP  z_uv,z0m,elm,elh,elh_rho,blending_option,cumulus,l_shallow_cth,zhpar,   &
!$OMP  ntdsc,weight_1dbl,weight_bltop,delta_smag,rneutml_sq,BL_diag,local_fa)
  !-----------------------------------------------------------------
  ! 2.1 Calculate asymptotic mixing lengths LAMBDAM and LAMBDAH
  !-----------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if (l_rp2) then
        lambdam = max ( lambda_min , par_mezcla_rp(rp_idx)*zh_local(i,j) )
      else
        lambdam = max ( lambda_min , lambda_fac*zh_local(i,j) )
      end if
      !-----------------------------------------------------------------
      ! Reduce mixing lengths above BL
      !-----------------------------------------------------------------
      if (k >= ntml_local(i,j)+2) then
        lambdam = lambda_min
      end if

      lambdah = lambdam
      lambdah_rho  = lambdah

      if ( local_fa == free_trop_layers ) then
        lambdam = max( lambdam, lambda_fac*turb_length(i,j,k) )
        lambdah = max( lambdah, lambda_fac*turb_length(i,j,k) )
        ! lambdah_rho does not need to be recalculated under
        ! local_fa option "free_trop_layers" as the full KH profile
        ! will be interpolated in bdy_expl2
      end if
      !-----------------------------------------------------------------------
      ! 2.2 Calculate mixing lengths ELH, ELM coincident with RI(K) and so
      !     at Z_TQ(K-1)
      !-----------------------------------------------------------------------
      !  Incorporate log profile corrections to the vertical finite
      !  differences into the definitions of ELM and ELH.
      !  Note that ELH_RHO is calculated (on rho levels) for direct inclusion
      !  in RHOKH and also (as elh) on theta levels for the unstable
      !  stability functions and inclusion in RHOKH before interpolation
      !  (under local_fa option "free_trop_layers").
      !  To save computing logarithms for all K, the values of ELM and ELH
      !  are unchanged for K > K_LOG_LAYR.

      if (k  <=  k_log_layr) then
        vkz   = vkman * ( z_uv(i,j,k) - z_uv(i,j,k-1) )
        f_log = log( ( z_uv(i,j,k) + z0m(i,j)   ) /                            &
                     ( z_uv(i,j,k-1) + z0m(i,j) ) )
        elm(i,j,k) = vkz / ( f_log + vkz/lambdam )
        elh(i,j,k) = vkz / ( f_log + vkz/lambdah )
        vkz   = vkman * ( z_tq(i,j,k) - z_tq(i,j,k-1) )
        f_log = log( ( z_tq(i,j,k) + z0m(i,j)   ) /                            &
                     ( z_tq(i,j,k-1) + z0m(i,j) ) )
        elh_rho(i,j,k) = vkz / ( f_log + vkz/lambdah_rho )
      else
        vkz = vkman * ( z_tq(i,j,k-1) + z0m(i,j) )
        elm(i,j,k) = vkz / (one + vkz/lambdam )
        elh(i,j,k) = vkz / (one + vkz/lambdah )
        vkz = vkman * ( z_uv(i,j,k) + z0m(i,j) )
        elh_rho(i,j,k) = vkz / (one + vkz/lambdah_rho )
      end if
    end do
  end do
!$OMP end do
!----------------------------------------------------------------
! 2.3 Blend mixing lengths between 1D and 3D Smagorinsky
!----------------------------------------------------------------
  if (blending_option /= off) then
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        zz = z_tq(i,j,k-1)  ! height of rhokm(k)
        ! turb_length is the greater of the local and non-local
        ! BL depths up to that bl top
        z_scale = max( zz, turb_length(i,j,k) )
        ! zht = interface between BL and FA
        zht = max( z_uv(i,j,ntml_nl(i,j)+1) , zh_local(i,j) )
        ! Relevant scale in cumulus layers can be cloud top height, zhpar
        if ( cumulus(i,j) .and. ( blending_option /= blend_cth_shcu_only .or.  &
                                  l_shallow_cth(i,j) ) ) then
          z_scale = max( z_scale, zhpar(i,j) )
          zht     = max( zht, zhpar(i,j) )
        end if
        ! BL top includes decoupled stratocu layer, if it exists
        if (ntdsc(i,j) > 0) zht = max( zht, z_uv(i,j,ntdsc(i,j)+1) )
        ! Need to restrict z_scale to dsc depth within a dsc layer
        ! (given by turb_length) and to distance from dsc top below the
        ! dsc layer
        if ( k-1 <= ntdsc(i,j) ) then
          z_scale = min( z_scale,                                              &
              max( turb_length(i,j,k), z_uv(i,j,ntdsc(i,j)+1)-zz ) )
        end if

        ! Finally calculate 1D BL weighting factor
        if ( blending_option == blend_except_cu .and.                          &
             cumulus(i,j) .and. ntdsc(i,j) == 0) then
          ! pure cumulus layer so revert to 1D BL scheme
          weight_1dbl(i,j,k) = one
        else

          if ( blending_option == blend_gridindep_fa .or.                      &
               blending_option == blend_cth_shcu_only ) then
            if (zz <= zht) then
              weight_1dbl(i,j,k) =                                             &
               one - tanh( beta_bl*z_scale/delta_smag(i,j)) *                  &
                 max( zero,                                                    &
                   min( one, (linear0-delta_smag(i,j)/z_scale)*rlinfac) )
              weight_bltop(i,j) = weight_1dbl(i,j,k)
            else ! above PBL
              ! Above the PBL top (at zht) increase weight to one smoothly
              ! between zht and zfa in order to default to 1D BL when not
              ! turbulent.  There is some arbitrariness here but:
              ! a) we want to use a physical height, to avoid grid dependence
              ! b) for shallow PBLs at high resolution it seems sensible to
              !    get well (a PBL depth) above the resolved PBL before
              !    reverting to 1D
              ! c) for deep PBLs we still want to revert to 1D reasonably
              !    quickly, hence within at most 1km of zht
              zfa=min( 2.0_r_bl*zht, zht+1000.0_r_bl )
              if (zz <= zfa ) then
                weight_1dbl(i,j,k) = one + one_half *                          &
                                     (weight_bltop(i,j) - one) *               &
                                     ( one + cos(pi*(zz-zht)/(zfa-zht)) )
              else
                weight_1dbl(i,j,k) = one
              end if
              if ( local_fa == free_trop_layers .and.                          &
                   ri(i,j,k) < ricrit(i,j) ) then
                ! Except in an elevated turbulent layer where we still use
                ! the standard blending weight
                z_scale = turb_length(i,j,k)
                weight_1dbl(i,j,k) =                                           &
                 one - tanh( beta_bl*z_scale/delta_smag(i,j)) *                &
                  max( zero,                                                   &
                   min( one, (linear0-delta_smag(i,j)/z_scale)*rlinfac) )

              end if
            end if ! test on zz < zht
          else
            zfa=zht+1000.0_r_bl
            if (zz <= zht) then
              beta=beta_bl
            else if (zz <= zfa) then
              beta = beta_bl*(zfa-zz)/(zfa-zht) +                              &
                     beta_fa*(zz-zht)/(zfa-zht)
            else
              beta=beta_fa
            end if
            weight_1dbl(i,j,k) =                                               &
             one - tanh( beta*z_scale/delta_smag(i,j)) * max( zero,            &
                 min( one, (linear0-delta_smag(i,j)/z_scale)*rlinfac) )
          end if
        end if

        elm(i,j,k) = elm(i,j,k)*weight_1dbl(i,j,k) +                           &
                     sqrt(rneutml_sq(i,j,k-1))*(one-weight_1dbl(i,j,k))
        elh(i,j,k) = elh(i,j,k)*weight_1dbl(i,j,k) +                           &
                     sqrt(rneutml_sq(i,j,k-1))*(one-weight_1dbl(i,j,k))
      end do
    end do
!$OMP end do
  end if  ! test on blending_option
!$OMP end PARALLEL
end do  ! loop over levels
!----------------------------------------------------------------
! 3.0 Calculate stability functions
!----------------------------------------------------------------
! 3.1 Set-up a BL weighting function, =1 near the ground (ie in the BL)
!                                     =0 in the free troposphere
! Rate and height at which transition occurs varys depending on choices
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none) private( i, j, k, z_scale, zpr)                   &
!$OMP SHARED( pdims, bl_levels, BL_weight, local_fa, z_tq, sg_orog_mixing,     &
!$OMP         sigma_h)
!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      BL_weight(i,j,k) = one
    end do
  end do
end do
!$OMP end do

if (local_fa == to_sharp_across_1km) then
  !---------------------------------------------------------
  ! Additional code to allow the local Ri scheme to use
  ! SHARPEST in the free atmosphere, ie above the BL top,
  ! regardless of the tail option selected above.
  ! Set Z_SCALE to 1km to mimic old value of BL_LEVELS,
  ! gives BL_weight~0 by 2km, ~0.95 at 500m
  !---------------------------------------------------------
  z_scale = 1000.0_r_bl
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        zpr = z_tq(i,j,k-1)/z_scale
        BL_weight(i,j,k) = one_half*(one - tanh(3.0_r_bl*(zpr-one) ) )
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (sg_orog_mixing /= off) then
  !-----------------------------------------------------------------
  ! Subgrid orographic height dependence for SBL tail (option 1)
  ! or orographic dependence of mixing lengths, lambdam,h (opt 2)
  ! Gives BL_weight~[1,0.95,0.5,0] at ZPR=[0,0.6,1,1.7]
  !----------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (sigma_h(i,j) > one ) then
          zpr = z_tq(i,j,k-1)/sigma_h(i,j)
          BL_weight(i,j,k) = one_half*( one - tanh(4.0_r_bl*(zpr-one) ) )
        end if
      end do
    end do
  end do
!$OMP end do NOWAIT
end if
!$OMP end PARALLEL
! ----------------------------------------------------------------
! 3.2 calculating stable stability function
! ----------------------------------------------------------------
! Load up 2D array FUNC with selected stability function for Ri>=0
!
!  SBL_OP                 Option
!
!  Long_tails             Long tails
!  Sharpest               SHARPEST function
!  Sharp_sea_long_land    SHARPEST over sea ; Long tails over land
!  Mes_tails              MESOSCALE model: Louis/SHARPEST blend
!  Louis_tails            Louis function
!  Sharp_sea_mes_land     SHARP over sea; Mes over land
!  Sharp_sea_Louis_land   SHARP over sea; Louis over land
! ----------------------------------------------------------------
do k = 2, bl_levels
  select case (sbl_op)

    !--------------------------------------------
    ! long TAILS
    !--------------------------------------------
  case (long_tails)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (ri(i,j,k) >= zero)                                                 &
           func(i,j)=one / ( one + g0 * ri(i,j,k) )
      end do
    end do

    !--------------------------------------------
    ! SHARP TAILS
    !--------------------------------------------
  case (sharpest)

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private( i, j )               &
!$OMP SHARED( pdims, ri, ritrans, func, g0, a_ri, b_ri, k)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (ri(i,j,k)  <  ritrans ) then
          func(i,j) = one - one_half * g0 * ri(i,j,k)
        else
          func(i,j) = one / ( a_ri + b_ri*ri(i,j,k) )
        end if
        func(i,j)=func(i,j)*func(i,j)
      end do
    end do
!$OMP end PARALLEL do

    !--------------------------------------------
    ! LEM TAILS (cut-off at Ric)
    !--------------------------------------------
  case (lem_stability)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( ri(i,j,k) >= zero .and. ri(i,j,k)< ric ) then
          ! here func is essentially giving fh and the LEM stable
          ! prandtl_number will take back out the linear Ri term for fm
          rifac = (one-ri(i,j,k)*ricinv)**4
          func(i,j) = rifac*(one-subg*ri(i,j,k))
        else if (ri(i,j,k) >= ric) then
          func(i,j) = zero
        end if
      end do
    end do

    !--------------------------------------------
    ! SHARP over sea; long tails over land
    !--------------------------------------------
  case (sharp_sea_long_land)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (flandg(i,j) < one_half) then
              ! SHARPEST over sea
          if (ri(i,j,k)  <  ritrans ) then
            func(i,j) = one - one_half * g0 * ri(i,j,k)
          else
            func(i,j) = one / ( a_ri + b_ri*ri(i,j,k) )
          end if
          func(i,j)=func(i,j)*func(i,j)
        else
              ! Long tails over land
          if (ri(i,j,k) >= zero)                                               &
            func(i,j)= one / ( one + g0 * ri(i,j,k) )
        end if
      end do
    end do

    !--------------------------------------------
    ! MESOSCALE MODEL TAILS
    !--------------------------------------------
  case (mes_tails)

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
            ! Louis function
        if (ri(i,j,k) >= zero) then
          fm = one / ( one + one_half * g0 * ri(i,j,k) )
          fm_louis = fm * fm
        else
          fm_louis = one
        end if
            ! code for SHARPEST
        if (ri(i,j,k)  < ritrans ) then
          fm = one - one_half * g0 * ri(i,j,k)
        else
          fm = one / ( a_ri + b_ri*ri(i,j,k) )
        end if
        fm_sharpest = fm * fm
            ! Linear weighting function giving Louis
            ! at z=0, SHARPEST above Z_SCALE
        z_scale = 200.0_r_bl
        if ( z_tq(i,j,k-1)  >=  z_scale ) then
          func(i,j) = fm_sharpest
        else
          func(i,j)= fm_louis *( one - z_tq(i,j,k-1)/z_scale )                 &
                     + fm_sharpest * z_tq(i,j,k-1)/z_scale
        end if
      end do
    end do

    !--------------------------------------------
    ! LOUIS TAILS
    !--------------------------------------------
  case (louis_tails)

        ! LOUIS function
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (ri(i,j,k) >= zero) then
          func(i,j)=one / ( one + one_half * g0 * ri(i,j,k) )
          func(i,j)=func(i,j)*func(i,j)
        end if
      end do
    end do

    !--------------------------------------------
    ! SHARP TAILS OVER SEA; MES TAILS OVER LAND
    !--------------------------------------------
  case (sharp_sea_mes_land)
        ! SHARP sea; MES land
    z_scale = 200.0_r_bl
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, g0, ri, ritrans, a_ri, b_ri, flandg, func, z_tq,          &
!$OMP         k, z_scale )                                                     &
!$OMP private( i, j, fm, fm_louis, fm_sharpest )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
            ! Louis function
        if (ri(i,j,k) >= zero) then
          fm = one / ( one + one_half * g0 * ri(i,j,k) )
          fm_louis = fm * fm
        else
          fm_louis = one
        end if
            ! code for SHARPEST family
        if (ri(i,j,k) < ritrans) then
          fm = one - one_half * g0 * ri(i,j,k)
        else
          fm = one / ( a_ri + b_ri*ri(i,j,k) )
        end if
        fm_sharpest = fm * fm

            ! Linear weighting function giving Louis at z=0,
            ! SHARPEST above Z_SCALE
        if (flandg(i,j) < one_half) then
              ! SHARPEST family over sea
          func(i,j) = fm_sharpest
        else
              ! MES land
          if ( z_tq(i,j,k-1)  >=  z_scale ) then
            func(i,j) = fm_sharpest
          else
            func(i,j) = fm_louis *( one - z_tq(i,j,k-1)/z_scale )              &
                       + fm_sharpest * z_tq(i,j,k-1)/z_scale
          end if

        end if  ! FLANDG(i,j) < one_half

      end do ! loop over i
    end do ! loop over j
!$OMP end PARALLEL do

    !--------------------------------------------
    ! SHARP TAILS OVER SEA; LOUIS TAILS OVER LAND
    !--------------------------------------------
  case (sharp_sea_louis_land)
        ! SHARP sea; Louis land
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if (flandg(i,j) < one_half) then
              ! SHARP sea
          if (ri(i,j,k)  < ritrans ) then
            fm = one - one_half * g0 * ri(i,j,k)
          else
            fm = one / ( a_ri + b_ri*ri(i,j,k) )
          end if
          func(i,j)=fm * fm
        else
              ! Louis land
          if (ri(i,j,k) >= zero) then
            fm_louis = one / ( one + one_half * g0 * ri(i,j,k) )
            func(i,j)= (one - WeightLouisToLong) * fm_louis * fm_louis +       &
                       WeightLouisToLong * one / ( one + g0 * ri(i,j,k) )
          end if   ! ri >= 0
        end if  ! FLANDG(i,j) < one_half

      end do ! loop over i
    end do ! loop over j

  end select ! SBL_OP

  !------------------------------------------------------------------
  ! Additional option to use SHARPEST in the free atmosphere, ie above
  ! the BL top, regardless of the tail option selected above.
  !------------------------------------------------------------------
  if (local_fa == to_sharp_across_1km) then

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, ri, ritrans, sharp, a_ri, b_ri, func, BL_weight,          &
!$OMP         k, g0 )   private( i, j )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        !----------------------------
        ! Calculate SHARPEST function
        !----------------------------
        if (ri(i,j,k)  < ritrans ) then
          sharp(i,j) = one - one_half * g0 * ri(i,j,k)
        else
          sharp(i,j) = one / ( a_ri + b_ri*ri(i,j,k) )
        end if
        sharp(i,j)=sharp(i,j)*sharp(i,j)

        func(i,j) = func(i,j) * BL_weight(i,j,k)                               &
                  + sharp(i,j)*( one - BL_weight(i,j,k) )

      end do
    end do
!$OMP end PARALLEL do

  end if
  !------------------------------------------------------------------
  ! Additional code to allow the local Ri scheme stable mixing to
  ! depend the size of subgrid orography.  Overwrites values of FUNC
  ! as calculated above
  !------------------------------------------------------------------
  if (sg_orog_mixing == extended_tail) then

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        !-------------------------------------------------------
        ! SBL tail dependent on subgrid orography
        !  - use SHARPEST function but with variable coefficient
        !    that reduces to sharpest both with height above
        !    orography and as orography gets smaller
        !-------------------------------------------------------
        if ( sigma_h(i,j) > 0.1_r_bl ) then
          ! Then additional near-surface orographic dependence
          g0_orog = g0 / ( one +                                               &
                           (sigma_h(i,j)/25.0_r_bl)*BL_weight(i,j,k) )

          if (ri(i,j,k) < one/g0_orog) then
            func(i,j) = one - one_half * g0_orog * ri(i,j,k)
          else
            func(i,j) = one / ( 2.0_r_bl * g0_orog * ri(i,j,k) )
          end if
          func(i,j) = func(i,j)*func(i,j)

        end if
      end do
    end do

  end if
  !---------------------------------------------------------------
  ! 3.3 Set stable Prandtl number (=KM/KH)
  !---------------------------------------------------------------
  if (sbl_op == lem_stability) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, prandtl_number, pr_n, ri, ric, subg, k )                  &
!$OMP private( i, j )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( ri(i,j,k) >= zero .and. ri(i,j,k) < ric) then
          prandtl_number(i,j) = pr_n/(one-subg*ri(i,j,k))
        else if (ri(i,j,k) >= ric) then
          prandtl_number(i,j) = pr_n/(one-subg*ric)
        else
          prandtl_number(i,j) = pr_n
        end if
      end do
    end do
!$OMP end PARALLEL do
  else if (Prandtl == LockMailhot2004) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( pdims, prandtl_number, pr_n, ri, k )                             &
!$OMP private( i, j )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        prandtl_number(i,j) = min( pr_max,                                     &
                                   pr_n*(one + 2.0_r_bl*ri(i,j,k)) )
      end do
    end do
!$OMP end PARALLEL do
  end if
!----------------------------------------------------------------
! 3.4 Calculate (values of) stability functions FH, FM.
!----------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( k, pdims,BL_diag,elm,elh,ri,func,prandtl_number,cbl_op,r_pr_n,   &
!$OMP         l_subfilter_vert,l_subfilter_horiz,fm_3d,fh_3d,rhokm,rhokh,      &
!$OMP         rho_wet_tq,dvdzm,l_mr_physics,local_fa,tke_loc,subb,subc,g0,dm,  &
!$OMP         dh  )       &
!$OMP PRIVATE( i, j, fm, fh, rtmri, rpr )
!$OMP do SCHEDULE(STATIC)
   do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      if (BL_diag%l_elm3d) BL_diag%elm3d(i,j,k)=elm(i,j,k)
      if (ri(i,j,k)  >=  zero) then
        !-----------------------------------------------------------
        ! Note that we choose to include the Pr dependence such that
        ! fm(Ri=0)=1 and decreases slower than func with increasing
        ! Ri, due to GW activity, rather than have fh decrease
        ! faster than func
        !---------------------------------------------------------
        fh = func(i,j) * r_pr_n
        fm = fh * prandtl_number(i,j)

      else   ! ri < 0
        if (cbl_op == neut_cbl) then
          ! Use neutral stability for unstable mixing
          fm = one
          fh = r_pr_n
        else if (cbl_op == lem_std .or. cbl_op == lem_conven                   &
                                   .or. cbl_op == lem_adjust) then
          fm = sqrt(one-subc*ri(i,j,k))
          fh = sqrt(one-subb*ri(i,j,k)) * r_pr_n
        else
          !           ! UM_std
          rtmri = (elm(i,j,k)/elh(i,j,k)) * sqrt ( -ri(i,j,k) )
          fm = one - ( g0*ri(i,j,k) / ( one + dm*rtmri ) )
          fh = (one - ( g0*ri(i,j,k) / ( one + dh*rtmri ) )) * r_pr_n
        end if
      end if

      !------------------------------------------------------------------
      ! 4.0 Calculate exchange coefficients RHO*KM(K), RHO*KH(K)
      !     both on TH-level K-1 at this stage (RHOKH will be interpolated
      !     onto uv-levels and then be multiplied by ELH in BDY_EXPL2 if
      !     local_fa is not "free_trop_layers")
      !------------------------------------------------------------------

      if (l_subfilter_vert .or. l_subfilter_horiz) then
        fm_3d(i,j,k)=fm
        fh_3d(i,j,k)=fh
      end if
      if (BL_diag%l_fm) BL_diag%fm(i,j,k)=fm
      if (BL_diag%l_fh) BL_diag%fh(i,j,k)=fh

      rhokm(i,j,k) = rho_wet_tq(i,j,k-1) * elm(i,j,k) * elm(i,j,k)             &
                                  * dvdzm(i,j,k) * fm
      if (l_mr_physics) then
          ! Note "RHO" here is always wet density (RHO_WET_TQ) so
          ! save multiplication of RHOKH to after interpolation
        rhokh(i,j,k) =                elm(i,j,k) * dvdzm(i,j,k) * fh
      else
        rhokh(i,j,k) = rho_wet_tq(i,j,k-1) * elm(i,j,k) * dvdzm(i,j,k)         &
                       * fh
      end if
      ! If using the FA mixing length profile it is simplest to
      ! interpolate the full KH profile, including elh (in bdy_expl2)
      if (local_fa == free_trop_layers)                                        &
                  rhokh(i,j,k) = rhokh(i,j,k) * elh(i,j,k)

      if (BL_diag%l_tke) then
        rpr = fh / max(fm, tiny(one) )
        tke_loc(i,j,k) = ( r_c_tke*elm(i,j,k)*dvdzm(i,j,k)*dvdzm(i,j,k)        &
                           *(rhokm(i,j,k)/rho_wet_tq(i,j,k-1))                 &
                           *max( one_half, min( 10.0_r_bl,                     &
                           one - ri(i,j,k)*rpr ) )                             &
                         )**two_thirds
      end if

    end do !i
  end do !j
!$OMP end do
!$OMP end PARALLEL
end do ! bl_levels

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ex_coef
end module ex_coef_mod
