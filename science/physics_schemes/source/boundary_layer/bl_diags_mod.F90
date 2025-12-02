! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  ---------------------------------------------------------------------
!  Sructure containing new boundary layer diagnostics.
!  This permits easier addition of new boundary layer
!  diagnostics without additional passing of arguments
!  though the boundary layer tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

!- ----------------------------------------------------------------------

module bl_diags_mod

use um_types, only: r_bl, integer_64

implicit none
save

type :: strnewbldiag

  ! Need to create a flag and a pointer
  ! Ordered by UM stash code
  ! Additional flags can be added for requesting calculation of specific
  ! diagnostics irrespective of whether they are required by STASH

  logical :: l_rhogamu =      .false. ! 130 countergradient taux
  logical :: l_rhogamv =      .false. ! 131 countergradient tauy
  logical :: l_rhogamt =      .false. ! 132 countergradient ftl
  logical :: l_rhogamq =      .false. ! 133 countergradient fqw
  logical :: l_elm =          .false. ! 134 mixing length
  logical :: l_tke_shr_prod = .false. ! 135 shear prod tke
  logical :: l_tke_boy_prod = .false. ! 136 buyo prod tke
  logical :: l_tke_dissp =    .false. ! 137 dissip tke
  logical :: l_sm =           .false. ! 138 diff coeff momentum
  logical :: l_sh =           .false. ! 139 diff coeff heat
  logical :: l_wb_ng =        .false. ! 140 buoy flux non-grad
  logical :: l_cf_trb =       .false. ! cf of tke scheme
  logical :: l_ql_trb =       .false. ! qcl of tke scheme
  logical :: l_sgm_trb =      .false. ! standard dev of qT distribution
  logical :: l_t_incr =       .false. ! 181 temp inc
  logical :: l_q_incr =       .false. ! 182 q inc
  logical :: l_qcl_incr =     .false. ! 183 qcl inc
  logical :: l_qcf_incr =     .false. ! 184 qcf inc
  logical :: l_u_incr =       .false. ! 185 u inc
  logical :: l_v_incr =       .false. ! 186 v inc
  logical :: l_w_incr =       .false. ! 187 w inc
  logical :: l_dtfric =       .false. ! 188 fric heating inc
  logical :: l_tl_incr =      .false. ! 189 TL inc
  logical :: l_qtl_incr =     .false. ! 190 qt inc
  logical :: l_cf_incr =      .false. ! 192 bcf inc
  logical :: l_cfl_incr =     .false. ! 193 cfl inc
  logical :: l_cff_incr =     .false. ! 194 cff inc
  logical :: l_ftl =          .false. ! 216 TL flux
  logical :: l_taux =         .false. ! 219 momentum x flux
  logical :: l_tauy =         .false. ! 220 momentum y flux
  logical :: l_fqw =          .false. ! 222 Qw flux
  logical :: l_zht =          .false. ! 304 turb mix height
  logical :: l_smltop =       .false. ! 356 surf mix layer top
  logical :: l_dsctop =       .false. ! 357 dec stratocu top
  logical :: l_zhlocal =      .false. ! 358 local zh top
  logical :: l_zhpar =        .false. ! 359 diag parcel top
  logical :: l_dscbase =      .false. ! 360 dec layer base
  logical :: l_cldbase =      .false. ! 361 scu cloud base
  logical :: l_weparm =       .false. ! 362 entrain rate surf mix layer
  logical :: l_weparm_dsc =   .false. ! 363 entrain rate dec scu
  logical :: l_dzh =          .false. ! 364 inv thickness
  logical :: l_oblen =        .false. ! 464 obukhov length
  logical :: l_ustar =        .false. ! 465 friction velocity
  logical :: l_wstar =        .false. ! 466 conv vel scale
  logical :: l_wbsurf =       .false. ! 467 surf buoy flux
  logical :: l_gradrich =     .false. ! 468 gradient ri
  logical :: l_dbdz =         .false. ! 469 vert buoy gradient
  logical :: l_dvdzm =        .false. ! 470 mod of wind shear
  logical :: l_rhokm =        .false. ! 471 momentum diffusivity
  logical :: l_rhokh =        .false. ! 472 heat diffusivity
  logical :: l_tke =          .false. ! 473 turb kin energy
  logical :: l_ostressx =     .false. ! 474 orog stress x
  logical :: l_ostressy =     .false. ! 475 orog stress y
  logical :: l_elm3d =        .false. ! 501 mix len momentum
  logical :: l_elh3d =        .false. ! 502 mix len heat
  logical :: l_rhokmloc =     .false. ! 503 local km
  logical :: l_rhokhloc =     .false. ! 504 local kh
  logical :: l_rhokmsurf =    .false. ! 505 surface driven km
  logical :: l_rhokhsurf =    .false. ! 506 surface driven kh
  logical :: l_rhokmsc =      .false. ! 507 cloud top driven km
  logical :: l_rhokhsc =      .false. ! 508 cloud top driven kh
  logical :: l_fh =           .false. ! 511 heat stability
  logical :: l_fm =           .false. ! 512 momentum stability
  logical :: l_weight1d =     .false. ! 513 weighting of 1D scheme
  logical :: l_slvar =        .false. ! 516 diagnostic of sL variance
  logical :: l_qwvar =        .false. ! 517 diagnostic of qw variance
  logical :: l_kblendtq =     .false. ! 581 coupling level (T/Q level)
  logical :: l_kblenduv =     .false. ! 582 coupling level (U/V level)
  logical :: l_slqw  =        .false. ! 518 diagnostic of sLqw covariance
  logical :: l_dissip =       .false. ! 519 diagnostic of dissipation rate
  logical :: l_grad_ftl =     .false. ! 800 explicit grad_ftl
  logical :: l_non_grad_ftl = .false. ! 801 explicit non_grad_ftl
  logical :: l_ftl_e =        .false. ! 802 explicit entrainment flux
  logical :: l_rhokhz_ex =    .false. ! 803 rhokhz from ex_flux_tq
  logical :: l_grad_t_adj =   .false. ! 804 grad_t_adj from ex_flux_tq

  ! Diagnostic calculation request switch - independent of STASH
  logical :: l_request_tke = .false.  ! Turbulent kinetic energy

  real(kind=r_bl), allocatable :: rhogamu(:, :, :)
  !                    130      counter gradient term of taux
  real(kind=r_bl), allocatable :: rhogamv(:, :, :)
  !                    131      counter gradient term of tauy
  real(kind=r_bl), allocatable :: rhogamt(:, :, :)
  !                    132      counter gradient term of ftl
  real(kind=r_bl), allocatable :: rhogamq(:, :, :)
  !                    133      counter gradient term of fqw
  real(kind=r_bl), allocatable :: elm(:, :, :)
  !                    134      mixing length
  real(kind=r_bl), allocatable :: tke_shr_prod(:, :, :)
  !                    135      production rate of TKE by shear
  real(kind=r_bl), allocatable :: tke_boy_prod(:, :, :)
  !                    136      production rate of TKE by buoyancy
  real(kind=r_bl), allocatable :: tke_dissp(:, :, :)
  !                    137      dissipation rate of TKE
  real(kind=r_bl), allocatable :: sm(:, :, :)
  !                    138      non-dimensional diffusion coef. for u, v
  real(kind=r_bl), allocatable :: sh(:, :, :)
  !                    139      non-dimensional diffusion coef. for t, q
  real(kind=r_bl), allocatable :: wb_ng(:, :, :)
  !                    140      non-gradient buoyancy flux
  real(kind=r_bl), allocatable :: cf_trb(:, :, :)
  !                    141      cloud fraction used in the TKE schemes
  real(kind=r_bl), allocatable :: ql_trb(:, :, :)
  !                    142      condensed water used in the TKE schemes
  real(kind=r_bl), allocatable :: sgm_trb(:, :, :)
  !                    143      standard deviation of the distribution
  !                             function in the TKE schemes
  real(kind=r_bl), allocatable :: t_incr(:,:,:)
  !                    181      temperature increment
  real(kind=r_bl), allocatable :: q_incr(:,:,:)
  !                    182      vapour increment
  real(kind=r_bl), allocatable :: qcl_incr(:,:,:)
  !                    183      liquid water increment
  real(kind=r_bl), allocatable :: qcf_incr(:,:,:)
  !                    184      ice water increment
  real(kind=r_bl), allocatable :: u_incr(:,:,:)
  !                    185      u wind increment
  real(kind=r_bl), allocatable :: v_incr(:,:,:)
  !                    186      v wind increment
  real(kind=r_bl), allocatable :: w_incr(:,:,:)
  !                    187      w wind increment
  real(kind=r_bl), allocatable :: dTfric(:, :, :)
  !                    188      Heating increment from turbulence dissipation
  real(kind=r_bl), allocatable :: cf_incr(:,:,:)
  !                    192      bulk cloud fraction increment
  real(kind=r_bl), allocatable :: cfl_incr(:,:,:)
  !                    193      liquid cloud fraction increment
  real(kind=r_bl), allocatable :: cff_incr(:,:,:)
  !                    194      frozen cloud fraction increment
  real(kind=r_bl), allocatable :: zht(:, :)
  !                    304      Turbulent mixing height
  real(kind=r_bl), allocatable :: smltop(:, :)
  !                    356      Top of surface mixed layer
  real(kind=r_bl), allocatable :: dsctop(:, :)
  !                    357      Top of decoupled stratocu layer
  real(kind=r_bl), allocatable :: zhlocal(:, :)
  !                    358      BL depth diagnosed from Ri>RiCrit
  real(kind=r_bl), allocatable :: zhpar(:, :)
  !                    359      Height of diagnosis parcel top
  real(kind=r_bl), allocatable :: dscbase(:, :)
  !                    360      Height of decoupled layer base
  real(kind=r_bl), allocatable :: cldbase(:, :)
  !                    361      Height of stratocumulus cloud base
  real(kind=r_bl), allocatable :: weparm(:, :)
  !                    362      Entrainment rate for SML
  real(kind=r_bl), allocatable :: weparm_dsc(:, :)
  !                    363      Entrainment rate for DSC
  real(kind=r_bl), allocatable :: dzh(:, :)
  !                    364      Inversion thickness
  real(kind=r_bl), allocatable :: oblen(:, :)
  !                    464      Surface Obukhov length
  real(kind=r_bl), allocatable :: ustar(:, :)
  !                    465      Friction velocity
  real(kind=r_bl), allocatable :: wstar(:, :)
  !                    466      Convective velocity scale
  real(kind=r_bl), allocatable :: wbsurf(:, :)
  !                    467      Surface buoyancy flux
  real(kind=r_bl), allocatable :: gradrich(:, :, :)
  !                    468      Gradient Richardson number
  real(kind=r_bl), allocatable :: dbdz(:, :, :)
  !                    469      Vertical buoyancy gradient
  real(kind=r_bl), allocatable :: dvdzm(:, :, :)
  !                    470      Modulus of wind shear
  real(kind=r_bl), allocatable :: rhokm(:, :, :)
  !                    471      BL Momentum diffusivity
  real(kind=r_bl), allocatable :: rhokh(:, :, :)
  !                    472      BL Thermal diffusivity
  real(kind=r_bl), allocatable :: tke(:, :, :)
  !                    473      Turbulent kinetic energy
  real(kind=r_bl), allocatable :: ostressx(:, :, :)
  !                    474      Orographic stress (x-component)
  real(kind=r_bl), allocatable :: ostressy(:, :, :)
  !                    475      Orographic stress (y-component)
  real(kind=r_bl), allocatable :: elm3d(:, :, :)
  !                    501      Mixing length for momentum
  real(kind=r_bl), allocatable :: elh3d(:, :, :)
  !                    502      Mixing length for heat and moisture
  real(kind=r_bl), allocatable :: rhokmloc(:, :, :)
  !                    503      Km diffusion coeff from local scheme
  real(kind=r_bl), allocatable :: rhokhloc(:, :, :)
  !                    504      Kh diffusion coeff from local scheme
  real(kind=r_bl), allocatable :: rhokmsurf(:, :, :)
  !                    505      Km diffusion coeff for surface-driven turb
  real(kind=r_bl), allocatable :: rhokhsurf(:, :, :)
  !                    506      Kh diffusion coeff for surface-driven turb
  real(kind=r_bl), allocatable :: rhokmsc(:, :, :)
  !                    507      Km diffusion coeff for  cloud-top-driven turb
  real(kind=r_bl), allocatable :: rhokhsc(:, :, :)
  !                    508      Kh diffusion coeff for  cloud-top-driven turb
  real(kind=r_bl), allocatable :: fh(:, :, :)
  !                    511      stability function for scalars
  real(kind=r_bl), allocatable :: fm(:, :, :)
  !                    512      stability function for momentum
  real(kind=r_bl), allocatable :: weight1d(:, :, :)
  !                    513      weighting applied to 1D BL scheme in Smag
  !                             blending
  real(kind=r_bl), allocatable :: slvar(:, :, :)
  !                    516      diagnostic of sl variance
  real(kind=r_bl), allocatable :: qwvar(:, :, :)
  !                    517      diagnostic of qw variance
  real(kind=r_bl), allocatable :: slqw(:, :, :)
  !                    518      diagnostic of sl-qw covariance
  real(kind=r_bl), allocatable :: dissip(:, :, :)
  !                    519      diagnostic of dissipation rate
  integer(kind=integer_64), allocatable :: kblendtq(:, :)
  !                    581      coupling level (T/Q level)
  integer(kind=integer_64), allocatable :: kblenduv(:, :)
  !                    582      coupling level (U/V level)
  real(kind=r_bl), allocatable :: grad_ftl(:,:,:)
  !                    800      explicit grad_ftl
  real(kind=r_bl), allocatable :: non_grad_ftl(:,:,:)
  !                    801      explicit non_grad_ftl
  real(kind=r_bl), allocatable :: ftl_e(:,:,:)
  !                    802      explicit entrainment flux
  real(kind=r_bl), allocatable :: rhokhz_ex(:,:,:)
  !                    803      rhokhz from ex_flux_tq
  real(kind=r_bl), allocatable :: grad_t_adj(:,:)
  !                    804      grad_t_adj from ex_flux_tq

end type strnewbldiag

type (Strnewbldiag) :: BL_diag

character(len=*), parameter, private :: ModuleName = 'BL_DIAGS_MOD'
! ----------------------------------------------------------------------
contains

! allocation of variables for the explicit part of BL code
! ordered by UM stash code
subroutine alloc_bl_expl(bl_diag, l_apply_diag)

use atm_fields_bounds_mod, only: pdims, pdims_s
use bl_option_mod, only: zero
use model_domain_mod, only: model_type, mt_single_column, mt_lfric
use nlsizes_namelist_mod, only: bl_levels
use stash_array_mod, only: sf
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

type (strnewbldiag), intent(in out) :: bl_diag
logical, intent(in) :: l_apply_diag

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

! Local variables
integer :: i, j, k

character(len=*), parameter :: RoutineName='ALLOC_BL_EXPL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set the logical switches for whether diagnostic is requested
select case (model_type)
case DEFAULT

  ! counter gradient term for u
  BL_diag%l_rhogamu      = (l_apply_diag .and. sf(130,3))
  ! counter gradient term for v
  BL_diag%l_rhogamv      = (l_apply_diag .and. sf(131,3))
  ! counter gradient term for t
  BL_diag%l_rhogamt      = (l_apply_diag .and. sf(132,3))
  ! counter gradient term for q
  BL_diag%l_rhogamq      = (l_apply_diag .and. sf(133,3))
  ! mixing length
  BL_diag%l_elm          = (l_apply_diag .and. sf(134,3))
  ! production rate of TKE by shear
  BL_diag%l_tke_shr_prod = (l_apply_diag .and. sf(135,3))
  ! production rate of TKE by buoyancy
  BL_diag%l_tke_boy_prod = (l_apply_diag .and. sf(136,3))
  ! dissipation rate of TKE
  BL_diag%l_tke_dissp    = (l_apply_diag .and. sf(137,3))
  ! non-dimensional diffusion coef. for u, v
  BL_diag%l_sm           = (l_apply_diag .and. sf(138,3))
  ! non-dimensional diffusion coef. for t, q
  BL_diag%l_sh           = (l_apply_diag .and. sf(139,3))
  ! non-gradient buoyancy flux
  BL_diag%l_wb_ng        = (l_apply_diag .and. sf(140,3))
  ! cloud fraction used in the TKE schemes
  BL_diag%l_cf_trb       = (l_apply_diag .and. sf(141,3))
  ! condensed water used in the TKE schemes
  BL_diag%l_ql_trb       = (l_apply_diag .and. sf(142,3))
  ! standard deviation of the distribution function in the TKE schemes
  BL_diag%l_sgm_trb      = (l_apply_diag .and. sf(143,3))
  ! Top of surface mixed layer (Ksurf profile)
  BL_diag%l_smltop       = (l_apply_diag .and. sf(356,3))
  ! Top of decoupled stratocu layer
  BL_diag%l_dsctop       = (l_apply_diag .and. sf(357,3))
  ! BL depth diagnosed from Ri>RiCrit
  BL_diag%l_zhlocal      = (l_apply_diag .and. sf(358,3))
  ! Height of diagnosis parcel top
  BL_diag%l_zhpar        = (l_apply_diag .and. sf(359,3))
  ! Decoupled stratocu base height
  BL_diag%l_zht          = (l_apply_diag .and. sf(304,3))
  ! Turbulent mixing height
  BL_diag%l_dscbase      = (l_apply_diag .and. sf(360,3))
  ! BL cloud base height
  BL_diag%l_cldbase      = (l_apply_diag .and. sf(361,3))
  ! Entrainment rate
  BL_diag%l_weparm       = (l_apply_diag .and. sf(362,3))
  ! Entrainment rate for decoupled stratocu
  BL_diag%l_weparm_dsc   = (l_apply_diag .and. sf(363,3))
  ! inversion thickness
  BL_diag%l_dzh          = (l_apply_diag .and. sf(364,3))
  ! Obukhov length, also required for gustiness diagnostic (463, 515)
  BL_diag%l_oblen        = (l_apply_diag .and.                                 &
                           (sf(464,3) .or. sf(463,3) .or. sf(515,3)))
  ! Friction velocity, also required for gustiness diagnostic (463, 515)
  BL_diag%l_ustar        = (l_apply_diag .and.                                 &
                           (sf(465,3) .or. sf(463,3) .or. sf(515,3)))
  ! Convective velocity scale
  BL_diag%l_wstar        = (l_apply_diag .and. sf(466,3))
  ! Surface buoyancy flux
  BL_diag%l_wbsurf       = (l_apply_diag .and. sf(467,3))
  ! Stratification
  BL_diag%l_dbdz         = (l_apply_diag .and. sf(469,3))
  ! Modulus of shear
  BL_diag%l_dvdzm        = (l_apply_diag .and. sf(470,3))
  ! Momentum diffusivity
  BL_diag%l_rhokm        = (l_apply_diag .and. sf(471,3))
  ! Thermal diffusivity
  BL_diag%l_rhokh        = (l_apply_diag .and. sf(472,3))
  ! x component of orographic stress
  BL_diag%l_ostressx     = (l_apply_diag .and. sf(474,3))
  ! y component of orographic stress
  BL_diag%l_ostressy     = (l_apply_diag .and. sf(475,3))
  ! local mixing length for momentum
  BL_diag%l_elh3d        = (l_apply_diag .and. sf(502,3))
  ! local momentum diffusion coefficient
  BL_diag%l_rhokmloc     = (l_apply_diag .and. sf(503,3))
  ! local scalar diffusion coefficient
  BL_diag%l_rhokhloc     = (l_apply_diag .and. sf(504,3))
  ! surface driven momentum diffusion coefficient
  BL_diag%l_rhokmsurf    = (l_apply_diag .and. sf(505,3))
  ! surface driven scalar diffusion coefficient
  BL_diag%l_rhokhsurf    = (l_apply_diag .and. sf(506,3))
  ! stratocu-top-driven momentum diffusion coefficient
  BL_diag%l_rhokmsc      = (l_apply_diag .and. sf(507,3))
  ! stratocu-top-driven scalar diffusion coefficient
  BL_diag%l_rhokhsc      = (l_apply_diag .and. sf(508,3))
  ! stability function for scalars
  BL_diag%l_fh           = (l_apply_diag .and. sf(511,3))
  ! stability function for momentum
  BL_diag%l_fm           = (l_apply_diag .and. sf(512,3))
  ! weighting applied to 1D BL scheme in Smag blending
  ! also needed for scale aware gust diag (515)
  BL_diag%l_weight1d     = (l_apply_diag .and. (sf(513,3) .or. sf(515,3)))
  BL_diag%l_slvar        = (l_apply_diag .and. sf(516,3))
  BL_diag%l_qwvar        = (l_apply_diag .and. sf(517,3))
  BL_diag%l_slqw         = (l_apply_diag .and. sf(518,3))
  BL_diag%l_dissip       = (l_apply_diag .and. sf(519,3))
  BL_diag%l_kblendtq     = (l_apply_diag .and. sf(581,3))
  BL_diag%l_kblenduv     = (l_apply_diag .and. sf(582,3))
  BL_diag%l_grad_ftl     = (l_apply_diag .and. sf(800,3))
  ! explicit grad_ftl
  BL_diag%l_non_grad_ftl = (l_apply_diag .and. sf(801,3))
  ! explicit non_grad_ftl
  BL_diag%l_ftl_e        = (l_apply_diag .and. sf(802,3))
  ! explicit entrainmnet flux
  BL_diag%l_rhokhz_ex     = (l_apply_diag .and. sf(803,3))
  ! rhokhz from ex_flux_tq
  BL_diag%l_grad_t_adj   = (l_apply_diag .and. sf(804,3))
  ! grad_t_adj from ex_flux_tq

case (mt_single_column)
  BL_diag%l_zht          = .true.
  BL_diag%l_tke_shr_prod = .true.
  BL_diag%l_tke_boy_prod = .true.
  BL_diag%l_tke_dissp    = .true.
  BL_diag%l_slvar        = .true.
  BL_diag%l_qwvar        = .true.
  BL_diag%l_slqw         = .true.

case (mt_lfric)
  !Empty

end select


!       Allocate space for those BL diagnostic arrays required and zero
!       the elements explicitly

if (BL_diag%l_rhogamu) then
  allocate(BL_diag%rhogamu(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhogamu(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhogamu(1,1,1))
end if

if (BL_diag%l_rhogamv) then
  allocate(BL_diag%rhogamv(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhogamv(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhogamv(1,1,1))
end if

if (BL_diag%l_rhogamt) then
  allocate(BL_diag%rhogamt(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhogamt(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhogamt(1,1,1))
end if

if (BL_diag%l_rhogamq) then
  allocate(BL_diag%rhogamq(pdims%i_start:pdims%i_end,                          &
                          pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhogamq(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhogamq(1,1,1))
end if

if (BL_diag%l_elm) then
  allocate(BL_diag%elm(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%elm(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%elm(1,1,1))
end if

if (BL_diag%l_tke_shr_prod) then
  allocate(BL_diag%tke_shr_prod(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%tke_shr_prod(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%tke_shr_prod(1,1,1))
end if

if (BL_diag%l_tke_boy_prod) then
  allocate(BL_diag%tke_boy_prod(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%tke_boy_prod(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%tke_boy_prod(1,1,1))
end if

if (BL_diag%l_tke_dissp) then
  allocate(BL_diag%tke_dissp(pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%tke_dissp(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%tke_dissp(1,1,1))
end if

if (BL_diag%l_sm) then
  allocate(BL_diag%sm(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%sm(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%sm(1,1,1))
end if

if (BL_diag%l_sh) then
  allocate(BL_diag%sh(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%sh(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%sh(1,1,1))
end if

if (BL_diag%l_wb_ng) then
  allocate(BL_diag%wb_ng(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%wb_ng(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%wb_ng(1,1,1))
end if

if (BL_diag%l_cf_trb) then
  allocate(BL_diag%cf_trb(pdims%i_start:pdims%i_end,                           &
                          pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%cf_trb(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%cf_trb(1,1,1))
end if

if (BL_diag%l_ql_trb) then
  allocate(BL_diag%ql_trb(pdims%i_start:pdims%i_end,                           &
                          pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%ql_trb(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%ql_trb(1,1,1))
end if

if (BL_diag%l_sgm_trb) then
  allocate(BL_diag%sgm_trb(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%sgm_trb(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%sgm_trb(1,1,1))
end if

if (BL_diag%l_smltop) then
  allocate(BL_diag%smltop(                                                     &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%smltop(1,1))
end if
if (BL_diag%l_dsctop) then
  allocate(BL_diag%dsctop(                                                     &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%dsctop(1,1))
end if
if (BL_diag%l_zhlocal) then
  allocate(BL_diag%zhlocal(                                                    &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%zhlocal(1,1))
end if
if (BL_diag%l_zhpar) then
  allocate(BL_diag%zhpar(                                                      &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%zhpar(1,1))
end if
if (BL_diag%l_zht) then
  allocate(BL_diag%zht(                                                        &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%zht(1,1))
end if
if (BL_diag%l_dscbase) then
  allocate(BL_diag%dscbase(                                                    &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%dscbase(1,1))
end if
if (BL_diag%l_cldbase) then
  allocate(BL_diag%cldbase(                                                    &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%cldbase(1,1))
end if
if (BL_diag%l_weparm) then
  allocate(BL_diag%weparm(                                                     &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%weparm(1,1))
end if
if (BL_diag%l_weparm_dsc) then
  allocate(BL_diag%weparm_dsc(                                                 &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%weparm_dsc(1,1))
end if
if (BL_diag%l_dzh) then
  allocate(BL_diag%dzh(                                                        &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%dzh(1,1))
end if

if (BL_diag%l_oblen) then
  allocate(BL_diag%oblen(                                                      &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%oblen(1,1))
end if
if (BL_diag%l_ustar) then
  allocate(BL_diag%ustar(                                                      &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%ustar(1,1))
end if
if (BL_diag%l_wstar) then
  allocate(BL_diag%wstar(                                                      &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  BL_diag%wstar(:,:) = zero
else
  allocate(BL_diag%wstar(1,1))
end if
if (BL_diag%l_wbsurf) then
  allocate(BL_diag%wbsurf(                                                     &
       pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  allocate(BL_diag%wbsurf(1,1))
end if
if (BL_diag%l_dbdz) then
  allocate(BL_diag%dbdz(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%dbdz(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%dbdz(1,1,1))
end if
if (BL_diag%l_dvdzm) then
  allocate(BL_diag%dvdzm(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%dvdzm(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%dvdzm(1,1,1))
end if
if (BL_diag%l_rhokm) then
  allocate(BL_diag%rhokm(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%rhokm(1,1,1))
end if
if (BL_diag%l_rhokh) then
  allocate(BL_diag%rhokh(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%rhokh(1,1,1))
end if
if (BL_diag%l_ostressx) then
  allocate(BL_diag%ostressx(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%ostressx(1,1,1))
end if
if (BL_diag%l_ostressy) then
  allocate(BL_diag%ostressy(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%ostressy(1,1,1))
end if
if (BL_diag%l_elh3d) then
  allocate(BL_diag%elh3d(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%elh3d(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%elh3d(1,1,1))
end if
if (BL_diag%l_rhokmloc) then
  allocate(BL_diag%rhokmloc(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%rhokmloc(1,1,1))
end if
if (BL_diag%l_rhokhloc) then
  allocate(BL_diag%rhokhloc(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%rhokhloc(1,1,1))
end if
if (BL_diag%l_rhokmsurf) then
  allocate(BL_diag%rhokmsurf(pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhokmsurf(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhokmsurf(1,1,1))
end if
if (BL_diag%l_rhokhsurf) then
  allocate(BL_diag%rhokhsurf(pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhokhsurf(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhokhsurf(1,1,1))
end if
if (BL_diag%l_rhokmsc) then
  allocate(BL_diag%rhokmsc(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhokmsc(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhokmsc(1,1,1))
end if
if (BL_diag%l_rhokhsc) then
  allocate(BL_diag%rhokhsc(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhokhsc(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhokhsc(1,1,1))
end if
if (BL_diag%l_fh) then
  allocate(BL_diag%fh(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%fh(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%fh(1,1,1))
end if
if (BL_diag%l_fm) then
  allocate(BL_diag%fm(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%fm(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%fm(1,1,1))
end if
if (BL_diag%l_weight1d) then
  allocate(BL_diag%weight1d(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
else
  allocate(BL_diag%weight1d(1,1,1))
end if
if (BL_diag%l_dissip) then
  allocate(BL_diag%dissip(pdims%i_start:pdims%i_end,                           &
                          pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%dissip(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%dissip(1,1,1))
end if
if (BL_diag%l_slvar) then
  allocate(BL_diag%slvar(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%slvar(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%slvar(1,1,1))
end if
if (BL_diag%l_qwvar) then
  allocate(BL_diag%qwvar(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%qwvar(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%qwvar(1,1,1))
end if
if (BL_diag%l_slqw) then
  allocate(BL_diag%slqw(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%slqw(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%slqw(1,1,1))
end if
if (BL_diag%l_kblendtq) then
  allocate(BL_diag%kblendtq(pdims%i_start:pdims%i_end,                         &
                         pdims%j_start:pdims%j_end))
  BL_diag%kblendtq(:,:) = 0
else
  allocate(BL_diag%kblendtq(1,1))
end if
if (BL_diag%l_kblenduv) then
  allocate(BL_diag%kblenduv(pdims_s%i_start:pdims_s%i_end,                     &
                         pdims_s%j_start:pdims_s%j_end))
  BL_diag%kblenduv(:,:) = 0
else
  allocate(BL_diag%kblenduv(1,1))
end if
if (BL_diag%l_grad_ftl) then
  allocate(BL_diag%grad_ftl(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%grad_ftl(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%grad_ftl(1,1,1))
end if
if (BL_diag%l_non_grad_ftl) then
  allocate(BL_diag%non_grad_ftl(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%non_grad_ftl(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%non_grad_ftl(1,1,1))
end if
if (BL_diag%l_ftl_e) then
  allocate(BL_diag%ftl_e(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%ftl_e(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%ftl_e(1,1,1))
end if
if (BL_diag%l_rhokhz_ex) then
  allocate(BL_diag%rhokhz_ex(pdims%i_start:pdims%i_end,                        &
                            pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%rhokhz_ex(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%rhokhz_ex(1,1,1))
end if
if (BL_diag%l_grad_t_adj) then
  allocate(BL_diag%grad_t_adj(pdims%i_start:pdims%i_end,                       &
                            pdims%j_start:pdims%j_end))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j)                                                             &
!$OMP SHARED(BL_diag,pdims)
  do j = pdims%j_start,pdims%j_end
    do i = pdims%i_start,pdims%i_end
      BL_diag%grad_t_adj(i,j) = zero
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%grad_t_adj(1,1))
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine alloc_bl_expl

! deallocation of variables from explicit part of BL code
! in reverse order from which they were allocated
subroutine dealloc_bl_expl(bl_diag)

implicit none

type (strnewbldiag), intent(in out) :: bl_diag

deallocate(BL_diag%weight1d)
deallocate(BL_diag%fm)
deallocate(BL_diag%fh)
deallocate(BL_diag%rhokhsc)
deallocate(BL_diag%rhokmsc)
deallocate(BL_diag%rhokhsurf)
deallocate(BL_diag%rhokmsurf)
deallocate(BL_diag%rhokhloc)
deallocate(BL_diag%rhokmloc)
deallocate(BL_diag%elh3d)
deallocate(BL_diag%ostressy)
deallocate(BL_diag%ostressx)
deallocate(BL_diag%rhokh)
deallocate(BL_diag%rhokm)
deallocate(BL_diag%dvdzm)
deallocate(BL_diag%dbdz)
deallocate(BL_diag%wbsurf)
deallocate(BL_diag%wstar)
deallocate(BL_diag%ustar)
deallocate(BL_diag%oblen)
deallocate(BL_diag%dzh)
deallocate(BL_diag%weparm_dsc)
deallocate(BL_diag%weparm)
deallocate(BL_diag%cldbase)
deallocate(BL_diag%dscbase)
deallocate(BL_diag%zht)
deallocate(BL_diag%zhpar)
deallocate(BL_diag%zhlocal)
deallocate(BL_diag%dsctop)
deallocate(BL_diag%smltop)
deallocate(BL_diag%sgm_trb)
deallocate(BL_diag%ql_trb)
deallocate(BL_diag%cf_trb)
deallocate(BL_diag%wb_ng)
deallocate(BL_diag%sh)
deallocate(BL_diag%sm)
deallocate(BL_diag%tke_dissp)
deallocate(BL_diag%tke_boy_prod)
deallocate(BL_diag%tke_shr_prod)
deallocate(BL_diag%elm)
deallocate(BL_diag%rhogamq)
deallocate(BL_diag%rhogamt)
deallocate(BL_diag%rhogamv)
deallocate(BL_diag%rhogamu)
deallocate(BL_diag%dissip)
deallocate(BL_diag%slvar)
deallocate(BL_diag%qwvar)
deallocate(BL_diag%slqw)
deallocate(BL_diag%kblendtq)
deallocate(BL_diag%kblenduv)
deallocate(BL_diag%grad_ftl)
deallocate(BL_diag%non_grad_ftl)
deallocate(BL_diag%ftl_e)
deallocate(BL_diag%rhokhz_ex)
deallocate(BL_diag%grad_t_adj)

return

end subroutine dealloc_bl_expl

! allocation of variables for the explicit part of BL code
! ordered by UM stash code
subroutine alloc_bl_imp(bl_diag, l_apply_diag, nscmdpkgs, l_scmdiags)

use atm_fields_bounds_mod, only: pdims
use bl_option_mod, only: fric_heating, off, zero
use model_domain_mod, only: model_type, mt_single_column
use nlsizes_namelist_mod, only: bl_levels
use s_scmop_mod, only: scmdiag_bl, scmdiag_incs
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim


implicit none

type (strnewbldiag), intent(in out) :: bl_diag
logical, intent(in) :: l_apply_diag

! Local variables
integer :: i, j, k

! SCM diagnostic switches, used to set bl_diag switches in this routine
integer, intent(in) :: nscmdpkgs
logical, intent(in) :: l_scmdiags(nscmdpkgs)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='ALLOC_BL_IMP'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set the logical switches for whether diagnostic is requested
BL_diag%l_ftl   = .true.
BL_diag%l_fqw   = .true.

if (BL_diag%l_dtfric) then
  allocate(BL_diag%dTfric(pdims%i_start:pdims%i_end,                           &
                          pdims%j_start:pdims%j_end,bl_levels))
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(BL_diag,bl_levels,pdims)
  do k = 1, bl_levels
    do j = pdims%j_start,pdims%j_end
      do i = pdims%i_start,pdims%i_end
        BL_diag%dTfric(i,j,k) = zero
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  allocate(BL_diag%dTfric(1,1,1))
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine alloc_bl_imp

! deallocation of variables from implicit part of BL code
! in reverse order from which they were allocated
subroutine dealloc_bl_imp(bl_diag)


implicit none

type (strnewbldiag), intent(in out) :: bl_diag

deallocate (BL_diag%cff_incr)
deallocate (BL_diag%cfl_incr)
deallocate (BL_diag%cf_incr)
deallocate (BL_diag%qcf_incr)
deallocate (BL_diag%qcl_incr)
deallocate (BL_diag%q_incr)
deallocate (BL_diag%t_incr)
deallocate(BL_diag%dTfric)

return

end subroutine dealloc_bl_imp

end module bl_diags_mod
