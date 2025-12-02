! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module ls_cld_mod

use ls_cld_c_mod, only: ls_cld_c

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LS_CLD_MOD'

contains
!  Large-scale Cloud Scheme.
! Subroutine Interface:
subroutine ls_cld(                                                             &
!      Pressure related fields
 p_theta_levels, rhcrit,                                                       &
!      Array dimensions
 levels, bl_levels,                                                            &
 rhc_row_length,rhc_rows,                                                      &
!      From convection diagnosis (only used if A05_4A)
 ntml, cumulus, l_mixing_ratio,                                                &
!      Prognostic Fields
 t, cf, q, qcf, qcl,                                                           &
!      Liquid and frozen ice cloud fractions
 cfl, cff,                                                                     &
 i_err)

use cv_run_mod,            only: l_param_conv
use vectlib_mod,           only: powr_v
use conversions_mod,       only: pi
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use cloud_inputs_mod,      only: cloud_fraction_method,                        &
 ice_fraction_method, i_eacf, overlap_ice_liquid, ctt_weight,                  &
 t_weight, qsat_fixed, sub_cld, smith_orig, cloud_top_temp,                    &
 min_liq_overlap, all_clouds, not_mixph

use qsat_mod, only: qsat_wat, qsat_wat_mix

use c_cldsgs_mod, only: qcfmin
use missing_data_mod, only: rmdi
implicit none
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme.

! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP No. 29


!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 levels,                                                                       &
!       No. of levels being processed.
   bl_levels,                                                                  &
!       No. of boundary layer levels
   rhc_row_length,rhc_rows

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 qcf(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Cloud ice content at processed levels (kg water per kg air).
   p_theta_levels(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,levels),                           &
!       pressure at all points (Pa).
   rhcrit(rhc_row_length,rhc_rows,levels)
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.

integer ::                                                                     &
 ntml(          tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       in Height of diagnosed BL top

logical ::                                                                     &
 l_mixing_ratio
!       in true if using mixing ratios

logical ::                                                                     &
 cumulus(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       in Logical indicator of convection

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 q(             tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                  (kg water per kg air).
   t(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

real(kind=real_umphys) ::                                                      &
                      !, intent(out)
 cf(            tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Cloud fraction at processed levels (decimal fraction).
   qcl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Cloud liquid water content at processed levels (kg per kg air).
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels)
!       Frozen cloud fraction at processed levels (decimal fraction).

!     Error Status:
integer :: i_err     !, intent(out)  0 if OK; 1 if bad arguments.

!  Local parameters and other physical constants------------------------
real(kind=real_umphys) :: rootwo       ! Sqrt(2.)
real(kind=real_umphys) :: subgrid ! Subgrid parameter in ice cloud calculation

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys) ::                                                      &
 phiqcf,                                                                       &
                      ! Arc-cosine term in Cloud ice fraction calc.
 cosqcf,                                                                       &
                      ! Cosine term in Cloud ice fraction calc.
 overlap_max,                                                                  &
                      ! Maximum possible overlap
 overlap_min,                                                                  &
                      ! Minimum possible overlap
 overlap_random,                                                               &
                      ! Random overlap
 temp0,                                                                        &
 temp1,                                                                        &
 temp2,                                                                        &
                      ! Temporaries for combinations of the
 qn_imp,                                                                       &
 qn_adj
!                       ! overlap parameters

!  (b) Others.
integer :: k,i,j       ! Loop counters: K - vertical level index.
!                                        I,J - horizontal field indices.

integer :: qc_points,                                                          &
                      ! No. points with non-zero cloud
        multrhc   ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  Local dynamic arrays-------------------------------------------------
!    6 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
 qcfrbs(        tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       qCF / bs
   qsl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
   qsl_out(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Saturated specific humidity for temp TL or T.
   qsl_ctt(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
   qsl_ctt_out(   tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Saturated specific humidity wrt liquid at cloud top temperature
   qn(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Cloud water normalised with BS.
   grid_qc(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Gridbox mean saturation excess at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   bs(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Maximum moisture fluctuation /6*sigma at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   ctt(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels)
!       Ice cloud top temperature (K) - as coded it is really TL

real(kind=real_umphys), allocatable :: cfl_max(:,:,:)
!      Maximum value of liquid cloud in a column

logical ::                                                                     &
 lqc(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       True for points with non-zero cloud
integer ::                                                                     &
 idx(tdims%j_len*tdims%i_len,2),                                               &
!       Index for points with non-zero cloud
   llwic(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Last Level With Ice Cloud
real(kind=real_umphys) :: rhcritx
!       Scalar copy of RHCRIT(I,J,K)

!       Variables for cache-blocking
integer            :: jj             !Block index

!       jblock is not a parameter. Value needs to be flexible at runtime.
integer            :: jblock


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_CLD'

!- End of Header

!Set value of blocking factor, jblock. May need to be changed on other
! platforms or for OMP threads > 4
jblock = 4

! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
i_err=0

if ( (rhc_row_length * rhc_rows)  >   1) then
  multrhc = 1
else
  multrhc = 0
end if

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (ice_fraction_method == min_liq_overlap) then
  allocate(cfl_max(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,levels))
else
  allocate(cfl_max(1,1,1))
end if

! Initialize cloud-top-temperature and last-level-with-ice-cloud arrays

! Levels_do1:
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED(llwic, ctt, levels, qcl, cfl, grid_qc, bs, p_theta_levels,       &
!$OMP  l_mixing_ratio, qcf, ice_fraction_method, ctt_weight, i_eacf,           &
!$OMP  rhc_row_length, rhc_rows, bl_levels, cloud_fraction_method, cf,         &
!$OMP  overlap_ice_liquid, cff, t_weight, sub_cld, q, t, jblock, cfl_max,      &
!$OMP  qsat_fixed, multrhc, cumulus, ntml, rhcrit, l_param_conv, tdims)        &
!$OMP  private(k, j, i, rhcritx, qc_points, rootwo, subgrid, qsl, qsl_ctt,     &
!$OMP  phiqcf, cosqcf, qn_imp, qn_adj, overlap_max, overlap_min,               &
!$OMP  overlap_random, temp0, temp1, temp2, qn, lqc, idx, qcfrbs, jj,          &
!$OMP  qsl_out, qsl_ctt_out)

!Cache-blocking applied to loop over j. It is used here to allow
!OpenMP parallelism to be over j (looping over k is order-dependent),
!but still have j and k in the correct order for good cache use.
!If jblock=rows=(tdims%j_len),
!then the inner "do j" loop does the most work, hence
!the loops over j and k are still the correct way around.
if (ice_fraction_method  ==  cloud_top_temp) then
!$OMP do SCHEDULE(DYNAMIC)
  do jj = tdims%j_start, tdims%j_end, jblock
    do j = jj, min((jj+jblock)-1,tdims%j_len)
      do i = tdims%i_start, tdims%i_end
        llwic(i,j)=0
      end do
    end do

    do k = levels, 1, -1
      do j = jj, min((jj+jblock)-1,tdims%j_len)
        do i = tdims%i_start, tdims%i_end

          if (llwic(i,j)  /=  k+1) then
            ctt(i,j,k)=t(i,j,k)
          else
            ctt(i,j,k)=ctt(i,j,k+1)
          end if
          if (qcf(i,j,k)  >   qcfmin) then
            llwic(i,j)=k
          end if
        end do
      end do
    end do
  end do
!$OMP end do
end if

!$OMP do SCHEDULE(DYNAMIC)
do k = 1, levels

  ! ----------------------------------------------------------------------
  ! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
  !    cloud water, sub-grid distribution and fraction arrays.
  !    This requires a preliminary calculation of the pressure.
  !    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
  ! ----------------------------------------------------------------------
  !
  !CDIR COLLAPSE
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      qcl(i,j,k) = 0.0
      cfl(i,j,k) = 0.0
      grid_qc(i,j,k) = rmdi
      bs(i,j,k) = rmdi
    end do ! i
  end do ! j

  if ( l_mixing_ratio ) then
    call qsat_wat_mix(qsl,t(:,:,k),p_theta_levels(:,:,k),                      &
          tdims%i_len,tdims%j_len)
  else
    call qsat_wat(qsl,t(:,:,k),p_theta_levels(:,:,k),                          &
          tdims%i_len,tdims%j_len)
  end if

  do j = tdims%j_start, tdims%j_end

    do i = tdims%i_start, tdims%i_end
      if (multrhc==1) then
        rhcritx = rhcrit(i,j,k)
      else
        rhcritx = rhcrit(1,1,k)
      end if

      ! Omit CUMULUS points below (and including) NTML+1

      if ( .not. l_param_conv .or. (l_param_conv .and.                         &
         (.not. cumulus(i,j) .or. ( cumulus(i,j)                               &
         .and. (k  >   ntml(i,j)+1) ))) ) then

        ! Rhcrit_if:
        if (rhcritx  <   1.0) then
          ! -----------------------------------------------------------------
          ! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
          !    if RHcrit is less than 1
          ! -----------------------------------------------------------------

          qn(i,j) = (q(i,j,k) / qsl(i,j) - 1.0) /                              &
                    (1.0 - rhcritx)

          ! -----------------------------------------------------------------
          ! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
          !    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
          ! -----------------------------------------------------------------

          lqc(i,j) = (qn(i,j)  >   -1.0)
        else
          ! -----------------------------------------------------------------
          ! 2.a Calculate QN = QW - QSL if RHcrit equals 1
          ! -----------------------------------------------------------------

          qn(i,j) = q(i,j,k) - qsl(i,j)

          ! -----------------------------------------------------------------
          ! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
          !     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
          ! -----------------------------------------------------------------

          lqc(i,j) = (qn(i,j)  >   0.0)
        end if ! Rhcrit_if
      else if (l_param_conv) then
        lqc(i,j) = .false.
      end if  ! Test on CUMULUS and NTML for A05_4A only
    end do ! i
  end do ! j

  ! ----------------------------------------------------------------------
  ! 4. Form index of points where non-zero liquid cloud fraction
  ! ----------------------------------------------------------------------

  qc_points=0

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (lqc(i,j)) then
        qc_points = qc_points + 1
        idx(qc_points,1) = i
        idx(qc_points,2) = j
      end if
    end do ! i
  end do ! j

  ! ----------------------------------------------------------------------
  ! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
  !                  water cloud fraction and determine temperature.
  ! ----------------------------------------------------------------------
  ! Qc_points_if:
  if (qc_points  >   0) then
    call ls_cld_c(p_theta_levels(1,1,k),rhcrit(1,1,k),qsl,qn,                  &
                  q(1,1,k),t(1,1,k),                                           &
                  qcl(1,1,k),cfl(1,1,k),grid_qc(1,1,k),bs(1,1,k),              &
                  idx,qc_points,rhc_row_length,rhc_rows,                       &
                  bl_levels,k, l_mixing_ratio)
  end if ! Qc_points_if

end do !k loop
!$OMP end do
  ! ----------------------------------------------------------------------
  ! 6. Calculate cloud fractions for ice clouds.
  !    Begin by calculating Qsat_wat(T,P), at Temp. T, for estimate of bs.
  ! ----------------------------------------------------------------------
if (ice_fraction_method == min_liq_overlap) then

  k = levels
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      cfl_max(i,j,k) = cfl(i,j,k)
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do k = levels-1, 1, -1
      do i = tdims%i_start, tdims%i_end
        if (cfl(i,j,k) > cfl_max(i,j,k+1)) then
          ! Cloud frac has increased over peak, so use this
          cfl_max(i,j,k) = cfl(i,j,k)
        else
          ! Cloud frac has decreased, so keep the same
          cfl_max(i,j,k) = cfl_max(i,j,k+1)
        end if
      end do
    end do
  end do
!$OMP end do
end if

rootwo = sqrt(2.0)

!$OMP do SCHEDULE(DYNAMIC)
do k = 1, levels

  if ( l_mixing_ratio ) then
    call qsat_wat_mix(qsl,t(:,:,k),p_theta_levels(:,:,k),                      &
          tdims%i_len,tdims%j_len)
  else
    call qsat_wat(qsl,t(:,:,k),p_theta_levels(:,:,k),                          &
          tdims%i_len,tdims%j_len)
  end if

  if (ice_fraction_method  ==  cloud_top_temp) then
    ! Use cloud top temperature and a fixed qsat to give QCFRBS

    if ( l_mixing_ratio ) then
      call qsat_wat_mix(qsl_ctt,ctt(:,:,k),p_theta_levels(:,:,k),              &
                            tdims%i_len,tdims%j_len)
    else
      call qsat_wat(qsl_ctt,ctt(:,:,k),p_theta_levels(:,:,k),                  &
                            tdims%i_len,tdims%j_len)
    end if

    call powr_v(                                                               &
       tdims%i_len*tdims%j_len,                                                &
       qsl_ctt,ctt_weight,qsl_ctt_out )
    qsl_ctt=qsl_ctt_out
    call powr_v(                                                               &
       tdims%i_len*tdims%j_len,                                                &
       qsl,t_weight,qsl_out )
    qsl=qsl_out

    subgrid = sub_cld ** (1.0-t_weight)                                        &
              / qsat_fixed ** (1.0-t_weight-ctt_weight)
  end if ! ice_fraction_method eq 2

  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if ( multrhc== 1) then
        rhcritx = rhcrit(i,j,k)
      else
        rhcritx = rhcrit(1,1,k)
      end if
      ! ----------------------------------------------------------------------
      ! 6a Calculate qCF/bs.
      ! ----------------------------------------------------------------------
      ! Rhcrit_if2:
      if (rhcritx  <   1.0) then

        if (ice_fraction_method  ==  smith_orig) then
          qcfrbs(i,j)=  qcf(i,j,k) / ((1.0-rhcritx) * qsl(i,j))
        else if (ice_fraction_method  ==  cloud_top_temp) then
          qcfrbs(i,j) = subgrid * qcf(i,j,k) / ((1.0-rhcritx)                  &
                   * qsl_ctt(i,j)*qsl(i,j))
        else if (ice_fraction_method == min_liq_overlap) then
          ! method described in appendix of Abel et al (2017, JAS)
          qcfrbs(i,j)=  max(1.0 - cfl_max(i,j,k), 0.05) * qcf(i,j,k) /         &
                        ((1.0-rhcritx) * qsl(i,j))
        else
          ! No ice cloud fraction method defined
        end if ! ice_fraction_method

        ! ----------------------------------------------------------------------
        ! 6b Calculate frozen cloud fraction from frozen cloud water content.
        ! ----------------------------------------------------------------------
        if (qcfrbs(i,j)  <=  0.0) then
          cff(i,j,k) = 0.0
        else if (0.0<qcfrbs(i,j) .and. (6.0*qcfrbs(i,j)) <= 1.0) then
          cff(i,j,k) = 0.5 * ((6.0 * qcfrbs(i,j))**(2.0/3.0))
        else if (1.0<(6.0*qcfrbs(i,j)) .and. qcfrbs(i,j) < 1.0) then
          phiqcf = acos(rootwo * 0.75 * (1.0 - qcfrbs(i,j)))
          cosqcf = cos((phiqcf + (4.0 * pi)) / 3.0)
          cff(i,j,k) = 1.0 - (4.0 * cosqcf * cosqcf)
        else if (qcfrbs(i,j)  >=  1.0) then
          cff(i,j,k) = 1.0
        end if
        if (i_eacf == all_clouds .or.                                          &
             (i_eacf == not_mixph .and. cfl(i,j,k) < 0.05) ) then
          ! Empirically adjusted cloud fraction
          ! Back out QN
          if (0.0< qcfrbs(i,j) .and. (6.0*qcfrbs(i,j)) <=  1.0) then
            qn_imp=sqrt(2.0*cff(i,j,k))-1.0
          else if (1.0<(6.0*qcfrbs(i,j)) .and. qcfrbs(i,j)<1.0) then
            qn_imp=1.0-sqrt((1.0-cff(i,j,k))*2.0)
          else
            qn_imp = 1.0
          end if

          ! Modify QN with EACF relationship
          if (k >  bl_levels) then
            qn_adj=(qn_imp+0.0955)/(1.0-0.0955)
          else
            qn_adj=(qn_imp+0.184)/(1.0-0.184)
          end if

          ! Recalculate ice cloud fraction with modified QN
          if (qcfrbs(i,j)  <=  0.0) then
            cff(i,j,k) = 0.0
          else if (qn_adj  <=  0.0) then
            cff(i,j,k) = 0.5 * (1.0 + qn_adj) * (1.0 + qn_adj)
          else if (qn_adj  <   1.0) then
            cff(i,j,k) = 1.0 - 0.5 * (1.0-qn_adj) * (1.0-qn_adj)
          else
            cff(i,j,k) = 1.0
          end if

        end if  ! i_eacf


      else ! RHcrit = 1, set cloud fraction to 1 or 0

        if (qcf(i,j,k)  >   0.0) then
          cff(i,j,k) = 1.0
        else
          cff(i,j,k) = 0.0
        end if

      end if
    end do ! i
  end do ! j

  ! ----------------------------------------------------------------------
  ! 6c Calculate combined cloud fraction.
  ! ----------------------------------------------------------------------

  if (cloud_fraction_method  ==  1) then

    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        !             Use minimum overlap condition
        cf(i,j,k) = min(cfl(i,j,k)+cff(i,j,k), 1.0)
      end do
    end do

  else if (cloud_fraction_method  ==  2) then

      !Hand-unrolled loop
    do j = tdims%j_start, tdims%j_len -                                        &
              mod(tdims%j_len,2), 2

      do i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=min(cfl(i,j,k),cff(i,j,k))
        overlap_min=max(cfl(i,j,k)+cff(i,j,k)-1.0,0.0)
        overlap_random=cfl(i,j,k)*cff(i,j,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,j,k)=cfl(i,j,k)+cff(i,j,k)                                        &
                -(temp0+temp1*overlap_ice_liquid                               &
                +temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,j,k)=min(cf(i,j,k),cfl(i,j,k)+cff(i,j,k))

      end do

      do i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=min(cfl(i,j+1,k),cff(i,j+1,k))
        overlap_min=max(cfl(i,j+1,k)+cff(i,j+1,k)-1.0,0.0)
        overlap_random=cfl(i,j+1,k)*cff(i,j+1,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,j+1,k)=cfl(i,j+1,k)+cff(i,j+1,k)                                  &
                -(temp0+temp1*overlap_ice_liquid                               &
                +temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,j+1,k)=min(cf(i,j+1,k),cfl(i,j+1,k)+cff(i,j+1,k))

      end do
    end do

      !Post-conditioning
    if (mod(tdims%j_len,2) == 1) then
      do i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=min( cfl(i,tdims%j_len,k),                                 &
                         cff(i,tdims%j_len,k) )
        overlap_min=max( cfl(i,tdims%j_len,k) +                                &
                         cff(i,tdims%j_len,k)-1.0,                             &
                         0.0 )
        overlap_random=  cfl(i,tdims%j_len,k) *                                &
                         cff(i,tdims%j_len,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,tdims%j_len,k)=                                                   &
          cfl(i,tdims%j_len,k) +                                               &
          cff(i,tdims%j_len,k) -                                               &
          (temp0+temp1*overlap_ice_liquid +                                    &
          temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,tdims%j_len,k)=                                                   &
          min(cf(i,tdims%j_len,k),                                             &
             cfl(i,tdims%j_len,k) +                                            &
             cff(i,tdims%j_len,k))

      end do
    end if !Post-conditioning

    ! CFF + CFL >= 1 implies CF = 1.
    ! Deal with this case separately to avoid roundoff issues:
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if (cfl(i,j,k)+cff(i,j,k) >= 1.0) cf(i,j,k) = 1.0
      end do
    end do

  else
    ! No total cloud fraction method defined


  end if ! cloud_fraction_method

end do ! Levels_do
!$OMP end do

!$OMP end PARALLEL

deallocate(cfl_max)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_cld
end module ls_cld_mod
