! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module kmkhz_9c_wtrac_mod

use um_types, only: r_bl

implicit none

! Description:
!  Contains water tracer subroutines called from kmkhz_9c, which ultimately
!  calculate the vertical entrainment flux of water tracer across the top of
!  the surface mixed layer (SML) or the decoupled stratocumulus layer (DSC).
!  These subroutines have been created to reduce the amount of water tracer
!  code in kmkhz_9c but should be kept consistent with the water code in
!  kmkhz_9c.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName='KMKHZ_9C_WTRAC_MOD'

contains

subroutine calc_dqw_inv_wtrac(bl_levels, k_inv_start, k_inv, dqw_meth,         &
                              disc_inv, rdz, z_tq, zh, inv_calc,               &
                              qw_lapse_zero, wtrac_bl, dqw_inv_wtrac)

! Description:
!  Calculate total water tracer change across (SML or DSC) inversion
!  (i.e. equivalents of dqw_sml and dqw_dsc in kmkhz_9c)
!

use atm_fields_bounds_mod,   only: pdims, tdims
use free_tracers_inputs_mod, only: n_wtrac, noniso_class_wtrac
use water_tracers_mod,       only: wtrac_info
use wtrac_bl_mod,            only: bl_wtrac_type
use bl_option_mod,           only: zero, one_half

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: bl_levels            ! No. of atmospheric levels for
                                            ! which boundary layer fluxes are
                                            !    calculated.

integer, intent(in) :: k_inv_start(pdims%i_start:pdims%i_end,                  &
                                   pdims%j_start:pdims%j_end)
                                ! Inversion level - starting value
integer, intent(in) :: k_inv(pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end)
                                ! Inversion level - updated value
integer, intent(in) :: disc_inv(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end)
                                ! Discontinuous inversion flag
integer, intent(in) :: dqw_meth(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end)
                                ! Method used in normal water calculation


real(kind=r_bl), intent(in) ::                                                 &
            rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                ! Reciprocal of distance between
                                ! full levels (m-1).  1/RDZ(,K) is
                                ! the vertical distance from level
                                ! K-1 to level K, except that for
                                ! K=1 it is the height of the
                                ! lowest atmospheric full level.
real(kind=r_bl), intent(in) ::                                                 &
           z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                ! Z_tq(*,K) is the height of the
                                ! k-th full level above the surface.
real(kind=r_bl), intent(in) ::                                                 &
             zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! Boundary layer height or cloud layer height

logical, intent(in) :: inv_calc(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end)
                                ! Flag to do calculations or not
                                ! (Always set to T for SML, and equal to array
                                !  dsc for DSC)
logical, intent(in) :: qw_lapse_zero(pdims%i_start:pdims%i_end,                &
                                    pdims%j_start:pdims%j_end)
                                ! Flag to indicate if the qw lapse rate has
                                ! been set to zero

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), intent(in out) :: wtrac_bl(n_wtrac)

real(kind=r_bl), intent(in out) ::                                             &
          dqw_inv_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                         n_wtrac)     ! Change in water tracer qw across
                                      ! inversion
! Local parameters

integer :: i, j, k, i_wt     ! Loop counters

real(kind=r_bl) :: dz_disc   ! height of ZH below Z_uv(NTML+2)
real(kind=r_bl) :: qw_lapse  ! Lapse rate of QW above inversion
                                    !  (kg/kg/m)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter ::  RoutineName = 'CALC_DQW_INV_WTRAC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(SHARED)                            &
!$OMP  private (i, j, k, dz_disc, qw_lapse)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      dqw_inv_wtrac(i,j,i_wt) = zero

      if (inv_calc(i,j)) then   ! do inversion calculation

        k = k_inv_start(i,j)       ! Use starting value
        !..by default, keep ZH at the half-level where it was diagnosed
        !..initially and use grid-level jumps

        dqw_inv_wtrac(i,j,i_wt) = wtrac_bl(i_wt)%qw(i,j,k+1)                   &
                                   - wtrac_bl(i_wt)%qw(i,j,k)

        if (disc_inv(i,j) == 1) then    ! Discontinuous inversion

          !-----------------------------------------------------------
          !..Calculate inversion discontinuous jumps of QW
          !-----------------------------------------------------------
          ! Allow for lapse rate above inversion, if known

          ! Note, k_inv can have been changed at this point, so reset k
          k = k_inv(i,j)
          dz_disc = z_tq(i,j,k+2) - zh(i,j)
          qw_lapse = zero

          if (wtrac_info(i_wt)%wt_class == noniso_class_wtrac) then
            ! Non-isotopic water tracers must follow exactly the same
            ! calculations as done for water in the code above.

            if ( k  <=  bl_levels-3 ) then
              if (qw_lapse_zero(i,j)) then
                ! The water lapse rate has been set to zero so set non-iso
                ! lapse rate to zero as well
                qw_lapse = zero
              else
                ! Calculate non-iso lapse rate and allow it to be positive
                ! to follow the exact same calculations as water
                qw_lapse =                                                     &
               ( wtrac_bl(i_wt)%qw(i,j,k+3) - wtrac_bl(i_wt)%qw(i,j,k+2) )     &
                 *rdz(i,j,k+3)
              end if
            end if

            !-----------------
            ! Next QW jump
            !-----------------

            if (dqw_meth(i,j) > 0) then

              dqw_inv_wtrac(i,j,i_wt) = wtrac_bl(i_wt)%qw(i,j,k+2)             &
                                        - wtrac_bl(i_wt)%qw(i,j,k)

              if (dqw_meth(i,j) == 2) then
                dqw_inv_wtrac(i,j,i_wt) = dqw_inv_wtrac(i,j,i_wt) -            &
                         one_half*dqw_inv_wtrac(i,j,i_wt)
              else if (dqw_meth(i,j) == 3) then
                dqw_inv_wtrac(i,j,i_wt) = dqw_inv_wtrac(i,j,i_wt) -            &
                         qw_lapse*dz_disc
              end if
            end if

          else  ! Other types of water tracers

            ! For 'normal' water tracer and water isotopes, which are
            ! realistic physical fields, the calculations should be based
            ! on the water tracer gradients.

            if ( k  <=  bl_levels-3 ) then
              qw_lapse = min( zero,                                            &
               ( wtrac_bl(i_wt)%qw(i,j,k+3) - wtrac_bl(i_wt)%qw(i,j,k+2) )     &
                 *rdz(i,j,k+3) )
            end if

            !-----------------
            ! Next QW jump
            !-----------------

            if (wtrac_bl(i_wt)%qw(i,j,k+2) < wtrac_bl(i_wt)%qw(i,j,k+1) .and.  &
                wtrac_bl(i_wt)%qw(i,j,k+1) < wtrac_bl(i_wt)%qw(i,j,k) ) then
              ! QW monotonically decreasing across inversion
              ! Only allow for QW lapse rate if both it and the
              ! 2 grid-level jump are negative (expected sign)
              dqw_inv_wtrac(i,j,i_wt) = wtrac_bl(i_wt)%qw(i,j,k+2)             &
                                      - wtrac_bl(i_wt)%qw(i,j,k)
              if ( dqw_inv_wtrac(i,j,i_wt) < zero ) then
                dqw_inv_wtrac(i,j,i_wt) = dqw_inv_wtrac(i,j,i_wt) -            &
                       max( one_half*dqw_inv_wtrac(i,j,i_wt), qw_lapse*dz_disc )
              end if
            else if (                                                          &
              wtrac_bl(i_wt)%qw(i,j,k+2) > wtrac_bl(i_wt)%qw(i,j,k+1) .and.    &
              wtrac_bl(i_wt)%qw(i,j,k+1) > wtrac_bl(i_wt)%qw(i,j,k) ) then
              ! QW monotonically increasing across inversion
              ! Suggests something unusual is going so not clear how
              ! to proceed, so currently leaving DQW as 2 level jump
              dqw_inv_wtrac(i,j,i_wt) = wtrac_bl(i_wt)%qw(i,j,k+2)             &
                                         - wtrac_bl(i_wt)%qw(i,j,k)
              ! else
              ! In this case, either:
              ! a) NTML has been lowered a level and it's possible that
              !    neither of the above two conditions are met for normal
              !    water (as k has changed since monotonic_inv was
              !    calculated).  Water tracer calculations will be the same
              !    as normal water. or
              ! b) qw is monotonic but the water tracers are not.
              !    This means the water tracer calculations will use
              !    the most appropriate gradient calculation given their
              !    vertical profile.
            end if
          end if ! Type of water tracer

        end if  ! disc_inv
      end if     ! Do inversion calculation?
    end do  ! i
  end do   !j
!$OMP end PARALLEL do
end do     ! i_wt


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine calc_dqw_inv_wtrac

! ----------------------------------------------------------------------

subroutine calc_fqw_inv_wtrac(bl_levels, k_inv, totqf_efl_meth1,               &
                              totqf_efl_meth2, t_frac, zh, zh_frac, zrzi,      &
                              z_uv_k, dzl, we_rho, w_ls, we_parm,              &
                              dqw_inv_wtrac, fq_nt_zh_wtrac, moisten,          &
                              inv_type, wtrac_bl)
!
! Description:
!  Calculate the vertical entrainment flux of water tracer across the top of
!  the surface mixed layer (SML) or the decoupled stratocumulus layer (DSC)
!  (i.e. water tracer equivalent of fqw_sml or fqw_dsc at entrainment level).
!
use atm_fields_bounds_mod,   only: pdims, tdims
use timestep_mod,            only: timestep
use free_tracers_inputs_mod, only: n_wtrac, noniso_class_wtrac
use water_tracers_mod,       only: wtrac_info
use wtrac_bl_mod,            only: bl_wtrac_type
use bl_option_mod,           only: zero, one
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: bl_levels
                               ! No. of atmospheric levels for which boundary
                               !  layer fluxes are calculated.

integer, intent(in) :: k_inv(pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end)
                               ! Inversion level - updated value

integer, intent(in) :: totqf_efl_meth1(pdims%i_start:pdims%i_end,              &
                                       pdims%j_start:pdims%j_end)
integer, intent(in) :: totqf_efl_meth2(pdims%i_start:pdims%i_end,              &
                                       pdims%j_start:pdims%j_end)
                               ! Indicators of methods used by normal water
                               !   calculations

real(kind=r_bl), intent(in) :: t_frac(pdims%i_start:pdims%i_end,               &
                                             pdims%j_start:pdims%j_end)
                                ! Fraction of timestep inversion is above
                                !   entr.t flux-level for SML or DSC layers

real(kind=r_bl), intent(in) :: zh(pdims%i_start:pdims%i_end,                   &
                                          pdims%j_start:pdims%j_end)
                                ! Boundary layer height or depth of DSC layer
real(kind=r_bl), intent(in) :: zh_frac(pdims%i_start:pdims%i_end,              &
                                              pdims%j_start:pdims%j_end)
                                ! (ZH-ZHALF)/DZ or (ZHSC-ZHALF)/DZ
real(kind=r_bl), intent(in) :: zrzi(pdims%i_start:pdims%i_end,                 &
                                           pdims%j_start:pdims%j_end)
                                ! (z-z_base)/(z_i-z_base)
real(kind=r_bl), intent(in) :: z_uv_k(pdims%i_start:pdims%i_end,               &
                                             pdims%j_start:pdims%j_end)
                                ! height of level ntml+1-1/2 or DSCDEPTH
real(kind=r_bl), intent(in) :: dzl(tdims%i_start:tdims%i_end,                  &
                                          tdims%j_start:tdims%j_end,bl_levels)
                                ! Layer depths (m).  DZL(,K) is the
                                !    distance from layer boundary K-1/2
                                !    to layer boundary K+1/2.  For K=1
                                !    the lower boundary is the surface.

real(kind=r_bl), intent(in) :: we_rho(pdims%i_start:pdims%i_end,               &
                                             pdims%j_start:pdims%j_end)
                                ! rho*entrainment rate
real(kind=r_bl), intent(in) :: w_ls(pdims%i_start:pdims%i_end,                 &
                                           pdims%j_start:pdims%j_end)
                                ! large-scale (subs) velocity at subgrid
                                ! inversion heights
real(kind=r_bl), intent(in) :: we_parm(pdims%i_start:pdims%i_end,              &
                                              pdims%j_start:pdims%j_end)
                                ! Parametrised entrainment rates (m/s)
                                !   for surf and DSC layers

real(kind=r_bl), intent(in) :: dqw_inv_wtrac                                   &
                                       (pdims%i_start:pdims%i_end,             &
                                        pdims%j_start:pdims%j_end,n_wtrac)
                                ! Water tracer QW changes across disc inv
real(kind=r_bl), intent(in) :: fq_nt_zh_wtrac                                  &
                                       (pdims%i_start:pdims%i_end,             &
                                        pdims%j_start:pdims%j_end,n_wtrac)
                                ! Water tracer FQ_NT at ZH or ZHSC

logical, intent(in) :: moisten(pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end)
                                ! Indicator of whether inversion grid-level
                                !   should moisten this timestep (or dry)

character(len=3), intent(in) :: inv_type  ! Type of inversion, 'SML' or 'DSC'

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), intent(in out) :: wtrac_bl(n_wtrac)

! Local variables

integer :: i,j,i_wt,k        ! Loop counters

real(kind=r_bl) :: ml_tend        ! mixed layer tendency (d/dt)
real(kind=r_bl) :: fa_tend        ! free atmospheric tendency (d/dt)
real(kind=r_bl) :: inv_tend       ! limit on inversion grid-level
                                         ! tendency (d/dt)
real(kind=r_bl) :: totqf_efl      ! Total water tracer QW flux at
                                         ! entrainment flux grid-level

logical :: moisten_wt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! Indicator of whether inversion grid-level
                                !   should moisten this timestep (or dry)
                                !   in terms of water tracers (rather than
                                !   actual water)


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter ::  RoutineName = 'CALC_FQW_INV_WTRAC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(SHARED)                            &
!$OMP  private (i, j, k, ml_tend, fa_tend, inv_tend, totqf_efl)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      if (inv_type == 'SML') then
        wtrac_bl(i_wt)%totqf_zh(i,j) = zero
      else if (inv_type == 'DSC') then
        wtrac_bl(i_wt)%totqf_zhsc(i,j) = zero
      end if

      k = k_inv(i,j)+1

      if ( t_frac(i,j)  >   zero ) then

        if (inv_type == 'SML') then
          ! Calculate total (turb+micro+subs) QW flux at subgrid
          ! inversion height
          wtrac_bl(i_wt)%totqf_zh(i,j) =                                       &
                                     - we_rho(i,j)*dqw_inv_wtrac(i,j,i_wt)     &
                                     + fq_nt_zh_wtrac(i,j,i_wt)
          ! Interpolate to entrainment flux-level below
          totqf_efl = wtrac_bl(i_wt)%fq_nt(i,j,1) + wtrac_bl(i_wt)%fqw(i,j,1)  &
                      +  zrzi(i,j) * ( wtrac_bl(i_wt)%totqf_zh(i,j)            &
                       - wtrac_bl(i_wt)%fq_nt(i,j,1)                           &
                       - wtrac_bl(i_wt)%fqw(i,j,1) )

          ml_tend = - (  wtrac_bl(i_wt)%totqf_zh(i,j)                          &
                      - wtrac_bl(i_wt)%fq_nt(i,j,1)                            &
                      - wtrac_bl(i_wt)%fqw(i,j,1) ) /zh(i,j)

        else if (inv_type == 'DSC') then
          ! Calculate total (turb+micro) QW flux at subgrid inversion
          wtrac_bl(i_wt)%totqf_zhsc(i,j) =                                     &
                            - we_rho(i,j)*dqw_inv_wtrac(i,j,i_wt)              &
                               + fq_nt_zh_wtrac(i,j,i_wt)
          ! Interpolate to entrainment flux-level
          totqf_efl = wtrac_bl(i_wt)%fq_nt_dscb(i,j) +                         &
              (wtrac_bl(i_wt)%totqf_zhsc(i,j)-wtrac_bl(i_wt)%fq_nt_dscb(i,j))  &
                              * zrzi(i,j)
          ! Note, here zh = dscdepth
          ml_tend = - ( wtrac_bl(i_wt)%totqf_zhsc(i,j)                         &
                         - wtrac_bl(i_wt)%fq_nt_dscb(i,j) )/zh(i,j)

        end if

        ! Need to ensure the total QW flux gradient in inversion
        ! grid-level is consistent with inversion rising or falling.
        ! If QW(K) is drier than mixed layer then inversion rising
        ! implies moistening in level K relative to mixed layer
        ! while falling would imply relative drying of level K.
        ! If QW(K) is moister than ML then want opposite tendencies.
        fa_tend = zero
        if ( k+1  <=  bl_levels )                                              &
          fa_tend = - ( wtrac_bl(i_wt)%fq_nt(i,j,k+2)                          &
                           - wtrac_bl(i_wt)%fq_nt(i,j,k+1) )                   &
                             / dzl(i,j,k+1)
        inv_tend =       zh_frac(i,j) * ml_tend                                &
                             + (one-zh_frac(i,j)) * fa_tend

        if (wtrac_info(i_wt)%wt_class == noniso_class_wtrac) then
          ! Non-isotopic water tracers must follow exactly the same
          ! calculations as done for water in the kmkhz_9c code.

          if ( moisten(i,j) ) then
            ! Ensure inversion level does moisten relative to ML
            if ( totqf_efl_meth1(i,j) == 1) then
              totqf_efl = wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)
            end if
            if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
              ! Ensure inversion level won't end up more moist than
              ! K_INV by end of timestep.
              ! Set INV_TEND to max allowable moistening rate, also
              ! allowing for change in ML_TEND arising from this change
              ! to TOTQF_EFL:
              if (inv_type == 'SML') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                           - wtrac_bl(i_wt)%qw(i,j,k))/timestep                &
                           + (wtrac_bl(i_wt)%fq_nt(i,j,1) +                    &
                              wtrac_bl(i_wt)%fqw(i,j,1))/z_uv_k(i,j)
              else if (inv_type == 'DSC') then
                ! Note, here z_uv_k = dscdepth
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                           - wtrac_bl(i_wt)%qw(i,j,k))/timestep                &
                           + wtrac_bl(i_wt)%fq_nt_dscb(i,j)/z_uv_k(i,j)
              end if
              if ( totqf_efl_meth2(i,j) == 1) then
                totqf_efl =                                                    &
                        (wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))    &
                         /(one+ dzl(i,j,k)/z_uv_k(i,j))
              end if
            end if
          else
            if (totqf_efl_meth1(i,j) == 1) then
              totqf_efl = wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)
            end if
            if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
              ! Ensure inversion level won't end up drier than
              ! K_INV by end of timestep.
              ! Set INV_TEND to max allowable drying rate:
              if (inv_type == 'SML') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                            - wtrac_bl(i_wt)%qw(i,j,k))/timestep               &
                            + (wtrac_bl(i_wt)%fq_nt(i,j,1 )                    &
                               + wtrac_bl(i_wt)%fqw(i,j,1))/z_uv_k(i,j)
              else if (inv_type == 'DSC') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                            - wtrac_bl(i_wt)%qw(i,j,k))/timestep               &
                            + wtrac_bl(i_wt)%fq_nt_dscb(i,j)/z_uv_k(i,j)
              end if

              if (totqf_efl_meth2(i,j) == 1) then
                totqf_efl =                                                    &
                       (wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))     &
                        /(one+ dzl(i,j,k)/z_uv_k(i,j))
              end if
            end if
          end if   ! moisten

        else

          ! For the 'normal' water tracer and water isotopes, which are
          ! realistic physical fields, the calculations are based
          ! on the water tracer gradients.

          if (we_parm(i,j)+w_ls(i,j) >=  zero) then
            ! inversion moving up so inversion will moisten/dry
            ! depending on relative QW in level below
            moisten_wt(i,j) =                                                  &
                (wtrac_bl(i_wt)%qw(i,j,k) <= wtrac_bl(i_wt)%qw(i,j,k-1))
          else
            ! inversion moving down so inversion will moisten/dry
            ! depending on relative QW in level above
            moisten_wt(i,j) =                                                  &
                (wtrac_bl(i_wt)%qw(i,j,k) <= wtrac_bl(i_wt)%qw(i,j,k+1))
          end if

          if ( moisten_wt(i,j) ) then
            ! Ensure inversion level does moisten relative to ML
            totqf_efl = max( totqf_efl,                                        &
                        wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
            if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
              ! Ensure inversion level won't end up more moist than
              ! K_INV by end of timestep.
              ! Set INV_TEND to max allowable moistening rate, also
              ! allowing for change in ML_TEND arising from this change
              ! to TOTQF_EFL:
              if (inv_type == 'SML') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                           - wtrac_bl(i_wt)%qw(i,j,k))/timestep                &
                           + (wtrac_bl(i_wt)%fq_nt(i,j,1) +                    &
                                 wtrac_bl(i_wt)%fqw(i,j,1))/z_uv_k(i,j)
              else if (inv_type == 'DSC') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                           - wtrac_bl(i_wt)%qw(i,j,k))/timestep                &
                           + wtrac_bl(i_wt)%fq_nt_dscb(i,j)/z_uv_k(i,j)
              end if
              totqf_efl = min( totqf_efl,                                      &
                        (wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))    &
                         /(one+ dzl(i,j,k)/z_uv_k(i,j))   )
            end if
          else
            totqf_efl = min( totqf_efl,                                        &
                           wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
            if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
              ! Ensure inversion level won't end up drier than
              ! K_INV by end of timestep.
              ! Set INV_TEND to max allowable drying rate:
              if (inv_type == 'SML') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                            - wtrac_bl(i_wt)%qw(i,j,k))/timestep               &
                            + (wtrac_bl(i_wt)%fq_nt(i,j,1 )                    &
                               + wtrac_bl(i_wt)%fqw(i,j,1))/z_uv_k(i,j)
              else if (inv_type == 'DSC') then
                inv_tend = (wtrac_bl(i_wt)%qw(i,j,k-1)                         &
                            - wtrac_bl(i_wt)%qw(i,j,k))/timestep               &
                            + wtrac_bl(i_wt)%fq_nt_dscb(i,j)/z_uv_k(i,j)
              end if

              totqf_efl = max( totqf_efl,                                      &
                         (wtrac_bl(i_wt)%fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))   &
                          /(one+ dzl(i,j,k)/z_uv_k(i,j))   )
            end if
          end if   ! moisten_wt

        end if  ! Type of water tracer

        wtrac_bl(i_wt)%fqw(i,j,k) = t_frac(i,j) *                              &
                       ( totqf_efl - wtrac_bl(i_wt)%fq_nt(i,j,k) )

      end if  ! t_frac

    end do
  end do
!$OMP end PARALLEL do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine calc_fqw_inv_wtrac

end module kmkhz_9c_wtrac_mod
