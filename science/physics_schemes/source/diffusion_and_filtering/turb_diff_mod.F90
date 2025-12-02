! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables for turbulent diffusion input
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: diffusion_and_filtering

module turb_diff_mod

! Subgrid turbulence scheme input.

use umPrintMgr, only: umprint, newline
use atmos_max_sizes, only: model_levels_max
use Control_Max_Sizes, only: max_121_rows, max_sponge_width,                   &
  max_updiff_levels
use missing_data_mod, only: imdi, rmdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

use um_types, only: real_umphys

implicit none

logical :: L_diff_active
logical :: L_subfilter_horiz = .false. ! activate horiz diffusion
logical :: L_subfilter_vert  = .false. ! activate vertical diffusion
logical :: L_blend_isotropic = .true.
                  ! Use same diffusion coefficient in all directions

integer :: turb_startlev_horiz = imdi  ! 1st lev for horiz subgrid turb
integer :: turb_endlev_horiz   = imdi  ! last lev for horiz subgrid turb
integer :: turb_startlev_vert  = imdi  ! 1st lev for vert subgrid turb
integer :: turb_endlev_vert    = imdi  ! last lev for vert subgrid turb
integer :: hyd_mix_opt         = imdi  ! mixing option for hydrometeors
integer, parameter :: liq_only = 1    ! only mix qcl
integer, parameter :: liq_and_ice = 2 ! mix qcl and qcf
integer, parameter :: all_hydro = 3   ! mix qcl, qcf, qr, qgraup, qcf2

real(kind=real_umphys) :: diff_factor = rmdi
real(kind=real_umphys) :: mix_factor  = rmdi

logical :: L_print_w        = .false.
logical :: L_print_div      = .false.
logical :: L_diag_print_ops = .false. ! diagnostic prints for ops
logical :: L_print_pe       = .false. ! print diagnostics on all pe's
logical :: L_print_shear    = .false. ! wind diagnostic prints
logical :: L_print_max_wind = .false. ! wind diagnostic prints
logical :: L_diag_L2norms   = .false. ! l2norm diagnostic prints
logical :: L_diag_L2helm    = .false. ! l2norm diagnostic prints from solver
logical :: L_diag_wind      = .false. ! wind diagnostic prints
logical :: L_diag_noise     = .false. ! w diagnostic prints
logical :: L_flush6         = .false. ! flush buffers on failure
logical :: L_diag_print     = .false. ! Print diagnostics
logical :: L_print_lapse    = .false. ! Print lapse_rate diagnostics
logical :: L_print_wmax     = .false. ! Print max w diagnostic
logical :: L_print_theta1   = .false. ! Print level 1 theta diagnostic

logical :: L_diffusion    = .false. ! Removed from NL, set in code
logical :: L_upper_ramp   = .false. ! ramp upper-level diffusion
logical :: L_pftheta      = .false. ! activate polar filter of theta
logical :: L_pfuv         = .false. ! activate polar filter of u,v
logical :: L_pfw          = .false. ! activate polar filter of w
logical :: L_pofil_new    = .false. ! activate new polar filter

logical :: L_filter      = .false. ! Activate polar filter or diffusion
                                   ! (not in NL)
logical :: L_filter_incs = .false. ! Activate polar filter of incs (not in NL)

logical :: L_pfexner  ! Activate polar filter of Exner pressure
                      ! Always equal to L_pofil_hadgem2
                      ! Removed from NL, set in code

logical :: L_pofil_hadgem2 = .false. ! Use HadGEM2 polar filter setting
logical :: L_diff_exner  = .false.   ! activate diffusion of Exner pressure
logical :: L_sponge      = .false.   ! activate lateral bndaries sponge zones
logical :: L_pfincs      = .false.   ! activate polar filter of incs
logical :: L_diff_thermo = .false.   ! horiz. diffusion of theta
logical :: L_diff_wind   = .false.   ! horiz. diffusion of u, v
logical :: L_diff_wind_ew_only = .false.  ! Omit horiz. diff. in NS dirn.
logical :: L_diff_w      = .false.   ! horiz. diffusion of w
logical :: L_diff_incs   = .false.   ! horiz. diffusion of increments
logical :: L_diff_auto   = .false.   ! UM calculates diffusion parameters
logical :: L_tardiff_q   = .false.   ! activate targeted diffusion q

logical :: L_diff_ctl                ! general diffusion control
                                     ! Removed from NL, set in code

! L_diff_ctl == L_diffusion .or. L_cdiffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print

logical :: L_cdiffusion = .false.    ! Removed from NL, set in code

integer :: hdiffopt = imdi ! Horizontal diffusion option
                         ! 0: off, 1: old, 2: conservative, 3: subgrid
integer :: hdiff_strlev_wind (model_levels_max)=imdi ! Start, end levels for
integer :: hdiff_endlev_wind (model_levels_max)=imdi !  horiz diffn of wind
real(kind=real_umphys)    :: hdiff_coeff_wind  (model_levels_max)=rmdi
                                                     ! Horiz diff coeffs wind
integer :: hdiff_order_wind  (model_levels_max)=imdi ! Horiz diff orders wind
integer :: hdiff_strlev_theta(model_levels_max)=imdi ! Start, end levels for
integer :: hdiff_endlev_theta(model_levels_max)=imdi !  horiz diffn of theta
real(kind=real_umphys)    :: hdiff_coeff_theta (model_levels_max)=rmdi
                                                     ! Horiz diff coeffs theta
integer :: hdiff_order_theta (model_levels_max)=imdi ! Horiz diff orders theta
integer :: hdiff_strlev_q    (model_levels_max)=imdi ! Start, end levels for
integer :: hdiff_endlev_q    (model_levels_max)=imdi !  horiz diffusion of q
real(kind=real_umphys)    :: hdiff_coeff_q     (model_levels_max)=rmdi
                                                     ! Horiz diffn coeffs q
integer :: hdiff_order_q     (model_levels_max)=imdi ! Horiz diffn orders q

integer :: pofil_opt = imdi ! Polar filter option (0:none, 1:combined, 2:old)

integer  :: diffusion_order_thermo(model_levels_max) = imdi
integer  :: diffusion_order_wind(model_levels_max) = imdi
integer  :: diffusion_order_q(model_levels_max) = imdi
integer  :: diffusion_order_w(model_levels_max) = imdi
real(kind=real_umphys) :: diffusion_coefficient_thermo(model_levels_max) = rmdi
real(kind=real_umphys) :: diffusion_coefficient_wind(model_levels_max) = rmdi
real(kind=real_umphys) :: diffusion_coefficient_q(model_levels_max) = rmdi
real(kind=real_umphys) :: diffusion_coefficient_w(model_levels_max) = rmdi

integer :: print_step       = imdi ! To control diagnostic printing interval
integer :: diag_interval    = imdi ! diagnostic printing sampling frequency
integer :: norm_lev_start   = imdi ! start level for norm diagnostics
integer :: norm_lev_end     = imdi ! end level for norm diagnostics
integer :: first_norm_print = imdi ! first timestep for norm printing
integer :: dom_w_in = 0     ! define start for block printing
integer :: dom_e_in = 0     ! define end for block printing
integer :: dom_s_in = 0     ! define start for block printing
integer :: dom_n_in = 0     ! define end for block printing
integer :: blockx_in = 0    ! define size for block printing
integer :: blocky_in = 0    ! define size for  block printing

integer :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
integer :: u_end(0:max_121_rows)   ! Sweep control on 121 filter
integer :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
integer :: v_end(0:max_121_rows)   ! Sweep control on 121 filter
integer :: u_sweeps(max_121_rows)  ! Sweep control on 121 filter
integer :: v_sweeps(max_121_rows)  ! Sweep control on 121 filter
integer :: max_sweeps = imdi       ! Max sweeps wanted for 121 filter
integer :: global_u_filter         ! Sweep control on 121 filter
integer :: global_v_filter         ! Sweep control on 121 filter
integer :: diff_order_thermo = imdi      ! diffusion order for theta
integer :: diff_order_wind   = imdi      ! diffusion order for winds
integer :: diff_timescale_thermo = imdi  ! diffusion timescale for theta
integer :: diff_timescale_wind   = imdi  ! diffusion timescale for wind
integer :: top_filt_start = imdi ! start level upper-level diffusion
integer :: top_filt_end   = imdi ! end level upper-level diffusion
integer :: sponge_ew = imdi      ! left/right boundaries sponge zone width
integer :: sponge_ns = imdi      ! north/south boundaries sponge zone width
integer :: sponge_power = imdi   ! sponge zone weighting order
integer :: tardiffq_test  = imdi ! test level test w targetted diffusion
integer :: tardiffq_start = imdi ! start level test w targetted diffusion
integer :: tardiffq_end   = imdi ! end level test w targetted diffusion

!  level - assume surfaces are horizontal
integer  :: horizontal_level = imdi

real(kind=real_umphys) :: scale_ratio = rmdi ! Pass control on 121 filter
real(kind=real_umphys) :: ref_lat_deg = rmdi ! Reference lat for auto diffusion
real(kind=real_umphys) :: top_diff = rmdi !  upper-level diffusion coefficient
real(kind=real_umphys) :: up_diff_scale = rmdi
                        ! upper-level diffusion ramping factor
real(kind=real_umphys) :: adjust_lapse_min = 0.0
                        ! min dtheta/dz in vertical adjustment

real(kind=real_umphys) :: up_diff(max_updiff_levels)
                        ! upper-level diffusion coeff
real(kind=real_umphys) :: sponge_wts_ew(max_sponge_width) ! sponge weights
real(kind=real_umphys) :: sponge_wts_ns(max_sponge_width) ! sponge weights

real(kind=real_umphys) :: diff_coeff_ref = rmdi
                        ! EW diffusion coefficient at polar cap
real(kind=real_umphys) :: diff_coeff_thermo = rmdi  ! NS theta diffusion coeff
real(kind=real_umphys) :: diff_coeff_wind = rmdi   ! NS u,v diffusion coeff
real(kind=real_umphys) :: diff_coeff_phi = rmdi
                        ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
real(kind=real_umphys) :: polar_cap = rmdi ! Apply 1-2-1 filter polewards
real(kind=real_umphys) :: tardiffq_factor = rmdi
                        ! targeted diffusion coefficient
real(kind=real_umphys) :: w_print_limit ! w Threshold for diagnostic printing
real(kind=real_umphys) :: w_conv_limit = rmdi
                        ! w Threshold for limiting convection

! Divergence damping control options:
logical :: L_divdamp = .false.  ! Obsolete option - removed from NL
real(kind=real_umphys) :: div_damp_coefficient(model_levels_max)

! Polar filter control options:
logical :: L_polar_filter = .false.
! T: use polar filter to filter increment
logical :: L_polar_filter_incs = .false.

! Moisture resetting control options:
logical :: l_qpos = .false.           ! logical to run qpos code

real(kind=real_umphys) :: qlimit = rmdi ! lowest allowed value of q

! Inputs for the Leonard terms
logical :: l_leonard_term = .false.
real(kind=real_umphys) :: leonard_kl = rmdi


! Inputs for enhanced Smagorinsky and TKE schemes.
logical :: L_fullstress      = .false. ! use the full 3D stress term


namelist/RUN_Diffusion/                                                        &
  hdiffopt,                                                                    &
  hdiff_strlev_wind, hdiff_endlev_wind,                                        &
  hdiff_coeff_wind, hdiff_order_wind,                                          &
  hdiff_strlev_theta, hdiff_endlev_theta,                                      &
  hdiff_coeff_theta, hdiff_order_theta,                                        &
  hdiff_strlev_q, hdiff_endlev_q,                                              &
  hdiff_coeff_q, hdiff_order_q,                                                &
  pofil_opt,                                                                   &
  diff_order_thermo, diff_timescale_thermo,                                    &
  diff_order_wind, diff_timescale_wind,                                        &
  horizontal_level,                                                            &
  L_polar_filter,L_polar_filter_incs,                                          &
  l_qpos, qlimit,                                                              &
  L_tardiff_q, w_conv_limit,                                                   &
  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,                &
  L_subfilter_horiz, L_subfilter_vert, hyd_mix_opt,                            &
  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,             &
  turb_startlev_vert, turb_endlev_vert,                                        &
  L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pofil_hadgem2,                         &
  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_wind_ew_only,                &
  L_diff_w,                                                                    &
  L_pofil_new, L_diff_auto,                                                    &
  diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,                         &
  max_sweeps, L_upper_ramp, up_diff_scale, top_diff,                           &
  top_filt_start, top_filt_end,                                                &
  L_sponge, sponge_power, sponge_ew, sponge_ns,                                &
  L_diag_print, L_diag_print_ops, L_print_pe,                                  &
  L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,                      &
  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,                   &
  L_diag_noise, L_diag_L2norms, L_diag_L2helm,                                 &
  norm_lev_start, norm_lev_end, first_norm_print,                              &
  print_step, diag_interval, w_print_limit, L_Flush6,                          &
  l_leonard_term, leonard_kl, L_fullstress

! DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='TURB_DIFF_MOD'

contains

subroutine check_run_diffusion()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the run_diffusion namelist.

use ereport_mod, only: ereport
use bl_option_mod, only: non_local_bl, blending_option, ng_stress, off,        &
                         i_bl_vn, i_bl_vn_0, i_bl_vn_1a, i_bl_vn_9c
use mym_option_mod, only: l_3dtke
use cv_run_mod, only: cldbase_opt_sh
use cv_param_mod, only: sh_wstar_closure
use chk_opts_mod, only: chk_var, def_src

implicit none

integer                       :: ErrorStatus   ! used for ereport
integer                       :: i,j
character (len=errormessagelength)    :: cmessage      ! used for ereport
character (len=*), parameter  :: RoutineName = 'CHECK_RUN_DIFFUSION'


real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

if (l_3dtke .and. i_bl_vn==i_bl_vn_1a ) then
  if (.not. l_subfilter_horiz) then
    cmessage =                                               newline//         &
    'Possible inconsistency: l_subfilter_horiz false even though 3D TKE chosen.'
    ErrorStatus = -100
    call ereport(RoutineName, ErrorStatus, cmessage)
  end if
end if

! set BL parameters dependent on diffusion options
if (i_bl_vn==i_bl_vn_9c .and. l_subfilter_vert                                 &
    .and. blending_option == off) then

  non_local_bl = off
  cmessage =                                                 newline//         &
  'Non_local_bl switched off since vertical Smagorinsky'//   newline//         &
  'chosen without blending.'
  ErrorStatus = -100
  call ereport(RoutineName, ErrorStatus, cmessage)

  ng_stress = off
  cmessage =                                                 newline//         &
  'Non-gradient stress parametrization switched off since'// newline//         &
  'vertical Smagorinsky chosen without blending.'
  ErrorStatus = -100
  call ereport(RoutineName, ErrorStatus, cmessage)

end if

! set convection parameters dependent on diffusion options
if (.not. l_subfilter_horiz) then
  cldbase_opt_sh = sh_wstar_closure
  cmessage =                                                 newline//         &
  'cldbase_opt_sh set to sh_wstar_closure since'//           newline//         &
  'Smagorinsky diffusion not chosen.'
  ErrorStatus = -100
  call ereport(RoutineName, ErrorStatus, cmessage)
end if

! check horizontal Smagorinsky options
if (l_subfilter_horiz) then
  call chk_var(hyd_mix_opt,'hyd_mix_opt',                                      &
       [0,liq_only,liq_and_ice,all_hydro])
end if

! UMUI-to-ROSE transition UM vn9.0.
! Under the UMUI regime, the variable hdiffopt was not in the name list but
! was an internal umui variable used to set a number of logicals.
! To facilitate transition to ROSE, hdiffopt was moved into the name list,
! while the dependent logicals were removed from the it.
! The value of hdiffopt is used here to set dependent logicals.
! Logical l_diffusion switches on the "old diffusion option", while
! l_cdiffusion switches on the "conservative diffusion option".
if (hdiffopt==1) l_diffusion = .true.
if (hdiffopt==2) l_cdiffusion = .true.

if (hdiffopt==1 .or. hdiffopt==2) then
  do i = 1, model_levels_max
    if (hdiff_strlev_wind(i) > 0) then
      do j = hdiff_strlev_wind(i), hdiff_endlev_wind(i)
        diffusion_coefficient_wind(j) = hdiff_coeff_wind(i)
        diffusion_order_wind      (j) = hdiff_order_wind(i)
      end do
    end if
    if (hdiff_strlev_theta(i) > 0) then
      do j = hdiff_strlev_theta(i), hdiff_endlev_theta(i)
        diffusion_coefficient_thermo(j) = hdiff_coeff_theta(i)
        diffusion_order_thermo      (j) = hdiff_order_theta(i)
      end do
    end if
    if (hdiff_strlev_q(i) > 0) then
      do j = hdiff_strlev_q(i), hdiff_endlev_q(i)
        diffusion_coefficient_q(j) = hdiff_coeff_q(i)
        diffusion_order_q      (j) = hdiff_order_q(i)
      end do
    end if
  end do
end if

if (hdiffopt /=0 .or. l_tardiff_q                                              &
     .or. l_diag_print) then
  L_diff_ctl = .true.
else
  L_diff_ctl = .false.
end if

! These two logicals were both in the NL pre-UM 9.0, but the umui
!  code gave them the same value. One of them has been removed from
!  the NL under ROSE rationalisation.
l_pfexner = l_pofil_hadgem2

! Check that Leonard term parameter within allowed range, if used.
if (l_leonard_term)  call chk_var( leonard_kl, 'leonard_kl', '[0.0:6.0]' )

! Check not trying to use Leonard terms without BL scheme; won't work as
! Leonard term calculation uses some arrays calculated in the BL scheme
if ( l_leonard_term .and. i_bl_vn == i_bl_vn_0 ) then
  cmessage =                                                 newline//         &
  'Cannot run with Leonard terms switched on'              //newline//         &
  '(l_leonard_term = .true.)'                              //newline//         &
  'but the boundary-layer scheme switched off'             //newline//         &
  '(i_bl_vn == i_bl_vn_0).'                                //newline//         &
  'You must switch on the BL scheme if you want to run'    //newline//         &
  'with Leonard terms, as the Leonard term calculations'   //newline//         &
  'depend on outputs from the BL code.'
  ErrorStatus = 100
  call ereport(RoutineName, ErrorStatus, cmessage)
end if


def_src = ''
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine check_run_diffusion

subroutine print_nlist_run_diffusion()
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_DIFFUSION'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_diffusion',                             &
    src='nl_cruntimc_h_run_diffusion')

write(lineBuffer,'(A,I0)')' hdiffopt = ',hdiffopt
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,'(A,I0)')' pofil_opt = ',pofil_opt
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diffusion_order_thermo = ',diffusion_order_thermo
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diffusion_order_wind = ',diffusion_order_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diffusion_order_q = ',diffusion_order_q
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)                                                            &
    ' diffusion_coefficient_thermo = ',diffusion_coefficient_thermo
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)                                                            &
    ' diffusion_coefficient_wind = ',diffusion_coefficient_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diffusion_coefficient_q = ',diffusion_coefficient_q
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_order_thermo = ',diff_order_thermo
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_timescale_thermo = ',diff_timescale_thermo
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_order_wind = ',diff_order_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_timescale_wind = ',diff_timescale_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' horizontal_level = ',horizontal_level
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_polar_filter = ',L_polar_filter
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_polar_filter_incs = ',L_polar_filter_incs
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' l_qpos = ',l_qpos
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' qlimit = ',qlimit
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_tardiff_q = ',L_tardiff_q
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' w_conv_limit = ',w_conv_limit
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' tardiffq_factor = ',tardiffq_factor
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' tardiffq_test = ',tardiffq_test
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' tardiffq_start = ',tardiffq_start
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' tardiffq_end = ',tardiffq_end
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_subfilter_horiz = ',L_subfilter_horiz
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_subfilter_vert = ',L_subfilter_vert
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,'(A,L1)')' L_fullstress = ',L_fullstress
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,'(A,I0)')' hyd_mix_opt = ',hyd_mix_opt
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_factor = ',diff_factor
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' mix_factor = ',mix_factor
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' turb_startlev_horiz = ',turb_startlev_horiz
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' turb_endlev_horiz = ',turb_endlev_horiz
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' turb_startlev_vert = ',turb_startlev_vert
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' turb_endlev_vert = ',turb_endlev_vert
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pftheta = ',L_pftheta
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pfuv = ',L_pfuv
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pfw = ',L_pfw
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pfincs = ',L_pfincs
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pofil_hadgem2 = ',L_pofil_hadgem2
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_incs = ',L_diff_incs
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_thermo = ',L_diff_thermo
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_wind = ',L_diff_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_wind_ew_only = ',L_diff_wind_ew_only
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_w = ',L_diff_w
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_pofil_new = ',L_pofil_new
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diff_auto = ',L_diff_auto
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diff_coeff_ref = ',diff_coeff_ref
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' polar_cap = ',polar_cap
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' scale_ratio = ',scale_ratio
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' ref_lat_deg = ',ref_lat_deg
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' max_sweeps = ',max_sweeps
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_upper_ramp = ',L_upper_ramp
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' up_diff_scale = ',up_diff_scale
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' top_diff = ',top_diff
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' top_filt_start = ',top_filt_start
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' top_filt_end = ',top_filt_end
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_sponge = ',L_sponge
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' sponge_power = ',sponge_power
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' sponge_ew = ',sponge_ew
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' sponge_ns = ',sponge_ns
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_print = ',L_diag_print
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_print_ops = ',L_diag_print_ops
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_pe = ',L_print_pe
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_w = ',L_print_w
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_wmax = ',L_print_wmax
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_lapse = ',L_print_lapse
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_theta1 = ',L_print_theta1
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_div = ',L_print_div
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_wind = ',L_diag_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_shear = ',L_print_shear
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_print_max_wind = ',L_print_max_wind
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_noise = ',L_diag_noise
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_L2norms = ',L_diag_L2norms
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_diag_L2helm = ',L_diag_L2helm
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' norm_lev_start = ',norm_lev_start
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' norm_lev_end = ',norm_lev_end
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' first_norm_print = ',first_norm_print
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' print_step = ',print_step
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' diag_interval = ',diag_interval
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' w_print_limit = ',w_print_limit
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' L_Flush6 = ',L_Flush6
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' l_leonard_term = ',l_leonard_term
call umPrint(lineBuffer,src='turb_diff_mod')
write(lineBuffer,*)' leonard_kl = ',leonard_kl
call umPrint(lineBuffer,src='turb_diff_mod')


call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='nl_cruntimc_h_run_diffusion')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_diffusion


end module turb_diff_mod
