!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of variables for JULES pftparm namelist

module jules_pftparm_init_mod

  implicit none

  private
  public :: jules_pftparm_init

contains

  !>@brief Initialise JULES pftparm variables used in JULES code
  !>@details Pass JULES pftparm variables from model config to pftparm_mod
  !> @param[in] config   The config of the model run
  subroutine jules_pftparm_init(config)
    use c_z0h_z0m,                only: z0h_z0m
    use config_mod,               only: config_type
    use constants_mod,            only: r_um, i_def
    use jules_pftparm_config_mod, only:                                        &
       c3_io_no, c3_io_yes,fsmc_mod_io_weight, fsmc_mod_io_average,            &
       orient_io_spherical, orient_io_horizontal
    use pftparm, only: a_wl, a_ws, act_jmax, act_vcmax, aef, albsnc_max,       &
       albsnc_min, albsnf_max, albsnf_maxl, albsnf_maxu, alnir, alnirl,        &
       alniru, alpar, alparl, alparu, alpha, alpha_elec, avg_ba, b_wl, c3,     &
       can_struct_a, catch0, ccleaf_max, ccleaf_min, ccwood_max, ccwood_min,   &
       ci_st, dcatch_dlai, deact_jmax, deact_vcmax, dfp_dcuo, dgl_dm, dgl_dt,  &
       dqcrit, ds_jmax, ds_vcmax, dust_veg_scj, dz0v_dh, emis_pft, eta_sl, f0, &
       fd, fef_bc, fef_c2h4, fef_c2h6, fef_c3h8, fef_ch4, fef_co, fef_co2,     &
       fef_dms, fef_hcho, fef_mecho, fef_nh3, fef_nox, fef_oc, fef_so2,        &
       fire_mort, fl_o3_ct, fsmc_mod, fsmc_of, fsmc_p0, g1_stomata, g_leaf_0,  &
       glmin, gpp_st, gsoil_f, hw_sw, ief, infil_f, jv25_ratio, kext, kn, knl, &
       kpar, lai_alb_lim, lma, mef, neff, nl0, nmass, nr, nr_nl, ns_nl, nsw,   &
       omega, omegal, omegau, omnir, omnirl, omniru, orient, psi_close,        &
       psi_open, q10_leaf, r_grow, rootd_ft, sigl, sox_a, sox_p50, sox_rp_min, &
       sug_g0, sug_grec, sug_yg, tef, tleaf_of, tlow, tupp, vint, vsl, z0v

    use jules_pftparm_nml_iterator_mod, only: jules_pftparm_nml_iterator_type
    use jules_pftparm_nml_mod,    only: jules_pftparm_nml_type
    use jules_surface_types_mod,  only: npft, brd_leaf, ndl_leaf, c3_grass,    &
       c4_grass, shrub

    use log_mod, only: log_event, log_scratch_space, log_level_error

    implicit none

    ! Model run working data set
    type(config_type), intent(in) :: config

    type(jules_pftparm_nml_iterator_type) :: iter
    type(jules_pftparm_nml_type), pointer :: jules_pftparm

    integer(kind=i_def) :: i, n

    character(len=*), parameter :: RoutineName='JULES_PFTPARM_INIT'

    n = 0
    call iter%initialise( config%jules_pftparm )
    do while ( iter%has_next() )
      n = n + 1
      jules_pftparm => iter%next()
      ! For now add mapping from instance to jules_surface_types
      select case ( trim( jules_pftparm%pft_name_io() ) )
      case ( 'brd_leaf' )
        i = brd_leaf
      case ( 'ndl_leaf' )
        i = ndl_leaf
      case ( 'c3_grass' )
        i = c3_grass
      case ( 'c4_grass' )
        i = c4_grass
      case ( 'shrub' )
        i = shrub
      case DEFAULT
        write(log_scratch_space,'(A)')                                         &
           'PFT name not recognised: ' // jules_pftparm%pft_name_io()
        call log_event(                                                        &
           RoutineName//': '//trim(log_scratch_space), log_level_error         &
           )
      end select

      ! Range of specified types (1:npft) checked by check_jules_surface_types
      if ( i < 1 ) then
        write(log_scratch_space,'(A)')                                         &
           'jules_pftparm and jules_surface_types inputs are inconsistent; '// &
           trim(jules_pftparm%pft_name_io()) //                                &
           ' is not specified in jules_surface_types'
        call log_event(                                                        &
           RoutineName//': '//trim(log_scratch_space), log_level_error         &
           )
      end if

      ! c3_io would make more sense as a logical
      ! (see MetOffice/jules/issues/106)
      select case ( jules_pftparm%c3_io() )
      case ( c3_io_no )
        c3(i) = 0
      case ( c3_io_yes )
        c3(i) = 1
      end select
      select case ( jules_pftparm%fsmc_mod_io() )
      case ( fsmc_mod_io_weight )
        fsmc_mod(i) = 0
      case ( fsmc_mod_io_average )
        fsmc_mod(i) = 1
      end select
      select case ( jules_pftparm%orient_io() )
      case ( orient_io_spherical )
        orient(i) = 0
      case ( orient_io_horizontal )
        orient(i) = 1
      end select

      a_wl(i) = real(jules_pftparm%a_wl_io(), r_um)
      a_ws(i) = real(jules_pftparm%a_ws_io(), r_um)
      act_jmax(i) = real(jules_pftparm%act_jmax_io(), r_um)
      act_vcmax(i) = real(jules_pftparm%act_vcmax_io(), r_um)
      aef(i) = real(jules_pftparm%aef_io(), r_um)
      albsnc_max(i) = real(jules_pftparm%albsnc_max_io(), r_um)
      albsnc_min(i) = real(jules_pftparm%albsnc_min_io(), r_um)
      albsnf_max(i) = real(jules_pftparm%albsnf_max_io(), r_um)
      albsnf_maxl(i) = real(jules_pftparm%albsnf_maxl_io(), r_um)
      albsnf_maxu(i) = real(jules_pftparm%albsnf_maxu_io(), r_um)
      alnir(i) = real(jules_pftparm%alnir_io(), r_um)
      alnirl(i) = real(jules_pftparm%alnirl_io(), r_um)
      alniru(i) = real(jules_pftparm%alniru_io(), r_um)
      alpar(i) = real(jules_pftparm%alpar_io(), r_um)
      alparl(i) = real(jules_pftparm%alparl_io(), r_um)
      alparu(i) = real(jules_pftparm%alparu_io(), r_um)
      alpha(i) = real(jules_pftparm%alpha_io(), r_um)
      alpha_elec(i) = real(jules_pftparm%alpha_elec_io(), r_um)
      avg_ba(i) = real(jules_pftparm%avg_ba_io(), r_um)
      b_wl(i) = real(jules_pftparm%b_wl_io(), r_um)
      can_struct_a(i) = real(jules_pftparm%can_struct_a_io(), r_um)
      catch0(i) = real(jules_pftparm%catch0_io(), r_um)
      ccleaf_max(i) = real(jules_pftparm%ccleaf_max_io(), r_um)
      ccleaf_min(i) = real(jules_pftparm%ccleaf_min_io(), r_um)
      ccwood_max(i) = real(jules_pftparm%ccwood_max_io(), r_um)
      ccwood_min(i) = real(jules_pftparm%ccwood_min_io(), r_um)
      ci_st(i) = real(jules_pftparm%ci_st_io(), r_um)
      dcatch_dlai(i) = real(jules_pftparm%dcatch_dlai_io(), r_um)
      deact_jmax(i) = real(jules_pftparm%deact_jmax_io(), r_um)
      deact_vcmax(i) = real(jules_pftparm%deact_vcmax_io(), r_um)
      dfp_dcuo(i) = real(jules_pftparm%dfp_dcuo_io(), r_um)
      dgl_dm(i) = real(jules_pftparm%dgl_dm_io(), r_um)
      dgl_dt(i) = real(jules_pftparm%dgl_dt_io(), r_um)
      dqcrit(i) = real(jules_pftparm%dqcrit_io(), r_um)
      ds_jmax(i) = real(jules_pftparm%ds_jmax_io(), r_um)
      ds_vcmax(i) = real(jules_pftparm%ds_vcmax_io(), r_um)
      dust_veg_scj(i) = real(jules_pftparm%dust_veg_scj_io(), r_um)
      dz0v_dh(i) = real(jules_pftparm%dz0v_dh_io(), r_um)
      emis_pft(i) = real(jules_pftparm%emis_pft_io(), r_um)
      eta_sl(i) = real(jules_pftparm%eta_sl_io(), r_um)
      f0(i) = real(jules_pftparm%f0_io(), r_um)
      fd(i) = real(jules_pftparm%fd_io(), r_um)
      fef_bc(i) = real(jules_pftparm%fef_bc_io(), r_um)
      fef_c2h4(i) = real(jules_pftparm%fef_c2h4_io(), r_um)
      fef_c2h6(i) = real(jules_pftparm%fef_c2h6_io(), r_um)
      fef_c3h8(i) = real(jules_pftparm%fef_c3h8_io(), r_um)
      fef_ch4(i) = real(jules_pftparm%fef_ch4_io(), r_um)
      fef_co(i) = real(jules_pftparm%fef_co_io(), r_um)
      fef_co2(i) = real(jules_pftparm%fef_co2_io(), r_um)
      fef_dms(i) = real(jules_pftparm%fef_dms_io(), r_um)
      fef_hcho(i) = real(jules_pftparm%fef_hcho_io(), r_um)
      fef_mecho(i) = real(jules_pftparm%fef_mecho_io(), r_um)
      fef_nh3(i) = real(jules_pftparm%fef_nh3_io(), r_um)
      fef_nox(i) = real(jules_pftparm%fef_nox_io(), r_um)
      fef_oc(i) = real(jules_pftparm%fef_oc_io(), r_um)
      fef_so2(i) = real(jules_pftparm%fef_so2_io(), r_um)
      fire_mort(i) = real(jules_pftparm%fire_mort_io(), r_um)
      fl_o3_ct(i) = real(jules_pftparm%fl_o3_ct_io(), r_um)
      fsmc_of(i) = real(jules_pftparm%fsmc_of_io(), r_um)
      fsmc_p0(i) = real(jules_pftparm%fsmc_p0_io(), r_um)
      g1_stomata(i) = real(jules_pftparm%g1_stomata_io(), r_um)
      g_leaf_0(i) = real(jules_pftparm%g_leaf_0_io(), r_um)
      glmin(i) = real(jules_pftparm%glmin_io(), r_um)
      gpp_st(i) = real(jules_pftparm%gpp_st_io(), r_um)
      gsoil_f(i) = real(jules_pftparm%gsoil_f_io(), r_um)
      hw_sw(i) = real(jules_pftparm%hw_sw_io(), r_um)
      ief(i) = real(jules_pftparm%ief_io(), r_um)
      infil_f(i) = real(jules_pftparm%infil_f_io(), r_um)
      jv25_ratio(i) = real(jules_pftparm%jv25_ratio_io(), r_um)
      kext(i) = real(jules_pftparm%kext_io(), r_um)
      kn(i) = real(jules_pftparm%kn_io(), r_um)
      knl(i) = real(jules_pftparm%knl_io(), r_um)
      kpar(i) = real(jules_pftparm%kpar_io(), r_um)
      lai_alb_lim(i) = real(jules_pftparm%lai_alb_lim_io(), r_um)
      lma(i) = real(jules_pftparm%lma_io(), r_um)
      mef(i) = real(jules_pftparm%mef_io(), r_um)
      neff(i) = real(jules_pftparm%neff_io(), r_um)
      nl0(i) = real(jules_pftparm%nl0_io(), r_um)
      nmass(i) = real(jules_pftparm%nmass_io(), r_um)
      nr(i) = real(jules_pftparm%nr_io(), r_um)
      nr_nl(i) = real(jules_pftparm%nr_nl_io(), r_um)
      ns_nl(i) = real(jules_pftparm%ns_nl_io(), r_um)
      nsw(i) = real(jules_pftparm%nsw_io(), r_um)
      omega(i) = real(jules_pftparm%omega_io(), r_um)
      omegal(i) = real(jules_pftparm%omegal_io(), r_um)
      omegau(i) = real(jules_pftparm%omegau_io(), r_um)
      omnir(i) = real(jules_pftparm%omnir_io(), r_um)
      omnirl(i) = real(jules_pftparm%omnirl_io(), r_um)
      omniru(i) = real(jules_pftparm%omniru_io(), r_um)
      psi_close(i) = real(jules_pftparm%psi_close_io(), r_um)
      psi_open(i) = real(jules_pftparm%psi_open_io(), r_um)
      q10_leaf(i) = real(jules_pftparm%q10_leaf_io(), r_um)
      r_grow(i) = real(jules_pftparm%r_grow_io(), r_um)
      rootd_ft(i) = real(jules_pftparm%rootd_ft_io(), r_um)
      sigl(i) = real(jules_pftparm%sigl_io(), r_um)
      sox_a(i) = real(jules_pftparm%sox_a_io(), r_um)
      sox_p50(i) = real(jules_pftparm%sox_p50_io(), r_um)
      sox_rp_min(i) = real(jules_pftparm%sox_rp_min_io(), r_um)
      sug_g0(i) = real(jules_pftparm%sug_g0_io(), r_um)
      sug_grec(i) = real(jules_pftparm%sug_grec_io(), r_um)
      sug_yg(i) = real(jules_pftparm%sug_yg_io(), r_um)
      tef(i) = real(jules_pftparm%tef_io(), r_um)
      tleaf_of(i) = real(jules_pftparm%tleaf_of_io(), r_um)
      tlow(i) = real(jules_pftparm%tlow_io(), r_um)
      tupp(i) = real(jules_pftparm%tupp_io(), r_um)
      vint(i) = real(jules_pftparm%vint_io(), r_um)
      vsl(i) = real(jules_pftparm%vsl_io(), r_um)
      z0v(i) = real(jules_pftparm%z0v_io(), r_um)
      z0h_z0m(i) = real(jules_pftparm%z0hm_pft_io(), r_um)
    end do

    if ( n /= npft ) then
      write(log_scratch_space,'(2(A,I0))')                                     &
         'Number of instances of jules_pftparm namelist (', n,                 &
         ') is not npft = ', npft
      call log_event(                                                          &
         RoutineName//': '//trim(log_scratch_space), log_level_error           &
         )
    end if

  end subroutine jules_pftparm_init

end module jules_pftparm_init_mod
