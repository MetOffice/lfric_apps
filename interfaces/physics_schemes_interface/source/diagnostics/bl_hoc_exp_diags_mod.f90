!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Processes diagnostics for bl_hoc_exp_alg

module bl_hoc_exp_diags_mod

  use constants_mod,       only: l_def
  use field_mod,           only: field_type
  use integer_field_mod,   only: integer_field_type
  use timing_mod,          only: start_timing, stop_timing, tik, LPROF
  use initialise_diagnostics_mod,     only : init_diag => init_diagnostic_field

  implicit none

  private

  type, public :: bl_hoc_diags_flag_type
    ! Logical indicating whether diagnostics are requested
    logical( l_def ) :: master_length
    logical( l_def ) :: rhogamu_bl
    logical( l_def ) :: rhogamv_bl
    logical( l_def ) :: rhogamt_bl
    logical( l_def ) :: rhogamq_bl
    logical( l_def ) :: tke_shr_prod
    logical( l_def ) :: tke_boy_prod
    logical( l_def ) :: tke_dissp
    logical( l_def ) :: sm25
    logical( l_def ) :: sh25
    logical( l_def ) :: dbdz
    logical( l_def ) :: dvdzm
  end type

  public :: initialise_diags_for_bl_hoc_exp
  public :: output_diags_for_bl_hoc_exp

contains

  !> @brief Initialise fields for locally-computed diagnostics
  !> @param[in,out] diags_flag    Logical switch of diagnostics
  !> @param[in,out] master_length Turbulent length scale
  !> @param[in,out] rhokm_bl      Momentum eddy diffusivity on BL levels
  !> @param[in,out] rhokh_bl      Heat eddy diffusivity on BL levels
  !> @param[in,out] rhogamu_bl    Counter Gradient Flux Term for U
  !> @param[in,out] rhogamv_bl    Counter Gradient Flux Term for V
  !> @param[in,out] rhogamt_bl    Counter Gradient Flux Term for tl
  !> @param[in,out] rhogamq_bl    Counter Gradient Flux Term for qt
  !> @param[in,out] tke_shr_prod  Production rate of TKE by shear
  !> @param[in,out] tke_boy_prod  Production rate of TKE by buoyancy
  !> @param[in,out] tke_dissp     Dissipation rate of TKE
  !> @param[in,out] sm25          Stability function for momentum
  !> @param[in,out] sh25          Stability function for scalar
  !> @param[in,out] dbdz          Vertical gradient of buoyancy
  !> @param[in,out] dvdzm         Modulus of wind shear
  subroutine initialise_diags_for_bl_hoc_exp( diags_flag,                      &
                                              master_length,                   &
                                              rhogamu_bl, rhogamv_bl,          &
                                              rhogamt_bl, rhogamq_bl,          &
                                              tke_shr_prod, tke_boy_prod,      &
                                              tke_dissp,                       &
                                              sm25, sh25, dbdz, dvdzm )

    implicit none

    type( bl_hoc_diags_flag_type ), intent(inout) :: diags_flag
    type( field_type ), intent(inout) :: master_length
    type( field_type ), intent(inout) :: rhogamu_bl
    type( field_type ), intent(inout) :: rhogamv_bl
    type( field_type ), intent(inout) :: rhogamt_bl
    type( field_type ), intent(inout) :: rhogamq_bl
    type( field_type ), intent(inout) :: tke_shr_prod
    type( field_type ), intent(inout) :: tke_boy_prod
    type( field_type ), intent(inout) :: tke_dissp
    type( field_type ), intent(inout) :: sm25
    type( field_type ), intent(inout) :: sh25
    type( field_type ), intent(inout) :: dbdz
    type( field_type ), intent(inout) :: dvdzm

    integer( tik ) :: id

    if ( LPROF ) call start_timing( id, 'diags.bl_exp' )

    diags_flag % master_length = init_diag(master_length, &
                                           'turbulence__master_length')
    diags_flag % rhogamu_bl    = init_diag(rhogamu_bl, 'turbulence__rhogamu')
    diags_flag % rhogamv_bl    = init_diag(rhogamv_bl, 'turbulence__rhogamv')
    diags_flag % rhogamt_bl    = init_diag(rhogamt_bl, 'turbulence__rhogamt')
    diags_flag % rhogamq_bl    = init_diag(rhogamq_bl, 'turbulence__rhogamq')
    diags_flag % tke_shr_prod  = init_diag(tke_shr_prod, &
                                           'turbulence__tke_shr_prod')
    diags_flag % tke_boy_prod  = init_diag(tke_boy_prod, &
                                           'turbulence__tke_boy_prod')
    diags_flag % tke_dissp     = init_diag(tke_dissp, 'turbulence__tke_dissp')
    diags_flag % sm25          = init_diag(sm25, 'turbulence__sm25')
    diags_flag % sh25          = init_diag(sh25, 'turbulence__sh25')
    diags_flag % dbdz          = init_diag(dbdz, 'turbulence__dbdz')
    diags_flag % dvdzm         = init_diag(dvdzm, 'turbulence__dvdzm')

    if ( LPROF ) call stop_timing( id, 'diags.bl_exp' )

  end subroutine initialise_diags_for_bl_hoc_exp

  !> @brief Output diagnostics from bl_hoc_exp_alg
  !> @param[in] diags_flag        Logical switch of diagnostics
  !> @param[in] ntml              Number of turbulently mixed levels
  !> @param[in] cumulus           Cumulus flag (true/false)
  !> @param[in] bl_type_ind       Diagnosed BL types
  !> @param[in] tke_bl            Turbulent kinetic energy (m2 s-2)
  !> @param[in] tsq_bl            Self covariant of tl'
  !> @param[in] qsq_bl            Self covariant of qw'
  !> @param[in] cov_bl            Correlation between tl' and qw'
  !> @param[in] master_length     Turbulent length scale
  !> @param[in] gradrinr          Gradient Richardson number in wth
  !> @param[in] rhokm_bl          Momentum eddy diffusivity on BL levels
  !> @param[in] rhokh_bl          Heat eddy diffusivity on BL levels
  !> @param[in] rhogamu_bl        Counter Gradient Flux Term for U
  !> @param[in] rhogamv_bl        Counter Gradient Flux Term for V
  !> @param[in] rhogamt_bl        Counter Gradient Flux Term for tl
  !> @param[in] rhogamq_bl        Counter Gradient Flux Term for qt
  !> @param[in] tke_shr_prod      Production rate of TKE by shear
  !> @param[in] tke_boy_prod      Production rate of TKE by buoyancy
  !> @param[in] tke_dissp         Dissipation rate of TKE
  !> @param[in] sm25              Stability function for momentum
  !> @param[in] sh25              Stability function for scalar
  !> @param[in] dbdz              Vertical gradient of buoyancy
  !> @param[in] dvdzm             Modulus of wind shear
  !> @param[in] dtrdz_tq_bl       dt/(rho*r*r*dz) in wth
  !> @param[in] rdz_tq_bl         1/dz in w3
  subroutine output_diags_for_bl_hoc_exp(diags_flag,                           &
                                         ntml, cumulus, bl_type_ind,           &
                                         tke_bl, tsq_bl, qsq_bl, cov_bl,       &
                                         master_length, gradrinr,              &
                                         rhokm_bl, rhokh_bl,                   &
                                         rhogamu_bl, rhogamv_bl,               &
                                         rhogamt_bl, rhogamq_bl,               &
                                         tke_shr_prod, tke_boy_prod, tke_dissp,&
                                         sm25, sh25, dbdz, dvdzm,              &
                                         dtrdz_tq_bl, rdz_tq_bl)

    implicit none

    type( bl_hoc_diags_flag_type ), intent(in) :: diags_flag
    ! Prognostic fields to output
    type( field_type ), intent(in)    :: tke_bl, tsq_bl, qsq_bl, cov_bl,       &
                                         master_length, gradrinr,              &
                                         rhokm_bl, rhokh_bl,                   &
                                         rhogamu_bl, rhogamv_bl,               &
                                         rhogamt_bl, rhogamq_bl,               &
                                         tke_shr_prod, tke_boy_prod,           &
                                         tke_dissp,                            &
                                         sm25, sh25, dbdz, dvdzm,              &
                                         dtrdz_tq_bl, rdz_tq_bl
    type(integer_field_type), intent(in) :: ntml, cumulus, bl_type_ind

    integer( tik ) :: id

    if ( LPROF ) call start_timing( id, 'diags.bl_exp' )

    ! Prognostic fields from turbulence collection
    call ntml%write_field('turbulence__ntml')
    call cumulus%write_field('turbulence__cumulus')
    call bl_type_ind%write_field('turbulence__bl_type_ind')
    call tke_bl%write_field('turbulence__tke')
    call tsq_bl%write_field('turbulence__tsq')
    call qsq_bl%write_field('turbulence__qsq')
    call cov_bl%write_field('turbulence__cov')
    call gradrinr%write_field('turbulence__gradrinr')
    call rhokm_bl%write_field('turbulence__rhokm')
    call rhokh_bl%write_field('turbulence__rhokh')
    call dtrdz_tq_bl%write_field('turbulence__dtrdz_tq')
    call rdz_tq_bl%write_field('turbulence__rdz_tq')

    if (diags_flag % master_length) then
      call master_length%write_field('turbulence__master_length')
    end if
    if (diags_flag % rhogamu_bl) then
      call rhogamu_bl%write_field('turbulence__rhogamu')
    end if
    if (diags_flag % rhogamv_bl) then
      call rhogamv_bl%write_field('turbulence__rhogamv')
    end if
    if (diags_flag % rhogamt_bl) then
      call rhogamt_bl%write_field('turbulence__rhogamt')
    end if
    if (diags_flag % rhogamq_bl) then
      call rhogamq_bl%write_field('turbulence__rhogamq')
    end if
    if (diags_flag % tke_shr_prod) then
      call tke_shr_prod%write_field('turbulence__tke_shr_prod')
    end if
    if (diags_flag % tke_boy_prod) then
      call tke_boy_prod%write_field('turbulence__tke_boy_prod')
    end if
    if (diags_flag % tke_dissp) then
      call tke_dissp%write_field('turbulence__tke_dissp')
    end if
    if (diags_flag % sm25) then
      call sm25%write_field('turbulence__sm25')
    end if
    if (diags_flag % sh25) then
      call sh25%write_field('turbulence__sh25')
    end if
    if (diags_flag % dbdz) then
      call dbdz%write_field('turbulence__dbdz')
    end if
    if (diags_flag % dvdzm) then
      call dvdzm%write_field('turbulence__dvdzm')
    end if

    if ( LPROF ) call stop_timing( id, 'diags.bl_exp' )

  end subroutine output_diags_for_bl_hoc_exp
end module bl_hoc_exp_diags_mod
