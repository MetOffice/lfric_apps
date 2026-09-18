!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Condensation and cloud condensation nuclei number concentrations
!>        from the GLOMAP-mode aerosol size distribution.
!>
!> @details Each GLOMAP mode is a lognormal distribution of dry particle
!>          diameter with geometric mean diameter drydp and fixed geometric
!>          standard deviation sigmag. The number of particles per unit
!>          volume with a dry diameter larger than a threshold dp0 is
!>          therefore available analytically as
!>
!>            0.5 * nd * ( 1 - erf( ln(dp0/drydp) / (sqrt(2)*ln(sigmag)) ) )
!>
!>          summed over the modes. This gives the condensation nuclei count
!>          (dp0 = 3 nm) and the cloud condensation nuclei counts at two
!>          activation thresholds (dp0 = 30 nm and dp0 = 50 nm).
!>
!>          Only the six modes that carry a dry modal diameter in LFRic are
!>          summed. The nucleation-soluble mode is omitted because LFRic has
!>          no drydp_nuc_sol field; its contribution is negligible at the
!>          30 nm and 50 nm thresholds but not at 3 nm, so the condensation
!>          nuclei count is a low estimate under prognostic UKCA.
!>
!>          The dust mass concentrations are the dust mass mixing ratios of
!>          the accumulation and coarse insoluble modes multiplied by the
!>          dry air density, which is the density the mixing ratios are
!>          defined against. They are computed here rather than as an XIOS
!>          expression so that every operand is on the aerosol mesh.
!>
!>          Every output field is only populated when it has been requested
!>          as a diagnostic; unrequested fields share the empty data array
!>          and are left untouched.

module glomap_ccn_diag_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_SCALAR,         &
                               GH_READ, GH_WRITE, &
                               CELL_COLUMN

  use fs_continuity_mod, only: WTHETA

  use kernel_mod,        only: kernel_type

  use empty_data_mod,    only: empty_real_data

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer

  type, public, extends(kernel_type) :: glomap_ccn_diag_kernel_type
    private
    type(arg_type) :: meta_args(24) = (/                &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! cn_number_conc
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! ccn_no_conc_30nm
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! ccn_no_conc_50nm
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! mconc_du_acc_ins
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA), & ! mconc_du_cor_ins
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! p_zero
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! one_over_kappa
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! theta_in_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! rho_in_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_ait_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_acc_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_cor_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_ait_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_acc_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! n_cor_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_ait_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_acc_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_cor_sol
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_ait_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_acc_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! drydp_cor_ins
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! acc_ins_du
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA)  & ! cor_ins_du
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: glomap_ccn_diag_code
  end type glomap_ccn_diag_kernel_type

  public :: glomap_ccn_diag_code

contains

!> @brief Sum the lognormal tail of each GLOMAP mode above three dry diameter
!>        thresholds to give condensation and cloud condensation nuclei
!>        counts, and convert the dust mass mixing ratios to concentrations.
!> @param[in]     nlayers              The number of layers
!> @param[in,out] cn_number_conc       Condensation nuclei number
!!                                      concentration, dry diameter > 3 nm
!> @param[in,out] ccn_number_conc_30nm Cloud condensation nuclei number
!!                                      concentration, dry diameter > 30 nm
!> @param[in,out] ccn_number_conc_50nm Cloud condensation nuclei number
!!                                      concentration, dry diameter > 50 nm
!> @param[in,out] mconc_du_acc_ins     Dust mass concentration in the
!!                                      accumulation insoluble mode
!> @param[in,out] mconc_du_cor_ins     Dust mass concentration in the
!!                                      coarse insoluble mode
!> @param[in]     p_zero               Reference surface pressure
!> @param[in]     one_over_kappa       Reciprocal of the ratio of the gas
!!                                      constant to the specific heat
!> @param[in]     theta_in_wth         Potential temperature field
!> @param[in]     exner_in_wth         Exner pressure in potential
!!                                      temperature space
!> @param[in]     rho_in_wth           Dry air density in potential
!!                                      temperature space
!> @param[in]     n_ait_sol            Aitken soluble mode number mixing ratio
!> @param[in]     n_acc_sol            Accumulation soluble mode number
!!                                      mixing ratio
!> @param[in]     n_cor_sol            Coarse soluble mode number mixing ratio
!> @param[in]     n_ait_ins            Aitken insoluble mode number
!!                                      mixing ratio
!> @param[in]     n_acc_ins            Accumulation insoluble mode number
!!                                      mixing ratio
!> @param[in]     n_cor_ins            Coarse insoluble mode number
!!                                      mixing ratio
!> @param[in]     drydp_ait_sol        Aitken soluble mode dry diameter
!> @param[in]     drydp_acc_sol        Accumulation soluble mode dry diameter
!> @param[in]     drydp_cor_sol        Coarse soluble mode dry diameter
!> @param[in]     drydp_ait_ins        Aitken insoluble mode dry diameter
!> @param[in]     drydp_acc_ins        Accumulation insoluble mode dry diameter
!> @param[in]     drydp_cor_ins        Coarse insoluble mode dry diameter
!> @param[in]     acc_ins_du           Accumulation insoluble mode dust mass
!!                                      mixing ratio
!> @param[in]     cor_ins_du           Coarse insoluble mode dust mass
!!                                      mixing ratio
!> @param[in]     ndf_wth              Number of degrees of freedom per cell
!!                                      for the potential temperature space
!> @param[in]     undf_wth             Number of unique degrees of freedom
!!                                      for the potential temperature space
!> @param[in]     map_wth              Dofmap for the cell at the base of the
!!                                      column for the potential temperature
!!                                      space
subroutine glomap_ccn_diag_code( nlayers,                                      &
                                 cn_number_conc,                               &
                                 ccn_number_conc_30nm,                         &
                                 ccn_number_conc_50nm,                         &
                                 mconc_du_acc_ins,                             &
                                 mconc_du_cor_ins,                             &
                                 p_zero,                                       &
                                 one_over_kappa,                               &
                                 theta_in_wth,                                 &
                                 exner_in_wth,                                 &
                                 rho_in_wth,                                   &
                                 n_ait_sol,                                    &
                                 n_acc_sol,                                    &
                                 n_cor_sol,                                    &
                                 n_ait_ins,                                    &
                                 n_acc_ins,                                    &
                                 n_cor_ins,                                    &
                                 drydp_ait_sol,                                &
                                 drydp_acc_sol,                                &
                                 drydp_cor_sol,                                &
                                 drydp_ait_ins,                                &
                                 drydp_acc_ins,                                &
                                 drydp_cor_ins,                                &
                                 acc_ins_du,                                   &
                                 cor_ins_du,                                   &
                                 ndf_wth, undf_wth, map_wth )

  use constants_mod,                   only: r_def, i_def
  use science_chemistry_constants_mod, only: boltzmann

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  ! Diagnostic outputs, which point at the shared empty data array when
  ! the diagnostic has not been requested
  real(kind=r_def), pointer, dimension(:), intent(inout) :: cn_number_conc
  real(kind=r_def), pointer, dimension(:), intent(inout) :: ccn_number_conc_30nm
  real(kind=r_def), pointer, dimension(:), intent(inout) :: ccn_number_conc_50nm
  real(kind=r_def), pointer, dimension(:), intent(inout) :: mconc_du_acc_ins
  real(kind=r_def), pointer, dimension(:), intent(inout) :: mconc_du_cor_ins

  real(kind=r_def), intent(in) :: p_zero
  real(kind=r_def), intent(in) :: one_over_kappa

  real(kind=r_def), intent(in), dimension(undf_wth) :: theta_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: exner_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: rho_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_ait_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_acc_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_cor_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_ait_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_acc_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: n_cor_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_ait_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_acc_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_cor_sol
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_ait_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_acc_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: drydp_cor_ins
  real(kind=r_def), intent(in), dimension(undf_wth) :: acc_ins_du
  real(kind=r_def), intent(in), dimension(undf_wth) :: cor_ins_du

  ! Internal variables

  ! The six GLOMAP modes that carry a dry modal diameter in LFRic, ordered
  ! aitken soluble, accumulation soluble, coarse soluble, aitken insoluble,
  ! accumulation insoluble, coarse insoluble
  integer(kind=i_def), parameter :: nmodes_diag = 6

  ! Geometric standard deviation of each mode, taken from the UKCA 7-mode
  ! setup i_sussbcocdu_7mode (ukca_mode_setup sigmag), with the leading
  ! nucleation-soluble entry dropped
  real(kind=r_def), parameter :: sigmag(nmodes_diag) =                        &
      (/ 1.59_r_def, 1.40_r_def, 2.00_r_def,                                  &
         1.59_r_def, 1.59_r_def, 2.00_r_def /)

  ! Dry diameter thresholds of the three diagnostics (m)
  real(kind=r_def), parameter :: dp0_cn   =  3.0e-9_r_def
  real(kind=r_def), parameter :: dp0_30nm = 30.0e-9_r_def
  real(kind=r_def), parameter :: dp0_50nm = 50.0e-9_r_def

  ! Smallest dry diameter admitted, well below any physical particle size.
  ! Guards the logarithm against unset diameters, which are left at zero
  ! until the first RADAER timestep on the climatology path.
  real(kind=r_def), parameter :: drydp_min = 1.0e-12_r_def

  ! Number of cubic centimetres in a cubic metre
  real(kind=r_def), parameter :: m3_to_cm3 = 1.0e+6_r_def

  ! Square root of two, written out for portability of constant
  ! expressions across compilers
  real(kind=r_def), parameter :: root_two = 1.4142135623730951_r_def

  ! Reciprocal of sqrt(2)*ln(sigmag) for each mode
  real(kind=r_def), dimension(nmodes_diag) :: recip_width

  ! Number concentration and dry modal diameter of each mode at one level
  real(kind=r_def), dimension(nmodes_diag) :: number_conc
  real(kind=r_def), dimension(nmodes_diag) :: drydp

  real(kind=r_def) :: exner_k        ! Exner pressure at this level
  real(kind=r_def) :: pressure       ! Pressure at this level (Pa)
  real(kind=r_def) :: temperature    ! Temperature at this level (K)
  real(kind=r_def) :: air_num_dens   ! Number density of air (cm-3)
  real(kind=r_def) :: tail_sum       ! Running sum over the modes

  ! Which diagnostics have real data behind them
  logical :: l_cn
  logical :: l_ccn_30nm
  logical :: l_ccn_50nm
  logical :: l_du_acc_ins
  logical :: l_du_cor_ins
  logical :: l_number_conc

  integer(kind=i_def) :: k, imode

  !---------------------------------------------------------------------------
  ! Determine which diagnostics have been requested
  !---------------------------------------------------------------------------

  l_cn         = .not. associated( cn_number_conc,       empty_real_data )
  l_ccn_30nm   = .not. associated( ccn_number_conc_30nm, empty_real_data )
  l_ccn_50nm   = .not. associated( ccn_number_conc_50nm, empty_real_data )
  l_du_acc_ins = .not. associated( mconc_du_acc_ins,     empty_real_data )
  l_du_cor_ins = .not. associated( mconc_du_cor_ins,     empty_real_data )

  l_number_conc = l_cn .or. l_ccn_30nm .or. l_ccn_50nm

  !---------------------------------------------------------------------------
  ! Lognormal width of each mode, which does not vary in the column
  !---------------------------------------------------------------------------

  if ( l_number_conc ) then
    do imode = 1, nmodes_diag
      recip_width(imode) = 1.0_r_def / ( root_two * log( sigmag(imode) ) )
    end do
  end if

  !---------------------------------------------------------------------------
  ! Column loop
  !---------------------------------------------------------------------------

  do k = 1, nlayers

    exner_k     = exner_in_wth(map_wth(1) + k)
    pressure    = p_zero * exner_k**one_over_kappa
    temperature = exner_k * theta_in_wth(map_wth(1) + k)

    !-------------------------------------------------------------------------
    ! Number concentrations above each dry diameter threshold
    !-------------------------------------------------------------------------

    if ( l_number_conc ) then

      ! Number density of air molecules from pressure and temperature
      air_num_dens = pressure / ( temperature * boltzmann * m3_to_cm3 )

      ! Number mixing ratios are per air molecule, so scaling by the air
      ! number density gives particles per cubic centimetre
      number_conc(1) = n_ait_sol(map_wth(1) + k) * air_num_dens
      number_conc(2) = n_acc_sol(map_wth(1) + k) * air_num_dens
      number_conc(3) = n_cor_sol(map_wth(1) + k) * air_num_dens
      number_conc(4) = n_ait_ins(map_wth(1) + k) * air_num_dens
      number_conc(5) = n_acc_ins(map_wth(1) + k) * air_num_dens
      number_conc(6) = n_cor_ins(map_wth(1) + k) * air_num_dens

      drydp(1) = max( drydp_ait_sol(map_wth(1) + k), drydp_min )
      drydp(2) = max( drydp_acc_sol(map_wth(1) + k), drydp_min )
      drydp(3) = max( drydp_cor_sol(map_wth(1) + k), drydp_min )
      drydp(4) = max( drydp_ait_ins(map_wth(1) + k), drydp_min )
      drydp(5) = max( drydp_acc_ins(map_wth(1) + k), drydp_min )
      drydp(6) = max( drydp_cor_ins(map_wth(1) + k), drydp_min )

      ! Condensation nuclei: dry diameter > 3 nm
      if ( l_cn ) then
        tail_sum = 0.0_r_def
        do imode = 1, nmodes_diag
          tail_sum = tail_sum + 0.5_r_def * number_conc(imode) *              &
              ( 1.0_r_def - erf( log( dp0_cn / drydp(imode) )                 &
                                 * recip_width(imode) ) )
        end do
        cn_number_conc(map_wth(1) + k) = tail_sum
      end if

      ! Cloud condensation nuclei: dry diameter > 30 nm
      if ( l_ccn_30nm ) then
        tail_sum = 0.0_r_def
        do imode = 1, nmodes_diag
          tail_sum = tail_sum + 0.5_r_def * number_conc(imode) *              &
              ( 1.0_r_def - erf( log( dp0_30nm / drydp(imode) )               &
                                 * recip_width(imode) ) )
        end do
        ccn_number_conc_30nm(map_wth(1) + k) = tail_sum
      end if

      ! Cloud condensation nuclei: dry diameter > 50 nm
      if ( l_ccn_50nm ) then
        tail_sum = 0.0_r_def
        do imode = 1, nmodes_diag
          tail_sum = tail_sum + 0.5_r_def * number_conc(imode) *              &
              ( 1.0_r_def - erf( log( dp0_50nm / drydp(imode) )               &
                                 * recip_width(imode) ) )
        end do
        ccn_number_conc_50nm(map_wth(1) + k) = tail_sum
      end if

    end if

    !-------------------------------------------------------------------------
    ! Dust mass concentrations: mass mixing ratio times the dry air density
    !-------------------------------------------------------------------------

    if ( l_du_acc_ins ) then
      mconc_du_acc_ins(map_wth(1) + k) = acc_ins_du(map_wth(1) + k) *         &
                                         rho_in_wth(map_wth(1) + k)
    end if

    if ( l_du_cor_ins ) then
      mconc_du_cor_ins(map_wth(1) + k) = cor_ins_du(map_wth(1) + k) *         &
                                         rho_in_wth(map_wth(1) + k)
    end if

  end do

  !---------------------------------------------------------------------------
  ! The zeroth level is redundant for the GLOMAP fields, so set it to the
  ! same as the first level. It appears in the diagnostic output but is not
  ! used in model evolution.
  !---------------------------------------------------------------------------

  if ( l_cn ) then
    cn_number_conc(map_wth(1)) = cn_number_conc(map_wth(1) + 1)
  end if
  if ( l_ccn_30nm ) then
    ccn_number_conc_30nm(map_wth(1)) = ccn_number_conc_30nm(map_wth(1) + 1)
  end if
  if ( l_ccn_50nm ) then
    ccn_number_conc_50nm(map_wth(1)) = ccn_number_conc_50nm(map_wth(1) + 1)
  end if
  if ( l_du_acc_ins ) then
    mconc_du_acc_ins(map_wth(1)) = mconc_du_acc_ins(map_wth(1) + 1)
  end if
  if ( l_du_cor_ins ) then
    mconc_du_cor_ins(map_wth(1)) = mconc_du_cor_ins(map_wth(1) + 1)
  end if

end subroutine glomap_ccn_diag_code

end module glomap_ccn_diag_kernel_mod
