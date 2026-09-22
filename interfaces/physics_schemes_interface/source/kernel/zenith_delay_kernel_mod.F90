!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Some of the content of this file has been produced with the assistance of
!> Met Office Claude Code Enterprise.
!> @brief Calculates the zenith total delay diagnostic.

module zenith_delay_kernel_mod

  use argument_mod,      only: arg_type,                  &
                               GH_FIELD, GH_SCALAR,       &
                               GH_READ, GH_WRITE,         &
                               GH_REAL, CELL_COLUMN,      &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: zenith_delay_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                              &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),     & ! exner_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! exner_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! mr_v
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! total_mass
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),     & ! height_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA), & ! height_wth
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! p_zero
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! recip_kappa
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! r
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! g
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! grcp
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! c_virtual
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! repsilon
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE,                        &
                  ANY_DISCONTINUOUS_SPACE_1)             & ! zenith_total_delay
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: zenith_delay_code
  end type zenith_delay_kernel_type

  public :: zenith_delay_code

contains

  !> @brief   Calculates the zenith total delay for one column.
  !> @details Port of the UM PWS diagnostic zenith_delay
  !>          (atmosphere/PWS_diagnostics/zenith_delay.F90). The Smith-
  !>          Weintraub refractivity N = nalpha*p/T + nbeta*p*q/(T**2 * ...)
  !>          is evaluated on every theta level (p in hPa) using a layer-mean
  !>          virtual temperature derived hydrostatically from the Exner
  !>          difference across the enclosing rho levels. N is assumed to
  !>          decay exponentially between adjacent theta levels and 1e-6*N is
  !>          integrated in height, with an analytic correction for the
  !>          atmosphere above the model top. Result is in metres.
  !>
  !>          UM reads a prognostic Exner value on an extra rho level above
  !>          the top theta level (model_levels+1) which LFRic does not hold.
  !>          UM defines the top theta-level Exner as the mean of the two
  !>          rho levels either side (calc_exner_at_theta.F90), so that
  !>          extra level is reconstructed here exactly by mirroring the top
  !>          W3 value about the top Wtheta value; the same mirror is used
  !>          for its height, as the UM routine itself does.
  !>
  !>          The 2D output is the last argument (as in pmsl_kernel_mod) so
  !>          that nlayers is taken from a 3D field.
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     exner_w3           Exner pressure on rho levels (W3)
  !> @param[in]     exner_wth          Exner pressure on theta levels (Wtheta)
  !> @param[in]     mr_v               Water vapour mixing ratio (kg/kg)
  !> @param[in]     total_mass         1 + sum of all mixing ratios (kg/kg)
  !> @param[in]     height_w3          Height of rho levels (m)
  !> @param[in]     height_wth         Height of theta levels (m)
  !> @param[in]     p_zero             Reference pressure (Pa)
  !> @param[in]     recip_kappa        1/kappa = cp/R
  !> @param[in]     r                  Gas constant for dry air (J/kg/K)
  !> @param[in]     g                  Gravitational acceleration (m/s**2)
  !> @param[in]     grcp               g/cp (K/m)
  !> @param[in]     c_virtual          1/epsilon - 1
  !> @param[in]     repsilon           epsilon = R/Rv
  !> @param[in,out] zenith_total_delay Zenith total delay (m)
  !> @param[in]     ndf_w3             Number of DOFs per cell for W3
  !> @param[in]     undf_w3            Number of total DOFs for W3
  !> @param[in]     map_w3             Dofmap for the cell at the base of the
  !>                                   column for W3
  !> @param[in]     ndf_wth            Number of DOFs per cell for Wtheta
  !> @param[in]     undf_wth           Number of total DOFs for Wtheta
  !> @param[in]     map_wth            Dofmap for the cell at the base of the
  !>                                   column for Wtheta
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d            Number of total DOFs for 2D fields
  !> @param[in]     map_2d             Dofmap for the cell at the base of the
  !>                                   column for 2D fields
  subroutine zenith_delay_code(nlayers,                         &
                               exner_w3,                        &
                               exner_wth,                       &
                               mr_v,                            &
                               total_mass,                      &
                               height_w3,                       &
                               height_wth,                      &
                               p_zero,                          &
                               recip_kappa,                     &
                               r,                               &
                               g,                               &
                               grcp,                            &
                               c_virtual,                       &
                               repsilon,                        &
                               zenith_total_delay,              &
                               ndf_w3, undf_w3, map_w3,         &
                               ndf_wth, undf_wth, map_wth,      &
                               ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in), dimension(undf_w3)  :: exner_w3
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: mr_v
    real(kind=r_def), intent(in), dimension(undf_wth) :: total_mass
    real(kind=r_def), intent(in), dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in) :: p_zero
    real(kind=r_def), intent(in) :: recip_kappa
    real(kind=r_def), intent(in) :: r
    real(kind=r_def), intent(in) :: g
    real(kind=r_def), intent(in) :: grcp
    real(kind=r_def), intent(in) :: c_virtual
    real(kind=r_def), intent(in) :: repsilon

    real(kind=r_def), intent(inout), dimension(undf_2d) :: zenith_total_delay

    ! Refractivity constants (UM zenith_delay.F90: nalpha, nbeta)
    real(kind=r_def), parameter :: nalpha = 77.6_r_def
    real(kind=r_def), parameter :: nbeta  = 3.73e5_r_def

    ! Internal variables
    integer(kind=i_def) :: k, kth, krho
    real(kind=r_def) :: h_rho_kp1, exner_rho_kp1
    real(kind=r_def) :: p_theta, q, temp1, temp2, tv, t
    real(kind=r_def) :: nwet, ndry, refrac, refrac_kp1, c
    real(kind=r_def) :: zen

    ! UM theta level k (1..model_levels) is Wtheta DOF map_wth(1)+k and UM
    ! rho level k is W3 DOF map_w3(1)+k-1.

    ! Correction for delay above the top of the model, from the pressure on
    ! the top theta level
    kth = map_wth(1) + nlayers
    p_theta = p_zero * exner_wth(kth)**recip_kappa
    zen = 1.0e-6_r_def * nalpha / 100.0_r_def * r * (p_theta / g)

    refrac     = 0.0_r_def
    refrac_kp1 = 0.0_r_def

    do k = nlayers, 1, -1

      kth  = map_wth(1) + k
      krho = map_w3(1) + k - 1

      ! Rho level above this theta level. For the top layer this is UM's
      ! model_levels+1, the mirror of the top rho level about the top theta
      ! level in both height and Exner.
      if (k == nlayers) then
        h_rho_kp1     = 2.0_r_def * height_wth(kth) - height_w3(krho)
        exner_rho_kp1 = 2.0_r_def * exner_wth(kth) - exner_w3(krho)
      else
        h_rho_kp1     = height_w3(krho + 1)
        exner_rho_kp1 = exner_w3(krho + 1)
      end if

      ! Compute mean layer virtual temperature
      temp1 = grcp * (h_rho_kp1 - height_w3(krho))
      temp2 = exner_w3(krho) - exner_rho_kp1
      tv    = exner_wth(kth) * temp1 / temp2

      ! Specific humidity and pressure on this theta level
      q       = mr_v(kth) / total_mass(kth)
      p_theta = p_zero * exner_wth(kth)**recip_kappa

      ! Compute wet refractivity
      t     = tv / (1.0_r_def + c_virtual * q)
      temp1 = 0.01_r_def * nbeta * p_theta * q
      temp2 = t * t * (repsilon + (1.0_r_def - repsilon) * q)
      nwet  = temp1 / temp2

      ! Compute dry refractivity for this level
      ndry = 0.01_r_def * nalpha * p_theta / t

      ! Total refractivity
      refrac_kp1 = refrac
      refrac     = ndry + nwet

      ! The top level only seeds the refractivity of the level above
      if (k == nlayers) cycle

      c = log(refrac_kp1 / refrac) /                                       &
          (height_wth(kth) - height_wth(kth + 1))

      ! Compute zenith delay for this level
      zen = zen + ( -1.0e-6_r_def * refrac *                               &
                    exp(c * height_wth(kth)) *                             &
                    ( exp(-1.0_r_def * c * height_wth(kth + 1)) -          &
                      exp(-1.0_r_def * c * height_wth(kth)) ) / c )

    end do

    zenith_total_delay(map_2d(1)) = zen

  end subroutine zenith_delay_code

end module zenith_delay_kernel_mod
