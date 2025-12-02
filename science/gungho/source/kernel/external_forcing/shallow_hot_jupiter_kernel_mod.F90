!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adds a shallow hot Jupiter forcing using the finite-difference
!!        representation of the fields.
!>
!> @details Kernel that adds the shallow hot Jupiter test based on Menou & Rauscher (2009),
!!          Atmospheric Circulation of Hot Jupiters: A Shallow Three-Dimensional Model,
!!          ApJ, 700, 887-897, 2009, DOI: 10.1088/0004-637X/700/1/887.
!!          Also performed in Mayne et al., (2014),
!!          Using the UM dynamical cores to reproduce idealised 3-D flows,
!!          Geoscientific Model Development, Volume 7, Issue 6, 2014, pp. 3059-3087,
!!          DOI: 10.5194/gmd-7-3059-2014.
!>
module shallow_hot_jupiter_kernel_mod

  use argument_mod,                     only: arg_type,                  &
                                              GH_FIELD, GH_REAL,         &
                                              GH_READ, GH_READWRITE,     &
                                              GH_SCALAR,                 &
                                              ANY_DISCONTINUOUS_SPACE_3, &
                                              GH_READ, CELL_COLUMN
  use constants_mod,                    only: r_def, i_def
  use sci_chi_transform_mod,            only: chi2llr
  use calc_exner_pointwise_mod,         only: calc_exner_pointwise
  use fs_continuity_mod,                only: Wtheta, Wchi
  use shallow_hot_jupiter_forcings_mod, only: shallow_hot_jupiter_newton_frequency, &
                                              shallow_hot_jupiter_equilibrium_theta
  use kernel_mod,                       only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: shallow_hot_jupiter_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                          &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,      Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                                 &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                                  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: shallow_hot_jupiter_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: shallow_hot_jupiter_code

contains

!> @brief Adds the shallow hot Jupiter test based on Menou & Rauscher (2009)
!!        and Mayne et al. (2014).
!> @param[in]     nlayers      The number of layers
!> @param[in,out] dtheta       Potential temperature increment data
!> @param[in]     theta        Potential temperature data
!> @param[in]     exner_in_wth The Exner pressure in Wtheta
!> @param[in]     height_wth   Height of Wtheta space levels above surface
!> @param[in]     chi_1        First component of the chi coordinate field
!> @param[in]     chi_2        Second component of the chi coordinate field
!> @param[in]     chi_3        Third component of the chi coordinate field
!> @param[in]     panel_id     A field giving the ID for mesh panels
!> @param[in]     kappa        Ratio of Rd and cp
!> @param[in]     dt           The model timestep length
!> @param[in]     ndf_wth      The number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wth     The number of unique degrees of freedom for Wtheta
!> @param[in]     map_wth      Dofmap for the cell at the base of the column for Wtheta
!> @param[in]     ndf_chi      The number of degrees of freedom per cell for Wchi
!> @param[in]     undf_chi     The number of unique degrees of freedom for Wchi
!> @param[in]     map_chi      Dofmap for the cell at the base of the column for Wchi
!> @param[in]     ndf_pid      Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid     Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid      Dofmap for the cell at the base of the column for panel_id
subroutine shallow_hot_jupiter_code(nlayers,                    &
                                    dtheta, theta,              &
                                    exner_in_wth, height_wth,   &
                                    chi_1, chi_2, chi_3,        &
                                    panel_id, kappa, dt,        &
                                    ndf_wth, undf_wth, map_wth, &
                                    ndf_chi, undf_chi, map_chi, &
                                    ndf_pid, undf_pid, map_pid  &
                                   )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid

  real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta
  real(kind=r_def), dimension(undf_wth), intent(in)    :: theta
  real(kind=r_def), dimension(undf_wth), intent(in)    :: exner_in_wth
  real(kind=r_def), dimension(undf_wth), intent(in)    :: height_wth
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def),                      intent(in)    :: kappa
  real(kind=r_def),                      intent(in)    :: dt

  integer(kind=i_def), dimension(ndf_wth),  intent(in) :: map_wth
  integer(kind=i_def), dimension(ndf_chi),  intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in) :: map_pid

  ! Internal variables
  integer(kind=i_def) :: k, df, location, ipanel

  real(kind=r_def)    :: theta_eq, exner
  real(kind=r_def)    :: lat, lon, radius

  real(kind=r_def) :: exner0 ! Lowest-level Exner value
  real(kind=r_def) :: sigma  ! exner/exner0**(1.0/kappa)

  real(kind=r_def) :: coords(3)
  real(kind=r_def), dimension(ndf_chi) :: chi_1_at_dof, chi_2_at_dof, chi_3_at_dof

  ! Local parameters
  real(kind=r_def), parameter :: Z_STRAT = 2.0e6_r_def ! Height (m) of Stratosphere

  ! Local variables
  integer(kind=i_def) :: strat_loc        ! The location of the stratosphere (index)
  real(kind=r_def)    :: newton_frequency ! Relaxation timescale
  real(kind=r_def)    :: sigma_strat      ! Sigma value of the stratosphere
  real(kind=r_def)    :: exner_strat      ! Exner value at the location of stratosphere
  real(kind=r_def)    :: layer_height     ! Height of the kth layer above the surface

  coords(:) = 0.0_r_def

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Calculate x, y and z at the centre of the lowest cell
  do df = 1, ndf_chi
    location = map_chi(df)
    chi_1_at_dof(df) = chi_1( location )
    chi_2_at_dof(df) = chi_2( location )
    chi_3_at_dof(df) = chi_3( location )
    coords(1) = coords(1) + chi_1( location )/ndf_chi
    coords(2) = coords(2) + chi_2( location )/ndf_chi
    coords(3) = coords(3) + chi_3( location )/ndf_chi
  end do

  call chi2llr(coords(1), coords(2), coords(3), ipanel, lon, lat, radius)

  exner0 = exner_in_wth(map_wth(1))

  ! Find location at the stratosphere
  strat_loc = minloc(abs(height_wth(map_wth(1):map_wth(1) + nlayers) - Z_STRAT), 1_i_def)
  exner_strat = exner_in_wth(map_wth(1) + strat_loc)

  ! Find sigma at the stratosphere
  sigma_strat = (exner_strat/exner0)**(1.0_r_def/kappa)

  ! Set the relaxation timescale
  newton_frequency = shallow_hot_jupiter_newton_frequency()

  do k = 0, nlayers

    exner = exner_in_wth(map_wth(1) + k)

    sigma = (exner/exner0)**(1.0_r_def/kappa)

    layer_height = height_wth(map_wth(1) + k)

    theta_eq = shallow_hot_jupiter_equilibrium_theta(exner, exner0, exner_strat, lat, lon, sigma, sigma_strat, &
                                                     strat_loc, kappa, layer_height, Z_STRAT)

    dtheta(map_wth(1) + k) = -newton_frequency &
        * (theta(map_wth(1) + k) - theta_eq) * dt

  end do

end subroutine shallow_hot_jupiter_code

end module shallow_hot_jupiter_kernel_mod
