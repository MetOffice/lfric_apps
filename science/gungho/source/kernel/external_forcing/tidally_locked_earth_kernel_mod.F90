!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adds a Tidally Locked Earth (TLE) temperature forcing to the
!!        finite-difference representation of the fields.
!>
!> @details Kernel that adds the TLE forcing to the temperature field.
!!          References:
!!          * Mayne, N. J., Baraffe, I., Acreman, D. M., Smith, C., Wood, N.,
!!          Amundsen, D. S., Thuburn, J., and Jackson, D. R.: Using the UM
!!          dynamical cores to reproduce idealised 3-D flows, Geosci. Model
!!          Dev., 7, 3059-3087, https://doi.org/10.5194/gmd-7-3059-2014, 2014.
!!          * Merlis, T. M., and Schneider, T. (2010), Atmospheric Dynamics of
!!          Earth-Like Tidally Locked Aquaplanets, J. Adv. Model. Earth Syst.,
!!          2, 13, doi:10.3894/JAMES.2010.2.13.
module tidally_locked_earth_kernel_mod

  use argument_mod,                      only: arg_type,                  &
                                               GH_FIELD, GH_REAL,         &
                                               GH_READ, GH_READWRITE,     &
                                               GH_SCALAR, ANY_SPACE_9,    &
                                               ANY_DISCONTINUOUS_SPACE_3, &
                                               GH_READ, CELL_COLUMN
  use constants_mod,                     only: r_def, i_def
  use sci_chi_transform_mod,             only: chi2llr
  use calc_exner_pointwise_mod,          only: calc_exner_pointwise
  use fs_continuity_mod,                 only: Wtheta
  use tidally_locked_earth_forcings_mod, only: &
    tidally_locked_earth_equilibrium_theta
  use held_suarez_forcings_mod,          only: held_suarez_newton_frequency
  use kernel_mod,                        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: tidally_locked_earth_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                          &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,      ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                                 &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                                  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tidally_locked_earth_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tidally_locked_earth_code

contains

!> @brief Adds the TLE forcing based on Merlis & Schneider (2010)
!!        and Mayne et al. (2014).
!> @param[in]     nlayers      The number of layers
!> @param[in,out] dtheta       Potential temperature increment data
!> @param[in]     theta        Potential temperature data
!> @param[in]     exner_in_wth The Exner pressure in Wtheta
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
subroutine tidally_locked_earth_code(nlayers,                     &
                                     dtheta, theta, exner_in_wth, &
                                     chi_1, chi_2, chi_3,         &
                                     panel_id, kappa, dt,         &
                                     ndf_wth, undf_wth, map_wth,  &
                                     ndf_chi, undf_chi, map_chi,  &
                                     ndf_pid, undf_pid, map_pid   &
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

  do k = 0, nlayers
    exner = exner_in_wth(map_wth(1) + k)

    sigma = (exner/exner0)**(1.0_r_def/kappa)

    theta_eq = tidally_locked_earth_equilibrium_theta(exner, lon, lat)

    dtheta(map_wth(1) + k) = -held_suarez_newton_frequency(sigma, lat)    &
       * (theta(map_wth(1) + k) - theta_eq) * dt
  end do

end subroutine tidally_locked_earth_code

end module tidally_locked_earth_kernel_mod
