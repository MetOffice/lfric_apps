!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains forcing terms for use in the Earth-Like kernels.
!>
!> @details Support functions for a kernel that adds the Earth-Like test
!!          based on Menou & Rauscher (2009),
!!          Atmospheric Circulation of Hot Jupiters: A Shallow Three-Dimensional Model,
!!          ApJ, 700, 887-897, 2009, DOI: 10.1088/0004-637X/700/1/887.
!!          Also performed in Mayne et al., (2014),
!!          Using the UM dynamical cores to reproduce idealised 3-D flows,
!!          Geoscientific Model Development, Volume 7, Issue 6, 2014, pp. 3059-3087,
!!          DOI: 10.5194/gmd-7-3059-2014.

module earth_like_forcings_mod

  use constants_mod,     only: r_def, i_def, pi
  use planet_config_mod, only: scaling_factor

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Local parameters
  !-------------------------------------------------------------------------------
  ! Earth-Like parameters
  real(kind=r_def), parameter :: TAU_NEWTON       = 1.296e6_r_def ! Heng et al (2011) use constant tau_rad
                                                                  ! of 15 days, therefore 1.296e6 seconds
  real(kind=r_def), parameter :: T_SURF           = 288.0_r_def   ! Surface temperature
  real(kind=r_def), parameter :: DT_EQ_POLE       = 60.0_r_def    ! Equator-Pole temperature diff (deltaT)_y
  real(kind=r_def), parameter :: STATIC_STABILITY = 2.0_r_def     ! Static stability temperature (delta \theta)_z
  real(kind=r_def), parameter :: GAMMA_TROP       = 6.5e-3_r_def  ! Lapse rate

  real(kind=r_def) :: t_vert    ! Stratospheric temperature profile
  real(kind=r_def) :: beta_trop ! Horizontal temperature gradient in troposphere

  public :: earth_like_newton_frequency
  public :: earth_like_equilibrium_theta

contains

!> @brief Function to calculate equilibrium theta profile for Earth-Like temperature forcing.
!> @param[in] exner         Exner pressure
!> @param[in] exner0        Exner pressure at the surface
!> @param[in] exner_strat   Exner pressure at the stratosphere
!> @param[in] lat           Latitude
!> @param[in] sigma         exner/exner0**(1.0/kappa)
!> @param[in] sigma_strat   Sigma value at the stratosphere
!> @param[in] strat_loc     Location of the stratosphere
!> @param[in] kappa         Ratio of Rd and cp
!> @param[in] layer_height  Height of the kth layer above the surface
!> @param[in] z_strat       Height (m) of stratosphere
!> @return    theta_eq      Equilibrium theta
function earth_like_equilibrium_theta(exner, exner0, exner_strat, lat, sigma, sigma_strat, strat_loc, &
                                      kappa, layer_height, z_strat) result(theta_eq)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: strat_loc
  real(kind=r_def), intent(in)    :: z_strat
  real(kind=r_def), intent(in)    :: exner, exner0, exner_strat
  real(kind=r_def), intent(in)    :: sigma, sigma_strat
  real(kind=r_def), intent(in)    :: lat, kappa, layer_height
  real(kind=r_def)                :: theta_eq ! Equilibrium theta

  if (layer_height <= z_strat) then

    ! In the stratosphere
    beta_trop = sin((pi*(sigma-sigma_strat))/(2.0_r_def*(1.0_r_def-sigma_strat)))

    t_vert = T_SURF-GAMMA_TROP*(z_strat+(layer_height-z_strat)*0.5_r_def) &
    + ((GAMMA_TROP*(layer_height-z_strat)*0.5_r_def)**2.0_r_def+STATIC_STABILITY**2.0_r_def)**0.5_r_def

  else

    ! Above stratosphere
    beta_trop = 0.0_r_def
    t_vert = T_SURF-GAMMA_TROP*z_strat+STATIC_STABILITY

  end if

  ! Recall using potential temperature
  ! Therefore, must convert the temperature to potential
  ! temperature using the Exner function (exner*theta=Temp)
  theta_eq = (t_vert + beta_trop * DT_EQ_POLE &
           * ((1.0_r_def/3.0_r_def)-sin(lat)*sin(lat)))/exner

end function earth_like_equilibrium_theta

!> @brief Function to calculate the Newton relaxation frequency for Earth-Like idealised test case.
function earth_like_newton_frequency() result(earth_like_frequency)

  implicit none

  real(kind=r_def) :: earth_like_frequency

  ! Set as a constant
  earth_like_frequency = 1.0_r_def/TAU_NEWTON

end function earth_like_newton_frequency

end module earth_like_forcings_mod
