!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Populate socrates arrays in the radiation_fields collection
!>        with the easyaerosol shortwave optical property variables.
!>        Also, dividing absorption and extinction by air density to convert
!>        to suitable units for socrates (from m-1 to m2 kg-1)

module easyaerosol_sw_kernel_mod
use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL, GH_READ,                       &
                             GH_SCALAR, GH_INTEGER,                            &
                             GH_WRITE, CELL_COLUMN,                            &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4

use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: easyaerosol_sw_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                                         &
       ! n_radaer_mode
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! n_radaer_step
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! n_sw_band
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! aer_sw_absorption
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! aer_sw_scattering
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! aer_sw_asymmetry
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! easy_absorption_sw
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),        &
       ! easy_extinction_sw
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),        &
       ! easy_asymmetry_sw
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),        &
       ! rho_wtheta
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),       &
       ! lit_fraction
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4)        &
       /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: easyaerosol_sw_code
end type

public :: easyaerosol_sw_code

contains

!> @brief         Easy Aerosol shortwave kernel required to populate socrates
!>                arrays with optical property values from the climatology.
!> @param[in]     nlayers             The number of layers
!> @param[in]     n_radaer_mode       Number of RADAER modes
!> @param[in]     n_radaer_step       Number of RADAER steps
!> @param[in]     n_sw_band           Number of bands in SW radiation spectrum
!> @param[in,out] aer_sw_absorption   aerosol SW absorption
!> @param[in,out] aer_sw_scattering   aerosol SW scattering
!> @param[in,out] aer_sw_asymmetry    aerosol SW asymmetry
!> @param[in]     easy_absorption_sw  easyaerosol SW absorption
!> @param[in]     easy_extinction_sw  easyaerosol SW scattering
!> @param[in]     easy_asymmetry_sw   easyaerosol SW asymmetry
!> @param[in]     rho_wtheta          Air density on theta levels
!> @param[in]     lit_fraction        Fraction of radiation timestep lit by SW
!> @param[in]     ndf_easy_sw         No. of DOFs per cell for easy_sw space
!> @param[in]     undf_easy_sw        No. unique of DOFs for easy_sw space
!> @param[in]     map_easy_sw         Dofmap for aer_sw space column base cell
!> @param[in]     ndf_aer_sw          No. of DOFs per cell for aer_sw space
!> @param[in]     undf_aer_sw         No. unique of DOFs for aer_sw space
!> @param[in]     map_aer_sw          Dofmap for aer_sw space column base cell
!> @param[in]     ndf_rho             No. DOFs per cell for 2D space
!> @param[in]     undf_rho            No. unique DOFs for rho_wtheta space
!> @param[in]     map_rho             Dofmap for rho_wtheta space
!> @param[in]     ndf_2d              No. DOFs per cell for 2D space
!> @param[in]     undf_2d             No. unique DOFs for 2D space
!> @param[in]     map_2d              Dofmap for 2D space column base cell

subroutine easyaerosol_sw_code( nlayers,                                       &
                                n_radaer_mode,                                 &
                                n_radaer_step,                                 &
                                n_sw_band,                                     &
                                aer_sw_absorption,                             &
                                aer_sw_scattering,                             &
                                aer_sw_asymmetry,                              &
                                easy_absorption_sw,                            &
                                easy_extinction_sw,                            &
                                easy_asymmetry_sw,                             &
                                rho_wtheta,                                    &
                                lit_fraction,                                  &
                                ndf_aer_sw,  undf_aer_sw,  map_aer_sw,         &
                                ndf_easy_sw, undf_easy_sw, map_easy_sw,        &
                                ndf_rho,     undf_rho,     map_rho,            &
                                ndf_2d,      undf_2d,      map_2d )

  use constants_mod,                     only: r_def, i_def, i_um

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: n_radaer_mode
  integer(kind=i_def), intent(in) :: n_radaer_step
  integer(kind=i_def), intent(in) :: n_sw_band

  integer(kind=i_def), intent(in) :: ndf_aer_sw
  integer(kind=i_def), intent(in) :: undf_aer_sw
  integer(kind=i_def), dimension(ndf_aer_sw), intent(in) :: map_aer_sw

  integer(kind=i_def), intent(in) :: ndf_easy_sw
  integer(kind=i_def), intent(in) :: undf_easy_sw
  integer(kind=i_def), dimension(ndf_easy_sw), intent(in) :: map_easy_sw

  integer(kind=i_def), intent(in) :: ndf_rho
  integer(kind=i_def), intent(in) :: undf_rho
  integer(kind=i_def), dimension(ndf_rho), intent(in) :: map_rho

  integer(kind=i_def), intent(in) :: ndf_2d
  integer(kind=i_def), intent(in) :: undf_2d
  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d

  real(kind=r_def), intent(inout), dimension(undf_aer_sw)  :: aer_sw_absorption
  real(kind=r_def), intent(inout), dimension(undf_aer_sw)  :: aer_sw_scattering
  real(kind=r_def), intent(inout), dimension(undf_aer_sw)  :: aer_sw_asymmetry
  real(kind=r_def), intent(in),    dimension(undf_easy_sw) :: easy_absorption_sw
  real(kind=r_def), intent(in),    dimension(undf_easy_sw) :: easy_extinction_sw
  real(kind=r_def), intent(in),    dimension(undf_easy_sw) :: easy_asymmetry_sw
  real(kind=r_def), intent(in),    dimension(undf_rho)     :: rho_wtheta
  real(kind=r_def), intent(in),    dimension(undf_2d)      :: lit_fraction


  ! Loop counters / array index
  integer(i_um) :: k, i_band, i_bnd_mode, i_easy, i_aer, i_rho

  !-----------------------------------------------------------------------
  ! SW aerosol optical properties are copied from
  ! the aerosol_fields arrays to radiation_fields arrays.
  !

  ! Only calculate SW on lit points
  ! If superstepping (n_radaer_step>1) then need to calculate on all points
  ! for use when the sun moves later
  !
  if (lit_fraction(map_2d(1)) > 0.0_r_def .or. n_radaer_step > 1) then
    i_bnd_mode = 0

    do i_band = 1, n_sw_band

      ! The socrates arrays will be populated with radaer outputs
      ! for modes 1 to n_radaer_mode so need to advance i_bnd_mode by
      ! n_radaer_mode + 1 each time this loop increments i_band
      !
      i_bnd_mode = i_bnd_mode + n_radaer_mode + 1

      ! The first model level is the surface where aerosols are zero
      k=0
      i_aer = map_aer_sw(1) + ((i_bnd_mode-1) * (nlayers+1)) + k
      aer_sw_absorption (i_aer) = 0.0_r_def
      aer_sw_scattering (i_aer) = 0.0_r_def
      aer_sw_asymmetry  (i_aer) = 0.0_r_def

      do k = 1, nlayers

        ! The easyaerosol arrays are arranged differently with only one mode
        ! and with the data ordered so that the waveband advances from one
        ! element to the next and data advances through the model levels in
        ! strides of n_sw_band

        i_easy = map_easy_sw(1)  + (k * n_sw_band) + i_band - 1
        i_aer  = map_aer_sw (1)  + ((i_bnd_mode-1) * (nlayers+1)) + k
        i_rho  = map_rho(1)      + k

        ! There is no need to divide absorption and scattering by aerosol
        ! mixing ratio here as mmr is unknown and was set to dummy values of
        ! 1.0 in easyaerosol_aer_kernel_mod.
        !
        aer_sw_absorption (i_aer) = easy_absorption_sw(i_easy)                 &
                                    / rho_wtheta(i_rho)
        aer_sw_scattering (i_aer) = (easy_extinction_sw(i_easy)                &
                                   - easy_absorption_sw(i_easy))               &
                                   / rho_wtheta(i_rho)
        aer_sw_asymmetry  (i_aer) = easy_asymmetry_sw (i_easy)

      end do
    end do

  else ! unlit points

    ! Dummy values to avoid problems in radiation code
    i_bnd_mode = 0
    do i_band = 1, n_sw_band

      i_bnd_mode = i_bnd_mode + n_radaer_mode + 1

      do k = 0, nlayers

        i_aer = map_aer_sw(1) + ((i_bnd_mode-1) * (nlayers+1)) + k

        aer_sw_absorption (i_aer) =  0.0_r_def
        aer_sw_scattering (i_aer) =  0.0_r_def
        aer_sw_asymmetry  (i_aer) =  0.0_r_def

      end do
    end do
  end if ! lit points

end subroutine easyaerosol_sw_code

end module easyaerosol_sw_kernel_mod
