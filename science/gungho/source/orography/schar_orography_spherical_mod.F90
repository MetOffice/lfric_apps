!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates Schar mountain orography profile in (lon,lat) coordinates.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of Schar mountain function from spherical
!>          polar coordinates: longitude (lambda) and latitude (phi).
!>          Reference: Wood et al. (2013), Section 7.1.
!>          Schar mountain parameters in (lon,lat) coordinates are:
!>          mountain_height - Height of Schar mountain function (m),
!>          half_width - Half-width of Schar mountain function (m),
!>          wavelength - Wavelength of cosine part of Schar mountain function (m),
!>          lambda_centre - Longitudinal centre of Schar mountain function (radian),
!>          phi_centre - Latitudinal centre of Schar mountain function (radian).
!-------------------------------------------------------------------------------
module schar_orography_spherical_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type

  implicit none

  private
  !> @brief Holds parameters and methods used to calculate Schar orography
  !>        profile in (lon,lat) coordinates.
  type, public, extends(analytic_orography_type) :: schar_spherical_type

    private
    ! Schar mountain function parameters in (lon,lat) coordinates
    real(kind=r_def) :: mountain_height
    real(kind=r_def) :: half_width
    real(kind=r_def) :: wavelength
    real(kind=r_def) :: lambda_centre
    real(kind=r_def) :: phi_centre

  contains

    procedure, public, pass(self) :: analytic_orography => schar_orography_spherical
    procedure                     :: schar_coordinate_spherical
    procedure                     :: write_schar_spherical_type

  end type schar_spherical_type

  ! Constructor for schar_spherical_type
  interface schar_spherical_type
    module procedure schar_spherical_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for Schar mountain function in (lon,lat) coordinates.
  !> @param[in] mountain_height Height of mountain function read from
  !>                            namelist (m)
  !> @param[in] half_width      Half-width of mountain function read from
  !>                            namelist (m)
  !> @param[in] wavelength      Wavelength of cosine part of  mountain function
  !>                            read from namelist (m)
  !> @param[in] lambda_centre   Longitudinal centre of mountain function read
  !>                            from namelist (m)
  !> @param[in] phi_centre      Longitudinal centre of mountain function read
  !>                            from namelist (m)
  !> @return    self            An object of type schar_spherical_type
  !=============================================================================
  type(schar_spherical_type) function schar_spherical_constructor(     &
                                                      mountain_height, &
                                                      half_width,      &
                                                      wavelength,      &
                                                      lambda_centre,   &
                                                      phi_centre )     &
                                                      result(self)

    use constants_mod, only : PI

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: mountain_height, &
                                    half_width,      &
                                    wavelength,      &
                                    lambda_centre,   &
                                    phi_centre

    ! Assign values
    self%mountain_height = mountain_height
    self%half_width      = half_width
    self%wavelength      = wavelength
    self%lambda_centre   = lambda_centre
    self%phi_centre      = phi_centre

    return
  end function schar_spherical_constructor

  !=============================================================================
  !> @brief Calculates Schar mountain function in (lon,lat) coordinates.
  !> @param[in] self      An object of type schar_spherical_type
  !> @param[in] chi_1     Longitude (lambda) (radian)
  !> @param[in] chi_2     Latitude (phi) (radian)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function schar_orography_spherical(self, chi_1, chi_2) result(chi_surf)

    implicit none

    ! Arguments
    class(schar_spherical_type), intent(in) :: self
    real(kind=r_def),            intent(in) :: chi_1, chi_2
    real(kind=r_def)                        :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(2)

    ! Calculate transformed/scaled function arguments
    call schar_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    ! Calculate Schar mountain surface height
    ! Reference: Wood et al. (2013), Section 7.1., Eq. 78
    chi_surf = self%mountain_height*exp(-chisurf_arg(1)**2)*(cos(chisurf_arg(2)))**2

    return
  end function schar_orography_spherical

  !=============================================================================
  !> @brief Transforms/scales coordinate for spherical Schar mountain function.
  !> @param[in]  self         An object of type schar_spherical_type
  !> @param[in]  chi_1        Longitude (lambda) (radian)
  !> @param[in]  chi_2        Latitude (phi) (radian)
  !> @param[out] chisurf_arg  Schar mountain function transformed/scaled
  !>                          arguments
  !=============================================================================
  subroutine schar_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    use constants_mod,     only : PI
    use planet_config_mod, only : scaled_radius

    implicit none

    ! Arguments
    class(schar_spherical_type), intent(in)  :: self
    real(kind=r_def),            intent(in)  :: chi_1, chi_2
    real(kind=r_def),            intent(out) :: chisurf_arg(2)
    ! Internal variables
    real(kind=r_def) :: chi_schar

    ! Initialise transformed/scaled function arguments
    chisurf_arg = 0.0_r_def

    ! Calculate transformed/scaled arguments
    ! Reference: Wood et al. (2013), Section 7., Eq. 76
    chi_schar = scaled_radius*acos( &
                  sin(chi_2)*sin(self%phi_centre) + &
                  cos(chi_2)*cos(self%phi_centre)*cos(chi_1 - self%lambda_centre) )
    ! Exponential function argument
    chisurf_arg(1) = chi_schar/self%half_width
    ! Cosine function argument
    chisurf_arg(2) = PI*chi_schar/self%wavelength

    return
  end subroutine schar_coordinate_spherical

  !=============================================================================
  !> @brief Writes out parameters of Schar mountain function in spherical
  !>        coordinates.
  !> @param[in] self An object of type schar_spherical_type
  !=============================================================================
  subroutine write_schar_spherical_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(schar_spherical_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'schar_sphere_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "Schar mountain parameters in (lon,lat) coordinates: "
    write(funit, fmtreal) "mountain_height = ", self%mountain_height
    write(funit, fmtreal) "half_width      = ", self%half_width
    write(funit, fmtreal) "wavelength      = ", self%wavelength
    write(funit, fmtreal) "lambda_centre   = ", self%lambda_centre
    write(funit, fmtreal) "phi_centre      = ", self%phi_centre
    close(funit)

    return
  end subroutine write_schar_spherical_type

end module schar_orography_spherical_mod

