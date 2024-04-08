!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Calculates double bell mountain orography profile in (lon,lat) coordinates
!>        on the sphere.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of bell case mountain function from
!>          spherical polar coordinates: longitude (lambda) and latitude (phi).
!>          Two ridge mountains of identical width, length and height.
!>          Reference:  Hughes & Jablonowski 2023.
!>          bell mountain parameters in (lon,lat) coordinates are:
!>          mountain_height - Height of bell mountain function (m),
!>          radius_lon - Longitudinal radius of bell mountain function (degrees),
!>          radius_lat - Latitudinal radius of bell mountain function (degrees),
!>          lambda_centre1 - Longitudinal centre of first bell mountain function (degrees),
!>          phi_centre1 - Latitudinal centre of first bell mountain function (degrees).
!>          lambda_centre2 - Longitudinal centre of second bell mountain function (degrees),
!>          phi_centre2 - Latitudinal centre of second bell mountain function (degrees).
!-------------------------------------------------------------------------------
module bell_orography_spherical_mod

  use constants_mod,          only : r_def, i_def, PI
  use analytic_orography_mod, only : analytic_orography_type

  implicit none

  private
  !> @brief Holds parameters and methods used to calculate bell orography
  !>        profile in (lon,lat) coordinates.
  type, public, extends(analytic_orography_type) :: bell_spherical_type

    private
    ! bell mountain function parameters in (lon,lat) coordinates
    real(kind=r_def) :: mountain_height
    real(kind=r_def) :: radius_lon
    real(kind=r_def) :: radius_lat
    real(kind=r_def) :: lambda_centre1
    real(kind=r_def) :: phi_centre1
    real(kind=r_def) :: lambda_centre2
    real(kind=r_def) :: phi_centre2

  contains

    procedure, public, pass(self) :: analytic_orography => bell_orography_spherical
    procedure                     :: bell_coordinate_spherical
    procedure                     :: write_bell_spherical_type

  end type bell_spherical_type

  ! Constructor for bell_spherical_type
  interface bell_spherical_type
    module procedure bell_spherical_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for bell mountain function in (lon,lat) coordinates.
  !> @param[in] mountain_height Height of mountain function read from
  !>                            namelist (m)
  !> @param[in] radius_lon      Longitudinal half-width of mountain function read from
  !>                            namelist (degrees)
  !> @param[in] radius_lat      Latitudinal half-width of mountain function read from
  !>                            namelist (degrees)
  !> @param[in] lambda_centre1  Longitudinal centre of first mountain function read
  !>                            from namelist (degrees)
  !> @param[in] lambda_centre2  Longitudinal centre of first mountain function read
  !>                            from namelist (degrees)
  !> @param[in] phi_centre1     Latitudinal centre of first mountain function read
  !>                            from namelist (degrees)
  !> @param[in] phi_centre2     Latitudinal centre of first mountain function read
  !>                            from namelist (degrees)
  !> @return    self            An object of type bell_spherical_type
  !=============================================================================
  type(bell_spherical_type) function bell_spherical_constructor(                 &
                                                                mountain_height, &
                                                                radius_lon,      &
                                                                radius_lat,      &
                                                                lambda_centre1,  &
                                                                lambda_centre2,  &
                                                                phi_centre1,     &
                                                                phi_centre2 )    &
                                                                result(self)

    use constants_mod, only : PI

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: mountain_height, &
                                    radius_lon,      &
                                    radius_lat,      &
                                    lambda_centre1,  &
                                    lambda_centre2,  &
                                    phi_centre1,     &
                                    phi_centre2

    ! Assign values, converting to radians
    self%mountain_height = mountain_height
    self%radius_lon      = radius_lon*PI/180.0_r_def
    self%radius_lat      = radius_lat*PI/180.0_r_def
    self%lambda_centre1  = lambda_centre1*PI/180.0_r_def
    self%lambda_centre2  = lambda_centre2*PI/180.0_r_def
    self%phi_centre1     = phi_centre1*PI/180.0_r_def
    self%phi_centre2     = phi_centre2*PI/180.0_r_def

    return
  end function bell_spherical_constructor

  !=============================================================================
  !> @brief Calculates bell mountain function in (lon,lat) coordinates.
  !> @param[in] self      An object of type bell_spherical_type
  !> @param[in] chi_1     Longitude (lambda) (radian)
  !> @param[in] chi_2     Latitude (phi) (radian)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function bell_orography_spherical(self, chi_1, chi_2) result(chi_surf)

    implicit none

    ! Arguments
    class(bell_spherical_type), intent(in) :: self
    real(kind=r_def),            intent(in) :: chi_1, chi_2
    real(kind=r_def)                        :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(4)

    ! Calculate transformed/scaled function arguments
    call bell_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    ! Calculate bell mountain surface height
    ! Reference:  Hughes & Jablonowski 2023
    if (chisurf_arg(1) < PI .and. chisurf_arg(2) < PI ) then
     chi_surf = self%mountain_height*(exp(-(chisurf_arg(1))**2 -(chisurf_arg(2))**6))
    else if (chisurf_arg(3) < PI .and. chisurf_arg(4) < PI) then
     chi_surf = self%mountain_height*(exp(-(chisurf_arg(3))**2 -(chisurf_arg(4))**6))
    else
     chi_surf = 0.0_r_def
    end if

    return
  end function bell_orography_spherical

  !=============================================================================
  !> @brief Transforms/scales coordinate for spherical bell mountain function.
  !> @param[in]  self         An object of type bell_spherical_type
  !> @param[in]  chi_1        Longitude (lambda) (radian)
  !> @param[in]  chi_2        Latitude (phi) (radian)
  !> @param[out] chisurf_arg  bell mountain function transformed/scaled
  !>                          arguments
  !=============================================================================
  subroutine bell_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    implicit none

    ! Arguments
    class(bell_spherical_type),  intent(in)  :: self
    real(kind=r_def),      intent(in)        :: chi_1, chi_2
    real(kind=r_def),      intent(out)       :: chisurf_arg(4)

    real(kind=r_def) :: d, c, d1, d2, l1, l2

    ! Initialise transformed/scaled function arguments
    chisurf_arg = 0.0_r_def

    c = self%radius_lon*(-log(0.1_r_def))**(-1.0_r_def/2.0_r_def)
    d = self%radius_lat*(-log(0.1_r_def))**(-1.0_r_def/6.0_r_def)
    d1 = (chi_1 + PI) - self%lambda_centre1
    l1 = min(d1, 2.0_r_def*PI)
    d2 = (chi_1 + PI) - self%lambda_centre2
    l2 = min(d2, 2.0_r_def*PI)
    ! Lon 1
    chisurf_arg(1) = l1/c
    ! Lat 1
    chisurf_arg(2) = (chi_2 - self%phi_centre1)/d
    ! Lon 2
    chisurf_arg(3) = l2/c
    ! Lat 2
    chisurf_arg(4) = (chi_2 - self%phi_centre2)/d


    return
  end subroutine bell_coordinate_spherical

  !=============================================================================
  !> @brief Writes out parameters of bell mountain function in spherical
  !>        coordinates.
  !> @param[in] self An object of type bell_spherical_type
  !=============================================================================
  subroutine write_bell_spherical_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(bell_spherical_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'bell_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "bell mountain parameters in (lon,lat) coordinates: "
    write(funit, fmtreal) "mountain_height     = ", self%mountain_height
    write(funit, fmtreal) "radius_lon          = ", self%radius_lon
    write(funit, fmtreal) "radius_lat          = ", self%radius_lat
    write(funit, fmtreal) "lambda_centre1      = ", self%lambda_centre1
    write(funit, fmtreal) "phi_centre1         = ", self%phi_centre1
    write(funit, fmtreal) "lambda_centre2      = ", self%lambda_centre2
    write(funit, fmtreal) "phi_centre2         = ", self%phi_centre2
    close(funit)

    return
  end subroutine write_bell_spherical_type

end module bell_orography_spherical_mod

