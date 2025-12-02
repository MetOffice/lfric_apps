!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Calculates DCMIP200 mountain orography profile in (lon,lat) coordinates.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of DCMIP200 case mountain function from
!>          spherical polar coordinates: longitude (lambda) and latitude (phi).
!>          Reference: Ulrich et al. (2012), Section 2.0.
!>          DCMIP200 mountain parameters in (lon,lat) coordinates are:
!>          mountain_height - Height of DCMIP200 mountain function (m),
!>          radius - Radius of DCMIP200 mountain function (radian),
!>          osc_half_width - Oscillation half-width of DCMIP200 mountain (radian),
!>          lambda_centre - Longitudinal centre of DCMIP200 mountain function (radian),
!>          phi_centre - Latitudinal centre of DCMIP200 mountain function (radian).
!-------------------------------------------------------------------------------
module dcmip200_orography_spherical_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type

  implicit none

  private
  !> @brief Holds parameters and methods used to calculate DCMIP200 orography
  !>        profile in (lon,lat) coordinates.
  type, public, extends(analytic_orography_type) :: dcmip200_spherical_type

    private
    ! DCMIP200 mountain function parameters in (lon,lat) coordinates
    real(kind=r_def) :: mountain_height
    real(kind=r_def) :: radius
    real(kind=r_def) :: osc_half_width
    real(kind=r_def) :: lambda_centre
    real(kind=r_def) :: phi_centre

  contains

    procedure, public, pass(self) :: analytic_orography => dcmip200_orography_spherical
    procedure                     :: dcmip200_coordinate_spherical
    procedure                     :: write_dcmip200_spherical_type

  end type dcmip200_spherical_type

  ! Constructor for dcmip200_spherical_type
  interface dcmip200_spherical_type
    module procedure dcmip200_spherical_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for DCMIP200 mountain function in (lon,lat) coordinates.
  !> @param[in] mountain_height Height of mountain function read from
  !>                            namelist (m)
  !> @param[in] radius          Half-width of DCMIP mountain function read from
  !>                            namelist (radian)
  !> @param[in] osc_half_width  Oscillation half-width of mountain function
  !>                            read from namelist (radian)
  !> @param[in] lambda_centre   Longitudinal centre of mountain function read
  !>                            from namelist (radian)
  !> @param[in] phi_centre      Latitudinal centre of mountain function read
  !>                            from namelist (radian)
  !> @return    self            An object of type dcmip200_spherical_type
  !=============================================================================
  type(dcmip200_spherical_type) function dcmip200_spherical_constructor(                 &
                                                                        mountain_height, &
                                                                        radius,          &
                                                                        osc_half_width,  &
                                                                        lambda_centre,   &
                                                                        phi_centre )     &
                                                                        result(self)

    use constants_mod, only : PI

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: mountain_height, &
                                    radius,          &
                                    osc_half_width,  &
                                    lambda_centre,   &
                                    phi_centre

    ! Assign values
    self%mountain_height = mountain_height
    self%radius          = radius
    self%osc_half_width  = osc_half_width
    self%lambda_centre   = lambda_centre
    self%phi_centre      = phi_centre

    return
  end function dcmip200_spherical_constructor

  !=============================================================================
  !> @brief Calculates DCMIP200 mountain function in (lon,lat) coordinates.
  !> @param[in] self      An object of type dcmip200_spherical_type
  !> @param[in] chi_1     Longitude (lambda) (radian)
  !> @param[in] chi_2     Latitude (phi) (radian)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function dcmip200_orography_spherical(self, chi_1, chi_2) result(chi_surf)

    use constants_mod,     only : PI

    implicit none

    ! Arguments
    class(dcmip200_spherical_type), intent(in) :: self
    real(kind=r_def),            intent(in) :: chi_1, chi_2
    real(kind=r_def)                        :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(2)

    ! Calculate transformed/scaled function arguments
    call dcmip200_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    ! Calculate DCMIP200 mountain surface height
    ! Reference: Ulrich et al. (2012), Section 2.0, Eq. 63
    if (chisurf_arg(1) < PI) then
      chi_surf = 0.5_r_def * self%mountain_height * &
                   (1.0_r_def + cos(chisurf_arg(1)))*(cos(chisurf_arg(2)))**2
    else
      chi_surf = 0.0_r_def
    end if

    return
  end function dcmip200_orography_spherical

  !=============================================================================
  !> @brief Transforms/scales coordinate for spherical DCMIP200 mountain function.
  !> @param[in]  self         An object of type dcmip200_spherical_type
  !> @param[in]  chi_1        Longitude (lambda) (radian)
  !> @param[in]  chi_2        Latitude (phi) (radian)
  !> @param[out] chisurf_arg  DCMIP200 mountain function transformed/scaled
  !>                          arguments
  !=============================================================================
  subroutine dcmip200_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    use constants_mod,     only : PI

    implicit none

    ! Arguments
    class(dcmip200_spherical_type),  intent(in)  :: self
    real(kind=r_def),      intent(in)  :: chi_1, chi_2
    real(kind=r_def),      intent(out) :: chisurf_arg(2)
    ! Internal variables
    real(kind=r_def) :: chi_dcmip200

    ! Initialise transformed/scaled function arguments
    chisurf_arg = 0.0_r_def

    ! Calculate transformed/scaled arguments
    ! Reference: Ulrich et al. (2012), Section 2.0, Eq. 64
    chi_dcmip200 = acos( sin(chi_2)*sin(self%phi_centre) + &
                     cos(chi_2)*cos(self%phi_centre)*cos(chi_1 - self%lambda_centre) )
    ! Cosine function argument
    chisurf_arg(1) = PI*chi_dcmip200/self%radius
    ! Cosine squared function argument
    chisurf_arg(2) = PI*chi_dcmip200/self%osc_half_width

    return
  end subroutine dcmip200_coordinate_spherical

  !=============================================================================
  !> @brief Writes out parameters of DCMIP200 mountain function in spherical
  !>        coordinates.
  !> @param[in] self An object of type dcmip200_spherical_type
  !=============================================================================
  subroutine write_dcmip200_spherical_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(dcmip200_spherical_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'dcmip200_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "DCMIP200 mountain parameters in (lon,lat) coordinates: "
    write(funit, fmtreal) "mountain_height     = ", self%mountain_height
    write(funit, fmtreal) "radius              = ", self%radius
    write(funit, fmtreal) "osc_half_width      = ", self%osc_half_width
    write(funit, fmtreal) "lambda_centre       = ", self%lambda_centre
    write(funit, fmtreal) "phi_centre          = ", self%phi_centre
    close(funit)

    return
  end subroutine write_dcmip200_spherical_type

end module dcmip200_orography_spherical_mod

