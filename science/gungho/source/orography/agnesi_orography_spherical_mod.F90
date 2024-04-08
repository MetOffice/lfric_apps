!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates Witch-of-Agnesi mountain orography profile in spherical
!>        coordinates.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of Witch-of-Agnesi mountain function from
!>          spherical polar coordinates: longitude (lambda) and latitude (phi).
!>          Reference: Wedi and Smolarkiewicz (2009), Section 4.1.
!>          Witch-of-Agnesi mountain parameters in (lon,lat) coordinates are:
!>          mountain_height - Height of Witch-of-Agnesi mountain function (m),
!>          half_width - Half-width of Witch-of-Agnesi mountain function (m),
!>          lambda_centre - Longitudinal centre of Witch-of-Agnesi mountain
!>                          function (radian),
!>          phi_centre - Latitudinal centre of Witch-of-Agnesi mountain
!>                       function (radian),
!>          lambda_focus - Longitudinal parameter of Witch-of-Agnesi mountain
!>                         function's centre (radian),
!>          phi_focus - Latitudinal parameter of Witch-of-Agnesi mountain
!>                     function's centre (radian).
!-------------------------------------------------------------------------------
module agnesi_orography_spherical_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type

  implicit none

  private
  !> @brief Holds parameters and methods used to calculate Witch-of-Agnesi
  !>        orography profile in (lon,lat) coordinates.
  type, public, extends(analytic_orography_type) :: agnesi_spherical_type

    private
    ! Witch-of-Agnesi mountain function parameters in (lon,lat) coordinates
    real(kind=r_def) :: mountain_height
    real(kind=r_def) :: half_width
    real(kind=r_def) :: lambda_centre
    real(kind=r_def) :: phi_centre
    real(kind=r_def) :: lambda_focus
    real(kind=r_def) :: phi_focus

  contains

    procedure, public, pass(self) :: analytic_orography => agnesi_orography_spherical
    procedure                     :: agnesi_coordinate_spherical
    procedure                     :: write_agnesi_spherical_type

  end type agnesi_spherical_type

  ! Constructor for agnesi_spherical_type
  interface agnesi_spherical_type
    module procedure agnesi_spherical_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for the Witch-of-Agnesi mountain function in spherical
  !>        coordinates.
  !>
  !> @param[in] mountain_height Height of mountain function read from
  !>                            namelist (m)
  !> @param[in] half_width      Half-width of mountain function read from
  !>                            namelist (m)
  !> @param[in] lambda_centre   Longitudinal centre of mountain function read
  !>                            from namelist (m)
  !> @param[in] phi_centre      Longitudinal centre of mountain function read
  !>                            from namelist (m)
  !> @param[in] lambda_focus    Longitudinal parameter of mountain function read
  !>                            from namelist (m)
  !> @param[in] phi_focus       Latitudinal parameter of mountain function read
  !>                            from namelist (m)
  !> @return    self            An object of type agnesi_spherical_type
  !=============================================================================
  type(agnesi_spherical_type) function agnesi_spherical_constructor(     &
                                                        mountain_height, &
                                                        half_width,      &
                                                        lambda_centre,   &
                                                        phi_centre,      &
                                                        lambda_focus,    &
                                                        phi_focus )      &
                                                        result(self)

    use constants_mod, only : PI

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: mountain_height, &
                                    half_width,      &
                                    lambda_centre,   &
                                    phi_centre,      &
                                    lambda_focus,    &
                                    phi_focus

    ! Assign values
    self%mountain_height = mountain_height
    self%half_width      = half_width
    self%lambda_centre   = lambda_centre
    self%phi_centre      = phi_centre
    self%lambda_focus    = lambda_focus
    self%phi_focus       = phi_focus

    return
  end function agnesi_spherical_constructor

  !=============================================================================
  !> @brief Calculates Witch-of-Agnesi mountain function in (lon,lat) coordinates.
  !>
  !> @param[in] self      An object of type agnesi_spherical_type
  !> @param[in] chi_1     Longitude (lambda) (radian)
  !> @param[in] chi_2     Latitude (phi) (radian)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function agnesi_orography_spherical(self, chi_1, chi_2) result(chi_surf)

    implicit none

    ! Arguments
    class(agnesi_spherical_type), intent(in) :: self
    real(kind=r_def),             intent(in) :: chi_1, chi_2
    real(kind=r_def)                         :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(2)

    ! Calculate transformed/scaled function arguments
    call agnesi_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    ! Calculate Witch-of-Agnesi mountain surface height
    ! Reference: Wedi and Smolarkiewicz (2009), Section 4.1., Eq. 9
    chi_surf = self%mountain_height/(1.0_r_def + sum(chisurf_arg**2))

    return
  end function agnesi_orography_spherical

  !=============================================================================
  !> @brief Transforms/scales coordinate for spherical Witch-of-Agnesi mountain
  !>        function.
  !>
  !> @param[in]  self         An object of type agnesi_spherical_type
  !> @param[in]  chi_1        Longitude (lambda) (radian)
  !> @param[in]  chi_2        Latitude (phi) (radian)
  !> @param[out] chisurf_arg  Witch-of-Agnesi mountain function
  !>                          transformed/scaled arguments
  !=============================================================================
  subroutine agnesi_coordinate_spherical(self, chi_1, chi_2, chisurf_arg)

    use planet_config_mod, only : scaled_radius

    implicit none

    ! Arguments
    class(agnesi_spherical_type), intent(in)  :: self
    real(kind=r_def),             intent(in)  :: chi_1, chi_2
    real(kind=r_def),             intent(out) :: chisurf_arg(2)
    ! Internal variables
    real(kind=r_def) :: half_width_agnesi(2), half_width_focus
    real(kind=r_def) :: sinp_foc, cosp_foc, sinp_cen, cosp_cen

    ! Initialise transformed/scaled function arguments
    chisurf_arg = 0.0_r_def

    ! Calculate transformed/scaled arguments
    ! Reference: Wedi and Smolarkiewicz (2009), Section 4.1., after Eq. 9
    sinp_foc = sin(self%phi_focus)
    cosp_foc = cos(self%phi_focus)
    sinp_cen = sin(self%phi_centre)
    cosp_cen = cos(self%phi_centre)
    ! Focal half-width
    half_width_focus = scaled_radius*acos( sinp_foc*sinp_cen +     &
                         cosp_foc*cosp_cen*cos(self%lambda_focus - &
                                               self%lambda_centre) )
    ! Half-width in longitudinal direction (lambda)
    half_width_agnesi(1) = self%half_width
    ! Half-width in latitudinal direction (phi)
    half_width_agnesi(2) = sqrt( abs(half_width_agnesi(1)**2 - &
                                     half_width_focus**2) )
    ! Non-scaled function arguments in lambda and phi directions
    chisurf_arg(1) = scaled_radius*acos( sinp_cen**2 + (cosp_cen**2)* &
                                         cos(chi_1 - self%lambda_centre) )
    chisurf_arg(2) = scaled_radius*acos( &
                         sin(chi_2)*sinp_cen + &
                         cos(chi_2)*cosp_cen )
    ! Scaled function arguments in lambda and phi directions
    chisurf_arg = chisurf_arg/half_width_agnesi

    return
  end subroutine agnesi_coordinate_spherical

  !=============================================================================
  !> @brief Writes out parameters of Witch-of-Agnesi mountain function in
  !>        (lon,lat) coordinates.
  !>
  !> @param[in] self An object of type agnesi_spherical_type
  !=============================================================================
  subroutine write_agnesi_spherical_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(agnesi_spherical_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'agnesi_sphere_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "Witch-of-Agnesi mountain parameters in (lon,lat) coordinates: "
    write(funit, fmtreal) "mountain_height = ", self%mountain_height
    write(funit, fmtreal) "half_width      = ", self%half_width
    write(funit, fmtreal) "lambda_centre   = ", self%lambda_centre
    write(funit, fmtreal) "phi_centre      = ", self%phi_centre
    write(funit, fmtreal) "lambda_focus    = ", self%lambda_focus
    write(funit, fmtreal) "phi_focus       = ", self%phi_focus
    close(funit)

    return
  end subroutine write_agnesi_spherical_type

end module agnesi_orography_spherical_mod

