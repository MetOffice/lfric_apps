!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates Schar mountain orography profile in Cartesian coordinates.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of Schar mountain function from Cartesian
!>          coordinates: x and y.
!>          Reference: Melvin et al. (2010), Section 4.4.
!>          Schar mountain parameters in Cartesian coordinates are:
!>          mountain_height - Height of Schar mountain function (m),
!>          half_width_x - Half-width of Schar mountain function in
!>                         x direction (m),
!>          half_width_y - Half-width of Schar mountain function in
!>                         y direction (m),
!>          wavelength - Wavelength of cosine part of Schar mountain
!>                       function (m),
!>          x_centre - x coordinate centre of Schar mountain function (m),
!>          y_centre - y coordinate centre of Schar mountain function (m),
!>          direction_cart - Direction of Schar mountain function
!>                           (x, y or both).
!-------------------------------------------------------------------------------
module schar_orography_cartesian_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type
  use orography_schar_cartesian_config_mod, only : direction_x, &
                                                   direction_y, &
                                                   direction_xy
  use log_mod,                only : log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR

  implicit none

  private
  !> @brief Holds parameters and methods used to calculate Schar orography
  !>        profile in Cartesian coordinates.
  type, public, extends(analytic_orography_type) :: schar_cartesian_type

    private
    ! Schar mountain function parameters in Cartesian coordinates
    real(kind=r_def)    :: mountain_height
    real(kind=r_def)    :: half_width_x
    real(kind=r_def)    :: half_width_y
    real(kind=r_def)    :: wavelength
    real(kind=r_def)    :: x_centre
    real(kind=r_def)    :: y_centre
    integer(kind=i_def) :: direction

  contains

    procedure, public, pass(self) :: analytic_orography => schar_orography_cartesian
    procedure                     :: schar_coordinate_cartesian
    procedure                     :: write_schar_cartesian_type

  end type schar_cartesian_type

  ! Constructor for schar_cartesian_type
  interface schar_cartesian_type
    module procedure schar_cartesian_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for Schar mountain function in Cartesian coordinates.
  !>
  !> @param[in] mountain_height Height of mountain function read from
  !>                            namelist (m)
  !> @param[in] half_width_x    Half-width of mountain function in x direction
  !>                            read from namelist (m)
  !> @param[in] half_width_y    Half-width of mountain function in y direction
  !>                            read from namelist (m)
  !> @param[in] wavelength      Wavelength of cosine part of  mountain function
  !>                            read from namelist (m)
  !> @param[in] x_centre        x coordinate centre of mountain function read
  !>                            from namelist (m)
  !> @param[in] y_centre        y coordinate centre of mountain function read
  !>                            from namelist (m)
  !> @param[in] direction       Direction of mountain function read from
  !>                            namelist (x, y or both)
  !> @return    self            An object of type schar_cartesian_type
  !=============================================================================
  type(schar_cartesian_type) function schar_cartesian_constructor(     &
                                                      mountain_height, &
                                                      half_width_x,    &
                                                      half_width_y,    &
                                                      wavelength,      &
                                                      x_centre,        &
                                                      y_centre,        &
                                                      direction )      &
                                                      result(self)

    implicit none

    ! Arguments
    real(kind=r_def),    intent(in) :: mountain_height, &
                                       half_width_x,    &
                                       half_width_y,    &
                                       wavelength,      &
                                       x_centre,        &
                                       y_centre
    integer(kind=i_def), intent(in) :: direction

    ! Assign values
    self%mountain_height = mountain_height
    self%half_width_x    = half_width_x
    self%half_width_y    = half_width_y
    self%wavelength      = wavelength
    self%x_centre        = x_centre
    self%y_centre        = y_centre
    self%direction       = direction

    return
  end function schar_cartesian_constructor

  !=============================================================================
  !> @brief Calculates Schar mountain function in Cartesian coordinates.
  !>
  !> @param[in] self      An object of type schar_cartesian_type
  !> @param[in] chi_1     x coordinate (m)
  !> @param[in] chi_2     y coordinate (m)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function schar_orography_cartesian(self, chi_1, chi_2) result(chi_surf)

    implicit none

    ! Arguments
    class(schar_cartesian_type), intent(in) :: self
    real(kind=r_def),            intent(in) :: chi_1, chi_2
    real(kind=r_def)                        :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(2)

    ! Calculate transformed/scaled function arguments
    call schar_coordinate_cartesian(self, chi_1, chi_2, chisurf_arg)

    ! Calculate Schar mountain surface height
    ! Reference: Melvin et al. (2010), Section 4.4., Eq. 70
    chi_surf = self%mountain_height*exp(-chisurf_arg(1)**2)*(cos(chisurf_arg(2)))**2

    return
  end function schar_orography_cartesian

  !=============================================================================
  !> @brief Transforms/scales coordinate for Cartesian Schar mountain function.
  !>
  !> @param[in]  self         An object of type schar_cartesian_type
  !> @param[in]  chi_1        x coordinate (m)
  !> @param[in]  chi_2        y coordinate (m)
  !> @param[out] chisurf_arg  Schar mountain function transformed/scaled
  !>                          arguments
  !=============================================================================
  subroutine schar_coordinate_cartesian(self, chi_1, chi_2, chisurf_arg)

    use constants_mod,                  only : PI
    use orography_helper_functions_mod, only : coord_transform_cart_biperiodic

    implicit none

    ! Arguments
    class(schar_cartesian_type), intent(in)  :: self
    real(kind=r_def),            intent(in)  :: chi_1, chi_2
    real(kind=r_def),            intent(out) :: chisurf_arg(2)
    ! Internal variables
    real(kind=r_def)               :: chi_1_per, chi_2_per

    ! Initialise transformed/scaled function arguments
    chisurf_arg = 0.0_r_def

    ! Tranform coordinates to biperiodic Cartesian domain
    call coord_transform_cart_biperiodic( chi_1,         &
                                          chi_2,         &
                                          self%x_centre, &
                                          self%y_centre, &
                                          chi_1_per,     &
                                          chi_2_per )

    ! Calculate transformed/scaled arguments depending on direction.
    ! chisurf_arg(1) is Exponential function argument,
    ! chisurf_arg(2) is Cosine function argument
    select case(self%direction)
      case (direction_x)
        ! Exponential function argument
        chisurf_arg(1) = chi_1_per/self%half_width_x
        ! Cosine function argument
        chisurf_arg(2) = PI*chi_1_per/self%wavelength
      case (direction_y)
        ! Exponential function argument
        chisurf_arg(1) = chi_2_per/self%half_width_y
        ! Cosine function argument
        chisurf_arg(2) = PI*chi_2_per/self%wavelength
      case (direction_xy)
        ! Exponential function argument
        chisurf_arg(1) = sqrt((chi_1_per/self%half_width_x)**2 + (chi_2_per/self%half_width_y)**2)
        ! Cosine function argument
        chisurf_arg(2) = PI*sqrt((chi_1_per**2+chi_2_per**2))/self%wavelength
      case default
        write(log_scratch_space,'(A)') "schar_coordinate_cartesian: "// &
              "No valid orography directions (x, y, or xy) selected. "
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select

    return
  end subroutine schar_coordinate_cartesian

  !=============================================================================
  !> @brief Writes out parameters of Schar mountain function in Cartesian
  !>        coordinates.
  !>
  !> @param[in] self An object of type schar_cartesian_type
  !=============================================================================
  subroutine write_schar_cartesian_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(schar_cartesian_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'schar_cart_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "Schar mountain parameters in Cartesian coordinates: "
    write(funit, fmtreal) "mountain_height = ", self%mountain_height
    write(funit, fmtreal) "half_width_x    = ", self%half_width_x
    write(funit, fmtreal) "half_width_y    = ", self%half_width_y
    write(funit, fmtreal) "wavelength      = ", self%wavelength
    write(funit, fmtreal) "x_centre        = ", self%x_centre
    write(funit, fmtreal) "y_centre        = ", self%y_centre
    write(funit, fmtint)  "direction       = ", self%direction
    close(funit)

    return
  end subroutine write_schar_cartesian_type

end module schar_orography_cartesian_mod
