!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines and implements bell-shaped hill orography profile in
!>        Cartesian coordinates.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of bell-shaped mountain function from
!>          Cartesian coordinates: x and y.
!>          Reference: Lock et al. (2012), Section 3b.
!>          Bell-shaped hill parameters in Cartesian coordinates are:
!>          mountain_height - Height of bell-shaped hill function (m),
!>          half_width_x - Half-width of bell-shaped hill function in
!>                         x direction (m),
!>          half_width_y - Half-width of bell-shaped hill function in
!>                         y direction (m),
!>          x_centre -  x coordinate centre of bell-shaped hill
!>                     function (m),
!>          y_centre -  y coordinate centre of bell-shaped hill
!>                    function (m),
!>          direction - Direction of bell-shaped hill function
!>                           (x, y, or xy).
!-------------------------------------------------------------------------------
module bell_orography_cartesian_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type
  use orography_bell_cartesian_config_mod, only : direction_x, &
                                                  direction_y, &
                                                  direction_xy
  use log_mod,                only : log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR

  implicit none

  private

  !> @brief Holds parameters and methods used to calculate bell-shaped
  !>        orography profile in Cartesian coordinates.
  type, public, extends(analytic_orography_type) :: bell_cartesian_type

    private
    ! Bell-shaped hill function parameters in Cartesian coordinates
    real(kind=r_def)    :: mountain_height
    real(kind=r_def)    :: half_width_x
    real(kind=r_def)    :: half_width_y
    real(kind=r_def)    :: x_centre
    real(kind=r_def)    :: y_centre
    integer(i_def)      :: direction

  contains

    procedure, public, pass(self) :: analytic_orography => bell_orography_cartesian
    procedure                     :: bell_coordinate_cartesian
    procedure                     :: write_bell_cartesian_type

  end type bell_cartesian_type

  ! Constructor for bell_cartesian_type
  interface bell_cartesian_type
    module procedure bell_cartesian_constructor
  end interface

contains

  !=============================================================================
  !> @brief Constructor for bell-shaped hill function in Cartesian
  !>        coordinates.
  !>
  !> @param[in] mountain_height Height of mountain function (m)
  !> @param[in] half_width_x    Half-width of mountain function in x direction (m)
  !> @param[in] half_width_y    Half-width of mountain function in y direction (m)
  !> @param[in] x_centre        x coordinate centre of mountain function (m)
  !> @param[in] y_centre        y coordinate centre of mountain function (m)
  !> @param[in] direction       Direction of mountain function read (x, y or both)
  !> @return    self            An object of type bell_cartesian_type
  !=============================================================================
  type(bell_cartesian_type) function bell_cartesian_constructor(     &
                                                        mountain_height, &
                                                        half_width_x,    &
                                                        half_width_y,    &
                                                        x_centre,        &
                                                        y_centre,        &
                                                        direction )      &
                                                        result(self)

    implicit none

    ! Arguments
    real(kind=r_def),    intent(in) :: mountain_height, &
                                       half_width_x,    &
                                       half_width_y,    &
                                       x_centre,        &
                                       y_centre
    integer(i_def),      intent(in) :: direction

    ! Assign values
    self%mountain_height = mountain_height
    self%half_width_x    = half_width_x
    self%half_width_y    = half_width_y
    self%x_centre        = x_centre
    self%y_centre        = y_centre
    self%direction       = direction

    return
  end function bell_cartesian_constructor

  !=============================================================================
  !> @brief Calculates bell-shaped hill function in Cartesian coordinates.
  !>
  !> @param[in] self      An object of type bell_cartesian_type
  !> @param[in] chi_1     x coordinate (m)
  !> @param[in] chi_2     y coordinate (m)
  !> @return    chi_surf  Surface height (m)
  !=============================================================================
  function bell_orography_cartesian(self, chi_1, chi_2) result(chi_surf)

    implicit none

    ! Arguments
    class(bell_cartesian_type),   intent(in) :: self
    real(kind=r_def),             intent(in) :: chi_1, chi_2
    real(kind=r_def)                         :: chi_surf
    ! Internal variables
    real(kind=r_def) :: chisurf_arg(2)

    ! Calculate transformed/scaled function arguments
    call bell_coordinate_cartesian(self, chi_1, chi_2, chisurf_arg)

    ! Calculate bell-shaped hill surface height
    ! Reference: Lock et al. (2012), Section 3b, Eq. 14
    chi_surf = self%mountain_height/(1.0_r_def + sum(chisurf_arg**2))**1.5_r_def

    return
  end function bell_orography_cartesian

  !=============================================================================
  !> @brief Transforms/scales coordinate for Cartesian bell-shaped hill
  !>       function.
  !>
  !> @param[in]  self         An object of type bell_cartesian_type
  !> @param[in]  chi_1        x coordinate (m)
  !> @param[in]  chi_2        y coordinate (m)
  !> @param[out] chisurf_arg  Bell-shaped hill function
  !>                          transformed/scaled arguments
  !=============================================================================
  subroutine bell_coordinate_cartesian(self, chi_1, chi_2, chisurf_arg)

    use orography_helper_functions_mod, only : coord_transform_cart_biperiodic

    implicit none

    ! Arguments
    class(bell_cartesian_type),   intent(in)  :: self
    real(kind=r_def),             intent(in)  :: chi_1, chi_2
    real(kind=r_def),             intent(out) :: chisurf_arg(2)
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
    ! chisurf_arg(1) is function argument in x direction,
    ! chisurf_arg(2) is function argument in y direction.
    select case(self%direction)
      case (direction_x)
        chisurf_arg(1) = chi_1_per/self%half_width_x
      case (direction_y)
        chisurf_arg(2) = chi_2_per/self%half_width_y
      case (direction_xy)
        chisurf_arg(1) = chi_1_per/self%half_width_x
        chisurf_arg(2) = chi_2_per/self%half_width_y
      case default
        write(log_scratch_space,'(A)') "bell_coordinate_cartesian: "// &
              "No valid orography directions (x, y, or xy) selected. "
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select

    return
  end subroutine bell_coordinate_cartesian

  !=============================================================================
  !> @brief Writes out parameters of bell-shaped hill function in
  !>        Cartesian coordinates.
  !>
  !> @param[in] self An object of type bell_cartesian_type
  !=============================================================================
  subroutine write_bell_cartesian_type(self)

    use constants_mod, only : str_short, str_max_filename

    implicit none

    ! Arguments
    class(bell_cartesian_type), intent(in) :: self
    ! Temporary write variables
    integer(kind=i_def),             parameter :: funit = 777
    character(len=str_max_filename), parameter :: fname = 'bell_cart_params.txt'
    character(len=str_short),        parameter :: fmtreal = '(A,ES15.3E3)', &
                                                  fmtint  = '(A,I15)'

    ! Temporary write
    open(funit, file = trim(fname), status = 'replace')
    write(funit,'(A)') &
          "Bell-shaped hill parameters in Cartesian coordinates: "
    write(funit, fmtreal) "mountain_height = ", self%mountain_height
    write(funit, fmtreal) "half_width_x    = ", self%half_width_x
    write(funit, fmtreal) "half_width_y    = ", self%half_width_y
    write(funit, fmtreal) "x_centre        = ", self%x_centre
    write(funit, fmtreal) "y_centre        = ", self%y_centre
    write(funit, fmtint)  "direction       = ", self%direction
    close(funit)

    return
  end subroutine write_bell_cartesian_type

end module bell_orography_cartesian_mod
