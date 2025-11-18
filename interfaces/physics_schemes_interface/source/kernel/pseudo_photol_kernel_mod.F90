!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets up a prescribed photolysis rate profile based on location and local time.

module pseudo_photol_kernel_mod

  use argument_mod,      only: arg_type, ANY_DISCONTINUOUS_SPACE_1, &
                               GH_FIELD, GH_REAL, GH_SCALAR,        &
                               GH_INTEGER, GH_READ, GH_READWRITE,   &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def, radians_to_degrees
  use fs_continuity_mod, only: WTHETA
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: pseudo_photol_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
         arg_type( GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),              & ! photol_rate_single
         arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),                           & ! current_time_hour
         arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),                           & ! current_time_minute
         arg_type( GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1 ), & ! latitude
         arg_type( GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1 )  & ! longitude
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pseudo_photol_code
  end type pseudo_photol_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: pseudo_photol_code

contains

!> @brief Uses a standard photolysis value and creates a 1-D profile by applying
!!        diurnal, vertical and geographical scaling.
!!        Uses current time and longitude to determine rough 'local time' with
!!        a view to applying the peak as observed in actual calculations.
!!        Values increase vertically with height. Uses COS(latitude) to ensure
!!        peak at equator (ignore N/S hemisphere seasons for now).
!>
!> @param[in]     nlayers             Number of layers
!> @param[in,out] photol_rate_single  Photolysis rates profile
!> @param[in]     current_hour        Current model hour
!> @param[in]     current_minutes     Current model minutes
!> @param[in]     latitude array      Latitude (radians)
!> @param[in]     longitude array     Longitude (radians)
!> @param[in]     ndf_wtheta          Number of degrees of freedom per cell for Wtheta space
!> @param[in]     undf_wtheta         Number of unique degrees of freedom for Wtheta space
!> @param[in]     map_wtheta          Dofmap for the cell at the base of the column for Wtheta space
!> @param[in]     ndf_2d              Number of degrees of freedom per cell for 2D fields
!> @param[in]     undf_2d             Number of unique degrees of freedom per cell for 2D fields
!> @param[in]     map_2d              Dofmap for cell for 2D fields

subroutine pseudo_photol_code(nlayers,                             &
                              photol_rate_single,                  &
                              current_hour,                        &
                              current_minute,                      &
                              latitude,                            &
                              longitude,                           &
                              ndf_wtheta, undf_wtheta, map_wtheta, &
                              ndf_2d, undf_2d, map_2d              &
                             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers   
  integer(kind=i_def), intent(in) :: ndf_wtheta
  integer(kind=i_def), intent(in) :: undf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d

  real(kind=r_def), intent(in out), dimension(undf_wtheta) :: photol_rate_single
  integer(kind=i_def), intent(in) :: current_hour
  integer(kind=i_def), intent(in) :: current_minute
  real(kind=r_def), intent(in), dimension(undf_2d) :: latitude
  real(kind=r_def), intent(in), dimension(undf_2d) :: longitude

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: tloc  ! local time : hour, rounded
  real(kind=r_def) :: lon_deg  ! Longitude - degrees
  real(kind=r_def) :: cos_lat  ! COS of Latitude
  real(kind=r_def) :: tgmt     ! Model time GMT : hour
  real(kind=r_def) :: tdiff    ! Local-GMT diff based on longitude
  real(kind=r_def) :: jrate

  ! Standard Rate (s^-1) i.e reference/ base rate used to generate the profile
  real(kind=r_def), parameter :: std_rate = 1.3_r_def

  ! Diurnal scaling factors for local time (from sample model values)
  real(kind=r_def), dimension(24), parameter :: diurn_scaling = [              &
  !   1:00        2:00        3:00        4:00        5:00        6:00
    0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def,    &
  !   7:00        8:00        9:00        10:00        11:00       12:00
    0.41_r_def, 0.78_r_def, 0.91_r_def, 1.00_r_def,  0.90_r_def, 0.64_r_def,   &
  !   13:00       14:00       15:00       16:00       17:00       18:00
    0.63_r_def, 0.62_r_def, 0.50_r_def, 0.23_r_def, 0.04_r_def, 0.00_r_def,    &
  !   19:00       20:00       21:00       22:00       23:00       24:00
    0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def, 0.00_r_def  ]

  ! Set up a single column with diurnal, vertical and geograpical scaling.
  ! All operations are done with LFRic KIND variables.

  ! Determine local time

  lon_deg = radians_to_degrees * longitude(map_2d(1))
  if ( lon_deg < 0.0_r_def )  lon_deg = lon_deg + 360.0_r_def

  tgmt = real(current_hour + (current_minute / 60.0_r_def), r_def)
  tdiff = 24.0_r_def * lon_deg / 360.0_r_def
  tloc = int( tgmt + tdiff, i_def )

  if (tloc > 24_i_def) tloc = tloc - 24_i_def
  if (tloc < 0_i_def)  tloc = tloc + 24_i_def
  tloc = max(tloc, 1_i_def)
  tloc = min(tloc, 24_i_def)

  ! Non-vertically-varying component - diurnal and geographical scaling
  cos_lat = cos(latitude(map_2d(1)))
  cos_lat = min(cos_lat, 1.0_r_def)
  cos_lat = max(cos_lat, 0.0_r_def)
  jrate = (std_rate/nlayers) * diurn_scaling(tloc) * cos_lat

  ! Apply vertical scaling
  do k = 0, nlayers
    photol_rate_single(map_wtheta(1)+k) = jrate * k
  end do

end subroutine pseudo_photol_code

end module pseudo_photol_kernel_mod
