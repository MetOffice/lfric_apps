!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-20.
!-----------------------------------------------------------------------------
!> @brief Clear-air turbulence potential including mountain waves.

module mountain_wave_pot_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_WRITE,         &
                           GH_REAL, GH_INTEGER,       &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2, &
                           STENCIL, REGION
  use constants_mod, only: r_def, i_def, rmdi, pi
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: mountain_wave_pot_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                   STENCIL(REGION)),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                   STENCIL(REGION)),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                   STENCIL(REGION)),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                   STENCIL(REGION)),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                   STENCIL(REGION)),           &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: mountain_wave_pot_code
  end type mountain_wave_pot_kernel_type

  public :: mountain_wave_pot_code

  !> Number of stress thresholds used to spread the stress downstream
  integer(kind=i_def), parameter :: nthresh = 3_i_def
  !> Stress thresholds, in increasing order (Pa)
  real(kind=r_def), parameter :: stress_thresh(nthresh) =                    &
                                 (/ 0.007_r_def, 0.065_r_def, 0.25_r_def /)
  !> Stress below which the turbulence probability is zero (Pa)
  real(kind=r_def), parameter :: stress_min = 0.0625_r_def
  !> Slope of the turbulence probability with stress (% Pa-1)
  real(kind=r_def), parameter :: prob_slope = 42.6666667_r_def
  !> Offset of the turbulence probability (%)
  real(kind=r_def), parameter :: prob_offset = 1.833333_r_def
  !> Largest turbulence probability (%)
  real(kind=r_def), parameter :: prob_max = 100.0_r_def
  !> Scaling of the mountain-wave probability, which is otherwise too
  !> strong compared with the clear-air turbulence probability
  real(kind=r_def), parameter :: prob_scaling = 0.1_r_def
  !> Half-width of a direction sector with no neighbour on one side
  real(kind=r_def), parameter :: default_half_width = 0.25_r_def * pi

contains

  !> @brief   Clear-air turbulence potential including mountain waves.
  !> @details The mountain-wave predictor follows Turner (1999), Met Office
  !>          Forecasting Research Technical Report 265. Starting from the
  !>          orographic stress in the cell, for each stress threshold in
  !>          turn, a cell below the threshold takes the largest stress
  !>          above the threshold of any neighbour whose wind blows towards
  !>          it. A neighbour's wind blows towards this cell if its
  !>          direction lies within the sector of directions closer to this
  !>          cell than to the cells either side, with the sector edges at
  !>          the midpoints between neighbouring cells. The stress is then
  !>          mapped to a turbulence probability. The output is the larger
  !>          of this and the clear-air turbulence probability, or the
  !>          clear-air turbulence probability alone where the stress is
  !>          missing. Missing data where the clear-air turbulence
  !>          probability is missing.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] turb_pot       Turbulence potential on pressure
  !>                               levels (%)
  !> @param[in]     stress_plev    Magnitude of the orographic stress on
  !>                               pressure levels (Pa), or missing data
  !> @param[in]     smap_s_size    Number of cells in the stencil
  !> @param[in]     smap_s         Stencil dofmap for the stress
  !> @param[in]     u_plev         Eastward wind on pressure levels (m s-1)
  !> @param[in]     smap_u_size    Number of cells in the stencil
  !> @param[in]     smap_u         Stencil dofmap for the eastward wind
  !> @param[in]     v_plev         Northward wind on pressure levels (m s-1)
  !> @param[in]     smap_v_size    Number of cells in the stencil
  !> @param[in]     smap_v         Stencil dofmap for the northward wind
  !> @param[in]     cat_prob       Clear-air turbulence probability on
  !>                               pressure levels (%)
  !> @param[in]     latitude       Latitude of the cell centre (radians)
  !> @param[in]     smap_ll_size   Number of cells in the stencil
  !> @param[in]     smap_ll        Stencil dofmap for latitude
  !> @param[in]     longitude      Longitude of the cell centre (radians)
  !> @param[in]     smap_lo_size   Number of cells in the stencil
  !> @param[in]     smap_lo        Stencil dofmap for longitude
  !> @param[in]     nplev          Number of pressure levels
  !> @param[in]     ndf_pl         Number of degrees of freedom per cell for
  !>                               the pressure-level fields
  !> @param[in]     undf_pl        Number of unique degrees of freedom for
  !>                               the pressure-level fields
  !> @param[in]     map_pl         Dofmap for the pressure-level fields
  !> @param[in]     ndf_2d         Number of degrees of freedom per cell for
  !>                               latitude and longitude
  !> @param[in]     undf_2d        Number of unique degrees of freedom for
  !>                               latitude and longitude
  !> @param[in]     map_2d         Dofmap for latitude and longitude
  subroutine mountain_wave_pot_code(nlayers, turb_pot,                   &
                                    stress_plev, smap_s_size, smap_s,    &
                                    u_plev, smap_u_size, smap_u,         &
                                    v_plev, smap_v_size, smap_v,         &
                                    cat_prob,                            &
                                    latitude, smap_ll_size, smap_ll,     &
                                    longitude, smap_lo_size, smap_lo,    &
                                    nplev,                               &
                                    ndf_pl, undf_pl, map_pl,             &
                                    ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d
    integer(kind=i_def), intent(in) :: smap_s_size, smap_u_size
    integer(kind=i_def), intent(in) :: smap_v_size
    integer(kind=i_def), intent(in) :: smap_ll_size, smap_lo_size
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_pl, smap_s_size) :: smap_s
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_pl, smap_u_size) :: smap_u
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_pl, smap_v_size) :: smap_v
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_2d, smap_ll_size) :: smap_ll
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_2d, smap_lo_size) :: smap_lo

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(inout), dimension(undf_pl) :: turb_pot
    real(kind=r_def), intent(in),    dimension(undf_pl) :: stress_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: u_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: v_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: cat_prob
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude
    real(kind=r_def), intent(in),    dimension(undf_2d) :: longitude
    integer(kind=i_def), intent(in) :: nplev

    ! Internal variables
    integer(kind=i_def) :: kp, m, j, ithresh, dof
    real(kind=r_def)    :: lat0, lat_n, dlon, rel, closest_plus, closest_minus
    real(kind=r_def)    :: stress, stress_nbr, wind_dir, offset, prob, cat
    ! Position of each stencil cell in the local tangent plane (radians)
    real(kind=r_def)    :: x_nbr(smap_s_size), y_nbr(smap_s_size)
    ! Direction of each stencil cell from the centre (radians)
    real(kind=r_def)    :: dir_nbr(smap_s_size)
    ! Angular half-widths of the sector of directions pointing from each
    ! neighbour towards this cell, on either side (radians)
    real(kind=r_def)    :: width_plus(smap_s_size), width_minus(smap_s_size)
    integer(kind=i_def) :: j_plus, j_minus

    ! Directions of the neighbours, cell 1 being this cell
    lat0 = latitude(smap_ll(1, 1))
    do m = 2, smap_s_size
      lat_n = latitude(smap_ll(1, m))
      dlon  = longitude(smap_lo(1, m)) - longitude(smap_lo(1, 1))
      x_nbr(m) = cos(lat_n) * sin(dlon)
      y_nbr(m) = cos(lat0) * sin(lat_n) - sin(lat0) * cos(lat_n) * cos(dlon)
      dir_nbr(m) = atan2(y_nbr(m), x_nbr(m))
    end do

    ! Sector edges. The direction from neighbour m back to this cell is
    ! opposite to dir_nbr(m); the neighbours of m are assumed to be laid
    ! out like those of this cell, so the sector half-widths are those of
    ! dir_nbr(m) among this cell's neighbours. Each edge is the direction
    ! of the midpoint between m and the next neighbour round.
    do m = 2, smap_s_size
      closest_plus  = pi
      closest_minus = -pi
      j_plus  = 0_i_def
      j_minus = 0_i_def
      do j = 2, smap_s_size
        if ( j == m ) cycle
        rel = wrap_angle(dir_nbr(j) - dir_nbr(m))
        if ( rel > 0.0_r_def .and. rel < closest_plus ) then
          closest_plus = rel
          j_plus = j
        else if ( rel <= 0.0_r_def .and. rel > closest_minus ) then
          closest_minus = rel
          j_minus = j
        end if
      end do
      width_plus(m)  = edge_width(m, j_plus)
      width_minus(m) = edge_width(m, j_minus)
    end do

    do kp = 1, nplev

      dof = map_pl(1) + kp - 1
      cat = cat_prob(dof)
      stress = stress_plev(smap_s(1, 1) + kp - 1)

      if ( cat == rmdi .or. stress == rmdi ) then
        ! No clear-air turbulence probability, or no mountain-wave
        ! predictor on this level
        turb_pot(dof) = cat
        cycle
      end if

      ! Spread the stress downstream, one threshold at a time
      do ithresh = 1, nthresh
        if ( stress >= stress_thresh(ithresh) ) cycle
        do m = 2, smap_s_size
          stress_nbr = stress_plev(smap_s(1, m) + kp - 1)
          if ( stress_nbr <= stress_thresh(ithresh) ) cycle
          if ( u_plev(smap_u(1, m) + kp - 1) == 0.0_r_def .and.              &
               v_plev(smap_v(1, m) + kp - 1) == 0.0_r_def ) cycle
          wind_dir = atan2(v_plev(smap_v(1, m) + kp - 1),                  &
                           u_plev(smap_u(1, m) + kp - 1))
          ! Wind direction relative to the direction from m to this cell
          offset = wrap_angle(wind_dir - dir_nbr(m) - pi)
          if ( offset < width_plus(m) .and. offset > -width_minus(m) ) then
            stress = max(stress, stress_nbr)
          end if
        end do
      end do

      if ( stress < stress_min ) then
        prob = 0.0_r_def
      else
        prob = min(prob_slope * stress + prob_offset, prob_max)
      end if

      turb_pot(dof) = max(prob_scaling * prob, cat)

    end do

  contains

    !> @brief  Half-width of the sector between two neighbours.
    !> @param[in] m Index of the neighbour
    !> @param[in] j Index of the next neighbour round, or 0 if none
    !> @return    Angle from neighbour m to the sector edge (radians)
    function edge_width(m, j) result(width)

      implicit none

      integer(kind=i_def), intent(in) :: m, j
      real(kind=r_def) :: width

      if ( j == 0_i_def ) then
        width = default_half_width
      else
        width = abs( wrap_angle( atan2(y_nbr(m) + y_nbr(j),                &
                                       x_nbr(m) + x_nbr(j)) - dir_nbr(m) ) )
      end if

    end function edge_width

  end subroutine mountain_wave_pot_code

  !> @brief  Wrap an angle into the range (-pi, pi].
  !> @param[in] angle Angle (radians)
  !> @return    Equivalent angle in (-pi, pi]
  function wrap_angle(angle) result(wrapped)

    implicit none

    real(kind=r_def), intent(in) :: angle
    real(kind=r_def) :: wrapped

    wrapped = angle - 2.0_r_def * pi * anint( angle / ( 2.0_r_def * pi ) )

  end function wrap_angle

end module mountain_wave_pot_kernel_mod
