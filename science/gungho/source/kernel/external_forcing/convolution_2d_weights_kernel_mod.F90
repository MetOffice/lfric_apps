!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes weights for convolution for spectral nudging.
!> @details This kernel computes weights for performing a convolution that is
!!          used for nudging the spectrum of a field towards the spectrum of
!!          a reference field.
!!          The filter takes the form of a sum of Legendre polynomials which
!!          help select certain wavenumbers to be nudged. The filter is
!!          enveloped by a Gaussian function to ensure that the weights go
!!          smoothly to zero.
!!          Only implemented for the lowest-order elements
module convolution_2d_weights_kernel_mod

  use argument_mod,              only: arg_type, func_type,                    &
                                       GH_FIELD, GH_SCALAR,                    &
                                       GH_REAL, GH_INTEGER,                    &
                                       GH_READ, GH_WRITE, GH_BASIS,            &
                                       CELL_COLUMN, STENCIL, REGION,           &
                                       GH_EVALUATOR,                           &
                                       ANY_DISCONTINUOUS_SPACE_3,              &
                                       ANY_DISCONTINUOUS_SPACE_7,              &
                                       ANY_DISCONTINUOUS_SPACE_9
  use constants_mod,             only: r_def, i_def, PI, EPS
  use coord_transform_mod,       only: central_angle
  use kernel_mod,                only: kernel_type
  use sci_chi_transform_mod,     only: chi2llr


  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: convolution_2d_weights_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_9), &
        arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7,  &
                                                            STENCIL(REGION)),  &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,  &
                                                            STENCIL(REGION)),  &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              &
    /)
    type(func_type) :: meta_funcs(1) = (/                                      &
        func_type(ANY_DISCONTINUOUS_SPACE_7, GH_BASIS)                         &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: convolution_2d_weights_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: convolution_2d_weights_code

contains

!> @brief Computes weights for performing convolution for spectral nudging.
!> @param[in]     nlayers           Number of layers
!> @param[in,out] weights           Weights to be computed: a multidata 2D field
!> @param[in]     chi1              First coordinate field
!> @param[in]     chi2              Second coordinate field
!> @param[in]     chi3              Third coordinate field
!> @param[in]     stencil_size_chi  Size of the stencil for the chi fields
!> @param[in]     stencil_map_chi   Stencil map for the coordinate fields
!> @param[in]     panel_id          Panel ID field
!> @param[in]     stencil_size_pid  Size of the stencil for the panel ID
!> @param[in]     stencil_map_pid   Stencil map for the panel ID field
!> @param[in]     kmax              Maximum wavenumber in filter
!> @param[in]     kmin              Minimum wavenumber in filter
!> @param[in]     sigma             Width of the Gaussian envelope to filter
!> @param[in]     geometry          Geometry of the mesh
!> @param[in]     topology          Topology of the mesh
!> @param[in]     coord_system      System used by coordinate fields
!> @param[in]     scaled_radius     Scaled radius of the planet
!> @param[in]     ndf_2d            Num DoFs per cell for weights field
!> @param[in]     undf_2d           Num DoFs in this partition for weights field
!> @param[in]     map_2d            DoFmap for weights field
!> @param[in]     ndf_chi           Num DoFs per cell for coordinate fields
!> @param[in]     undf_chi          Num DoFs in this partition for chi fields
!> @param[in]     map_chi           DoFmap for coordinate fields
!> @param[in]     basis_chi         Basis funcs for chi fields evaluated at W3
!> @param[in]     ndf_pid           Num DoFs per cell for panel ID
!> @param[in]     undf_pid          Num DoFs in this partition for panel ID
!> @param[in]     map_pid           DoFmap for panel ID field
subroutine convolution_2d_weights_code(nlayers,                                &
                                       weights,                                &
                                       chi1,                                   &
                                       chi2,                                   &
                                       chi3,                                   &
                                       stencil_size_chi,                       &
                                       stencil_map_chi,                        &
                                       panel_id,                               &
                                       stencil_size_pid,                       &
                                       stencil_map_pid,                        &
                                       kmax,                                   &
                                       kmin,                                   &
                                       sigma,                                  &
                                       geometry,                               &
                                       topology,                               &
                                       coord_system,                           &
                                       scaled_radius,                          &
                                       ndf_2d,                                 &
                                       undf_2d,                                &
                                       map_2d,                                 &
                                       ndf_chi,                                &
                                       undf_chi,                               &
                                       map_chi,                                &
                                       basis_chi,                              &
                                       ndf_pid,                                &
                                       undf_pid,                               &
                                       map_pid)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_2d, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in)    :: undf_2d, undf_chi, undf_pid
  integer(kind=i_def), intent(in)    :: stencil_size_chi, stencil_size_pid
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: stencil_map_chi(ndf_chi, stencil_size_chi)
  integer(kind=i_def), intent(in)    :: stencil_map_pid(ndf_pid, stencil_size_pid)
  real(kind=r_def),    intent(in)    :: basis_chi(1, ndf_chi, ndf_2d)
  integer(kind=i_def), intent(in)    :: kmax, kmin
  real(kind=r_def),    intent(in)    :: sigma
  integer(kind=i_def), intent(in)    :: geometry, topology, coord_system
  real(kind=r_def),    intent(in)    :: scaled_radius
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: chi1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi3(undf_chi)
  real(kind=r_def),    intent(inout) :: weights(undf_2d)

  ! Internal variables
  integer(kind=i_def) :: i, l, df, idx_w
  integer(kind=i_def) :: panel_c, panel_i
  real(kind=r_def)    :: chi1_c, chi2_c, chi3_c
  real(kind=r_def)    :: chi1_i, chi2_i, chi3_i
  real(kind=r_def)    :: lon_c, lat_c, lon_i, lat_i
  real(kind=r_def)    :: radius
  real(kind=r_def)    :: dist, cos_dist, sum_weights
  real(kind=r_def)    :: P_0, P_1, P_l, P_lm1, P_lm2, l_real

  ! Convolution is a 2D operation, so we can just act in the lowest layer
  ! Get coordinates for the centre of this column
  chi1_c = 0.0_r_def
  chi2_c = 0.0_r_def
  chi3_c = 0.0_r_def
  do df = 1, ndf_chi
    chi1_c = chi1_c + basis_chi(1, df, 1) * chi1(map_chi(df))
    chi2_c = chi2_c + basis_chi(1, df, 1) * chi2(map_chi(df))
    chi3_c = chi3_c + basis_chi(1, df, 1) * chi3(map_chi(df))
  end do

  ! Convert to longitude and latitude
  panel_c = INT(panel_id(map_pid(1)), i_def)
  call chi2llr(                                                                &
      chi1_c, chi2_c, chi3_c, panel_c,                                         &
      geometry, topology, coord_system, scaled_radius,                         &
      lon_c, lat_c, radius                                                     &
  )

  ! Loop through cells in the stencil
  do i = 1, stencil_size_chi
    ! Get coordinates for this cell
    chi1_i = 0.0_r_def
    chi2_i = 0.0_r_def
    chi3_i = 0.0_r_def
    do df = 1, ndf_chi
      chi1_i = chi1_i + basis_chi(1, df, 1) * chi1(stencil_map_chi(df, i))
      chi2_i = chi2_i + basis_chi(1, df, 1) * chi2(stencil_map_chi(df, i))
      chi3_i = chi3_i + basis_chi(1, df, 1) * chi3(stencil_map_chi(df, i))
    end do

    panel_i = INT(panel_id(stencil_map_pid(1, i)), i_def)
    call chi2llr(                                                              &
        chi1_i, chi2_i, chi3_i, panel_i,                                       &
        geometry, topology, coord_system, scaled_radius,                       &
        lon_i, lat_i, radius                                                   &
    )

    ! Use central angle to get distance on sphere
    if (ABS(lon_c - lon_i) < 10.0_r_def*EPSILON(1.0_r_def)                     &
        .and. ABS(lat_c - lat_i) < 10.0_r_def*EPSILON(1.0_r_def)) then
      dist = 0.0_r_def
    else
      ! Only call this for points that are not identical, which gets around
      ! potential issues in the central angle routine
      call central_angle(lon_c, lat_c, lon_i, lat_i, dist)
    end if

    ! Convert distance to cos(dist) for use in spherical harmonics
    cos_dist = COS(dist)

    ! Set weight for this cell in the stencil
    idx_w = map_2d(1) + i - 1  ! Multidata index

    ! Convolving function: enveloped Legendre polynomial -----------------------
    weights(idx_w) = 0.0_r_def

    ! First Legendre polynomials
    P_0 = 1.0_r_def
    P_1 = cos_dist

    if (kmin == 0) then
      weights(idx_w) = weights(idx_w) + 1.0_r_def/(4.0_r_def*PI) * P_0
    end if
    if (kmin <= 1 .and. kmax >= 1) then
      weights(idx_w) = weights(idx_w) + 3.0_r_def/(4.0_r_def*PI) * P_1
    end if

    ! Obtain further Legendre polynomials using recurrence relation
    P_lm2 = P_0
    P_lm1 = P_1
    if (kmax >= 2) then
      do l = 2, kmax
        l_real = REAL(l, r_def)
        ! Determine next Legendre polynomial P_l
        P_l = (                                                                &
            (2.0_r_def*l_real - 1.0_r_def) * cos_dist * P_lm1                  &
            - (l_real - 1.0_r_def) * P_lm2                                     &
        ) / l_real
        ! Update convolution weight based on Legendre polynomial:
        ! (2l+1)/(4pi)) P_l
        weights(idx_w) = weights(idx_w) + (                                    &
            (2.0_r_def*l_real + 1.0_r_def ) / (4.0_r_def*PI) * P_l             &
        )
        ! Update previous two polynomials
        P_lm2 = P_lm1
        P_lm1 = P_l
      end do
    end if

    ! Apply envelope to weights, to ensure they go smoothly to zero at the
    ! edge of the stencil, which can help avoid amplifying some wavenumbers
    weights(idx_w) = weights(idx_w) * EXP(-0.5_r_def * (dist / sigma)**2)
  end do

  ! Normalise the weights
  idx_w = map_2d(1)
  sum_weights = SUM(weights(idx_w : idx_w+stencil_size_chi-1))
  if (ABS(sum_weights) > 0.0_r_def) then
    weights(idx_w : idx_w+stencil_size_chi-1) =                                &
        weights(idx_w : idx_w+stencil_size_chi-1) / sum_weights
  end if

end subroutine convolution_2d_weights_code

end module convolution_2d_weights_kernel_mod
