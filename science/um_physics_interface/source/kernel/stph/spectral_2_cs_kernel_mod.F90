!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Update spectral coefficients and perform spectra to cubesphere transformation
!> @details This kernel updates the spectral coefficient with the random component for the
!!          stph forcing pattern, then applies the phase shifting dependent
!!          on the stph vertical level, and finally it peforms the spectral
!!          to cubedsphere transformation for each level where stph is active.
module spectral_2_cs_kernel_mod
  ! TO DO after PSyclone ticket 1312
  ! at https://github.com/stfc/PSyclone/issues/1312
  ! Once GH_ARRAY and NRANKS
  ! uncomment lines below and removed the next "use argument_mod" call.
  ! use argument_mod,      only: arg_type, GH_FIELD,        &
  !                              GH_SCALAR, GH_ARRAY,       &
  !                              GH_WRITE, GH_READ,         &
  !                              GH_INTEGER, GH_REAL,       &
  !                              ANY_DISCONTINUOUS_SPACE_1, &
  !                              ANY_DISCONTINUOUS_SPACE_2, &
  !                              NRANKS, CELL_COLUMN

  use argument_mod,      only: arg_type, GH_FIELD,        &
                               GH_SCALAR,                 &
                               GH_WRITE, GH_READ,         &
                               GH_INTEGER, GH_REAL,       &
                               ANY_DISCONTINUOUS_SPACE_1, &
                               ANY_DISCONTINUOUS_SPACE_2, &
                               ANY_DISCONTINUOUS_SPACE_3, &
                               CELL_COLUMN

  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  ! Spectral coefficients (to be updated everytimestep)

  ! TO DO after PSyclone ticket 1312
  ! at https://github.com/stfc/PSyclone/issues/1312
  ! Once GH_ARRAY and NRANKS.
  ! Uncomment lines below and removed the next "type ... :: spectral_2_cs_kernel_type" call.

  ! !> Metadata describing the kernel to PSyclone
  ! !>
  ! type, public, extends(kernel_type) :: spectral_2_cs_kernel_type
  !   private
  !   !type(arg_type) :: meta_args(9) = (/                                   &
  !        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! fp
  !        arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),  & ! longitude
  !        arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3),  & ! Pnm_star
  !        arg_type(GH_ARRAY, GH_REAL, GH_READ, NRANKS*2),                   & ! stph_spectral_coeffc
  !        arg_type(GH_ARRAY, GH_REAL, GH_READ, NRANKS*2),                   & ! stph_spectral_coeffs
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_level_bottom
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_level_top
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_n_min
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! stph_n_max
  !        /)
  !        integer :: operates_on = CELL_COLUMN

  ! contains
  !   procedure, nopass ::  spectral_2_cs_code
  ! end type spectral_2_cs_kernel_type

  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spectral_2_cs_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! fp
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),  & ! longitude
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3),  & ! Pnm_star
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_level_top
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! stph_n_min
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! stph_n_max
         /)
         integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass ::  spectral_2_cs_code
  end type spectral_2_cs_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: spectral_2_cs_code
contains

  !> @brief Update spectral coefficients and perform spectral to cubedsphere transformation
  !> @param[in]      nlayers                The number of layers
  !> @param[in,out]  fp                     Forcing Pattern
  !> @param[in]      longitude              2D array with longitudes
  !> @param[in]      Pnm_star               Spherical harmonic coefficients
  !> @param[in]      nranks_array           No. Ranks (shape) of stph_spectral_coeff arrays
  !> @param[in]      dims_array             Dimension of stphh_spectral_coeff arrays
  !> @param[in]      coeffc_phase           Array with real (cosine) spectral coefficients
  !> @param[in]      coeffs_phase           Array with imaginary (sine) spectral coefficients
  !> @param[in]      stph_level_bottom      Bottom level where stph is applied
  !> @param[in]      stph_level_top         Top level where stph is applied
  !> @param[in]      stph_n_min             stph minimum wavenumber
  !> @param[in]      stph_n_max             stph maximum wavenumber
  !> @param[in]      ndf_hgt                Number of degrees of freedom per cell for wtheta
  !> @param[in]      undf_hgt               Number of total degrees of freedom for wtheta
  !> @param[in]      map_hgt                Dofmap for the cell at the base of the column for wthera
  !> @param[in]      ndf_2d                 Number of degrees of freedom per cell for 2d space
  !> @param[in]      undf_2d                Number of unique degrees of freedom for  2d space
  !> @param[in]      map_2d                 Dofmap for the cell at the base of the column for  2d space
  !> @param[in]      ndf_sp                 Number of degrees of freedom per cell for spectral space
  !> @param[in]      undf_sp                Number of unique degrees of freedom for spectral space
  !> @param[in]      map_sp                 Dofmap for the cell at the base of the column for spectral space

  subroutine  spectral_2_cs_code(nlayers,              &
                                 fp,                   &
                                 longitude,            &
                                 Pnm_star,             &
                                 nranks_array,         &
                                 dims_array,           &
                                 coeffc_phase,         &
                                 coeffs_phase,         &
                                 stph_level_bottom,    &
                                 stph_level_top,       &
                                 stph_n_min,           &
                                 stph_n_max,           &
                                 ndf_hgt,              &
                                 undf_hgt,             &
                                 map_hgt,              &
                                 ndf_2d,               &
                                 undf_2d,              &
                                 map_2d,               &
                                 ndf_sp,               &
                                 undf_sp,              &
                                 map_sp                &
                                 )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_hgt, ndf_2d, ndf_sp
    integer(kind=i_def), intent(in) :: undf_hgt, undf_2d, undf_sp
    integer(kind=i_def), intent(in) :: nranks_array
    integer(kind=i_def), intent(in), dimension(ndf_hgt)  :: map_hgt
    integer(kind=i_def), intent(in), dimension(ndf_sp)  ::  map_sp
    integer(kind=i_def), intent(in), dimension(ndf_2d)  ::  map_2d
    integer(kind=i_def), intent(in), dimension(nranks_array) :: dims_array
    !fields
    real(kind=r_def), intent(inout), dimension(undf_hgt) :: fp

    real(kind=r_def), intent(in), dimension(undf_2d)  :: longitude
    real(kind=r_def), intent(in), dimension(undf_sp)  :: Pnm_star
    real(kind=r_def), intent(in), dimension(dims_array(1), dims_array(2)) :: &
         coeffc_phase, coeffs_phase

    ! stph scalars
    integer(kind=i_def), intent(in) :: stph_level_bottom
    integer(kind=i_def), intent(in) :: stph_level_top
    integer(kind=i_def), intent(in) :: stph_n_min
    integer(kind=i_def), intent(in) :: stph_n_max

    ! Spectral coefficients for vertical phasing (local to the timestep)
    real(kind=r_def) :: cos_fac
    real(kind=r_def) :: sin_fac

    ! Integers for iteration
    integer(kind=i_def) :: k, m, n, n_row

    ! Initialize n_row to stph_n_min-1
    n_row = 0
    do n = 1,stph_n_min-1
      n_row = n_row + n
    end do
    ! Apply vertical scaling for each spectral coeff.
    do n = stph_n_min, stph_n_max
      n_row = n_row + n
      do m = 0, n

        cos_fac = Pnm_star(map_sp(1) + n_row+m)*cos(m*longitude(map_2d(1)))
        sin_fac = Pnm_star(map_sp(1) + n_row+m)*sin(m*longitude(map_2d(1)))

        ! Compute the inverse transformation for each stph level
        do k = stph_level_bottom, stph_level_top

          fp(map_hgt(1) + k)  = fp(map_hgt(1) + k) +  &
                  coeffc_phase(k,n_row+m)*cos_fac +   &
                  coeffs_phase(k,n_row+m)*sin_fac
        end do
      end do
    end do

  end subroutine  spectral_2_cs_code

end module spectral_2_cs_kernel_mod
