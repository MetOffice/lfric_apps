!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to calculate PMSL

module pmsl_kernel_mod

  use argument_mod,         only: arg_type,            &
                                  GH_FIELD, GH_SCALAR, &
                                  GH_READ, GH_WRITE, GH_INTEGER, &
                                  GH_REAL, CELL_COLUMN, ANY_DISCONTINUOUS_SPACE_1
  use fs_continuity_mod,    only: WTHETA, W3
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use planet_constants_mod, only: g, r, recip_kappa, lapse, p_zero, power

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                        & ! exner_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                    & ! exner_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                    & ! theta_wth
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                        & ! height_w3
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                    & ! height_wth
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                          & ! levelupper
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! pmsl
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! t_at_mean_sea_level
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_code
  end type pmsl_kernel_type

  public :: pmsl_code

contains

  !> @brief Calculate pressure at mean sea level
  !> @details Basic formation is based on the old UM diagnostic which is
  !>          described in UM Documentation Paper 80.
  !>          https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_080.pdf
  !>          This kernel calculates the unsmoothed version,
  !>          i.e. with just a basic extrapolation below orography
  !> @param[in]     nlayers     The number of layers
  !> @param[in]     exner_w3    Exner pressure in w3 space
  !> @param[in]     exner_wth   Exner pressure in theta space
  !> @param[in]     theta_wth   Potential temperature
  !> @param[in]     height_w3   Height of w3 levels above mean sea level
  !> @param[in]     height_wth  Height of wth levels above mean sea level
  !> @param[in]     levelupper  Level above boundary layer to use for PMSL calculation
  !> @param[in,out] pmsl        Pressure at mean sea level
  !> @param[in,out] t_at_mean_sea_level Temperature at mean sea level
  !> @param[in]     ndf_w3      Number of degrees of freedom per cell for wrho
  !> @param[in]     undf_w3     Number of total degrees of freedom for wrho
  !> @param[in]     map_w3      Dofmap for the cell at the base of the column for wrho
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth    Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column for wtheta
  !> @param[in]     ndf_2d      Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d     Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d      Dofmap for the cell at the base of the column for 2d fields

  subroutine pmsl_code(nlayers,            &
                       exner_w3,           &
                       exner_wth,          &
                       theta_wth,          &
                       height_w3,          &
                       height_wth,         &
                       levelupper,         &
                       pmsl,               &
                       t_at_mean_sea_level,&
                       ndf_w3,             &
                       undf_w3,            &
                       map_w3,             &
                       ndf_wth,            &
                       undf_wth,           &
                       map_wth,            &
                       ndf_2d,             &
                       undf_2d,            &
                       map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d
    integer(kind=i_def), intent(in) :: undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_w3)  :: exner_w3
    real(kind=r_def),    intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def),    intent(in), dimension(undf_wth) :: theta_wth
    real(kind=r_def),    intent(in), dimension(undf_w3)  :: height_w3
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height_wth
    integer(kind=i_def), intent(in) :: levelupper
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: pmsl
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: t_at_mean_sea_level

    ! Internal variables
    real(kind=r_def) :: pressure
    real(kind=r_def) :: t_ref_level_1

    ! Calculate pressure from exner_w3
    ! Used same code structure as in spectral_gwd_kernel_mod

    pressure = p_zero*(exner_w3(map_w3(1)))**recip_kappa

    ! Calculate T at reference level and mean sea level

    t_ref_level_1 = theta_wth(map_wth(1)+levelupper) * exner_wth(map_wth(1)+levelupper) &
                  + lapse * ( height_wth(map_wth(1)+levelupper)   &
                            - height_w3(map_w3(1)) )

    t_at_mean_sea_level(map_2d(1)) = t_ref_level_1                           &
                                   + lapse * (height_w3(map_w3(1)) )

    ! Calculate PMSL from the previously calculated variables above

    pmsl(map_2d(1)) = pressure * (t_at_mean_sea_level(map_2d(1)) / t_ref_level_1)**power

  end subroutine pmsl_code

end module pmsl_kernel_mod
