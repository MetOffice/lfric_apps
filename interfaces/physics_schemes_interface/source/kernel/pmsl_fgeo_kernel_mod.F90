!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to calculate geostrophic forcing of smoothed PMSL

module pmsl_fgeo_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  STENCIL, CROSS2D
  use fs_continuity_mod,    only: WTHETA, W2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use planet_constants_mod, only: kappa, p_zero, grcp

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_fgeo_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA, STENCIL(CROSS2D)),  & ! exner
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! theta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA, STENCIL(CROSS2D)),  & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! pmsl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! t_msl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                        & ! dx_at_w2
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! f_vg
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! f_ug
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! exner_m
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_fgeo_code
  end type pmsl_fgeo_kernel_type

  public :: pmsl_fgeo_code

contains

  !> @brief Calculate geostrophic forcing functions
  !> @details Basic formation is based on the old UM diagnostic which is
  !>          described in UM Documentation Paper 80.
  !>          https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_080.pdf
  !>          This kernel calculates the geostrophic forcing functions
  !>          u & v that will be used in the calculation of the
  !>          smoothed pmsl
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     exner         exner pressure in theta space
  !> @param[in]     smap_wt_size  Size of the stencil map in each direction
  !> @param[in]     wt_len        Max size of the stencil map in any direction
  !> @param[in]     smap_wt       Stencil map
  !> @param[in]     theta         potential temperature
  !> @param[in]     height_wth    Height of wth levels above mean sea level
  !> @param[in]     smap_wth_size Size of the stencil map in each direction
  !> @param[in]     wth_len       Max size of the stencil map in any direction
  !> @param[in]     smap_wth      Stencil map
  !> @param[in]     pmsl          unsmoothed pressure at mean sea level
  !> @param[in]     t_at_mean_sea_level temperature at mean sea level
  !> @param[in]     dx_at_w2      cell sizes at w2 dofs
  !> @param[in,out] f_vg          'meridional' geostrophic forcing function
  !> @param[in,out] f_ug          'zonal' geostrophic forcing function
  !> @param[in,out] exner_m       exner pressure at mean sea level
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth      Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for wtheta
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for w2 fields
  !> @param[in]     undf_w2       Number of total degrees of freedom for w2 fields
  !> @param[in]     map_w2        Dofmap for the cell at the base of the column for w2 fields

  subroutine pmsl_fgeo_code(nlayers,                          &
                            exner,                            &
                            smap_wt_size, wt_len, smap_wt,    &
                            theta,                            &
                            height_wth,                       &
                            smap_wth_size, wth_len, smap_wth, &
                            pmsl,                             &
                            t_at_mean_sea_level,              &
                            dx_at_w2,                         &
                            f_vg,                             &
                            f_ug,                             &
                            exner_m,                          &
                            ndf_wth, undf_wth, map_wth,       &
                            ndf_2d, undf_2d, map_2d,          &
                            ndf_w2, undf_w2, map_w2)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d
    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2

    integer(kind=i_def), intent(in) :: wt_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_wt_size
    integer(kind=i_def), dimension(ndf_wth,wt_len,4), intent(in) :: smap_wt

    integer(kind=i_def), intent(in) :: wth_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_wth_size
    integer(kind=i_def), dimension(ndf_wth,wth_len,4), intent(in) :: smap_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_wth) :: exner
    real(kind=r_def),    intent(in), dimension(undf_wth) :: theta
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(in), dimension(undf_2d) :: pmsl
    real(kind=r_def),    intent(in), dimension(undf_2d) :: t_at_mean_sea_level
    real(kind=r_def),    intent(in), dimension(undf_w2) :: dx_at_w2
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: f_vg
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: f_ug
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: exner_m

    ! Internal variables
    real(kind=r_def) :: theta_m, idx, idy, fac
    integer(kind=i_def) :: xp1, xm1, yp1, ym1

    ! Calculate exner and theta at mean sea level
    exner_m(map_2d(1)) = (pmsl(map_2d(1)) / p_zero)**kappa
    theta_m = t_at_mean_sea_level(map_2d(1)) / exner_m(map_2d(1))

    ! Calculate inverse dx values - distance between centres of adjacent cells
    idx = 1.0_r_def / (dx_at_w2(map_w2(1))+dx_at_w2(map_w2(3)))
    idy = 1.0_r_def / (dx_at_w2(map_w2(2))+dx_at_w2(map_w2(4)))

    ! Calculate which cell in the x branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    xp1 = smap_wt_size(3)
    xm1 = smap_wt_size(1)
    ! If one of the stencil points doesn't exist, we are doing a difference
    ! over half the dx distance, hence multiply the result by 2 to get an
    ! estimate of the total distance
    if (xp1 == 1 .or. xm1 == 1) then
      fac = 2.0_r_def
    else
      fac = 1.0_r_def
    end if

    ! Calculate geostrophic wind function f*vg
    f_vg(map_2d(1)) = fac * idx * ( grcp * (height_wth(smap_wth(1,xp1,3))      &
                                          - height_wth(smap_wth(1,xm1,1)) )    &
                                  + theta(map_wth(1)) *                        &
                       (exner(smap_wt(1,xp1,3)) - exner(smap_wt(1,xm1,1))) )   &
                       / theta_m

    ! Calculate which cell in the y branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    yp1 = smap_wt_size(4)
    ym1 = smap_wt_size(2)
    ! If one of the stencil points doesn't exist, we are doing a difference
    ! over half the dx distance, hence multiply the result by 2 to get an
    ! estimate of the total distance
    if (yp1 == 1 .or. ym1 == 1) then
      fac = 2.0_r_def
    else
      fac = 1.0_r_def
    end if

    ! Calculate geostrophic wind function f*ug
    ! For a true geostrophic wind, there should be a minus sign here, but
    ! that will cancel later on so we omit it
    f_ug(map_2d(1)) = fac * idy * ( grcp * (height_wth(smap_wth(1,yp1,4))      &
                                          - height_wth(smap_wth(1,ym1,2)) )    &
                                  + theta(map_wth(1)) *                        &
                       (exner(smap_wt(1,yp1,4)) - exner(smap_wt(1,ym1,2))) )   &
                       / theta_m

  end subroutine pmsl_fgeo_code

end module pmsl_fgeo_kernel_mod
