!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to solve the PMSL equation

module pmsl_solve_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  STENCIL, CROSS2D
  use fs_continuity_mod,    only: WTHETA, W2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_solve_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                        & ! dx_at_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! f_func
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), & ! exner_in
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /) ! exner_out
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_solve_code
  end type pmsl_solve_kernel_type

  public :: pmsl_solve_code

contains

  !> @brief Solve equation for PMSL
  !> @details Basic formation is based on the old UM diagnostic which is
  !>          described in UM Documentation Paper 80.
  !>          https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_080.pdf
  !>          This kernel iteratively solves the PMSL equation to
  !>          obtain the smoothed value
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dx_at_w2      cell sizes at w2 dofs
  !> @param[in]     f_func        total forcing function to relax against PMSL
  !> @param[in]     exner_in      Initial guess for exner at MSL
  !> @param[in]     smap_2d_size  Size of the stencil map in each direction
  !> @param[in]     sm_len        Max size of the stencil map in any direction
  !> @param[in]     smap_2d       Stencil map
  !> @param[in]     height_wth    Height of wth levels above mean sea level
  !> @param[in,out] exner_out     Next guess for exner at MSL
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for w2 fields
  !> @param[in]     undf_w2       Number of total degrees of freedom for w2 fields
  !> @param[in]     map_w2        Dofmap for the cell at the base of the column for w2 fields
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth      Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for wtheta
  subroutine pmsl_solve_code(nlayers,                         &
                             dx_at_w2,                        &
                             f_func,                          &
                             exner_in,                        &
                             smap_2d_size, sm_len, smap_2d,   &
                             height_wth,                      &
                             exner_out,                       &
                             ndf_w2, undf_w2, map_w2,         &
                             ndf_2d, undf_2d, map_2d,         &
                             ndf_wth, undf_wth, map_wth)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    integer(kind=i_def), intent(in) :: sm_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_2d_size
    integer(kind=i_def), dimension(ndf_2d,sm_len,4), intent(in) :: smap_2d

    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2

    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_w2) :: dx_at_w2
    real(kind=r_def),    intent(in), dimension(undf_2d) :: f_func
    real(kind=r_def),    intent(in), dimension(undf_2d) :: exner_in
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: exner_out

    ! Internal variables
    real(kind=r_def) :: dx2, dy2, pre_factor
    real(kind=r_def), parameter :: pmsl_smooth_height=500.0_r_def
    integer(kind=i_def) :: xp1, xm1, yp1, ym1

    ! Calculate which cell in the x branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    xp1 = smap_2d_size(3)
    xm1 = smap_2d_size(1)
    ! Calculate which cell in the y branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    yp1 = smap_2d_size(4)
    ym1 = smap_2d_size(2)

    ! Only calculated above a certain height or if stencil point exists
    if (height_wth(map_wth(1)) > pmsl_smooth_height .and. &
         xp1 == 2 .and. xm1 == 2 .and. yp1 ==2  .and. ym1 == 2) then

      ! Calculate inverse dx values - distance between centres of adjacent cells
      dx2 = ((dx_at_w2(map_w2(1))+dx_at_w2(map_w2(3)))/2.0_r_def)**2_i_def
      dy2 = ((dx_at_w2(map_w2(2))+dx_at_w2(map_w2(4)))/2.0_r_def)**2_i_def

      ! Factor obtained from re-arrangement of Poisson equation to put
      ! centre exner term on the LHS
      pre_factor = (dx2*dy2) / (2.0_r_def*dx2+2.0_r_def*dy2)

      ! Calculate average exner value
      exner_out(map_2d(1)) = pre_factor * (                 &
                               (1.0_r_def/dx2) *            &
                              (exner_in(smap_2d(1,xm1,1)) + &
                               exner_in(smap_2d(1,xp1,3)))  &
                             + (1.0_r_def/dy2) *            &
                              (exner_in(smap_2d(1,ym1,2)) + &
                               exner_in(smap_2d(1,yp1,4)))  &
                             - f_func(map_2d(1)) )

    else

      exner_out(map_2d(1)) = exner_in(smap_2d(1,1,1))

    end if

  end subroutine pmsl_solve_code

end module pmsl_solve_kernel_mod
