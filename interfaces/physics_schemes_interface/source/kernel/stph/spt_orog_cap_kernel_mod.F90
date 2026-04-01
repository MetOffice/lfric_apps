!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Remove SPT increments over points with high stddev of ororgraphy and forcing pattern values.
module spt_orog_cap_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,        &
                               GH_WRITE, GH_REAL,         &
                               GH_SCALAR, GH_INTEGER,     &
                               GH_READ,                   &
                               ANY_DISCONTINUOUS_SPACE_1, &
                               CELL_COLUMN

  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def, l_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spt_orog_cap_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & !dX
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & !fp_spt
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),& !sd_orog
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_top
         arg_type(GH_SCALAR, GH_REAL, GH_READ),         & ! spt_orog_forcing_pattern_thresh
         arg_type(GH_SCALAR, GH_REAL, GH_READ)          & ! spt_stddev_orog_thres
         /)
         integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: spt_orog_cap_code
  end type spt_orog_cap_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public spt_orog_cap_code
contains

  !> @brief Apply the orography capping for SPT increments
  !> @details SPT increments are removed in columns where
  !>          where the standard deviation of the orography and
  !>          the forcing pattern are both too high.
  !> @param[in]     nlayers     The number of layers
  !> @param[in,out] dX          Field with SPT perturbations
  !> @param[in]     fp_spt      Forcing pattern for orography
  !> @param[in]     sd_orog     Std dev of orography
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth    Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column
  !> @param[in]     ndf_2d      Number of degrees of freedom per cell for density space
  !> @param[in]     undf_2d     Number of unique degrees of freedom for density space
  !> @param[in]     map_2d      Dofmap for the cell at the base of the column for density space
  !> @param[in]     spt_level_bottom      Bottom level of the stochastic scheme
  !> @param[in]     spt_level_top         Top level of the stochastic scheme
  !> @param[in]     spt_orog_forcing_pattern_thresh        spt_orog_forcing_pattern_thresh in stochastic_physics_config_mod
  !> @param[in]     spt_stddev_orog_thres                  spt_stddev_orog_thres in stochastic_physics_config_mod

  subroutine spt_orog_cap_code(nlayers,                         &
                               dX,                              &
                               fp_spt,                          &
                               sd_orog,                         &
                               spt_level_bottom,                &
                               spt_level_top,                   &
                               spt_orog_forcing_pattern_thresh, &
                               spt_stddev_orog_thres,           &
                               ndf_wth,                         &
                               undf_wth,                        &
                               map_wth,                         &
                               ndf_2d,                          &
                               undf_2d,                         &
                               map_2d                           &
                               )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d)  ::  map_2d
    integer(kind=i_def), intent(in) :: spt_level_bottom
    integer(kind=i_def), intent(in) :: spt_level_top
    real(kind=r_def), intent(in) :: spt_orog_forcing_pattern_thresh
    real(kind=r_def), intent(in) :: spt_stddev_orog_thres

    ! Fields perturbations + tendencies
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dX
    real(kind=r_def), intent(in), dimension(undf_wth) :: fp_spt
    real(kind=r_def), intent(in), dimension(undf_2d)  :: sd_orog

    integer(kind=i_def) :: k
    logical(kind=l_def) :: orog_mask

    !!!!  Apply orographic capping
    orog_mask = .false.
    do k = spt_level_bottom, spt_level_top
      ! Set mask if it is over land points
      ! where std dev of orog is higher than threshold
      ! and forcing pattern higher to threshold in abs value
      if ((sd_orog(map_2d(1)) > spt_stddev_orog_thres) .and.                     &
          (abs(fp_spt(map_wth(1) + k)) > spt_orog_forcing_pattern_thresh )) then
        orog_mask = .true.
        exit
      end if
    end do

    if (orog_mask) then
      do k = spt_level_bottom, spt_level_top
        dX(map_wth(1) + k) = 0.0_r_def
      end do
    end if

  end subroutine spt_orog_cap_code

end module spt_orog_cap_kernel_mod
