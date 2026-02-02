!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute SPT moisture increments that comply with vertical conservation
module spt_moisture_conservation_kernel_mod

    use argument_mod,      only: arg_type, GH_FIELD, &
                                 GH_WRITE, GH_REAL,  &
                                 GH_SCALAR, GH_INTEGER, &
                                 GH_READ, CELL_COLUMN
    use fs_continuity_mod, only: Wtheta
    use constants_mod,     only: r_def, i_def
    use kernel_mod,        only: kernel_type

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------
    !> Metadata describing the kernel to PSyclone
    !>
    type, public, extends(kernel_type) :: spt_moisture_conservation_kernel_type
      private
      type(arg_type) :: meta_args(6) = (/                 &
           arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & !dmv
           arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & !mv
           arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & !dz_wth
           arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & !rho_wth
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & !spt_level_bottom
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ)       & !spt_level_top
           /)
           integer :: operates_on = CELL_COLUMN

    contains
      procedure, nopass :: spt_moisture_conservation_code
    end type spt_moisture_conservation_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
    public spt_moisture_conservation_code
  contains

    !> @brief Apply the SPT moisture conservation in the column
    !> @details Applies the SPT moisture conservation scheme to the
    !>          water vapour in the column
    !> @param[in]     nlayers     The number of layers
    !> @param[in,out] dmv         Forcing field for WV mixing ratio
    !> @param[in]     mv          WV mixing ratio
    !> @param[in]     dz_wth      Delta z at wtheta levels
    !> @param[in]     rho_wth     Dry density interpolated to the wtheta space
    !> @param[in]     ndf_wth     Number of DOFs per cell for potential temperature space
    !> @param[in]     undf_wth    Number of unique DOFs for potential temperature space
    !> @param[in]     map_wth     dofmap for the cell at the base of the column for potential temperature space
    !> @param[in]     spt_level_bottom      Bottom level of the stochastic scheme
    !> @param[in]     spt_level_top         Top level of the stochastic scheme

    subroutine spt_moisture_conservation_code(nlayers,          &
                                              dmv,              &
                                              mv,               &
                                              rho_wth,          &
                                              dz_wth,           &
                                              spt_level_bottom, &
                                              spt_level_top,    &
                                              ndf_wth,          &
                                              undf_wth,         &
                                              map_wth           &
                                              )

      implicit none

      !Arguments
      integer(kind=i_def), intent(in) :: nlayers
      integer(kind=i_def), intent(in) :: ndf_wth
      integer(kind=i_def), intent(in) :: undf_wth
      integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
      integer(kind=i_def), intent(in) :: spt_level_bottom
      integer(kind=i_def), intent(in) :: spt_level_top

      ! Fields
      real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv
      real(kind=r_def), intent(in),    dimension(undf_wth) :: mv
      real(kind=r_def), intent(in),    dimension(undf_wth) :: dz_wth
      real(kind=r_def), intent(in),    dimension(undf_wth) :: rho_wth

      integer(kind=i_def) :: k
      real (kind=r_def) :: q_bef, q_aft, alpha

      q_bef = 0.0_r_def
      q_aft = 0.0_r_def

      do k = spt_level_bottom, spt_level_top
        ! Compute total vertical tendency before SPT forcing
        ! Q_bef = integral ( q x rho x delta_Z )
        q_bef = q_bef + rho_wth(map_wth(1)+k) * mv(map_wth(1)+k) * dz_wth(map_wth(1)+k)

        ! Compute total vertical tendency After SPT forcing
        q_aft = q_aft + rho_wth(map_wth(1)+k) *             &
                (mv(map_wth(1) + k)+ dmv(map_wth(1) + k)) * &
                dz_wth(map_wth(1)+k)
      end do

      ! Get ratio between both vertical estimates
      alpha= q_bef / q_aft

      ! Rescale perturbations
      do k = spt_level_bottom, spt_level_top
        dmv(map_wth(1) + k) = dmv(map_wth(1) + k) * alpha +                        &
                              mv(map_wth(1) + k) * ( alpha -1.0_r_def)
      end do

    end subroutine spt_moisture_conservation_code

  end module spt_moisture_conservation_kernel_mod
