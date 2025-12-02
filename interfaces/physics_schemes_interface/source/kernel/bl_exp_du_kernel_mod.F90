!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate explicit estimate of turbulent momentum diffusion

module bl_exp_du_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type, func_type,                    &
                                      GH_FIELD, GH_READ, CELL_COLUMN,         &
                                      ANY_SPACE_1, ANY_SPACE_2,               &
                                      ANY_SPACE_3, ANY_SPACE_4,               &
                                      ANY_SPACE_5, ANY_SPACE_6,               &
                                      GH_REAL, GH_WRITE, GH_SCALAR, GH_INTEGER
  use constants_mod,            only: r_def, i_def, r_bl
  use fs_continuity_mod,        only: W1, W2
  use kernel_mod,               only: kernel_type
  use nlsizes_namelist_mod,     only: bl_levels
  use jules_surface_config_mod, only: formdrag, formdrag_dist_drag

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: bl_exp_du_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                     &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           &! nfaces
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! tau w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_2), &! tau_land 2d w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_2), &! tau_ssi 2d w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_3), &! rhokm wth/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_4), &! rdz wth/fd1
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_1), &! u_physics w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_5), &! surf_interp mult w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_3), &! ngstress wth/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_1), &! fd_tau w3/w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_6)  &! sea_current mix w3/w2
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_exp_du_code
  end type bl_exp_du_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public bl_exp_du_code

contains

  !> @brief Code based on the same underlying science as the UM BL scheme
  !>        but designed for cell faces
  !> @param[in]     nlayers        The number of layers in a column
  !> @param[in]     nfaces         The number of faces to calculate
  !> @param[in,out] tau            Turbulent stress
  !> @param[in,out] tau_land       Wind stress over land
  !> @param[in,out] tau_ssi        Wind stress over sea and sea-ice
  !> @param[in]     rhokm          Momentum eddy diffusivity
  !> @param[in]     rdz            1/dz
  !> @param[in]     u_physics      Wind in native space at time n
  !> @param[in]     surf_interp    Surface variables which need interpolating
  !> @param[in]     ngstress       NG stress function
  !> @param[in]     fd_tau         Stress from turbulent form-drag
  !> @param[in]     sea_current    Ocean surface current
  !> @param[in]     ndf_half       Number of DOFs per cell for half levels
  !> @param[in]     undf_half      Number of unique DOFs for half levels
  !> @param[in]     map_half       Dofmap for the cell at the base of the column
  !> @param[in]     ndf_2d         Number of DOFs per cell for the 2D space
  !> @param[in]     undf_2d        Number of unique DOFs for 2D space
  !> @param[in]     map_2d         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_full       Number of DOFs per cell for full levels
  !> @param[in]     undf_full      Number of unique DOFs for full levels
  !> @param[in]     map_full       Dofmap for the cell at the base of the column
  !> @param[in]     ndf_rdz        Number of DOFs per cell for rdz space
  !> @param[in]     undf_rdz       Number of unique DOFs for rdz space
  !> @param[in]     map_rdz        Dofmap for the cell at the base of the column
  !> @param[in]     ndf_surf       Number of DOFs per cell for surface space
  !> @param[in]     undf_surf      Number of unique DOFs for surface space
  !> @param[in]     map_surf       Dofmap for the cell at the base of the column
  !> @param[in]     ndf_curr       Number of DOFs per cell for current space
  !> @param[in]     undf_curr      Number of unique DOFs for current space
  !> @param[in]     map_curr       Dofmap for the cell at the base of the column
  subroutine bl_exp_du_code(nlayers,       &
                            nfaces,        &
                            tau,           &
                            tau_land,      &
                            tau_ssi,       &
                            rhokm,         &
                            rdz,           &
                            u_physics,     &
                            surf_interp,   &
                            ngstress,      &
                            fd_tau,        &
                            sea_current,   &
                            ndf_half,      &
                            undf_half,     &
                            map_half,      &
                            ndf_2d,        &
                            undf_2d,       &
                            map_2d,        &
                            ndf_full,      &
                            undf_full,     &
                            map_full,      &
                            ndf_rdz,       &
                            undf_rdz,      &
                            map_rdz,       &
                            ndf_surf,      &
                            undf_surf,     &
                            map_surf,      &
                            ndf_curr,      &
                            undf_curr,     &
                            map_curr)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ex_flux_uv_mod, only: ex_flux_uv
    use atm_fields_bounds_mod, only: pdims

    implicit none

    ! Arguments
    integer, intent(in) :: nlayers, nfaces

    integer(kind=i_def), intent(in) :: ndf_half, undf_half
    integer(kind=i_def), intent(in) :: map_half(ndf_half)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: ndf_curr, undf_curr
    integer(kind=i_def), intent(in) :: map_curr(ndf_curr)
    integer(kind=i_def), intent(in) :: ndf_rdz, undf_rdz
    integer(kind=i_def), intent(in) :: map_rdz(ndf_rdz)
    integer(kind=i_def), intent(in) :: ndf_full, undf_full
    integer(kind=i_def), intent(in) :: map_full(ndf_full)
    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf)

    real(kind=r_def), dimension(undf_half), intent(inout) :: tau
    real(kind=r_def), dimension(undf_2d), intent(inout) :: tau_land, tau_ssi

    real(kind=r_def), dimension(undf_half), intent(in) :: u_physics, &
         fd_tau
    real(kind=r_def), dimension(undf_full), intent(in) :: rhokm, ngstress
    real(kind=r_def), dimension(undf_rdz),  intent(in) :: rdz
    real(kind=r_def), dimension(undf_surf), intent(in) :: surf_interp
    real(kind=r_def), dimension(undf_curr), intent(in) :: sea_current

    ! Internal variables
    integer(kind=i_def) :: df, k
    real(kind=r_bl) :: rhokm_land, rhokm_ssi, fland, flandfac, fseafac
    real(kind=r_bl) :: zh_nonloc(1)
    real(kind=r_bl), dimension(0:bl_levels-1) :: tau_grad, tau_non_grad, u_sp, &
         rdz_sp, rhokm_sp, ngstress_sp, tau_sp, fd_tau_sp

    ! loop over all faces of the cell
    do df = 1,nfaces

      ! Only calculate face if it's not already been done
      if (tau(map_half(df)) == 0.0_r_def) then

        !================================================================
        ! In the UM this happens in bdy_expl3
        !================================================================
        fland      = surf_interp(map_surf(df) + 0)
        rhokm_land = surf_interp(map_surf(df) + 1)
        rhokm_ssi  = surf_interp(map_surf(df) + 2)
        flandfac   = surf_interp(map_surf(df) + 3)
        fseafac    = surf_interp(map_surf(df) + 4)
        zh_nonloc  = surf_interp(map_surf(df) + 5)

        do k = 0, bl_levels-1
          u_sp(k) = u_physics(map_half(df)+k)
          rdz_sp(k) = rdz(map_rdz(df)+k)
          rhokm_sp(k) = rhokm(map_full(df)+k)
          ngstress_sp(k) = ngstress(map_full(df)+k)
          fd_tau_sp(k) = fd_tau(map_half(df)+k)
          tau_sp(k) = tau(map_half(df)+k)
        end do

        tau_land(map_2d(df)) = rhokm_land * u_sp(0) * flandfac

        ! If sea surface current (sea_current) has not been obtained
        ! from coupling fields or from an ancillary file then it will simply
        ! contain uniform zeros.
        tau_ssi(map_2d(df)) = rhokm_ssi *                                      &
                 (u_sp(0) - sea_current(map_curr(df))) * fseafac
        tau_sp(0) = fland * tau_land(map_2d(df))                               &
                        + (1.0_r_bl - fland) * tau_ssi(map_2d(df))

        if (formdrag == formdrag_dist_drag) then
          if (fland > 0.0_r_bl) then
            tau_land(map_2d(df)) = tau_land(map_2d(df))                        &
                                    + fd_tau(map_half(df)) / fland
          end if
        end if

        call ex_flux_uv(pdims, pdims, pdims, bl_levels,                        &
                        u_sp(0:bl_levels-1),                                   &
                        zh_nonloc,                                             &
                        rdz_sp(1:bl_levels-1),                                 &
                        rhokm_sp(0:bl_levels-1),                               &
                        ngstress_sp(1:bl_levels-1),                            &
                        fd_tau_sp(0:bl_levels-1),                              &
                        tau_sp(0:bl_levels-1),                                 &
                        tau_grad(0:bl_levels-1), tau_non_grad(0:bl_levels-1))

        do k = 0, bl_levels-1
          tau(map_half(df)+k) = tau_sp(k)
        end do

      end if ! this face needs calculating

    end do ! loop over df

  end subroutine bl_exp_du_code

end module bl_exp_du_kernel_mod
