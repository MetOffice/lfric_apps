!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate implicit solution for turbulent momentum diffusion

module bl_imp_du_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type, func_type,                    &
                                      GH_FIELD, GH_READ, CELL_COLUMN,         &
                                      ANY_SPACE_1, ANY_SPACE_2,               &
                                      GH_INTEGER, GH_REAL, GH_SCALAR, GH_WRITE
  use constants_mod,            only: r_def, i_def, r_bl
  use extrusion_config_mod,     only: planet_radius
  use fs_continuity_mod,        only: W1, W2, WTheta
  use kernel_mod,               only: kernel_type
  use nlsizes_namelist_mod,     only: bl_levels
  use timestepping_config_mod,  only: outer_iterations
  use blayer_config_mod,        only: fric_heating, bl_mix_w

  implicit none

  private

  real(kind=r_bl), parameter :: sqrt2 = sqrt(2.0_r_bl)

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: bl_imp_du_kernel_type
    private
    type(arg_type) :: meta_args(21) = (/                     &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           &! outer
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),          &! du_bl
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),          &! dissip
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),          &! tau
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! wind10m
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! wind10m_neut
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_1), &! tau_land
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! tau_ssi
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! pseudotau
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! rhokm
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W1),          &! rdz
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! dtrdz
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! wetrho
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! u_physics
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! u_physics_star
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2), &! surf_interp
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! du_conv
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTheta),      &! dw_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! dA
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W1),          &! height_w1
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)           &! height_w2
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_imp_du_code
  end type bl_imp_du_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public bl_imp_du_code

contains

  !> @brief Code based on the same underlying science as the UM BL scheme
  !>        but designed for cell faces
  !> @param[in]     nlayers        Number of layers
  !> @param[in]     outer          Outer loop counter
  !> @param[in,out] du_bl          BL wind increment
  !> @param[in,out] dissip         Molecular dissipation rate
  !> @param[in,out] tau            Turbulent stress
  !> @param[in,out] wind10m        10m wind
  !> @param[in,out] wind10m_neut   neutral 10m wind
  !> @param[in]     tau_land       Wind stress over land
  !> @param[in,out] tau_ssi        Wind stress over sea and sea-ice
  !> @param[in,out] pseudotau      pseudo surface wind stress
  !> @param[in]     rhokm          Momentum eddy diffusivity mapped to cell face
  !> @param[in]     rdz            1/dz mapped to cell faces
  !> @param[in]     dtrdz          dt/(r*r*dz) mapped to cell faces
  !> @param[in]     wetrho         Wet density
  !> @param[in]     u_physics      Wind in native space at time n
  !> @param[in]     u_physics_star Wind in native space after advection
  !> @param[in]     surf_interp    Surface variables which need interpolating
  !> @param[in]     du_conv        Wind increment from convection
  !> @param[in]     dw_bl          Vertical wind increment from explicit BL
  !> @param[in]     dA             Area of faces
  !> @param[in]     height_w1      Height of cell top/bottom above surface
  !> @param[in]     height_w2      Height of cell centre above surface
  !> @param[in]     ndf_w2         Number of DOFs per cell for w2 space
  !> @param[in]     undf_w2        Number of unique DOFs for w2 space
  !> @param[in]     map_w2         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w2_2d      Number of DOFs per cell for w2 surface space
  !> @param[in]     undf_w2_2d     Number of unique DOFs for w2 surface space
  !> @param[in]     map_w2_2d      Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w1         Number of DOFs per cell for w1 space
  !> @param[in]     undf_w1        Number of unique DOFs for w1 space
  !> @param[in]     map_w1         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w2_surf    Number of DOFs per cell for w2 surface space
  !> @param[in]     undf_w2_surf   Number of unique DOFs for w2 surface space
  !> @param[in]     map_w2_surf    Dofmap for the cell at the base of the column
  !> @param[in]     ndf_wth        No of DOFs per cell for wtheta space
  !> @param[in]     undf_wth       No of unique DOFs for wtheta space
  !> @param[in]     map_wth        DOFmap for cell at base of wtheta column
  subroutine bl_imp_du_code(nlayers,       &
                            outer,         &
                            du_bl,         &
                            dissip,        &
                            tau,           &
                            wind10m,       &
                            wind10m_neut,  &
                            tau_land,      &
                            tau_ssi,       &
                            pseudotau,     &
                            rhokm,         &
                            rdz,           &
                            dtrdz,         &
                            wetrho,        &
                            u_physics,     &
                            u_physics_star,&
                            surf_interp,   &
                            du_conv,       &
                            dw_bl,         &
                            dA,            &
                            height_w1,     &
                            height_w2,     &
                            ndf_w2,        &
                            undf_w2,       &
                            map_w2,        &
                            ndf_w2_2d,     &
                            undf_w2_2d,    &
                            map_w2_2d,     &
                            ndf_w1,        &
                            undf_w1,       &
                            map_w1,        &
                            ndf_w2_surf,   &
                            undf_w2_surf,  &
                            map_w2_surf,   &
                            ndf_wth,       &
                            undf_wth,      &
                            map_wth)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use bl_option_mod, only: puns, pstb

    implicit none

    ! Arguments
    integer, intent(in) :: nlayers, outer

    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in) :: ndf_w2_2d, undf_w2_2d
    integer(kind=i_def), intent(in) :: map_w2_2d(ndf_w2_2d)
    integer(kind=i_def), intent(in) :: ndf_w2_surf, undf_w2_surf
    integer(kind=i_def), intent(in) :: map_w2_surf(ndf_w2_surf)
    integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
    integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)

    real(kind=r_def), dimension(undf_w2),  intent(inout) :: du_bl, dissip, tau
    real(kind=r_def), dimension(undf_w2_2d), intent(inout) :: wind10m,         &
         wind10m_neut, tau_ssi, pseudotau

    real(kind=r_def), dimension(undf_w2),  intent(in) :: rhokm, dtrdz,         &
         u_physics, u_physics_star, du_conv, dA, height_w2, wetrho
    real(kind=r_def), dimension(undf_w2_2d), intent(in) :: tau_land
    real(kind=r_def), dimension(undf_w2_surf), intent(in) :: surf_interp
    real(kind=r_def), dimension(undf_w1),   intent(in) :: height_w1, rdz
    real(kind=r_def), dimension(undf_wth),  intent(in) :: dw_bl

    ! Internal variables
    integer(kind=i_def) :: k, df, k_blend
    real(kind=r_bl) :: pnonl, i1, e1, e2, gamma1, gamma2, cdr10m, cdr10m_neut,&
         du_1, cq_cm_1, tau_land_star, tau_ssi_star, fb_surf, fland,           &
         tau_land_loc, tau_ssi_loc, cd10m_neut
    real(kind=r_bl), dimension(0:bl_levels-1) :: tau_star, cq_cm, du_star,     &
         du_nt, tau_loc, dtr_rhodz, rhokm_sp, du_bl_sp, rdz_sp

    !================================================================
    ! In the UM this happens in imp_solver - predictor section
    !================================================================
    ! loop over all faces of the cell
    do df = 1,4

      ! Only calculate face if it's not already been done
      if (du_bl(map_w2(df)) == 0.0_r_def .and. &
          du_bl(map_w2(df)+1) == 0.0_r_def) then

        fland   = surf_interp(map_w2_surf(df) + 0)
        fb_surf = surf_interp(map_w2_surf(df) + 6)
        k_blend = int(surf_interp(map_w2_surf(df) + 7), i_def)
        cdr10m  = surf_interp(map_w2_surf(df) + 8)
        cd10m_neut  = surf_interp(map_w2_surf(df) + 9)
        cdr10m_neut = surf_interp(map_w2_surf(df) + 10)

        if (fb_surf > 0.0_r_bl) then
          pnonl = puns
        else
          pnonl = pstb
        end if

        ! Set implicit weights for predictor section
        i1 = (1.0_r_bl + 1.0_r_bl/sqrt2)*(1.0_r_bl + pnonl)
        e1 = (1.0_r_bl + 1.0_r_bl/sqrt2)*(pnonl + 1.0_r_bl/sqrt2 +          &
                                  sqrt(pnonl*(sqrt2 - 1.0_r_bl) + 0.5_r_bl) )
        e2 = (1.0_r_bl + 1.0_r_bl/sqrt2)*(pnonl + 1.0_r_bl/sqrt2 -          &
                                  sqrt(pnonl*(sqrt2 - 1.0_r_bl) + 0.5_r_bl) )
        gamma1 = i1
        gamma2 = i1 - e1

        ! Take local copies so we don't over-write input values
        ! on each outer iteration
        tau_loc(0:bl_levels-1) = tau(map_w2(df):map_w2(df)+bl_levels-1)
        tau_ssi_loc = tau_ssi(map_w2_2d(df))
        tau_land_loc = tau_land(map_w2_2d(df))

        !================================================================
        ! In the UM this happens in bdy_impl3 - predictor section
        !================================================================

        do k = 0, bl_levels-1
          du_nt(k) = u_physics_star(map_w2(df) + k) -                          &
                     u_physics(map_w2(df) + k) +                               &
                     du_conv(map_w2(df) + k)
          dtr_rhodz(k) = dtrdz(map_w2(df) + k) / wetrho(map_w2(df) + k)
          rhokm_sp(k) = rhokm(map_w2(df) + k)
          du_bl_sp(k) = du_bl(map_w2(df) + k)
          rdz_sp(k) = rdz(map_w1(df) + k)
        end do

        call matrix_sweep(du_bl_sp, du_1, cq_cm, cq_cm_1, df, gamma1, gamma2,  &
                          height_w1, dtr_rhodz, tau_loc, du_nt,rhokm_sp,rdz_sp,&
                          k_blend, ndf_w1, undf_w1, map_w1)

        !================================================================
        ! In Jules this happens in im_sf_pt2 - predictor section
        !================================================================
        if (fland > 0.0_r_bl) then
          tau_land_star = ( gamma2 * tau_land_loc +                      &
                            gamma1 * rhokm_sp(0) * du_1 ) /              &
                          (1.0_r_bl + gamma1 * rhokm_sp(0) * cq_cm_1)
        else
          tau_land_star = 0.0_r_bl
        end if
        if (fland < 1.0_r_bl) then
          tau_ssi_star = ( gamma2 * tau_ssi_loc +                        &
                           gamma1 * rhokm_sp(0) * du_1 ) /               &
                         (1.0_r_bl + gamma1 * rhokm_sp(0) * cq_cm_1 )
        else
          tau_ssi_star = 0.0_r_bl
        end if
        tau_star(0) = fland * tau_land_star + (1.0_r_bl - fland) * tau_ssi_star

        !================================================================
        ! In the UM this happens in bdy_impl4 - predictor section
        !================================================================
        du_bl_sp(0) = du_bl_sp(0) - cq_cm(0) * tau_star(0)

        do k = 1, bl_levels-1
          du_bl_sp(k) = du_bl_sp(k) -                      &
                                  cq_cm(k) * du_bl_sp(k-1)
          tau_star(k) = gamma2 * tau_loc(k) +              &
                        gamma1 * rhokm_sp(k) * rdz_sp(k) * &
                        ( du_bl_sp(k) - du_bl_sp(k-1) )
        end do

        do k = 0, bl_levels-1
          du_star(k) = du_bl_sp(k)
        end do

        !================================================================
        ! In the UM this happens in bdy_impl3 - corrector section
        !================================================================

        gamma2 = i1 - e2

        ! Update explicit fluxes using predictor X* value as needed by the
        ! 2nd stage of the scheme.
        do k = 1, bl_levels-1
          tau_loc(k) = tau_loc(k) + rhokm_sp(k) *          &
                           ( du_bl_sp(k) - du_bl_sp(k-1) ) &
                           * rdz_sp(k)
        end do

        call matrix_sweep(du_bl_sp, du_1, cq_cm, cq_cm_1, df, gamma1, gamma2,  &
                          height_w1, dtr_rhodz, tau_loc, du_nt,rhokm_sp,rdz_sp,&
                          k_blend, ndf_w1, undf_w1, map_w1)

        !================================================================
        ! In Jules this happens in im_sf_pt2 - corrector section
        !================================================================
        if (fland > 0.0_r_bl) then
          tau_land_loc = ( gamma2 * ( tau_land_loc +                     &
                                      rhokm_sp(0) * du_star(0) )         &
                                    + gamma1 * rhokm_sp(0) * du_1 ) /    &
                         (1.0_r_bl + gamma1 * rhokm_sp(0) * cq_cm_1 )
        else
          tau_land_loc = 0.0_r_bl
        end if
        if (fland < 1.0_r_bl) then
          tau_ssi_loc = ( gamma2 * ( tau_ssi_loc +                       &
                                     rhokm_sp(0) * du_star(0) )          &
                                   + gamma1 * rhokm_sp(0) * du_1 ) /     &
                        (1.0_r_bl + gamma1 * rhokm_sp(0) * cq_cm_1 )
        else
          tau_ssi_loc = 0.0_r_bl
        end if
        tau_loc(0) = fland * tau_land_loc + (1.0_r_bl - fland) * tau_ssi_loc

        !================================================================
        ! In the UM this happens in bdy_impl4 - corrector section
        !================================================================
        du_bl_sp(0) = du_bl_sp(0) - cq_cm(0) * tau_loc(0)
        tau_loc(0) = tau_loc(0) + tau_star(0)

        do k = 1, bl_levels-1
          du_bl_sp(k) = du_bl_sp(k) -                      &
                                  cq_cm(k) * du_bl_sp(k-1)
        end do

        do k = 0, bl_levels-1
          du_bl_sp(k) = du_bl_sp(k) + du_star(k)
        end do

        ! Calculate the molecular dissipation rate for frictional heating

        ! first calculate tau
        do k = 1, bl_levels-1
          tau_loc(k) = tau_star(k) + gamma2 * tau_loc(k) + &
                       gamma1 * rhokm_sp(k) * rdz_sp(k) *  &
                       ( du_bl_sp(k) - du_bl_sp(k-1) )
        end do

        !================================================================
        ! In the UM this happens in imp_solver - end
        !================================================================
        if (fric_heating) then
          do k = 0, bl_levels-2
            dissip(map_w2(df) + k) = tau_loc(k+1) * rdz_sp(k+1) *    &
                     ( u_physics(map_w2(df) + k+1) + du_bl_sp(k+1) - &
                         u_physics(map_w2(df) + k) - du_bl_sp(k) )
          end do

          ! Add on dissipation in the surface layer
          dissip(map_w2(df)) = ( tau_loc(0) *                                  &
                             (u_physics(map_w2(df)) + du_bl_sp(0) ) +          &
                             dissip(map_w2(df)) *                              &
                          (height_w2(map_w2(df)+1) - height_w2(map_w2(df))) )/ &
                          (height_w2(map_w2(df)+1) - height_w1(map_w1(df)))

          ! And nothing on the top level
          dissip(map_w2(df)+bl_levels-1) = 0.0_r_bl
        end if

        ! Diagnostic calculations at final iteration only
        if (outer == outer_iterations) then
          tau_ssi(map_w2_2d(df))  = (tau_ssi_loc + tau_ssi_star)               &
                                    * dA(map_w2(df))
          tau_land_loc            = (tau_land_loc + tau_land_star)             &
                                    * dA(map_w2(df))
          wind10m(map_w2_2d(df))  = (u_physics(map_w2(df)) + du_bl_sp(0))      &
                                    * cdr10m * dA(map_w2(df))
          wind10m_neut(map_w2_2d(df)) =                                        &
                                   (u_physics(map_w2(df)) + du_bl_sp(0))       &
                                   * cdr10m_neut * dA(map_w2(df))
          ! tau in bottom level is grid-bix mean surface stress for diagnostics
          tau(map_w2(df)) = fland * tau_land_loc +                             &
                            (1.0_r_bl - fland) * tau_ssi(map_w2_2d(df))
          ! save end-of-timestep stress profile for diagnostics
          do k = 1, bl_levels-1
            tau(map_w2(df) + k) = tau_loc(k)*dA(map_w2(df))
          end do
          pseudotau(map_w2_2d(df)) = tau(map_w2(df)) / cd10m_neut
        end if

        ! final increment calculation
        do k = 0, bl_levels-1
          du_bl(map_w2(df) + k) = du_bl_sp(k) -                      &
               ( u_physics_star(map_w2(df) + k) - u_physics(map_w2(df) + k) )
        end do

        ! Above BL levels it's just the convection increment
        do k = bl_levels, nlayers-1
          du_bl(map_w2(df) + k) = du_conv(map_w2(df) + k)
        end do

      end if ! this face needs calculating

    end do ! loop over df

    if (bl_mix_w) then
      ! Copy dw_bl increment into du_bl
      do k = 1, bl_levels
        du_bl(map_w2(5)+k) = dw_bl(map_wth(1)+k)
      end do
    end if

  end subroutine bl_imp_du_code

  subroutine matrix_sweep(du_bl_sp, du_1, cq_cm, cq_cm_1, df, gamma1, gamma2,  &
                          height_w1, dtrdz, tau, du_nt,rhokm_sp,rdz_sp,k_blend,&
                          ndf_w1, undf_w1, map_w1)

    implicit none

    !Passed in
    integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
    integer(kind=i_def), intent(in) :: map_w1(ndf_w1)

    real(kind=r_def), dimension(undf_w1), intent(in) :: height_w1
    real(kind=r_bl), dimension(0:bl_levels-1), intent(in) :: du_nt, tau, dtrdz,&
         rhokm_sp, rdz_sp

    integer(kind=i_def), intent(in) :: df, k_blend
    real(kind=r_bl), intent(in) :: gamma1, gamma2

    ! Passed out
    real(kind=r_bl), dimension(0:bl_levels-1),  intent(inout) :: du_bl_sp
    real(kind=r_bl), dimension(0:bl_levels-1), intent(out) :: cq_cm
    real(kind=r_bl), intent(out) :: du_1, cq_cm_1

    ! Internal variables
    integer(kind=i_def) :: k
    real(kind=r_bl) :: r_sq, rr_sq, rbm, am, cu_prod

    ! Topmost level
    k = bl_levels-1
    r_sq = (height_w1(map_w1(df) + k) + planet_radius)**2
    du_bl_sp(k) = -dtrdz(k) * r_sq * tau(k)
    ! addition of non-turbulent increments
    du_bl_sp(k) = gamma2 * ( du_bl_sp(k) + du_nt(k) )

    cq_cm(k) = -dtrdz(k) * gamma1 * r_sq * rhokm_sp(k)               &
               * rdz_sp(k)

    rbm = 1.0_r_bl / (1.0_r_bl - cq_cm(k) )
    du_bl_sp(k) = rbm * du_bl_sp(k)
    cq_cm(k) = rbm * cq_cm(k)

    ! Downward sweep
    do k = bl_levels-2,1,-1
      r_sq = (height_w1(map_w1(df) + k) + planet_radius)**2
      rr_sq = (height_w1(map_w1(df) + k+1) + planet_radius)**2

      du_bl_sp(k) = dtrdz(k) * ( rr_sq * tau(k+1) - r_sq * tau(k) )
      ! addition of non-turbulent increments
      du_bl_sp(k) = gamma2 * ( du_bl_sp(k) + du_nt(k) )

      am = -dtrdz(k) * gamma1 * rr_sq * rhokm_sp(k+1)                &
           * rdz_sp(k+1)

      cq_cm(k) = -dtrdz(k) * gamma1 * r_sq * rhokm_sp(k)             &
                 * rdz_sp(k)

      rbm = 1.0_r_bl / (1.0_r_bl - cq_cm(k) - am * (1.0_r_bl + cq_cm(k+1) ) )

      du_bl_sp(k) = du_bl_sp(k) -                          &
                              am * du_bl_sp(k+1)

      du_bl_sp(k) = rbm * du_bl_sp(k)
      cq_cm(k) = rbm * cq_cm(k)
    end do

    ! Surface level
    k = 0
    r_sq = (height_w1(map_w1(df) + k) + planet_radius)**2
    rr_sq = (height_w1(map_w1(df) + k+1) + planet_radius)**2

    du_bl_sp(k) = dtrdz(k) * rr_sq * tau(k+1)
    ! addition of non-turbulent increments
    du_bl_sp(k) = gamma2 * ( du_bl_sp(k) + du_nt(k) )

    am = -dtrdz(k) * gamma1 * rr_sq * rhokm_sp(k+1)                  &
         * rdz_sp(k+1)

    rbm = 1.0_r_bl / (1.0_r_bl - am * (1.0_r_bl + cq_cm(k+1) ) )

    du_bl_sp(k) = du_bl_sp(k) - am * du_bl_sp(k+1)

    du_bl_sp(k) = rbm * du_bl_sp(k)
    cq_cm(k) = rbm * r_sq * dtrdz(k)

    ! Calculate the coefficients for calculation of implicit u momentum
    ! fluxes, i.e. Upward sweep of matrix (to k_blend)
    ! NOTE: Gives coeffs for coupling to the bottom model level
    ! when k_blend = 0 (i.e. blend_height_opt == blend_level1)
    k = k_blend
    du_1 = du_bl_sp(k)
    cu_prod = cq_cm(k)
    do k = k_blend-1,0,-1
      du_1 = du_1 + ( (-1) ** k_blend ) * du_bl_sp(k) * cu_prod
      cu_prod = cu_prod * cq_cm(k)
    end do
    k = 0
    cq_cm_1 = ( (-1) ** (k_blend + k) ) * cu_prod

  end subroutine matrix_sweep

end module bl_imp_du_kernel_mod
