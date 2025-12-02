!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Turbulent mixing of a tracer
module tracer_mix_kernel_mod

  use argument_mod,           only : arg_type,                   &
                                     GH_FIELD, GH_REAL,          &
                                     GH_INTEGER,                 &
                                     GH_READ,                    &
                                     GH_READWRITE, DOMAIN,       &
                                     ANY_DISCONTINUOUS_SPACE_1,  &
                                     ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: tracer_mix_kernel_type
    private
    type(arg_type) :: meta_args(16) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                 &! tracer
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   W3),                       &! height_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   WTHETA),                   &! height_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   W3),                       &! rdz_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1),&! level_ent
         arg_type(GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1),&! level_ent_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! ent_we_lim
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! ent_t_frac
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! ent_zrzi
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! ent_we_lim_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! ent_t_frac_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2) &! ent_zrzi_dsc
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: tracer_mix_code
  end type

  public :: tracer_mix_code

contains

  !> @details Mix an arbitrary tracer due to turbulent processes
  !>           using the UM tr_mix routine
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     tracer                 Tracer to be mixed
  !> @param[in]     height_w3              Height of density space above surface
  !> @param[in]     height_wth             Height of theta space above surface
  !> @param[in]     rdz_w3                 Inverse Layer depths at w3 points
  !> @param[in]     rhokh_bl               Heat eddy diffusivity on BL levels
  !> @param[in]     dtrdz_tq_bl            dt/(rho*r*r*dz) in wth
  !> @param[in]     zh_nonloc              Depth of non-local BL scheme
  !> @param[in]     zhsc_2d                Height of decoupled layer top
  !> @param[in]     level_ent              Level of surface mixed layer inversion
  !> @param[in]     level_ent_dsc          Level of decoupled stratocumulus inversion
  !> @param[in]     ent_we_lim             Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
  !> @param[in]     ent_t_frac             Fraction of time surface ML inversion is above level
  !> @param[in]     ent_zrzi               Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in]     ent_we_lim_dsc         Rho * entrainment rate at DSC inversion (kg m-2 s-1)
  !> @param[in]     ent_t_frac_dsc         Fraction of time DSC inversion is above level
  !> @param[in]     ent_zrzi_dsc           Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in]     ndf_wth                Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth               Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth                Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3                 Number of DOFs per cell for density space
  !> @param[in]     undf_w3                Number of unique DOFs for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Number of unique DOFs for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_ent                Number of DOFs per cell for entrainment levels
  !> @param[in]     undf_ent               Number of total DOFs for entrainment levels
  !> @param[in]     map_ent                Dofmap for cell for entrainment levels
  subroutine tracer_mix_code(nlayers, seg_len,                      &
                             tracer,                                &
                             height_w3,                             &
                             height_wth,                            &
                             rdz_w3,                                &
                             rhokh_bl,                              &
                             dtrdz_tq_bl,                           &
                             zh_nonloc,                             &
                             zhsc_2d,                               &
                             level_ent,                             &
                             level_ent_dsc,                         &
                             ent_we_lim,                            &
                             ent_t_frac,                            &
                             ent_zrzi,                              &
                             ent_we_lim_dsc,                        &
                             ent_t_frac_dsc,                        &
                             ent_zrzi_dsc,                          &
                             ndf_wth, undf_wth, map_wth,            &
                             ndf_w3, undf_w3, map_w3,               &
                             ndf_2d, undf_2d, map_2d,               &
                             ndf_ent, undf_ent, map_ent)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use atm_fields_bounds_mod, only: pdims
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: planet_radius

    ! subroutines used
    use tr_mix_mod, only: tr_mix

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth, seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3, seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d, seg_len)
    integer(kind=i_def), intent(in) :: ndf_ent, undf_ent
    integer(kind=i_def), intent(in) :: map_ent(ndf_ent, seg_len)

    real(kind=r_def), dimension(undf_wth), intent(inout):: tracer

    real(kind=r_def), dimension(undf_w3),  intent(in)   :: height_w3, rdz_w3,  &
                                                           rhokh_bl
    real(kind=r_def), dimension(undf_wth), intent(in)   :: height_wth,         &
                                                           dtrdz_tq_bl
    real(kind=r_def), dimension(undf_2d), intent(in) :: zhsc_2d, zh_nonloc

    integer(kind=i_def), dimension(undf_2d), intent(in) :: level_ent
    integer(kind=i_def), dimension(undf_2d), intent(in) :: level_ent_dsc
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_we_lim
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_t_frac
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_zrzi
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_we_lim_dsc
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_t_frac_dsc
    real(kind=r_def), dimension(undf_ent),  intent(in)  :: ent_zrzi_dsc

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    integer(i_def) :: k, i

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1,nlayers) :: z_rho, r_rho_levels

    ! profile field on boundary layer levels
    real(r_um), dimension(seg_len,1,bl_levels) :: dtrdz_charney_grid, &
         tracer_mixed, tracer_flux

    ! profile fields from level 2 upwards
    real(r_um), dimension(seg_len,1,2:bl_levels) :: rhokh_mix_bl

    ! profile fields from level 0 upwards
    real(r_um), dimension(seg_len,1,0:nlayers) :: r_theta_levels

    ! single level real fields
    real(r_um), dimension(seg_len,1) :: zhnl, zhsc, surf_dep_flux, zeroes, &
         rhokh_mix_surf

    real(r_um), dimension(seg_len,1,3) :: t_frac, t_frac_dsc, we_lim, &
         we_lim_dsc, zrzi, zrzi_dsc

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: kent, kent_dsc

    real(r_um), dimension(bl_levels) :: alpha_tr

    alpha_tr = 1.0_r_um
    zeroes = 0.0_r_um

    do i = 1, seg_len
      do k = 0, nlayers
        ! height of theta levels from centre of planet
        r_theta_levels(i,1,k) = height_wth(map_wth(1,i) + k) + planet_radius
      end do
      do k = 1, nlayers
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
        ! height of levels above surface
        z_rho(i,1,k) = r_rho_levels(i,1,k)-r_theta_levels(i,1,0)
      end do
    end do

    ! surface and interior scalar diffusivity
    do i = 1, seg_len
      rhokh_mix_surf(i,1) = rhokh_bl(map_w3(1,i))
      do k = 2, bl_levels
        rhokh_mix_bl(i,1,k) = rhokh_bl(map_w3(1,i) + k-1) * rdz_w3(map_w3(1,i) + k-1)
      end do
    end do

    ! mixed layer depths
    do i = 1, seg_len
      zhnl(i,1) = zh_nonloc(map_2d(1,i))
      zhsc(i,1) = zhsc_2d(map_2d(1,i))
    end do

    ! tracer to mix
    do i = 1, seg_len
      do k=1,bl_levels
        dtrdz_charney_grid(i,1,k) = dtrdz_tq_bl(map_wth(1,i) + k)
        tracer_mixed(i,1,k) = tracer(map_wth(1,i) + k)
      end do
    end do

    ! Entrainment mixing parameters
    do i = 1, seg_len
      kent(i,1) = level_ent(map_2d(1,i))
      kent_dsc(i,1) = level_ent_dsc(map_2d(1,i))
      do k = 1, 3
        we_lim(i,1,k) = ent_we_lim(map_ent(1,i) + k - 1)
        t_frac(i,1,k) = ent_t_frac(map_ent(1,i) + k - 1)
        zrzi(i,1,k) = ent_zrzi(map_ent(1,i) + k - 1)
        we_lim_dsc(i,1,k) = ent_we_lim_dsc(map_ent(1,i) + k - 1)
        t_frac_dsc(i,1,k) = ent_t_frac_dsc(map_ent(1,i) + k - 1)
        zrzi_dsc(i,1,k) = ent_zrzi_dsc(map_ent(1,i) + k - 1)
      end do
    end do

    call  tr_mix (                                                           &
         ! IN fields
         r_theta_levels, r_rho_levels, pdims, bl_levels, alpha_tr,           &
         rhokh_mix_bl, rhokh_mix_surf, dtrdz_charney_grid, zeroes, zeroes,   &
         kent, we_lim, t_frac, zrzi,                                         &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                         &
         zhnl, zhsc, z_rho,                                                  &
         ! INOUT / OUT fields
         tracer_mixed, tracer_flux, surf_dep_flux                            &
         )

    ! copy back mixed variable, and update 0th level
    do k = 1, bl_levels
      do i = 1, seg_len
        tracer(map_wth(1,i)+k) = tracer_mixed(i,1,k)
      end do
    end do
    do i = 1, seg_len
      tracer(map_wth(1,i)+0) = tracer(map_wth(1,i)+1)
    end do

  end subroutine tracer_mix_code

end module tracer_mix_kernel_mod
