!-------------------------------------------------------------------------------
!(c) Crown copyright 2023 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the SKEB2 numerical dissipation
module skeb_biharm_diss_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,          &
                               GH_REAL, GH_WRITE, GH_READ,  &
                               CELL_COLUMN, GH_INTEGER,     &
                               GH_SCALAR, STENCIL, CROSS,   &
                               GH_LOGICAL
  use fs_continuity_mod, only: W3, Wtheta, W1, W2
  use constants_mod,     only: r_def, i_def, l_def
  use kernel_mod,        only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: skeb_biharm_diss_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                      &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                & ! ndisp
    arg_type(GH_FIELD, GH_REAL, GH_READ, W1),                 & ! vorticity
    arg_type(GH_FIELD, GH_REAL, GH_READ, W3, STENCIL(CROSS)), & ! divergence
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                & ! dx_in_w2
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                & ! skeb_level_bottom
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                & ! skeb_level_top
    arg_type(GH_SCALAR, GH_REAL, GH_READ ),                   & ! dt
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                & ! norm_xi
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                & ! norm_div
    arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                 & ! norm_xi_flag
    arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                  & ! norm_div_flag

    /)
    integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: skeb_biharm_diss_code
  end type skeb_biharm_diss_kernel_type

  public skeb_biharm_diss_code
contains

  !> @brief Calculate the numerical dissipation for skeb
  !> @param[in]    nlayers     The number of layers
  !> @param[inout] ndisp       Numerical dissipation
  !> @param[in]    vorticity   2D Vorticity
  !> @param[in]    divergence  2D Divergence
  !> @param[in]    dx_in_w2    delta-x and w2 dofs
  !> @param[in]    skeb_level_bottom Bottom SKEB level
  !> @param[in]    skeb_level_top    Top SKEB level
  !> @param[inout] norm_xi     Norm of the vorticity gradient
  !> @param[inout] norm_div    Norm of the divergence gradient
  !> @param[in]    ndf_w3      Number of DOFs per cell for density space
  !> @param[in]    undf_w3     Number of unique DOFs  for density space
  !> @param[in]    map_w3      dofmap for the cell at the base of the column for density space
  !> @param[in]    ndf_w1      Number of DOFs per cell for w1 space
  !> @param[in]    undf_w1     Number of unique DOFs  for w1 space
  !> @param[in]    map_w1      dofmap for the cell at the base of the column for w1 space
  !> @param[in]    ndf_w2      Number of DOFs per cell for w2 space
  !> @param[in]    undf_w2     Number of unique DOFs  for w2 space
  !> @param[in]    map_w2      dofmap for the cell at the base of the column for w2 space
  !> @param[in]    norm_xi_flag  Control whether norm_xi calculation is needed
  !> @param[in]    norm_div_flag  Control whether norm_div calculation is needed

  subroutine skeb_biharm_diss_code(nlayers,           &
                                   ndisp,             &
                                   vorticity,         &
                                   divergence,        &
                                   map_w3_sten_size,  &
                                   map_w3_sten,       &
                                   dx_at_w2,          &
                                   skeb_level_bottom, &
                                   skeb_level_top,    &
                                   dt,                &
                                   norm_xi,           &
                                   norm_div,          &
                                   norm_xi_flag,      &
                                   norm_div_flag,     &
                                   ndf_w3,            &
                                   undf_w3,           &
                                   map_w3,            &
                                   ndf_w1,            &
                                   undf_w1,           &
                                   map_w1,            &
                                   ndf_w2,            &
                                   undf_w2,           &
                                   map_w2)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, ndf_w2, ndf_w1
    integer(kind=i_def), intent(in) :: undf_w3, undf_w2, undf_w1
    integer(kind=i_def), intent(in) :: map_w3_sten_size
    integer(kind=i_def), intent(in),  dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in),  dimension(ndf_w2)   :: map_w2
    integer(kind=i_def), intent(in),  dimension(ndf_w3,map_w3_sten_size) :: map_w3_sten
    integer(kind=i_def), intent(in),  dimension(ndf_w1)   :: map_w1
    logical(kind=l_def), intent(in) :: norm_xi_flag
    logical(kind=l_def), intent(in) :: norm_div_flag

    ! Fields
    real(kind=r_def),    intent(in),    dimension(undf_w1)  :: vorticity
    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: divergence
    real(kind=r_def),    intent(in),    dimension(undf_w2)  :: dx_at_w2
    real(kind=r_def),    intent(inout), dimension(undf_w3)  :: ndisp
    real(kind=r_def), pointer, intent(inout) :: norm_xi(:), norm_div(:)

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top
    real(kind=r_def), intent(in) :: dt

    integer(kind=i_def) :: k, df
    real(kind=r_def) :: amp_k, biharmonic_x_xi, biharmonic_y_xi, &
         biharmonic_x_div, biharmonic_y_div

    amp_k = 3.0_r_def / (128.0_r_def * dt)

    ! K is the biharmonic viscosity tensor. K = 3/128 x delta_i^4 / delta_t
    ! The tensorial decomposition of Dnum is:
    ! Dnum= Kx {(d(vor)/dx)**2 + (d(div)/dx)**2 ) +
    !       Ky {(d(vor)/dy)**2 + (d(div)/dy)**2 ) +
    !
    ! Derivates are averaged between the boxes and  dx^2 from derivates are
    ! cross out. The final expresion is thus
    ! Dnum_x=[ d(vor)**2 + d(div)**2 ) ] * (dx)**2
    !
    ! d(dvor)/dx is at w2 so it does need averaging between
    ! faces to place it into w3 points
    ! d(div)/dx is at w2 so it does averaging between cells

    do k = skeb_level_bottom, skeb_level_top

      biharmonic_x_xi = 0.5_r_def*( ((vorticity(map_w1(7)+k-1) -   &
                                      vorticity(map_w1(8)+k-1)) *  &
                                    dx_at_w2(map_w2(3)+k-1) )**2 + &
                                    ((vorticity(map_w1(6)+k-1) -   &
                                      vorticity(map_w1(5)+k-1)) *  &
                                    dx_at_w2(map_w2(3)+k-1) )**2 )
      biharmonic_x_div = 0.5_r_def*( ((divergence(map_w3_sten(1,4)+k-1) -   &
                                       divergence(map_w3_sten(1,1)+k-1)) *  &
                                     dx_at_w2(map_w2(3)+k-1) )**2 +         &
                                     ((divergence(map_w3_sten(1,1)+k-1) -   &
                                       divergence(map_w3_sten(1,2)+k-1)) *  &
                                     dx_at_w2(map_w2(1)+k-1) )**2 )

      biharmonic_y_xi = 0.5_r_def*( ((vorticity(map_w1(7)+k-1) -    &
                                      vorticity(map_w1(6)+k-1)) *   &
                                    dx_at_w2(map_w2(4)+k-1) )**2 +  &
                                    ((vorticity(map_w1(8)+k-1) -    &
                                      vorticity(map_w1(5)+k-1)) *   &
                                    dx_at_w2(map_w2(4)+k-1) )**2 )
      biharmonic_y_div = 0.5_r_def*( ((divergence(map_w3_sten(1,5)+k-1) -  &
                                       divergence(map_w3_sten(1,1)+k-1)) * &
                                     dx_at_w2(map_w2(4)+k-1) )**2 +        &
                                     ((divergence(map_w3_sten(1,1)+k-1) -  &
                                       divergence(map_w3_sten(1,3)+k-1)) * &
                                     dx_at_w2(map_w2(2)+k-1) )**2 )

      ndisp(map_w3(1)+k-1) = (biharmonic_x_div + biharmonic_y_div + &
                              biharmonic_x_xi + biharmonic_y_xi) * amp_K

      if (norm_xi_flag) then
        norm_xi(map_w3(1)+k-1) = biharmonic_x_xi + biharmonic_y_xi
      end if
      if (norm_div_flag) then
        norm_div(map_w3(1)+k-1) = biharmonic_x_div + biharmonic_y_div
      end if

    end do

  end subroutine skeb_biharm_diss_code

end module skeb_biharm_diss_kernel_mod
