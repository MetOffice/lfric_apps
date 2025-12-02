!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Apply 3D viscosity mu * (d2dx2 + d2dy2 + d2dz2) to the components
!>        of the momentum equation for lowest order elements.
!>
module momentum_viscosity_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_WRITE,         &
                                GH_SCALAR,                 &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: momentum_viscosity_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                 &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE,  W2),                  &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2, STENCIL(CROSS)),  &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2, STENCIL(CROSS)),  &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: momentum_viscosity_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: momentum_viscosity_code

contains

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] u_inc Diffusion increment for wind field
!! @param[in] u_n Input wind field
!! @param[in] map_w2_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2 Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] dx_at_w2 Grid length at the w2 dofs
!! @param[in] map_dx_stencil_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_dx_stencil Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] viscosity_mu Viscosity constant
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  Number of unique degrees of freedom  for wind_space
!! @param[in] cell_map_w2 Array holding the dofmap for the cell at the base of the column for w2
subroutine momentum_viscosity_code(nlayers,                               &
                                   u_inc, u_n,                            &
                                   map_w2_size, map_w2,                   &
                                   dx_at_w2,                              &
                                   map_dx_stencil_size, map_dx_stencil,   &
                                   viscosity_mu,                          &
                                   ndf_w2, undf_w2, cell_map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2_size, map_dx_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_size), intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_w2,map_dx_stencil_size), intent(in)  :: map_dx_stencil
  integer(kind=i_def), dimension(ndf_w2),             intent(in)  :: cell_map_w2

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n, dx_at_w2

  real(kind=r_def), intent(in) :: viscosity_mu

  ! Internal variables
  integer(kind=i_def)                      :: k, km, kp, df
  real(kind=r_def)                         :: d2dx, d2dy, d2dz
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2, idz2
  real(kind=r_def), dimension(0:nlayers-1,4) :: idx2_w2, idy2_w2, idz2_w2

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! If the full stencil isn't available, we must be at the domain edge.
  ! The increment is already 0, so we just exit the routine.
  if (map_w2_size < 5_i_def) then
    return
  end if

  ! Loop over horizontal faces of the cell for horizontal diffusion
  do df = 1,4

    ! Only calculate face if it's not already been done
    if (u_inc(cell_map_w2(df)+1) == 0.0_r_def) then

      ! Compute horizontal grid spacing
      if (df == 1 .or. df == 3) then
        ! For Dofs 1 & 3, dx is given by the input field,
        ! whilst dy needs to be computed from the 4 neighbouring dy points
        do k = 0, nlayers - 1
          idx2_w2(k,df) = 1.0_r_def/(dx_at_w2(map_dx_stencil(df,1)+k))**2
          idy2_w2(k,df) = (4.0_r_def/(dx_at_w2(map_dx_stencil(2,1)+k)+dx_at_w2(map_dx_stencil(4,1)+k)+ &
                                      dx_at_w2(map_dx_stencil(2,df+1)+k)+dx_at_w2(map_dx_stencil(4,df+1)+k)))**2
          idz2_w2(k,df) = (2.0_r_def/(dx_at_w2(map_dx_stencil(5,1)+k)+dx_at_w2(map_dx_stencil(5,df+1)+k)))**2
        end do
      else
        ! For Dofs 2 & 4, dy is given by the input field,
        ! whilst dx needs to be computed from the 4 neighbouring dx points
        do k = 0, nlayers - 1
          idx2_w2(k,df) = (4.0_r_def/(dx_at_w2(map_dx_stencil(1,1)+k)+dx_at_w2(map_dx_stencil(3,1)+k)+ &
                                      dx_at_w2(map_dx_stencil(1,df+1)+k)+dx_at_w2(map_dx_stencil(3,df+1)+k)))**2
          idy2_w2(k,df) = 1.0_r_def/(dx_at_w2(map_dx_stencil(df,1)+k))**2
          idz2_w2(k,df) = (2.0_r_def/(dx_at_w2(map_dx_stencil(5,1)+k)+dx_at_w2(map_dx_stencil(5,df+1)+k)))**2
        end do
      end if

      ! Horizontal Velocity diffusion
      k = 0
      km = k
      kp = k+1

      d2dx = (u_n(map_w2(df,2) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,4) + k) )*idx2_w2(k,df)
      d2dy = (u_n(map_w2(df,3) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,5) + k) )*idy2_w2(k,df)
      d2dz = (u_n(map_w2(df,1) + kp) - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,1) + km))*idz2_w2(k,df)
      u_inc(cell_map_w2(df)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

      do k = 1, nlayers-2
        km = k - 1
        kp = k + 1

        d2dx = (u_n(map_w2(df,2) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,4) + k) )*idx2_w2(k,df)
        d2dy = (u_n(map_w2(df,3) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,5) + k) )*idy2_w2(k,df)
        d2dz = (u_n(map_w2(df,1) + kp) - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,1) + km))*idz2_w2(k,df)
        u_inc(cell_map_w2(df)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
      end do

      k = nlayers-1
      km = k - 1
      kp = k

      d2dx = (u_n(map_w2(df,2) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,4) + k) )*idx2_w2(k,df)
      d2dy = (u_n(map_w2(df,3) + k)  - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,5) + k) )*idy2_w2(k,df)
      d2dz = (u_n(map_w2(df,1) + kp) - 2.0_r_def*u_n(map_w2(df,1) + k) + u_n(map_w2(df,1) + km))*idz2_w2(k,df)
      u_inc(cell_map_w2(df)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    end if

  end do

  ! Vertical Velocity diffusion
  do k = 1, nlayers-1
    km = k - 1

    idx2(k) = (2.0_r_def/(dx_at_w2(map_dx_stencil(1,1)+k)+dx_at_w2(map_dx_stencil(3,1)+k)))**2
    idy2(k) = (2.0_r_def/(dx_at_w2(map_dx_stencil(2,1)+k)+dx_at_w2(map_dx_stencil(4,1)+k)))**2
    idz2(k) = 1.0_r_def/(dx_at_w2(map_dx_stencil(5,1)+k))**2

    d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
    u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  end do

end subroutine momentum_viscosity_code

end module momentum_viscosity_kernel_mod
