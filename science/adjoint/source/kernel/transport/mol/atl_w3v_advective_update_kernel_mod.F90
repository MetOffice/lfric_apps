!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint vertical advective update through fitting a high order
!!        upwind reconstruction.

module atl_w3v_advective_update_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_OPERATOR,       &
                                GH_INC, GH_READ,   &
                                ANY_DISCONTINUOUS_SPACE_1, &
                                ANY_W2, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private
  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: atl_w3v_advective_update_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                      &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, ANY_W2),                    &
         arg_type(GH_FIELD,    GH_REAL, GH_INC,  ANY_W2),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ, W3, W3)                     &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: atl_w3v_advective_update_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: atl_w3v_advective_update_code
  contains

  !> @brief Computes the adjoint vertical advective update for a tracer in W3.
  !> @param[in]     cell                Horizontal cell index
  !> @param[in]     nlayers             Number of layers
  !> @param[in]     advective_increment Advective update field to compute
  !> @param[in]     tracer              Pointwise tracer field to advect stored on cell faces
  !> @param[in]     ls_wind             Linear wind field
  !> @param[in,out] wind                Perturbation wind field
  !> @param[in]     ncell_3d            Total number of cells
  !> @param[in]     m3_inv              Inverse mass matrix for W3 space
  !> @param[in]     ndf_w3              Number of degrees of freedom per cell
  !> @param[in]     undf_w3             Number of unique degrees of freedom for the
  !!                                    advective_update field
  !> @param[in]     map_w3              Dofmap for the cell at the base of the column
  !> @param[in]     ndf_md              Number of degrees of freedom per cell
  !> @param[in]     undf_md             Number of unique degrees of freedom for the
  !!                                    tracer field
  !> @param[in]     map_md              Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w2              Number of degrees of freedom per cell for the wind fields
  !> @param[in]     undf_w2             Number of unique degrees of freedom for the wind fields
  !> @param[in]     map_w2              Dofmap for the cell at the base of the column for the wind fields
  subroutine atl_w3v_advective_update_code( cell,                &
                                            nlayers,             &
                                            advective_increment, &
                                            tracer,              &
                                            ls_wind,             &
                                            wind,                &
                                            ncell_3d,            &
                                            m3_inv,              &
                                            ndf_w3,              &
                                            undf_w3,             &
                                            map_w3,              &
                                            ndf_md,              &
                                            undf_md,             &
                                            map_md,              &
                                            ndf_w2,              &
                                            undf_w2,             &
                                            map_w2               &
                                           )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                     :: nlayers, cell, ncell_3d
    integer(kind=i_def), intent(in)                     :: ndf_w3, ndf_w2, ndf_md
    integer(kind=i_def), intent(in)                     :: undf_w3, undf_w2, undf_md
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_md),  intent(in) :: map_md

    real(kind=r_def), dimension(undf_w3), intent(in)    :: advective_increment
    real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
    real(kind=r_def), dimension(undf_w2), intent(inout) :: wind
    real(kind=r_def), dimension(undf_md), intent(in)    :: tracer

    real(kind=r_def), dimension(ncell_3d, ndf_w3, ndf_w3), intent(in) :: m3_inv

    ! Internal variables
    integer(kind=i_def) :: k, ik, df, offset
    real(kind=r_def)    :: ls_w, w, dtdz, t_U, t_D

    if ( ndf_w2 == 2 ) then
      ! W2v space, reconstruction has ndata=2
      df = 1
      offset = 0
    else
      ! W2 space, reconstruction has ndata=6
      df = 5
      offset = 4*nlayers
    end if

    ! Vertical advective update
    do k = nlayers - 1, 0, -1
      ls_w =  0.5_r_def*( ls_wind(map_w2(df) + k) + ls_wind(map_w2(df) + k + 1) )

      if ( ls_w > 0.0_r_def .and. k > 0 ) then
        ! Overwrite t_D from cell below
        t_D = tracer(map_md(1) + offset + nlayers + k - 1)
      else
        ! Values from this cell
        t_D = tracer(map_md(1) + offset + k)
      end if

      if ( ls_w <= 0.0_r_def .and. k < nlayers-1 ) then
        ! Overwrite t_U from cell above
        t_U = tracer(map_md(1) + offset + k + 1)
      else
        ! Values from this cell
        t_U = tracer(map_md(1) + offset + nlayers + k)
      end if

      dtdz = t_U - t_D
      ik = 1 + k + (cell-1)*nlayers
      w = dtdz * advective_increment(map_w3(1) + k) * m3_inv(ik,1,1)
      wind(k + map_w2(df)) = wind(k + map_w2(df)) + 0.5_r_def * w
      wind(k + map_w2(df) + 1) = wind(k + map_w2(df) + 1) + 0.5_r_def * w
    end do

  end subroutine atl_w3v_advective_update_code

end module atl_w3v_advective_update_kernel_mod
