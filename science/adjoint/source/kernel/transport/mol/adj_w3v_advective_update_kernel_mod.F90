!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of vertical advective update through fitting
!!        a high order upwind reconstruction.
module adj_w3v_advective_update_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_OPERATOR,       &
                                GH_READWRITE, GH_READ, &
                                ANY_DISCONTINUOUS_SPACE_1, &
                                ANY_W2, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def, r_tran
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private
  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: adj_w3v_advective_update_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                           &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_W2),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, W3)                     &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_w3v_advective_update_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: adj_w3v_advective_update_code
  contains

  !> @brief Computes the adjoint horizontal advective update for a tracer in W3.
  !> @param[in]     cell                Horizontal cell index
  !> @param[in]     nlayers             Number of layers
  !> @param[in]     advective_increment Advective update field to compute
  !> @param[in,out] tracer              Pointwise tracer field to advect stored on cell faces
  !> @param[in]     wind                Wind field
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
  subroutine adj_w3v_advective_update_code( cell,                &
                                            nlayers,             &
                                            advective_increment, &
                                            tracer,              &
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

    real(kind=r_tran), dimension(undf_w3), intent(in)    :: advective_increment
    real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind
    real(kind=r_tran), dimension(undf_md), intent(inout) :: tracer

    real(kind=r_def), dimension(ncell_3d, ndf_w3, ndf_w3), intent(in) :: m3_inv

    ! Internal variables
    integer(kind=i_def) :: k, ik, df, offset
    real(kind=r_tran)    :: w, dtdz

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
      w =  0.5_r_tran*( wind(map_w2(df) + k) + wind(map_w2(df) + k + 1) )
      ik = 1 + k + (cell-1)*nlayers
      dtdz = advective_increment(map_w3(1) + k) * w * real(m3_inv(ik,1,1), kind=r_tran)

      if (w <= 0.0_r_tran .and. k < nlayers - 1) then
        tracer(map_md(1) + offset + k + 1) = tracer(map_md(1) + offset + k + 1) + dtdz
      else
        tracer(map_md(1) + offset + nlayers + k) = tracer(map_md(1) + offset + nlayers + k) + dtdz
      end if
      if (w > 0.0_r_tran .and. k > 0) then
        tracer(map_md(1) + offset + nlayers + k - 1) = tracer(map_md(1) + offset + nlayers + k - 1) - dtdz
      else
        tracer(map_md(1) + offset + k) = tracer(map_md(1) + offset + k) - dtdz
      end if
    end do

  end subroutine adj_w3v_advective_update_code

end module adj_w3v_advective_update_kernel_mod
