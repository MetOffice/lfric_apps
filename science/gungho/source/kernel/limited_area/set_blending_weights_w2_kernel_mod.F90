!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Populates the w2 weights field with the blending weights obtained from the
!>        configuration.
!> @details Uses the onion layer values as array index for blending_weights array
!>          to populate a field of blending weights.
!>          The onion layers are on the discontinuous W3 space.  For the W2 dofs,
!>          average values are taken using the weights from either side of the face.
module set_blending_weights_w2_kernel_mod

  use argument_mod,              only : arg_type,             &
                                        GH_SCALAR, GH_FIELD,  &
                                        GH_READ, GH_REAL,     &
                                        GH_INC,               &
                                        GH_INTEGER, GH_BASIS, &
                                        CELL_COLUMN,          &
                                        STENCIL, CROSS
  use fs_continuity_mod,         only : W3, W2
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type
  use reference_element_mod,     only : T, B

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: set_blending_weights_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                              &
         arg_type(GH_FIELD,   GH_REAL, GH_INC, W2),                  & ! weights_field
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3, STENCIL(CROSS)), & ! onion_layers
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                   & ! depth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_blending_weights_w2_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: set_blending_weights_w2_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in,out] weights_field The output field
!> @param[in] onion_layers The input field
!> @param[in] depth     Depth of the weight array
!> @param[in] ndf_out   Number of degrees of freedom for weights_field
!> @param[in] undf_out  Total number of degrees of freedom for weights_field
!> @param[in] map_out   Dofmap for the cell at the base of the column for weights_field
!> @param[in] ndf_in    Number of degrees of freedom for onion_layers
!> @param[in] undf_in   Total number of degrees of freedom for onion_layers
!> @param[in] map_in    Dofmap for the cell at the base of the column for onion_layers
subroutine set_blending_weights_w2_code( nlayers,       &
                                         weights_field, &
                                         onion_layers,  &
                                         stencil_size,  &
                                         stencil_map,   &
                                         depth,     &
                                         ndf_out,   &
                                         undf_out,  &
                                         map_out,   &
                                         ndf_in,    &
                                         undf_in,   &
                                         map_in)

  use boundaries_config_mod,        only : blending_weights

  implicit none

  ! Arguments
  integer(kind=i_def),                                 intent(in)    :: nlayers
  integer(kind=i_def),                                 intent(in)    :: ndf_out, undf_out
  integer(kind=i_def),                                 intent(in)    :: ndf_in, undf_in
  real(kind=r_def), dimension(undf_out),               intent(inout) :: weights_field
  real(kind=r_def), dimension(undf_in),                intent(in)    :: onion_layers
  integer(kind=i_def),                                 intent(in)    :: stencil_size
  integer(kind=i_def), dimension(ndf_in,stencil_size), intent(in)    :: stencil_map
  integer(kind=i_def),                                 intent(in)    :: depth
  integer(kind=i_def), dimension(ndf_out),             intent(in)    :: map_out
  integer(kind=i_def), dimension(ndf_in),              intent(in)    :: map_in

  ! Internal variables
  integer(kind=i_def) :: k, df
  integer(kind=i_def) :: index
  integer(kind=i_def) :: onion_layer, cell_next_layer

  onion_layer = int(onion_layers(stencil_map(1,1)))

  ! W2 weights should not go right up to the inner region
  if (onion_layer > 0_i_def)then
    index = depth - onion_layer + 1

    ! Vertical dofs first - these are all set to the blending weight

    do k=0,nlayers-1
      do df=B,T
        weights_field(map_out(df)+k) = blending_weights(index)
      end do
    end do

    ! Next the horizontal dofs

    do k=0,nlayers-1
      do df=1,stencil_size-1
        cell_next_layer = int(onion_layers(stencil_map(1,df+1)))
        ! If the onion_layer > 1, then we write to all dofs
        ! If the onion_layer == 1, then only write to dofs with
        ! neighbouring cells that are within the blending region
        if ( cell_next_layer > 0_i_def .or. onion_layer > 1_i_def ) then
          ! W2 weights are averages of neighouring W3 weights
          weights_field(map_out(df)+k) = min( 1.0_r_def, weights_field(map_out(df)+k) &
                                                   + 0.5_r_def*blending_weights(index) )
        end if
      end do
    end do

  end if

end subroutine set_blending_weights_w2_code

end module set_blending_weights_w2_kernel_mod
