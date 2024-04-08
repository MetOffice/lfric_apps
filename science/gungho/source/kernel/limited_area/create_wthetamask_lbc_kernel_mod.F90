!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the mask to define the lateral boundary condition (LBC)
!!        region for Wtheta space.
!> @details The coordinates of the Wtheta dofs are compared with the
!!          coordinates of the boundary. The mask takes the value 1 in the
!!          LBC region and 0 otherwise. This mask is used to overwrite data
!!          in the LBC region by the boundary condition data from the driving
!!          model.
module create_wthetamask_lbc_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL, GH_BASIS,   &
                                        CELL_COLUMN, GH_EVALUATOR
  use constants_mod,             only : r_def, i_def, l_def
  use fs_continuity_mod,         only : Wtheta, Wchi
  use kernel_mod,                only : kernel_type
  use base_mesh_config_mod,      only : geometry,            &
                                        geometry_spherical
  use finite_element_config_mod, only : coord_system, &
                                        coord_system_xyz

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_wthetamask_lbc_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                  &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wchi)    &
         /)
    type(func_type) :: meta_funcs(1) = (/                 &
         func_type(Wchi, GH_BASIS)                        &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: create_wthetamask_lbc_code
  end type

    !-------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------
    public :: create_wthetamask_lbc_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in] boundary_s   Southern outer boundary coordinates
!> @param[in] boundary_n   Northern outer boundary coordinates
!> @param[in] boundary_e   Eastern  outer boundary coordinates
!> @param[in] boundary_w   Western  outer boundary coordinates
!> @param[in] rim_s        Southern inner boundary coordinates
!> @param[in] rim_n        Northern inner boundary coordinates
!> @param[in] rim_e        Eastern  inner boundary coordinates
!> @param[in] rim_w        Western  inner boundary coordinates
!> @param[in,out] theta_mask The mask to define the LBC region for Wtheta
!> @param[in] chi_1        X component of the chi coordinate field
!> @param[in] chi_2        Y component of the chi coordinate field
!> @param[in] chi_3        Z component of the chi coordinate field
!> @param[in] ndf_wtheta   Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wtheta  Total number of degrees of freedom for Wtheta
!> @param[in] map_wtheta   Dofmap for the cell at the base of the column for Wtheta
!> @param[in] ndf_chi      Number of degrees of freedom per cell for chi
!> @param[in] undf_chi     Number of degrees of freedom for chi
!> @param[in] map_chi      Dofmap for the cell at the base of the column for chi
!> @param[in] chi_basis    Basis functions evaluated at gaussian quadrature point
subroutine create_wthetamask_lbc_code( nlayers,     &
                                       boundary_s,  &
                                       boundary_n,  &
                                       boundary_e,  &
                                       boundary_w,  &
                                       rim_s,       &
                                       rim_n,       &
                                       rim_e,       &
                                       rim_w,       &
                                       theta_mask,  &
                                       chi_1,       &
                                       chi_2,       &
                                       chi_3,       &
                                       ndf_wtheta,  &
                                       undf_wtheta, &
                                       map_wtheta,  &
                                       ndf_chi,     &
                                       undf_chi,    &
                                       map_chi,     &
                                       chi_basis )

  implicit none

  ! Arguments
  integer(kind=i_def),                              intent(in) :: nlayers,     &
                                                                  ndf_wtheta,  &
                                                                  ndf_chi,     &
                                                                  undf_wtheta, &
                                                                  undf_chi
  integer(kind=i_def), dimension(ndf_wtheta),       intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi),          intent(in) :: map_chi
  real(kind=r_def),                                 intent(in) :: boundary_s,  &
                                                                  boundary_n,  &
                                                                  boundary_e,  &
                                                                  boundary_w,  &
                                                                  rim_s,       &
                                                                  rim_n,       &
                                                                  rim_e,       &
                                                                  rim_w
  real(kind=r_def), dimension(undf_wtheta),         intent(inout) :: theta_mask
  real(kind=r_def), dimension(1,ndf_chi,ndf_wtheta),intent(in)    :: chi_basis
  real(kind=r_def), dimension(undf_chi),            intent(in)    :: chi_1,    &
                                                                     chi_2,    &
                                                                     chi_3

  ! Internal variables
  integer(kind=i_def)                    :: k, df1, df2
  real(kind=r_def), dimension(ndf_chi)   :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: x(3),           &
                                            query_value_ns, &
                                            query_value_ew, &
                                            mask_value
  logical(kind=l_def)                    :: innerbox, outerbox

  ! Define for the bottom vertical layer and copy to other layers
  k = 0

  ! Get the chi values in the kth layer
  ! chi_e (u) = chi ( f(u) + k )
  do df2 = 1, ndf_chi
    chi_1_e(df2) = chi_1( map_chi(df2) + k)
    chi_2_e(df2) = chi_2( map_chi(df2) + k)
    chi_3_e(df2) = chi_3( map_chi(df2) + k)
  end do

  ! Get the x values (x,y,z) for each cell centre, by summing
  ! over W0 to Wtheta basis functions
  ! x(v) = sum_u chi_e(u) * g(u,v)
  do df1 = 1, ndf_wtheta
    x(:) = 0.0_r_def
    do df2 = 1,ndf_chi
      x(1) = x(1) + chi_1_e(df2)*chi_basis(1,df2,df1)
      x(2) = x(2) + chi_2_e(df2)*chi_basis(1,df2,df1)
      x(3) = x(3) + chi_3_e(df2)*chi_basis(1,df2,df1)
    end do

    ! Change of coordinates for spherical geometry
    if ( geometry == geometry_spherical .and. &
         coord_system == coord_system_xyz ) then
      ! in alpha beta space - and only pick the face where x(2)>0
      if ( x(2) > 0.0_r_def ) then
        query_value_ns = atan2(x(3),x(2))
        query_value_ew = atan2(x(1),x(2))
      else
        query_value_ns = 10.0_r_def
        query_value_ew = 10.0_r_def
      endif
    else
      ! in cartesian space
      query_value_ns = x(2)
      query_value_ew = x(1)
    end if

    ! Define the LBC mask by querying whether the dof is in the inner box
    ! and the outerbox.

    ! 000000000000000000000000
    ! 000-----outerbox-----000
    ! 000|1111111111111111|000
    ! 000|1111111111111111|000
    ! 000|11--innerbox--11|000
    ! 000|11|0000000000|11|000
    ! 000|11|0000000000|11|000
    ! 000|11------------11|000
    ! 000|1111111111111111|000
    ! 000|1111111111111111|000
    ! 000------------------000
    ! 000000000000000000000000

    if ( query_value_ns > boundary_s .and. &
         query_value_ns < boundary_n .and. &
         query_value_ew > boundary_w .and. &
         query_value_ew < boundary_e ) then
      outerbox = .true.
    else
      outerbox = .false.
    endif

    if ( query_value_ns > rim_s .and. &
         query_value_ns < rim_n .and. &
         query_value_ew > rim_w .and. &
         query_value_ew < rim_e ) then
      innerbox = .true.
    else
      innerbox = .false.
    endif

    mask_value = 1.0_r_def

    ! Zero outside the outerbox
    if ( .not. outerbox) then
      mask_value = 0.0_r_def
    endif

    ! Zero in the inner box
    if ( innerbox ) then
      mask_value = 0.0_r_def
    endif

    ! Loop over vertical layers
    do k = 0, nlayers-1
      theta_mask(map_wtheta(df1) + k) = mask_value
    enddo

  end do

end subroutine create_wthetamask_lbc_code

end module create_wthetamask_lbc_kernel_mod
