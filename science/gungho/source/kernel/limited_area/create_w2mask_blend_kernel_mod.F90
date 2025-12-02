!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the W2 mask for the lateral boundary condition (LBC)
!!        blending region.
!> @details The mask is defined using the W2 dofs. The coordinates of the dofs
!!          are compared with the coordinates of 3 boundaries: outer, rim and
!!          blend. From outer to blend, the mask takes the value 1. From blend
!!          to rim, the mask takes a value that ramps linearly from 1 to 0.
!!          The mask is used to blend the data in the LBC region with the
!!          boundary conditions from the driving model.
module create_w2mask_blend_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_REAL, GH_INTEGER, &
                                        GH_READ, GH_WRITE,   &
                                        GH_BASIS,            &
                                        CELL_COLUMN, GH_EVALUATOR
  use constants_mod,             only : r_def, i_def, l_def
  use fs_continuity_mod,         only : W2, Wchi
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

  type, public, extends(kernel_type) :: create_w2mask_blend_kernel_type
    private
    type(arg_type) :: meta_args(16) = (/                  &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),       &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),       &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2),  &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi) &
         /)
    type(func_type) :: meta_funcs(1) = (/                 &
         func_type(Wchi, GH_BASIS)                        &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: create_w2mask_blend_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_w2mask_blend_code

contains

!> @param[in] nlayers       Number of layers
!> @param[in] outer_s       Southern outer LBC boundary coordinates
!> @param[in] outer_n       Northern outer LBC boundary coordinates
!> @param[in] outer_e       Eastern outer LBC boundary coordinates
!> @param[in] outer_w       Western outer LBC boundary coordinates
!> @param[in] rim_s         Southern inner LBC boundary coordinates
!> @param[in] rim_n         Northern inner LBC boundary coordinates
!> @param[in] rim_e         Eastern inner LBC boundary coordinates
!> @param[in] rim_w         Western inner LBC boundary coordinates
!> @param[in] blend_s       Southern blending boundary coordinates
!> @param[in] blend_n       Northern blending boundary coordinates
!> @param[in] blend_e       Eastern blending boundary coordinates
!> @param[in] blend_w       Western blending boundary coordinates
!> @param[in] rim_width_ns  If ge 0 apply LBC in north south
!> @param[in] rim_width_ew  If ge 0 apply LBC in east west
!> @param[in] ndf_w2        Number of degrees of freedom per cell
!> @param[in] undf_w2       Total number of degrees of freedom
!> @param[in] map_w2        Dofmap for the cell at the base of the column
!> @param[in,out] wind_mask Blending mask for W2 fields
!> @param[in] ndf_chi       Number of degrees of freedom per cell for chi
!> @param[in] undf_chi      Number of degrees of freedom for chi
!> @param[in] map_chi       Dofmap for the cell at the base of the column for chi
!> @param[in] chi_1         X component of the chi coordinate field
!> @param[in] chi_2         Y component of the chi coordinate field
!> @param[in] chi_3         Z component of the chi coordinate field
!> @param[in] chi_basis     Basis functions evaluated at gaussian quadrature points
subroutine create_w2mask_blend_code( nlayers,      &
                                     outer_s,      &
                                     outer_n,      &
                                     outer_e,      &
                                     outer_w,      &
                                     rim_s,        &
                                     rim_n,        &
                                     rim_e,        &
                                     rim_w,        &
                                     blend_s,      &
                                     blend_n,      &
                                     blend_e,      &
                                     blend_w,      &
                                     rim_width_ns, &
                                     rim_width_ew, &
                                     wind_mask,    &
                                     chi_1,        &
                                     chi_2,        &
                                     chi_3,        &
                                     ndf_w2,       &
                                     undf_w2,      &
                                     map_w2,       &
                                     ndf_chi,      &
                                     undf_chi,     &
                                     map_chi,      &
                                     chi_basis )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers,      &
                                                         ndf_w2,       &
                                                         ndf_chi,      &
                                                         undf_w2,      &
                                                         undf_chi,     &
                                                         rim_width_ns, &
                                                         rim_width_ew
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def),                        intent(in) :: outer_s,    &
                                                         outer_n,    &
                                                         outer_e,    &
                                                         outer_w,    &
                                                         rim_s,      &
                                                         rim_n,      &
                                                         rim_e,      &
                                                         rim_w,      &
                                                         blend_s,    &
                                                         blend_n,    &
                                                         blend_e,    &
                                                         blend_w
  real(kind=r_def), dimension(undf_w2),          intent(inout) :: wind_mask
  real(kind=r_def), dimension(1,ndf_chi,ndf_w2), intent(in)    :: chi_basis
  real(kind=r_def), dimension(undf_chi),         intent(in)    :: chi_1, &
                                                                  chi_2, &
                                                                  chi_3

  ! Internal variables
  integer(kind=i_def)                    :: k, df1, df2
  real(kind=r_def), dimension(ndf_chi)   :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: x(3)
  real(kind=r_def), parameter            :: tol = 1.0e-15_r_def
  real(kind=r_def)                       :: query_value_ns, query_value_ew
  real(kind=r_def)                       :: mask_value
  real(kind=r_def)                       :: weight_n, weight_s, weight_e, &
                                            weight_w, weight
  logical(kind=l_def)                    :: ew_innerbox, ns_innerbox
  logical(kind=l_def)                    :: ew_outerbox, ns_outerbox
  logical(kind=l_def)                    :: ew_blendbox, ns_blendbox

  ! Define for the bottom vertical layer and copy to other layers
  k = 0

  ! Get the chi values in the kth layer
  ! chi_e (u) = chi ( f(u) + k )
  do df2 = 1, ndf_chi
    chi_1_e(df2) = chi_1( map_chi(df2) + k)
    chi_2_e(df2) = chi_2( map_chi(df2) + k)
    chi_3_e(df2) = chi_3( map_chi(df2) + k)
  enddo

  ! Get the x values (x,y,z) for each cell w2 dof, by summing
  ! over W0 to W2 basis functions
  ! x(v) = sum_u chi_e(u) * g(u,v)
  do df1 = 1, ndf_w2
    x(:) = 0.0_r_def
    do df2 = 1,ndf_chi
      x(1) = x(1) + chi_1_e(df2) * chi_basis(1,df2,df1)
      x(2) = x(2) + chi_2_e(df2) * chi_basis(1,df2,df1)
      x(3) = x(3) + chi_3_e(df2) * chi_basis(1,df2,df1)
    enddo

    ! Loop over vertical layers
    do k = 0, nlayers-1
     wind_mask(map_w2(df1) + k) = 0.0_r_def
    enddo

    ! Change of coordinates for spherical geometry
    if ( geometry == geometry_spherical .and. &
         coord_system == coord_system_xyz ) then
      ! In alpha beta space - and only pick the face where x(2)>0
      if ( x(2) > 0.0_r_def ) then
        query_value_ns = atan2(x(3),x(2))
        query_value_ew = atan2(x(1),x(2))
      else
        ! Set to an arbitrary value greater than pi/4
        query_value_ns = 10.0_r_def
        query_value_ew = 10.0_r_def
      endif
    else
      ! In cartesian space
      query_value_ns = x(2)
      query_value_ew = x(1)
    endif

    ! Define the blending region

    ! 0   0   0   0   0   0    0 0 0  0
    ! --------outerbox--------------| 0
    ! 1   1   1   1   1   1    1 1 1| 0
    ! 1   1   1   1   1   1    1 1 1| 0
    ! 1   1   1   1   1   1    1 1 1| 0
    ! --------blendbox-------|      |
    ! 0.8 0.8 0.8 0.8 0.8 0.8| 1 1 1| 0
    ! 0.5 0.5 0.5 0.5 0.5 0.8| 1 1 1| 0
    ! 0.1 0.1 0.1 0.1 0.5 0.8| 1 1 1| 0
    ! -innerbox-- 0.1 0.5 0.8| 1 1 1| 0
    ! 0   0   0 | 0.1 0.5 0.8| 1 1 1| 0
    ! 0   0   0 | 0.1 0.5 0.8| 1 1 1| 0

    ! To allow for Dirichlet boundaries in e.g. NS and periodic boundaries
    ! in e.g. EW, define the outerbox and inner box for the NS and EW
    ! directions separately. This is necessary for W2 space because the dofs
    ! are on the boundaries.

    ! We set the dofs on the outer box boundary to 1.

    ! Is the dof in (and on, hence >= ) the Outerbox boundary?
    if ( query_value_ns >= (outer_s - tol) .and. &
         query_value_ns <= (outer_n + tol) ) then
      ns_outerbox = .true.
    else
      ns_outerbox = .false.
    endif

    if ( query_value_ew >= (outer_w - tol).and. &
         query_value_ew <= (outer_e + tol) ) then
      ew_outerbox = .true.
    else
      ew_outerbox = .false.
    endif

    ! Is the dof in ( but not on, hence > ) the Innerbox boundary?
    if ( query_value_ns > (rim_s + tol).and. &
         query_value_ns < (rim_n - tol) ) then
      ns_innerbox = .true.
    else
      ns_innerbox = .false.
    endif
    if ( query_value_ew > (rim_w + tol).and. &
         query_value_ew < (rim_e - tol) ) then
      ew_innerbox = .true.
    else
      ew_innerbox = .false.
    endif

    ! Is the dof in ( but not on, hence > ) the Blendbox boundary?
    if ( query_value_ns > (blend_s + tol).and. &
         query_value_ns < (blend_n - tol) ) then
      ns_blendbox = .true.
    else
      ns_blendbox = .false.
    endif
    if ( query_value_ew > (blend_w + tol).and. &
         query_value_ew < (blend_e - tol) ) then
      ew_blendbox = .true.
    else
      ew_blendbox = .false.
    endif

    ! Define the mask_value depending on whether the dof is
    ! in the innerbox and the outerbox.

    ! Step 1. Set all the dofs to 0
    mask_value = 0.0_r_def

    ! Step 2. Set all the dofs outside (and on) the innerbox to 1.
    if ( rim_width_ns >= 0_i_def ) then
      if (.not. ns_innerbox) then
        mask_value = 1.0_r_def
      endif
    endif
    if ( rim_width_ew >= 0_i_def ) then
      if (.not. ew_innerbox) then
        mask_value = 1.0_r_def
      endif
    endif

    ! Step 3. Apply blending weights
    weight_n = 0.0_r_def
    weight_s = 0.0_r_def
    weight_e = 0.0_r_def
    weight_w = 0.0_r_def

    if ( rim_width_ns >= 0_i_def ) then
      if (ns_blendbox) then
        if ( (blend_n -rim_n) /= 0.0_r_def ) then
          weight_n = (query_value_ns - rim_n)/ (blend_n - rim_n)
        endif
        if ( (rim_s - blend_s) /= 0.0_r_def) then
          weight_s = (rim_s - query_value_ns)/ (rim_s - blend_s)
        endif
      endif
    endif
    if ( rim_width_ew >= 0_i_def ) then
      if (ew_blendbox) then
        if ( (blend_e -rim_e) /= 0.0_r_def ) then
          weight_e = (query_value_ew - rim_e)/ (blend_e - rim_e)
        endif
        if ( (rim_w - blend_w) /= 0.0_r_def) then
          weight_w = (rim_w - query_value_ew)/ (rim_w - blend_w)
        endif
      endif
    endif
    if (ns_blendbox .or. ew_blendbox) then
      weight = max( weight_n, weight_e, weight_s, weight_w)
      mask_value = mask_value * weight
    endif

    ! Step 4. Set all the dofs outside (and on) the blendbox to 1.
    if ( rim_width_ns >= 0_i_def ) then
      if (.not. ns_blendbox) then
        mask_value = 1.0_r_def
      endif
    endif
    if ( rim_width_ew >= 0_i_def ) then
      if (.not. ew_blendbox) then
        mask_value = 1.0_r_def
      endif
    endif

    ! Step 5. Set all the dofs outside the outerbox to 0.
    if ( rim_width_ns >= 0_i_def ) then
      if (.not. ns_outerbox) then
        mask_value = 0.0_r_def
      endif
    endif
    if ( rim_width_ew >= 0_i_def ) then
      if (.not. ew_outerbox) then
        mask_value = 0.0_r_def
      endif
    endif

    ! Loop over vertical layers
    do k = 0, nlayers-1
      wind_mask(map_w2(df1) + k) = mask_value
    enddo

  enddo

end subroutine create_w2mask_blend_code

end module create_w2mask_blend_kernel_mod
