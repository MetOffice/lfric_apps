!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the coefficients of the Helmholz operator for lowest-order
!!        quadrilateral elements with analytically eliminated theta.
!> @details The Helmholtz operator can be written as a sequence of
!!          coefficients (helm_X) multiplied by the Exner pressure increment
!!          at a location, e.g. for a tridiagonal vertical operator this would
!!          give lhs(k) = helm_u(k) * Exner(k+1) + helm_c(k) * Exner(k) + helm_d(k) * Exner(k-1).
!!          The Helmholtz operator in gungho with analytically eliminated theta, at lowest order, comprises
!!          a 7-pt stencil (tri-diagonal in the vertical + the 4 horizontal
!!          neighbours). The same interface as the non-eliminated version (9-pt
!!          stencil) is used and the unused coefficients are set to zero.
!!          This kernel computes the values of the coefficients
!!          from the constituent element matrices.
!!          For more details, see the solver section of
!!          https://code.metoffice.gov.uk/trac/lfric/wiki/GhaspSupport/Documentation
module elim_helmholtz_operator_kernel_mod

  use argument_mod,      only : arg_type,                   &
                                GH_FIELD, GH_OPERATOR,      &
                                GH_REAL, GH_READ, GH_WRITE, &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : i_def, r_solver
  use fs_continuity_mod, only : W2, W3, W2v
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: elim_helmholtz_operator_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                &
         arg_type(GH_FIELD*9,  GH_REAL, GH_WRITE, W3),                 & ! Helmholtz operator
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2, STENCIL(CROSS)), & ! hb_lumped_inv
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),                 & ! u_normalisation
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W3),             & ! div_star
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),             & ! M3_exner
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W2),             & ! Q32
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),             & ! M3^-1
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2)                  & ! W2 mask
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: elim_helmholtz_operator_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: elim_helmholtz_operator_code

contains

!> @brief Compute the coefficients of the Helmholtz operator, stored in a
!!        sequence of W3 fields: helm_X where X is the compass direction
!!        or vertical level that component applies to.
!> @param[in]     stencil_size    Number of cells in the horizontal stencil
!> @param[in]     cell_stencil    Stencil of horizontal cell indices
!> @param[in]     nlayers         Number of layers
!> @param[in,out] helm_c          Diagonal entry to Helmholtz matrix
!> @param[in,out] helm_n          North (j+1) entry to Helmholtz matrix
!> @param[in,out] helm_e          East (i+1) entry to Helmholtz matrix
!> @param[in,out] helm_s          South (j-1) entry to Helmholtz matrix
!> @param[in,out] helm_w          West (j-1) entry to Helmholtz matrix
!> @param[in,out] helm_u          Upper (k+1) entry to Helmholtz matrix
!> @param[in,out] helm_uu         2nd Upper (k+2) entry to Helmholtz matrix
!> @param[in,out] helm_d          Lower (k-1) entry to Helmholtz matrix
!> @param[in,out] helm_dd         2nd Lower (k-2) entry to Helmholtz matrix
!> @param[in]     hb_lumped_inv   Lumped inverse of the HB (mass matrix + buoyancy)
!!                                term
!> @param[in]     smap_size_w2    Number of cells in W2 stencil dofmap
!> @param[in]     smap_w2         Stencil dofmap for the W2 space
!> @param[in]     u_normalisation Normalisation used for the velocity equation
!> @param[in]     ncell_3d_1      Total number of cells for divergence matrix
!> @param[in]     div_star        Weighted transpose of the divergence operator
!> @param[in]     ncell_3d_2      Total number of cells for m3_exner_star matrix
!> @param[in]     m3_exner_star   Weighted W3 mass matrix
!> @param[in]     ncell_3d_3      Total number of cells for q32 matrix
!> @param[in]     q32             Weighted projection operator from W2 to W3
!> @param[in]     w2_mask         LAM mask for W2 space
!> @param[in]     ndf_w3          Number of degrees of freedom per cell for the pressure space
!> @param[in]     undf_w3         Unique number of degrees of freedom  for the pressure space
!> @param[in]     map_w3          Dofmap for the cell at the base of the column for the pressure space
!> @param[in]     ndf_w2          Number of degrees of freedom per cell for the velocity space
!> @param[in]     undf_w2         Unique number of degrees of freedom  for the velocity space
!> @param[in]     map_w2          Dofmap for the cell at the base of the column for the velocity space
subroutine elim_helmholtz_operator_code(stencil_size,                     &
                                        cell_stencil,                     &
                                        nlayers,                          &
                                        helm_c,                           &
                                        helm_n, helm_e, helm_s, helm_w,   &
                                        helm_u, helm_uu, helm_d, helm_dd, &
                                        hb_lumped_inv,                    &
                                        smap_size_w2, smap_w2,            &
                                        u_normalisation,                  &
                                        ncell_3d_1,                       &
                                        div_star,                         &
                                        ncell_3d_2,                       &
                                        m3_exner_star,                    &
                                        ncell_3d_3,                       &
                                        q32,                              &
                                        w2_mask,                          &
                                        ndf_w3, undf_w3, map_w3,          &
                                        ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                                  intent(in) :: nlayers
  integer(kind=i_def),                                  intent(in) :: stencil_size
  integer(kind=i_def), dimension(stencil_size),         intent(in) :: cell_stencil
  integer(kind=i_def),                                  intent(in) :: ncell_3d_1, ncell_3d_2, &
                                                                      ncell_3d_3
  integer(kind=i_def),                                  intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                                  intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                                  intent(in) :: smap_size_w2
  integer(kind=i_def), dimension(ndf_w3),               intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2),               intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w2, smap_size_w2), intent(in) :: smap_w2

  ! Fields
  real(kind=r_solver), dimension(undf_w3), intent(inout) :: helm_c,                         &
                                                            helm_n, helm_e, helm_s, helm_w, &
                                                            helm_u, helm_uu, helm_d, helm_dd
  real(kind=r_solver), dimension(undf_w2), intent(in)    :: hb_lumped_inv,   &
                                                            u_normalisation, &
                                                            w2_mask

  ! Operators
  real(kind=r_solver), dimension(ncell_3d_1, ndf_w2, ndf_w3), intent(in) :: div_star
  real(kind=r_solver), dimension(ncell_3d_2, ndf_w3, ndf_w3), intent(in) :: m3_exner_star
  real(kind=r_solver), dimension(ncell_3d_3, ndf_w3, ndf_w2), intent(in) :: q32

  ! Internal variables
  integer(kind=i_def) :: k, ik, kk, df, e, stencil_ik

  ! Integer mappings for neighbours in W2 spaces
  ! ( x, y )
  !      N
  !   |--4--|
  ! W 1     3 E
  !   |--2--|
  !      S
  !
  ! ( x, z )
  !     UU
  !   |--8--|
  !   |     |
  !   |  U  |
  !   |--6--|
  ! W 1     3 E
  !   |--5--|
  !   |  D  |
  !   |     |
  !   |--7--|
  !     DD

  ! Number of horizontal neighbours,
  ! this code is specific for quadrilateral meshes
  ! so we can hardwire this here instead of passing in.
  integer(kind=i_def), parameter :: nfaces_h = 4

  ! Number of fields in the Helmholtz operator,
  ! (this is the number of helm_X fields passed in)
  ! and is hardwired to 9.
  integer(kind=i_def), parameter :: helm_operator_size = 9

  ! The Compass points and up/down need to match the location
  ! of W2 DoF's.
  integer(kind=i_def), parameter :: centre   = 0
  integer(kind=i_def), parameter :: west     = 1
  integer(kind=i_def), parameter :: south    = 2
  integer(kind=i_def), parameter :: east     = 3
  integer(kind=i_def), parameter :: north    = 4
  integer(kind=i_def), parameter :: down     = 5
  integer(kind=i_def), parameter :: up       = 6

  ! The Compass points needed to match the location
  ! of W2 DoF's in neighbouring cells.
  integer(kind=i_def), dimension(nfaces_h) :: adjacent_face
  integer(kind=i_def)                      :: d1, d2
  integer(kind=i_def), parameter           :: ndf_w2v = 2

  real(kind=r_solver), dimension(ndf_w2, ndf_w3, 0:helm_operator_size-1) :: a_op
  real(kind=r_solver), dimension(ndf_w3, ndf_w2)                         :: b_op
  real(kind=r_solver), dimension(ndf_w3, ndf_w3)                         :: c_op

  ! The mixed system used by the approximate Schur complement with analytically
  ! eliminated theta and discretely eliminated rho is:
  ! (u     Exner)
  ! | I   -Nu*Hb^-1*div_star |
  ! | Q32  M3^-1*M3_exner    |

  ! This can be more compactly written as:
  ! | I A |
  ! | B C |

  ! With:                       : fs mapping
  ! I \equiv Identity           : W2 -> W2
  ! A \equiv -Nu*Hb^-1*div_star : W3 -> W2
  ! B \equiv  Q32               : W2 -> W3
  ! C \equiv  M3^-1*M3_exner    : W3 -> W3
  !
  ! I*u + A*p = Ru
  ! B*u + C*p = Rp
  ! u = Ru - A*p
  ! (C - B*A)*p = Rp + B*Ru
  !
  ! Given this the Helmholtz operator can be written as:
  ! [C - B*A]*Exner' = RHS
  !
  ! In the code these operators are:
  ! A = a_op
  ! B = b_op
  ! C = c_op
  ! And so we seek to rewrite this operator as the coefficients to
  ! each Exner value in the stencil.

  ! Compute direction faces in neighbouring cells that
  ! adjoin each direction of the central cell, i.e
  !
  !          |-----|
  !          |  N  |
  !          |E   W|
  !          |  S  |
  !    |-----|-----|-----|
  !    |  E  |  N  |  N  |
  !    |S   N|E   W|E   W|
  !    |  W  |  S  |  S  |
  !    |-----|-----|-----|
  !          |  E  |
  !          |N   S|
  !          |  W  |
  !          |-----|
  !
  ! Then adjacent face would be
  ! (/ N, E, E, S /)
  ! Since this kernel is only valid for lowest order we can use
  ! the W2 dofmap (one DoF per face) to compute this.
  adjacent_face(:) = -1
  do d1 = 1, nfaces_h
    do d2 = 1, nfaces_h
      if ( smap_w2(d2,1+d1) == smap_w2(d1,1) ) adjacent_face(d1) = d2
    end do
  end do

  do k = 0, nlayers - 1
    ik = 1 + k + (cell_stencil(1)-1)*nlayers

    ! Compute A for all cells in the stencil.
    do df = 1, ndf_w2
      ! A maps from W3 points to W2 points,
      ! first index is the location of the local edge of the cell,
      ! last index is the location of the cell in the stencil.
      ! Horizontal stencil:
      !
      !                |--------------|
      !                |   A(N,1,N)   |
      !                |              |
      !                |   A(S,1,N)   |
      ! |--------------|--------------|--------------|
      ! |A(W,1,W)      |              |A(W,1,E)      |
      ! |              |              |              |
      ! |      A(E,1,W)|              |      A(E,1,E)|
      ! |--------------|--------------|--------------|
      !                |   A(S,1,N)   |
      ! y              |              |
      ! |              |   A(S,1,S)   |
      ! 0--x           |--------------|

      do e = 1,stencil_size
        stencil_ik = 1 + k + (cell_stencil(e)-1)*nlayers
        a_op(df,:,e-1) = -u_normalisation(smap_w2(df,e)+k) &
                         *w2_mask(smap_w2(df,e)+k)         &
                         *hb_lumped_inv(smap_w2(df,e)+k)   &
                         *div_star(stencil_ik,df,:)
      end do

      ! Vertical stencil:
      !
      !     |--------------|
      !     |   A(U,1,U)   |
      ! k+1 |              |
      !     |   A(D,1,U)   |
      !     |--------------|
      !     |   A(U,1,C)   |
      ! k   |              |
      !     |   A(D,1,C)   |
      !     |--------------|
      !     |   A(U,1,D)   |
      ! k-1 |              |
      !     |   A(D,1,D)   |
      !     |--------------|
      kk = -1
      if ( k > 0 ) then
        a_op(df,:,down) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                          *hb_lumped_inv(map_w2(df)+k+kk)*div_star(ik+kk,df,:)
      else
        a_op(df,:,down) = 0.0_r_solver
      end if
      kk = 1
      if ( k < nlayers-1 ) then
        a_op(df,:,up) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                        *hb_lumped_inv(map_w2(df)+k+kk)*div_star(ik+kk,df,:)
      else
        a_op(df,:,up) = 0.0_r_solver
      end if
    end do
    ! Apply boundary conditions to A, these are terms that would multiply DoFs on
    ! the boundaries.
    if ( k == 0 )           a_op(down,:,centre:north) = 0.0_r_solver
    if ( k == 1 )           a_op(down,:,down)         = 0.0_r_solver
    if ( k == nlayers - 2 ) a_op(up,:,up)             = 0.0_r_solver
    if ( k == nlayers - 1 ) a_op(up,:,centre:north)   = 0.0_r_solver

    ! Compute B for all cells in the stencil,
    ! B maps from W2 points to W3 points and we only need it for the central
    ! cell.
    b_op = q32(ik,:,:)

    ! Compute C for all cells in the stencil,
    ! C maps from W3 points to W3 points and we only need it for the central
    ! cell.
    c_op = m3_exner_star(ik,:,:)

    ! Now compute the coefficients:

    ! -B*A,
    ! helm_e is the term that maps from Exner in the east cell to the central cell,
    ! therefore it is the combination of A for the east cell that maps to the
    ! common face of the central cell and the east cell (given by adjacent_face(east)):
    ! a(adjacent_face(east),1,east) multiplied the component of B that maps a
    ! W2 DoF on the east face to the W3 point for this cell: b(1,east).
    ! Most other terms are computed similarly.
    ! helm_c is all terms that map from the central cell to the W2 points
    ! of this cell and then back again and hence has six entries.
    helm_e(map_w3(1)+k)  = -b_op(1,east)*a_op(adjacent_face(east),1,east)
    helm_w(map_w3(1)+k)  = -b_op(1,west) *a_op(adjacent_face(west),1,west)
    helm_n(map_w3(1)+k)  = -b_op(1,north)*a_op(adjacent_face(north),1,north)
    helm_s(map_w3(1)+k)  = -b_op(1,south)*a_op(adjacent_face(south),1,south)
    helm_u(map_w3(1)+k)  = -b_op(1,up)   *a_op(down,1,up)
    helm_d(map_w3(1)+k)  = -b_op(1,down) *a_op(up,1,down)
    helm_uu(map_w3(1)+k) = 0.0_r_solver
    helm_dd(map_w3(1)+k) = 0.0_r_solver
    helm_c(map_w3(1)+k)  = - b_op(1,east) *a_op(east,1,centre)  &
                           - b_op(1,west) *a_op(west,1,centre)  &
                           - b_op(1,north)*a_op(north,1,centre) &
                           - b_op(1,south)*a_op(south,1,centre) &
                           - b_op(1,up)   *a_op(up,1,centre)    &
                           - b_op(1,down) *a_op(down,1,centre)
    ! C
    ! C maps from W3 to W3 and simply maps from the current cell to itself.
    helm_c(map_w3(1)+k) = helm_c(map_w3(1)+k) + c_op(1,1)

  end do

end subroutine elim_helmholtz_operator_code

end module elim_helmholtz_operator_kernel_mod
