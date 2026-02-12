!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Rayleigh Friction
!>
module rayleigh_friction_kernel_mod

  use argument_mod,             only: arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_READ, GH_INC,   &
                                      GH_SCALAR,         &
                                      CELL_COLUMN
  use constants_mod,            only: r_def, i_def, PI
  use fs_continuity_mod,        only: W2, Wtheta, W3
  use kernel_mod,               only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: rayleigh_friction_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                &
         arg_type(GH_FIELD,  GH_REAL, GH_INC,  W2),     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2),     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2),     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2),     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W3),     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),         &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),         &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),         &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),         &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)          &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: rayleigh_friction_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: rayleigh_friction_code

contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in,out] du u increment data
!! @param[in] u u data
!! @param[in] w2_rmultiplicity Reciprocal of multiplicity for W2
!! @param[in] exner The exner pressure
!! @param[in] exner_in_wth The exner pressure in Wtheta
!! @param[in] kappa Ratio of Rd and cp
!! @param[in] dt The model timestep length
!! @param[in] ndf_w2 The number of degrees of freedom per cell for W2
!! @param[in] undf_w2 The number of unique degrees of freedom for W2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the
!>            base of the column for W2
!! @param[in] ndf_w3 The number of degrees of freedom per cell for W3
!! @param[in] undf_w3 The number of unique degrees of freedom for W3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the
!>            base of the column for W3
!! @param[in] ndf_wth The number of degrees of freedom per cell for Wtheta
!! @param[in] undf_wth The number of unique degrees of freedom for Wtheta
!! @param[in] map_wth Integer array holding the dofmap for the cell at the
!>            base of the column for Wtheta
subroutine rayleigh_friction_code( nlayers,                                    &
                                   du, u, u_ref,                               &
                                   w2_rmultiplicity,                           &
                                   exner, exner_in_wth,                        &
                                   kappa, p_zero, p_cutoff, tau, dt,           &
                                   ndf_w2, undf_w2, map_w2,                    &
                                   ndf_w3, undf_w3, map_w3,                    &
                                   ndf_wth, undf_wth, map_wth                  &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: du
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u, u_ref, w2_rmultiplicity
  real(kind=r_def), dimension(undf_wth), intent(in)    :: exner_in_wth
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: exner
  real(kind=r_def),                      intent(in)    :: kappa, p_zero, p_cutoff, tau
  real(kind=r_def),                      intent(in)    :: dt

  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  ! Internal variables
  integer(kind=i_def) :: k, df, w2_idx
  real(kind=r_def)    :: kR(nlayers)
  real(kind=r_def)    :: exner_top, p_top
  real(kind=r_def)    :: p(nlayers), p_switch(nlayers)

  exner_top = exner_in_wth(map_wth(1)+nlayers)
  p_top = exner_top ** (1.0_r_def / kappa) * p_zero
  p(:) = exner(map_w3(1) : map_w3(1)+nlayers-1) ** (1.0_r_def / kappa) * p_zero
  ! Should only kick in when pressure is lower than p_cutoff
  p_switch(:) = 0.5_r_def + 0.5_r_def*SIGN(1.0_r_def, p_cutoff - p(:))
  kR(:) = p_switch / tau * (SIN(0.5_r_def*PI*(LOG(p_cutoff/p(:)) / LOG(p_cutoff/p_top))))**2

  do df = 1, 4
    do k = 1, nlayers
      w2_idx = map_w2(df) + k - 1
      ! Backward Euler discretisation
      du(w2_idx) = du(w2_idx) + w2_rmultiplicity(w2_idx) * (                   &
        (1.0_r_def / (1.0_r_def + kR(k) * dt) - 1.0_r_def) * u(w2_idx)         &
        + kR(k) * dt / (1.0_r_def + kR(k) * dt) * u_ref(w2_idx)                &
      )
    end do
  end do

end subroutine rayleigh_friction_code

end module rayleigh_friction_kernel_mod
