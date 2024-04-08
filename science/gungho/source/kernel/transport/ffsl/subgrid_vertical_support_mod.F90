!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing vertical Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_vertical_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN

implicit none

private

! Edge reconstructions
public :: second_order_vertical_edge
public :: second_order_vertical_gradient
public :: fourth_order_vertical_edge
! Field reconstructions
public :: vertical_nirvana_recon
public :: vertical_ppm_recon
! Monotonic limiters
public :: fourth_order_vertical_mono
public :: vertical_nirvana_mono_strict
public :: vertical_nirvana_mono_relax
public :: vertical_nirvana_positive
public :: vertical_ppm_mono_strict
public :: vertical_ppm_mono_relax
public :: vertical_ppm_positive

contains

  ! ========================================================================== !
  ! EDGE RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a second-order interpolation.
  !> @details Uses a second-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge rho value.
  !!
  !> @param[in]   rho        Density values of two cells which have the ordering
  !!                         | 1 | 2 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 |
  !!                         with edges  0   1   2
  !> @param[out]  edge_value The interpolated edge value at edge_to_do
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_edge(rho, dz, edge_to_do, edge_value)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(in)  :: rho(2)
    real(kind=r_tran),   intent(in)  :: dz(2)
    integer(kind=i_def), intent(in)  :: edge_to_do
    real(kind=r_tran),   intent(out) :: edge_value

    ! Internal Variables
    real(kind=r_tran) :: z(0:2), edge_height
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get edge height to interpolate rho
    edge_height = z(edge_to_do)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*rho(1)
    cmass(2) = cmass(1) + dz(2)*rho(2)

    ! Calculate derivative of the quadratic at z = edge_height
    edge_value =   ( 2.0_r_tran*edge_height - z(2) ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                 + ( 2.0_r_tran*edge_height - z(1) ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge gradient, taking into account the height
  !!        between layers, using a second-order method.
  !> @details Uses a second-order method to find the vertical cell edge
  !!          gradient of rho. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated twice
  !!          to give the cell edge gradient.
  !!
  !> @param[in]   rho           Density values of two cells which have the ordering
  !!                            | 1 | 2 |
  !> @param[in]   dz            Height of each layer, with index the same as rho
  !> @param[out]  edge_gradient The gradient at the edge
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_gradient(rho, dz, edge_gradient)

    implicit none

    ! Arguments
    real(kind=r_tran), intent(in)  :: rho(1:2)
    real(kind=r_tran), intent(in)  :: dz(1:2)
    real(kind=r_tran), intent(out) :: edge_gradient

    ! Internal Variables
    real(kind=r_tran) :: z(0:2)
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*rho(1)
    cmass(2) = cmass(1) + dz(2)*rho(2)

    ! Calculate second derivative of the quadratic
    edge_gradient =   ( 2.0_r_tran ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                    + ( 2.0_r_tran ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_gradient

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!
  !> @param[in]   rho        Density values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below The edge value located below layer k
  !!                         (layer k corresponds to cell 3 index above)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge(rho, dz, edge_to_do, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: rho(1:4)
    real(kind=r_tran),    intent(in)    :: dz(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: z(0:4), dzs(1:4), dzsum, edge_height
    real(kind=r_tran) :: dmass(1:4)
    real(kind=r_tran) :: cmass(0:4)
    real(kind=r_tran) :: poly_mass(1:4)
    real(kind=r_tran) :: dl_dz(1:4)

    integer(kind=i_def) :: i

    ! Get scaling value
    dzsum = sum(dz)

    ! Get scaled dz
    dzs = dz / dzsum

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 4
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate rho to
    edge_height = z(edge_to_do)

    ! Get mass scaled by height
    dmass = rho * dzs

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    do i = 1, 4
      cmass(i) = cmass(i-1) + dmass(i)
    end do

    ! Get cumulative mass divided by denominator of polynomial
    poly_mass(1) = cmass(1)/((z(1))*(z(1)-z(2))*(z(1)-z(3))*(z(1)-z(4)))
    poly_mass(2) = cmass(2)/((z(2))*(z(2)-z(1))*(z(2)-z(3))*(z(2)-z(4)))
    poly_mass(3) = cmass(3)/((z(3))*(z(3)-z(1))*(z(3)-z(2))*(z(3)-z(4)))
    poly_mass(4) = cmass(4)/((z(4))*(z(4)-z(1))*(z(4)-z(2))*(z(4)-z(3)))

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz    = 4.0_r_tran*edge_height**3
    dl_dz(1) = dl_dz(1) - 3.0_r_tran*(z(2)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(2)*z(3) + z(2)*z(4))*edge_height - z(2)*z(3)*z(4)
    dl_dz(2) = dl_dz(2) - 3.0_r_tran*(z(1)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(1)*z(3) + z(1)*z(4))*edge_height - z(1)*z(3)*z(4)
    dl_dz(3) = dl_dz(3) - 3.0_r_tran*(z(1)+z(2)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(4) + z(1)*z(2) + z(1)*z(4))*edge_height - z(1)*z(2)*z(4)
    dl_dz(4) = dl_dz(4) - 3.0_r_tran*(z(1)+z(2)+z(3))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(3) + z(1)*z(2) + z(1)*z(3))*edge_height - z(1)*z(2)*z(3)

    ! Calculate value of edge below layer k
    edge_below = sum( poly_mass * dl_dz )

  end subroutine fourth_order_vertical_edge

  ! ========================================================================== !
  ! FIELD RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief  Returns the vertical Nirvana reconstruction. This can be used to
  !!         compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction is a third-order reconstruction at the cell's
  !!         vertical edges, and is based on a quadratic subgrid reconstruction.
  !!
  !! @param[out]  recon       The Nirvana reconstruction
  !! @param[in]   dep         The fractional departure distance for the reconstruction point.
  !!                          For dep>0 the recon lives above the cell
  !!                          For dep<0 the recon lives below the cell
  !! @param[in]   field       Field value of the cell
  !! @param[in]   grad_below  Estimate of gradient at bottom edge
  !! @param[in]   grad_above  Estimate of gradient at top edge
  !! @param[in]   dz          Height of cell
  !----------------------------------------------------------------------------
  subroutine vertical_nirvana_recon(recon,      &
                                    dep,        &
                                    field,      &
                                    grad_below, &
                                    grad_above, &
                                    dz)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: recon
    real(kind=r_tran),   intent(in)  :: dep
    real(kind=r_tran),   intent(in)  :: field
    real(kind=r_tran),   intent(in)  :: dz
    real(kind=r_tran),   intent(in)  :: grad_below
    real(kind=r_tran),   intent(in)  :: grad_above

    ! Reconstruction weights
    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! Compute reconstruction weights
    if (dep >= 0.0_r_tran) then
      cp = 2.0_r_tran*sixth - 0.5_r_tran*dep + sixth*dep**2
      cc = 1.0_r_tran
      cm = sixth - sixth*dep**2
    else
      cp = -sixth + sixth*dep**2
      cc = 1.0_r_tran
      cm = -2.0_r_tran*sixth - 0.5_r_tran*dep - sixth*dep**2
    end if

    ! Apply weights to gradients and field
    recon = cm*grad_below*dz + cc*field + cp*grad_above*dz

  end subroutine vertical_nirvana_recon

  !----------------------------------------------------------------------------
  !> @brief  Returns the vertical PPM reconstruction (also used for the reversible
  !!         Nirvana reconstruction). This can be used to compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction is a third-order reconstruction at the cell's
  !!         vertical edges, and is based on a quadratic subgrid reconstruction that
  !!         uses the field interpolated to cell edges.
  !!
  !! @param[out]  recon       The PPM reconstruction
  !! @param[in]   dep         The fractional departure distance for the reconstruction point.
  !!                          For dep>0 the recon lives above the cell
  !!                          For dep<0 the recon lives below the cell
  !! @param[in]   field       Field values in the cell
  !! @param[in]   edge_below  Estimate of edge value below the cell
  !! @param[in]   edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_recon(recon,      &
                                dep,        &
                                field,      &
                                edge_below, &
                                edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: recon
    real(kind=r_tran),   intent(in)  :: dep
    real(kind=r_tran),   intent(in)  :: field
    real(kind=r_tran),   intent(in)  :: edge_below
    real(kind=r_tran),   intent(in)  :: edge_above

    ! Reconstruction weights
    real(kind=r_tran) :: cm, cc, cp

    ! Compute reconstruction weights
    if (dep >= 0.0_r_tran) then
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = - dep + dep**2
    else
      cp = dep + dep**2
      cc = -3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = 1.0_r_tran + 2.0_r_tran*dep + dep**2
    end if

    ! Apply weights to field and field edge values
    recon = cm*edge_below + cc*field + cp*edge_above

  end subroutine vertical_ppm_recon

  ! ========================================================================== !
  ! MONOTONIC LIMITERS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief Applies monotonicity to a 4th-order edge reconstruction.
  !> @param[in]    rho        Density values of four cells which have the ordering
  !!                          | 1 | 2 | 3 | 4 |
  !> @param[in]    edge_to_do Tells routine which edge to do based on
  !!                          cells       | 1 | 2 | 3 | 4 |
  !!                          with edges  0   1   2   3   4
  !> @param[inout] edge_below The edge value located below layer k
  !!                          (layer k corresponds to cell 3 index above)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_mono(rho,        &
                                        edge_to_do, &
                                        edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: rho(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(inout) :: edge_below

    real(kind=r_tran) :: t1, tmin, tmax

    ! Strict Monotonicity
    if ( edge_to_do > 0_i_def .AND. edge_to_do < 4_i_def) then
      t1 = ( edge_below - rho(edge_to_do) )*( rho(edge_to_do+1) - edge_below )
      if ( t1 < 0.0_r_tran ) then
        tmin = min(rho(edge_to_do+1),rho(edge_to_do))
        tmax = max(rho(edge_to_do+1),rho(edge_to_do))
        edge_below = min( tmax, max(edge_below,tmin) )
      end if
    else if ( edge_to_do == 0_i_def ) then
      edge_below = min( max( rho(2), rho(1) ), max( edge_below, min( rho(2), rho(1) ) ) )
    else if ( edge_to_do == 4_i_def ) then
      edge_below = min( max( rho(4), rho(3) ), max( edge_below, min( rho(4), rho(3) ) ) )
    end if

  end subroutine fourth_order_vertical_mono

  !----------------------------------------------------------------------------
  !> @brief  Applies a strict monotonic limiter to a Nirvana reconstruction.
  !!
  !! @param[inout] recon       The Nirvana reconstruction
  !! @param[in]    field       Field values of three cells which have the ordering
  !!                           | 1 | 2 | 3 |. Cells 1 and 3 are only used for monotonicity
  !! @param[in]    grad_below  Estimate of gradient at z = 0 of cell 2, i.e.
  !!                           at the edge between cells 1 and 2
  !! @param[in]    grad_above  Estimate of gradient at z = 1 of cell 2, i.e.
  !!                           at the edge between cells 2 and 3
  !! @param[in]    dz          Height of cell 2
  !----------------------------------------------------------------------------
  subroutine vertical_nirvana_mono_strict(recon,      &
                                          field,      &
                                          grad_below, &
                                          grad_above, &
                                          dz)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: field(1:3)
    real(kind=r_tran),   intent(in)    :: dz
    real(kind=r_tran),   intent(in)    :: grad_below
    real(kind=r_tran),   intent(in)    :: grad_above

    ! Monotonicity variables
    real(kind=r_tran) :: p0, p1, pmin0, pmin1, pmax0, pmax1, t1

    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! Strict monotonicity
    t1 = (grad_below*dz)/(grad_above*dz-grad_below*dz + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      recon = field(2)
    else
      p0 = field(2) - grad_below*dz *2.0_r_tran*sixth - grad_above*dz * sixth
      p1 = field(2) + grad_above*dz *2.0_r_tran*sixth + grad_below*dz * sixth
      pmin0 = min( field(1), field(2))
      pmax0 = max( field(1), field(2))
      pmin1 = min( field(3), field(2))
      pmax1 = max( field(3), field(2))
      if ( p0 .gt. pmax0 .OR. p0 .lt. pmin0 .OR. p1 .gt. pmax1 .OR. p1 .lt. pmin1) then
        recon = field(2)
      end if
    end if

  end subroutine vertical_nirvana_mono_strict

  !----------------------------------------------------------------------------
  !> @brief  Applies a relaxed monotonic limiter to a Nirvana reconstruction.
  !!
  !! @param[inout] recon       The Nirvana reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point.
  !!                           For dep>0 the recon lives between cells 2 and 3
  !!                           For dep<0 the recon lives between cells 1 and 2
  !! @param[in]    field       Field values of three cells which have the ordering
  !!                           | 1 | 2 | 3 |. Cells 1 and 3 are only used for monotonicity
  !! @param[in]    grad_below  Estimate of gradient at z = 0 of cell 2, i.e.
  !!                           at the edge between cells 1 and 2
  !! @param[in]    grad_above  Estimate of gradient at z = 1 of cell 2, i.e.
  !!                           at the edge between cells 2 and 3
  !! @param[in]    dz          Height of cell 2
  !----------------------------------------------------------------------------
  subroutine vertical_nirvana_mono_relax(recon,      &
                                         dep,        &
                                         field,      &
                                         grad_below, &
                                         grad_above, &
                                         dz)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field(1:3)
    real(kind=r_tran),   intent(in)    :: dz
    real(kind=r_tran),   intent(in)    :: grad_below
    real(kind=r_tran),   intent(in)    :: grad_above

    ! Internal variables
    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: edge_below, edge_above, t1, t2, t3
    real(kind=r_tran) :: p0, p1, pmin0, pmax0, pmin1, pmax1

    ! Relaxed monotonicity
    t1 = -0.5_r_tran*( field(2)-field(1) ) / &
         ( (field(1)-2.0_r_tran*field(2)+field(3))*0.5_r_tran + EPS_R_TRAN )
    ! Compute edge values of the parabolic subgrid reconstruction
    p0 = field(2) - grad_below*dz / 3.0_r_tran - grad_above*dz / 6.0_r_tran
    p1 = field(2) + grad_above*dz / 3.0_r_tran + grad_below*dz / 6.0_r_tran
    pmin0 = min( field(1), field(2))
    pmax0 = max( field(1), field(2))
    pmin1 = min( field(3), field(2))
    pmax1 = max( field(3), field(2))
    ! Check if stationary point lies within the cell, or if edge values exceed
    ! neighbouring field values
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran .OR. &
           p0 .gt. pmax0 .OR. p0 .lt. pmin0 .OR. p1 .gt. pmax1 .OR. p1 .lt. pmin1) then
      ! Linear interpolation for edge values
      edge_below = 0.5_r_tran*(field(1)+field(2))
      edge_above = 0.5_r_tran*(field(2)+field(3))
      t2 = (edge_above-field(2)) * (field(2)-edge_below)
      t3 = abs(field(2)-edge_below) - abs(edge_above-field(2))
      if (t2 < 0.0_r_tran) then
        ! Revert to constant reconstruction
        recon = field(2)
      else if (t3 < 0.0_r_tran .and. dep >= 0.0_r_tran) then
        ! Ensure subgrid reconstruction is bounded by edge values
        cp = 0.0_r_tran
        cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
        cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
        recon = cm*edge_below + cc*field(2) + cp*edge_above
      else if (t3 < 0.0_r_tran) then
        cp = 0.0_r_tran
        cc = dep**2
        cm = 1.0_r_tran - dep**2
        recon = cm*edge_below + cc*field(2) + cp*edge_above
      else if (dep >= 0.0_r_tran) then
        cp = 1.0_r_tran - dep**2
        cc = dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field(2) + cp*edge_above
      else
        cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
        cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field(2) + cp*edge_above
      end if
    end if

  end subroutine vertical_nirvana_mono_relax

  !----------------------------------------------------------------------------
  !> @brief  Applies a positive limiter to a Nirvana reconstruction.
  !!
  !! @param[inout] recon       The Nirvana reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point
  !! @param[in]    field       Field value of the cell upwind of the reconstruction
  !! @param[in]    grad_below  Estimate of gradient above cell
  !! @param[in]    grad_above  Estimate of gradient below cell
  !! @param[in]    dz          Height of cell
  !----------------------------------------------------------------------------
  subroutine vertical_nirvana_positive(recon,      &
                                       dep,        &
                                       field,      &
                                       grad_below, &
                                       grad_above, &
                                       dz)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: dz
    real(kind=r_tran),   intent(in)    :: grad_below
    real(kind=r_tran),   intent(in)    :: grad_above

    ! Internal variables
    real(kind=r_tran) :: p0, p1, t1, t2

    ! Positive definite limiter finds the stationary point of parabolic subgrid reconstruction
    t1 = (grad_below*dz)/(grad_above*dz-grad_below*dz + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If stationary point lies within the grid cell and makes the subgrid reconstruction
      ! negative, we revert to constant reconstruction
      t2 = (field - (grad_below*dz)/2.0_r_tran - (grad_above*dz - grad_below*dz)/3.0_r_tran) + &
           (grad_below*dz) * t1 + (grad_above*dz - grad_below*dz)  * t1 * t1
      if ( t2 < 0.0_r_tran ) then
        recon = field
      end if
    else
      ! If the end points of the parabolic subgrid reconstruction are negative
      ! we revert to constant reconstruction
      p0 = field - grad_below*dz / 3.0_r_tran - grad_above*dz / 6.0_r_tran
      p1 = field + grad_above*dz / 3.0_r_tran + grad_below*dz / 6.0_r_tran
      if ( p0 .lt. 0.0_r_tran .OR. p1 .lt. 0.0_r_tran) then
        recon = field
      end if
    end if

  end subroutine vertical_nirvana_positive

  !----------------------------------------------------------------------------
  !> @brief  Applies a strict limiter to a PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_mono_strict(recon,      &
                                      field,      &
                                      edge_below, &
                                      edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Monotonicity variable
    real(kind=r_tran) :: t1

    ! Strict monotonicity
    t1 = (2.0_r_tran*edge_below + edge_above - 3.0_r_tran*field) &
         / (3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If subgrid reconstruction has extrema in the cell then revert to constant reconstruction
      recon = field
    end if

  end subroutine vertical_ppm_mono_strict

  !----------------------------------------------------------------------------
  !> @brief  Applies a relaxed limiter to a vertical PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point.
  !!                           For dep>0 the recon lives above the cell
  !!                           For dep<0 the recon lives below the cell
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_mono_relax(recon,      &
                                     dep,        &
                                     field,      &
                                     edge_below, &
                                     edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Internal variables
    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: t1, t2, t3

    ! Relaxed monotonicity
    t1 = (2.0_r_tran*edge_below + edge_above - 3.0_r_tran*field) &
         / (3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If subgrid reconstruction has extrema in the cell then check smoothness of field
      t2 = (edge_above - field)*(field - edge_below)
      t3 = abs(field - edge_below) - abs(edge_above - field)
      if ( t2 < 0.0_r_tran ) then
        ! Revert to constant reconstruction
        recon = field
      else if ( t3 < 0.0_r_tran .and. dep >= 0.0_r_tran ) then
        ! Ensure subgrid reconstruction is bounded by edge values
        cp = 0.0_r_tran
        cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
        cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
        recon = cm*edge_below + cc*field + cp*edge_above
      else if ( t3 < 0.0_r_tran ) then
        cp = 0.0_r_tran
        cc = dep**2
        cm = 1.0_r_tran - dep**2
        recon = cm*edge_below + cc*field + cp*edge_above
      else if (dep >= 0.0_r_tran) then
        cp = 1.0_r_tran - dep**2
        cc = dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field + cp*edge_above
      else
        cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
        cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field + cp*edge_above
      end if
    end if

  end subroutine vertical_ppm_mono_relax

  !----------------------------------------------------------------------------
  !> @brief  Applies a positive limiter to a vertical PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point.
  !!                           For dep>0 the recon lives above the cell
  !!                           For dep<0 the recon lives below the cell
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_positive(recon,      &
                                   dep,        &
                                   field,      &
                                   edge_below, &
                                   edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Internal variables
    real(kind=r_tran) :: t1, t2, aa, bb

    ! Positive definite limiter
    aa = -4.0_r_tran*edge_below - 2.0_r_tran*edge_above + 6.0_r_tran*field
    bb = 3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field
    ! Find stationary point of parabolic subgrid reconstruction
    t1 = -0.5_r_tran*aa/(bb+EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If stationary point lies within the grid cell and makes the subgrid reconstruction
      ! negative, we revert to constant reconstruction
      t2 = edge_below + aa * t1 + bb  * t1 * t1
      if ( t2 < 0.0_r_tran ) then
        recon = field
      end if
    end if

  end subroutine vertical_ppm_positive

end module subgrid_vertical_support_mod
