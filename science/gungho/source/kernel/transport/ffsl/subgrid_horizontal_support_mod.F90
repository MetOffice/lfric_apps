!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing horizontal Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_horizontal_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use transport_enumerated_types_mod, only: horizontal_monotone_strict,   &
                                          horizontal_monotone_relaxed,  &
                                          horizontal_monotone_positive, &
                                          horizontal_monotone_none

implicit none

private

public :: fourth_order_horizontal_edge
public :: horizontal_ppm_recon
public :: horizontal_nirvana_recon
! Routines from SPT
public :: horizontal_ppm_recon_spt_edges
public :: ppm_density_at_any_edge
public :: horizontal_nirvana_coeffs_general
public :: horizontal_nirvana_recon_spt_edges
public :: horizontal_nirvana_case
public :: horizontal_ppm_case
public :: calculate_parabola_coeffs
public :: parabola_mono

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the estimated density at the edge of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of rho. The function is passed four density values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | 4 | and returns the estimated density value between
  !!         cells 2 and 3. The cells are assumed to be uniform in spacing.
  !!         Monotonicity options are provided.
  !!
  !! @param[in]   density            Has dof map of the form | 1 | 2 | 3 | 4 |
  !! @param[in]   monotone           Monotone option to ensures no over/undershoots
  !! @return      density_at_edge    Interpolated density value at edge between
  !!                                 cells 2 and 3.
  !----------------------------------------------------------------------------
  function fourth_order_horizontal_edge(density,monotone) result(density_at_edge)

    implicit none

    real(kind=r_tran),   intent(in) :: density(1:4)
    integer(kind=i_def), intent(in) :: monotone

    real(kind=r_tran) :: density_at_edge
    real(kind=r_tran) :: t1, tmax, tmin

    ! As the cell widths are assumed to be constant the edge value reduces to that given in
    ! Colella and Woodward, JCP, 54, 1984, equation (1.9)
    density_at_edge = (7.0_r_tran/12.0_r_tran) * (density(2)+density(3)) &
                     -(1.0_r_tran/12.0_r_tran) * (density(1)+density(4))

    if ( monotone == horizontal_monotone_strict .OR. &
         monotone == horizontal_monotone_relaxed ) then
      ! Limit edge value to be bounded by neighbouring values
      t1 = ( density_at_edge - density(2) )*( density(3) - density_at_edge )
      if ( t1 < 0.0_r_tran ) then
         tmin = min(density(3),density(2))
         tmax = max(density(3),density(2))
         density_at_edge = min( tmax, max(density_at_edge,tmin) )
      end if
    else if ( monotone == horizontal_monotone_positive ) then
      ! Positivity - ensure edge value is positive
      density_at_edge = max(density_at_edge,0.0_r_tran)
    end if

  end function fourth_order_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal PPM reconstruction. This can be used to
  !!         compute the flux as:
  !!         flux = u * reconstruction
  !!         The dofmap for the field values is of the form
  !!         | 1 | 2 | 3 | 4 | 5 | where the reconstructions are being
  !!         computed for the edges of cell 3. The reconstruction is third-order,
  !!         and is based on the quadratic subgrid reconstruction of PPM.
  !!
  !! @param[out]  recon       The PPM reconstruction
  !! @param[in]   dep         The fractional departure distance for the reconstruction point.
  !!                          For dep>0 the recon lives between cells 3 and 4
  !!                          For dep<0 the recon lives between cells 2 and 3
  !! @param[in]   field       Field values in the 5 cells with ordering
  !!                          | 1 | 2 | 3 | 4 | 5 |
  !! @param[in]   monotone    Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine horizontal_ppm_recon(recon,dep,field,monotone)

    implicit none

    real(kind=r_tran),    intent(out) :: recon
    real(kind=r_tran),    intent(in)  :: dep
    real(kind=r_tran),    intent(in)  :: field(1:5)
    integer(kind=i_def),  intent(in)  :: monotone

    real(kind=r_tran) :: edge_left
    real(kind=r_tran) :: edge_right
    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: t1, t2, t3, aa, bb

    ! Get PPM edge values
    edge_left = fourth_order_horizontal_edge(field(1:4),monotone)
    edge_right = fourth_order_horizontal_edge(field(2:5),monotone)

    ! Compute reconstruction weights
    if (dep >= 0.0_r_tran) then
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = - dep + dep**2
    else
      cp = dep + dep**2.0_r_tran
      cc = -3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = 1.0_r_tran + 2.0_r_tran*dep + dep**2
    end if

    ! Apply weights to field and field edge values
    recon = cm*edge_left + cc*field(3) + cp*edge_right

    ! Apply monotonicity if needed
    if ( monotone == horizontal_monotone_strict ) then
      ! Strict monotonicity
      t1 = (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field(3)) &
           / (3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field(3) + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If subgrid reconstruction has extrema in the cell then revert to constant reconstruction
        recon = field(3)
      end if
    else if ( monotone == horizontal_monotone_relaxed ) then
      ! Relaxed monotonicity
      t1 = (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field(3)) &
           / (3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field(3) + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If subgrid reconstruction has extrema in the cell then check smoothness of field
        t2 = (edge_right - field(3))*(field(3) - edge_left)
        t3 = abs(field(3) - edge_left) - abs(edge_right - field(3))
        if ( t2 < 0.0_r_tran ) then
          ! Revert to constant reconstruction
          recon = field(3)
        else
          ! Ensure subgrid reconstruction is bounded by edge values
          if ( t3 < 0.0_r_tran ) then
            if (dep >= 0.0_r_tran) then
              cp = 0.0_r_tran
              cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
              cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
            else
              cp = 0.0_r_tran
              cc = dep**2
              cm = 1.0_r_tran - dep**2
            end if
          else
            if (dep >= 0.0_r_tran) then
              cp = 1.0_r_tran - dep**2
              cc = dep**2.0_r_tran
              cm = 0.0_r_tran
            else
              cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
              cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
              cm = 0.0_r_tran
            end if
          end if
          recon = cm*edge_left + cc*field(3) + cp*edge_right
        end if
      end if
    else if ( monotone == horizontal_monotone_positive ) then
      ! Positive definite limiter
      aa = -4.0_r_tran*edge_left - 2.0_r_tran*edge_right + 6.0_r_tran*field(3)
      bb = 3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field(3)
      ! Find stationary point of parabolic subgrid reconstruction
      t1 = -0.5_r_tran*aa / (bb + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If stationary point lies within the grid cell and makes the subgrid reconstruction
        ! negative, we revert to constant reconstruction
        t2 = edge_left + aa * t1 + bb * t1 * t1
        if ( t2 .le. 0.0_r_tran ) then
          recon = field(3)
        end if
      end if
    end if

  end subroutine horizontal_ppm_recon

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal Nirvana reconstruction at a cell edge.
  !!         This can be used to compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction is third-order, and is based on a quadratic
  !!         subgrid reconstruction
  !!
  !! @param[out]  recon     The Nirvana reconstruction
  !! @param[in]   dep       The fractional departure distance for the reconstruction point.
  !!                        For dep>0 the recon lives between cells 2 and 3
  !!                        For dep<0 the recon lives between cells 1 and 2
  !! @param[in]   field     Field values of three cells which have the ordering
  !!                        | 1 | 2 | 3 |
  !! @param[in]   monotone  Monotone option to ensures no over/undershoots
  subroutine horizontal_nirvana_recon(recon, dep, field, monotone)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: recon
    real(kind=r_tran),   intent(in)  :: dep
    real(kind=r_tran),   intent(in)  :: field(1:3)
    integer(kind=i_def), intent(in)  :: monotone

    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: p0, p1, pmin0, pmin1, pmax0, pmax1
    real(kind=r_tran) :: q0, q1, t1, t2, t3
    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! Compute reconstruction weights based on sign of fractional departure distance
    if (dep >= 0.0_r_tran) then
      cp = 2.0_r_tran*sixth - 0.5_r_tran*dep + sixth*dep**2
      cc = 5.0_r_tran*sixth + 0.5_r_tran*dep - 2.0_r_tran*sixth*dep**2
      cm = -sixth + sixth*dep**2
    else
      cp = -sixth + sixth*dep**2
      cc = 5.0_r_tran*sixth - 0.5_r_tran*dep - 2.0_r_tran*sixth*dep**2
      cm = 2.0_r_tran*sixth + 0.5_r_tran*dep + sixth*dep**2
    end if

    ! Apply reconstruction weights to the field
    recon = cm*field(1) + cc*field(2) + cp*field(3)

    ! Apply monotonicity if needed
    if ( monotone == horizontal_monotone_strict ) then
      ! Strict monotonicity
      ! Find stationary point of the parabolic subgrid reconstruction
      t1 = -0.5_r_tran*( field(2)-field(1) ) / &
           ( (field(1)-2.0_r_tran*field(2)+field(3))/2.0_r_tran + EPS_R_TRAN )
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If subgrid reconstruction has extrema in the cell then revert to constant reconstruction
        recon = field(2)
      else
        ! Compute edge values of the parabolic subgrid reconstruction
        p0 = (-field(3) + 5.0_r_tran * field(2) + 2.0_r_tran * field(1)) * sixth
        p1 = (2.0_r_tran*field(3) + 5.0_r_tran * field(2) -  field(1)) * sixth
        pmin0 = min( field(1), field(2))
        pmax0 = max( field(1), field(2))
        pmin1 = min( field(3), field(2))
        pmax1 = max( field(3), field(2))
        if ( p0 .gt. pmax0 .OR. p0 .lt. pmin0 .OR. p1 .gt. pmax1 .OR. p1 .lt. pmin1) then
          ! If edge values exceed neighbouring field values then revert to constant reconstruction
          recon = field(2)
        end if
      end if
    else if ( monotone == horizontal_monotone_relaxed ) then
      ! Relaxed monotonicity
      ! Find stationary point of the parabolic subgrid reconstruction
      t1 = -0.5_r_tran*( field(2)-field(1) ) / &
           ( (field(1)-2.0_r_tran*field(2)+field(3))*0.5_r_tran + EPS_R_TRAN )
      ! Compute edge values of the parabolic subgrid reconstruction
      p0 = (-field(3) + 5.0_r_tran * field(2) + 2.0_r_tran * field(1)) / 6.0_r_tran
      p1 = (2.0_r_tran*field(3) + 5.0_r_tran * field(2) -  field(1)) / 6.0_r_tran
      pmin0 = min( field(1), field(2))
      pmax0 = max( field(1), field(2))
      pmin1 = min( field(3), field(2))
      pmax1 = max( field(3), field(2))
      ! Check if stationary point lies within the cell, or if edge values exceed
      ! neighbouring field values
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran .OR. &
             p0 .ge. pmax0 .OR. p0 .le. pmin0 .OR. p1 .ge. pmax1 .OR. p1 .le. pmin1) then
        ! Check smoothness of field
        q0 = 0.5_r_tran*(field(1)+field(2))
        q1 = 0.5_r_tran*(field(2)+field(3))
        t2 = (q1-field(2)) * (field(2)-q0)
        t3 = abs(field(2)-q0) - abs(q1-field(2))
        if (t2 < 0.0_r_tran) then
          ! Revert to constant reconstruction
          recon = field(2)
        else
          ! Ensure subgrid reconstruction is bounded by edge values
          if (t3 < 0.0_r_tran) then
            if (dep >= 0.0_r_tran) then
              cp = 0.0_r_tran
              cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
              cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
            else
              cp = 0.0_r_tran
              cc = dep**2
              cm = 1.0_r_tran - dep**2
            end if
          else
            if (dep >= 0.0_r_tran) then
              cp = 1.0_r_tran -  dep**2
              cc = dep**2
              cm = 0.0_r_tran
            else
              cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
              cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
              cm = 0.0_r_tran
            end if
          end if
          recon = cm*q0 + cc*field(2) + cp*q1
        end if
      end if
    else if ( monotone == horizontal_monotone_positive ) then
      ! Positive definite limiter finds the stationary point of parabolic subgrid reconstruction
      t1 = -0.5_r_tran*( field(2)-field(1) ) / &
           ( (field(1)-2.0_r_tran*field(2)+field(3))/2.0_r_tran + EPS_R_TRAN )
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If stationary point lies within the grid cell and makes the subgrid reconstruction
        ! negative, we revert to constant reconstruction
        t2 = (-field(3) + 5.0_r_tran*field(2) + 2.0_r_tran*field(1))/6.0_r_tran + &
             (field(2) - field(1)) * t1 + 0.5_r_tran*(field(3) - 2.0_r_tran*field(2) + field(1))  * t1 * t1
        if ( t2 < 0.0_r_tran ) then
          recon = field(2)
        end if
      else
        ! If the end points of the parabolic subgrid reconstruction are negative
        ! we revert to constant reconstruction
        p0 = (-field(3) + 5.0_r_tran * field(2) + 2.0_r_tran * field(1)) / 6.0_r_tran
        p1 = (2.0_r_tran*field(3) + 5.0_r_tran * field(2) -  field(1)) / 6.0_r_tran
        if ( p0 .le. 0.0_r_tran .OR. p1 .le. 0.0_r_tran ) then
          recon = field(2)
        end if
      end if
    end if

  end subroutine horizontal_nirvana_recon

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal PPM reconstruction. This is similar to the
  !!         subroutine "horizontal_ppm_recon" with special treatment of panel edges
  !!         and larger stencil around the cell (7 Cells instead of 5)
  !!
  !! @param[out]  recon       The PPM reconstruction
  !! @param[in]   dep         The fractional departure distance for the reconstruction point.
  !! @param[in]   field       Field values in the 7 cells with ordering
  !!                          | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  !! @param[in]   ipanel      Panel IDs for the 7 cells
  !! @param[in]   monotone    Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine horizontal_ppm_recon_spt_edges(recon,dep,field,ipanel,monotone)

    implicit none

    real(kind=r_tran),    intent(out) :: recon
    real(kind=r_tran),    intent(in)  :: dep
    real(kind=r_tran),    intent(in)  :: field(1:7)
    integer(kind=i_def),  intent(in)  :: ipanel(1:7)
    integer(kind=i_def),  intent(in)  :: monotone
    !Local variables
    real(kind=r_tran)   :: density_left, density_right, s
    real(kind=r_tran)   :: coeffs(1:3)
    integer(kind=i_def) :: rec_case, cell
    real(kind=r_tran), parameter  :: half = 0.5_r_tran
    real(kind=r_tran), parameter  :: third = 1.0_r_tran/3.0_r_tran
    real(kind=r_tran), parameter  :: one = 1.0_r_tran
    real(kind=r_tran), parameter  :: two = 2.0_r_tran

    ! Get reconstruction case

    call horizontal_ppm_case(ipanel,rec_case)

    ! Get PPM edge values depending on which case

    select case (rec_case)
    case (1) !all points on same panel
      density_left  = ppm_density_at_any_edge(field(2:5),3_i_def)
      density_right = ppm_density_at_any_edge(field(3:6),3_i_def)
    case (2) !2-3-4-5 on same panel and 6 on different panel (edge 5/6)
      density_left  = ppm_density_at_any_edge(field(2:5),3_i_def)
      density_right = ppm_density_at_any_edge(field(2:5),4_i_def)
    case (3) !1-2-3-4 on same panel and 5-6 on different panel (edge 4/5)
      density_left  = ppm_density_at_any_edge(field(1:4),4_i_def)
      density_right = ppm_density_at_any_edge(field(1:4),5_i_def)
    case (4) !3-4-5-6 on same panel and 2 on different panel (edge 2/3)
      density_left  = ppm_density_at_any_edge(field(3:6),2_i_def)
      density_right = ppm_density_at_any_edge(field(3:6),3_i_def)
    case (5) !4-5-6-7 on same panel and 2-3 on different panel (edge 3/4)
      density_left  = ppm_density_at_any_edge(field(4:7),1_i_def)
      density_right = ppm_density_at_any_edge(field(4:7),2_i_def)
    end select

    ! Compute parabola coeffs
    cell = 4_i_def
    call calculate_parabola_coeffs(density_left,density_right,field(cell),coeffs)

    !If monotone correct coeffs
    if (monotone /= horizontal_monotone_none ) then
       call parabola_mono(coeffs,field(3),field(4),field(5),monotone)
    end if

    !Compute reconstruction from dep and coeffs
    s = abs(dep)
    if ( dep >= 0.0_r_tran) then
      recon = ( third*(s**two) - s + one )*coeffs(3) +  &
                (one - half*s)*coeffs(2) + coeffs(1)
    else
      recon = third*(s**two)*coeffs(3) + half*s*coeffs(2) + coeffs(1)
    end if

  end subroutine horizontal_ppm_recon_spt_edges

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal PPM reconstruction for any cell edge (1:5) from
  !!         the 4 cells data density(1:4)
  !! @param[in]  density         Field values for 4 cells | 1 | 2 | 3 | 4 |
  !! @param[in]  edge            The edge=(1:5) where the density is required
  !! @param[out] density_at_edge Density at the edge specified
  !----------------------------------------------------------------------------
  function ppm_density_at_any_edge(density,edge) result(density_at_edge)

    implicit none

    real(kind=r_tran),   intent(in) :: density(1:4)
    integer(kind=i_def), intent(in) :: edge

    real(kind=r_tran)   :: density_at_edge
    integer(kind=i_def) :: i
    real(kind=r_tran), parameter, dimension(4,5) :: ppm_w =  reshape(  &
            (/25.0_r_tran, -23.0_r_tran, 13.0_r_tran, -3.0_r_tran,     &
               3.0_r_tran,  13.0_r_tran, -5.0_r_tran,  1.0_r_tran,     &
              -1.0_r_tran,   7.0_r_tran,  7.0_r_tran, -1.0_r_tran,     &
               1.0_r_tran,  -5.0_r_tran, 13.0_r_tran,  3.0_r_tran,     &
              -3.0_r_tran,  13.0_r_tran,-23.0_r_tran, 25.0_r_tran/)/12.0_r_tran, shape(ppm_w))

    density_at_edge = 0.0_r_tran
    do i = 1,4
      density_at_edge = density_at_edge + ppm_w(i,edge)*density(i)
    end do

  end function ppm_density_at_any_edge

  !----------------------------------------------------------------------------------------
  !> @brief Identify which case for ppm reconstruction
  !! This subroutine identify if it leaves the reconstruction centred around
  !! cell 4 (default=case 1) away from panel edges, or shift the stencil left
  !! (by either 1 or 2 cells) or shift right (by either 1 or 2 cells)
  !! case 1 = all points (2:6) on same panel
  !! case 2 = points (2:5) on same panel and 6 on different one   (shift left by 1 cell)
  !! case 3 = points (1:4) on same panel and 5/6 on different one (shift left by 2 cells)
  !! case 4 = points (3:6) on same panel and 2 on different one (shift right by 1 cell)
  !! case 5 = points (4:7) on same panel and 1/2 on different one (shift right by 2 cells)
  !! @param[in]   ipanel  Panel IDs
  !! @param[out]  case    Identifies which cases 1:5
  !----------------------------------------------------------------------------------------
  subroutine horizontal_ppm_case(ipanel,case)

    implicit none

    integer(kind=i_def),  intent(in)   :: ipanel(1:7)
    integer(kind=i_def),  intent(out)  :: case

    integer(kind=i_def) :: test(1:4)
    integer(kind=i_def) :: test1, test2

    test(1) = ipanel(4) - ipanel(3)
    test(2) = ipanel(4) - ipanel(2)
    test(3) = ipanel(4) - ipanel(5)
    test(4) = ipanel(4) - ipanel(6)
    test1 = test(1)+test(2)+test(3)+test(4)

    if ( test1 == 0_i_def ) then !all points on same panel
      case =  1_i_def
    else !Not all points on same panel
      test2 = test(1)+test(2)
      if (test2 == 0_i_def ) then ! edge on right (5 or 6 on different panel)
          if ( test(3) == 0_i_def ) then !edge between 5/6
             case =  2_i_def
          else  !edge between 4/5
             case =  3_i_def
          end if
      else ! 2 or 3 on different panel / edge on left
         if ( test(1) == 0_i_def ) then  !edge between 2/3
            case =  4_i_def
         else  !edge between 3/4
            case =  5_i_def
         end if
      end if
    end if
  end subroutine horizontal_ppm_case

  !----------------------------------------------------------------------------
  !> @brief Identify which case for nirvana reconstruction
  !! This subroutine identify if it leaves the reconstruction centred around
  !! cell 3 (default=case 1) away from panel edges, or shift the stencil left
  !! (case 2) or shift right (case 3)
  !! case 1 = all points (2:4) on same panel
  !! case 2 = points (1:3) on same panel and 4/5 on different one (shift left)
  !! case 3 = points (3:5) on same panel and 1/2 on different one (shift right)
  !! @param[in]   ipanel  Panel IDs
  !! @param[out]  case    Identifies which cases 1:3
  !----------------------------------------------------------------------------
  subroutine horizontal_nirvana_case(ipanel,case)

    implicit none

    integer(kind=i_def),  intent(in)   :: ipanel(1:5)
    integer(kind=i_def),  intent(out)  :: case
    integer(kind=i_def) :: test1, test2, test3

    test1 = ipanel(3) - ipanel(2)
    test2 = ipanel(3) - ipanel(4)
    test3 = test1+test2

    if ( test3 == 0_i_def ) then !the 3 poins (2:4) on same panel
      case =  1_i_def
    else ! Not all 3 points on same panel
      if ( test1 == 0_i_def ) then ! point 4 is on different panel
        case =  2_i_def
      else ! point 2 is on different panel
        case =  3_i_def
      end if
    end if
  end subroutine horizontal_nirvana_case

!----------------------------------------------------------------------------
!> @brief  Returns the horizontal Nirvana reconstruction parabola coefficients
!!         for any cell(i),i=1,3 from the array density(1:3).
!! @param[in]  rho        Field values for the 3 cells | 1 | 2 | 3 |
!! @param[in]  cell       The cell=(1:3) where the coefficients are required
!! @param[out] coeffs     Parabola coefficients for the required cell
!----------------------------------------------------------------------------
  subroutine horizontal_nirvana_coeffs_general(coeffs,rho,cell)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: coeffs(1:3)
    real(kind=r_tran),   intent(in)  :: rho(1:3)
    integer(kind=i_def), intent(in)  :: cell

    ! Internal variables
    real(kind=r_tran), parameter, dimension(3,3) :: c1_w =  reshape(     &
                               (/11.0_r_tran, -7.0_r_tran,  2.0_r_tran,  &
                                  2.0_r_tran,  5.0_r_tran, -1.0_r_tran,  &
                                -1.0_r_tran,  5.0_r_tran,  2.0_r_tran/)/6.0_r_tran, shape(c1_w))
    real(kind=r_tran), parameter, dimension(3,3) :: c2_w =  reshape(    &
                               (/-2.0_r_tran,  3.0_r_tran, -1.0_r_tran, &
                                 -1.0_r_tran,  1.0_r_tran,  0.0_r_tran, &
                                  0.0_r_tran, -1.0_r_tran,  1.0_r_tran/), shape(c2_w))
    real(kind=r_tran), parameter, dimension(3) :: c3_w = (/0.5_r_tran, -1.0_r_tran, 0.5_r_tran/)
    integer(kind=i_def) :: i

    ! Initialise coefficients to be zero
    coeffs(:) = 0.0_r_tran

    do i = 1, 3
      coeffs(1) = coeffs(1) + c1_w(i,cell)*rho(i)
      coeffs(2) = coeffs(2) + c2_w(i,cell)*rho(i)
      coeffs(3) = coeffs(3) + c3_w(i)*rho(i)
    end do

  end subroutine horizontal_nirvana_coeffs_general

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal Nirvana reconstruction at a cell edge with
  !!         a special treatment of cells adjacent to panel edges.
  !! @param[out]  recon     The Nirvana reconstruction
  !! @param[in]   dep       The fractional departure distance for the reconstruction point.
  !! @param[in]   field     Field values of three cells which have the ordering
  !!                        | 1 | 2 | 3 | 4 | 5 |
  !! @param[in]   ipanel    Panel IDs for the 5 cells
  !! @param[in]   monotone  Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine horizontal_nirvana_recon_spt_edges(recon, dep, field, ipanel, monotone)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: recon
    real(kind=r_tran),   intent(in)  :: dep
    real(kind=r_tran),   intent(in)  :: field(1:5)
    integer(kind=i_def), intent(in)  :: ipanel(1:5)
    integer(kind=i_def), intent(in)  :: monotone

    ! Internal variables
    real(kind=r_tran)   :: coeffs(1:3)
    integer(kind=i_def) :: rec_case
    real(kind=r_tran)   :: s
    real(kind=r_tran), parameter  :: half = 0.5_r_tran
    real(kind=r_tran), parameter  :: third = 1.0_r_tran/3.0_r_tran
    real(kind=r_tran), parameter  :: one = 1.0_r_tran
    real(kind=r_tran), parameter  :: two = 2.0_r_tran

    ! Identify which reconstruction case (rec_case)

    call horizontal_nirvana_case(ipanel,rec_case)

    ! Compute parabola coefficients

    select case (rec_case)
    case (1) !points (2:4) on same panel
      call horizontal_nirvana_coeffs_general(coeffs,field(2:4),2_i_def)
    case (2) !point (1:3) on same panel
      call horizontal_nirvana_coeffs_general(coeffs,field(1:3),3_i_def)
    case (3) !point (3:5) on same panel
      call horizontal_nirvana_coeffs_general(coeffs,field(3:5),1_i_def)
    end select

    !If monotone correct coeffs
    if (monotone /= horizontal_monotone_none ) then
       call parabola_mono(coeffs,field(2),field(3),field(4),monotone)
    end if

    !Compute reconstruction from dep and coeffs
    s = abs(dep)
    if ( dep >= 0.0_r_tran) then
      recon = ( third*(s**two) - s + one )*coeffs(3) +  &
                (one - half*s)*coeffs(2) + coeffs(1)
    else
       recon = third*(s**two)*coeffs(3) + half*s*coeffs(2) + coeffs(1)
    end if

  end subroutine horizontal_nirvana_recon_spt_edges

  !----------------------------------------------------------------------------------
  !> @brief check if the parabola is monotone otherwise change it to a monotone one
  !! @param[inout]  coeffs          The cell parabola coefficients
  !! @param[in]     density_left    Density at the left of the cell
  !! @param[in]     density         Density of the cell
  !! @param[in]     density_right   Density at the right of the cell
  !! @param[in]     monotone        Monotone option
  !----------------------------------------------------------------------------------
  subroutine parabola_mono(coeffs,density_left,density,density_right,monotone)

    implicit none

    real(kind=r_tran),   intent(inout) :: coeffs(1:3)
    real(kind=r_tran),   intent(in)    :: density
    real(kind=r_tran),   intent(in)    :: density_left, density_right
    integer(kind=i_def), intent(in)    :: monotone
    real(kind=r_tran) :: t, test1, test2, test3, test4
    real(kind=r_tran) :: val_left, val_right, val_mid

    t = - coeffs(2)/(2.0_r_tran*coeffs(3) + EPS_R_TRAN)
    test1 = (t+EPS_R_TRAN)*((1.0_r_tran + EPS_R_TRAN)-t)
    val_left  = coeffs(1)
    val_right = coeffs(1) + coeffs(2) + coeffs(3)

    ! If test1 > 0 then there is a turning point inside the interval [0,1],
    ! Therefore, we need to modify the parabola to a monotone reconstruction
    ! We have 2 options:
    ! If (monotone=horizontal_monotone_strict) we reduce the reconstruction to constant.
    !     With this option we also reduce to linear if the end-values (val_left,val_right)
    !     of the parabola are outside the immidiate neighbours
    ! If (monotone/=horizontal_monotone_strict) then we reduces to a monotone parbola

    if ( monotone == horizontal_monotone_strict ) then
      test2 = (val_left - density_left)*(density - val_left)
      test3 = (val_right - density)*(density_right - val_right)
      if ( (test1 > 0.0_r_tran) .or.  &
           (test2 < 0.0_r_tran) .or.  &
           (test3 < 0.0_r_tran)         ) then
         coeffs(1) = density
         coeffs(2) = 0.0_r_tran
         coeffs(3) = 0.0_r_tran
      end if
    else
      if ( test1 > 0.0_r_tran ) then
         val_mid   = 0.5_r_tran*(val_left + val_right)
         test4 = (density - val_left)*(val_mid - density)
         if ( test4 > 0.0_r_tran ) then !Left flat monotone parabola
           coeffs(1) = val_left
           coeffs(2) = 0.0_r_tran
           coeffs(3) = 3.0_r_tran*density -  3.0_r_tran*val_left
         else !Right flat monotone parabola
           coeffs(1) =  3.0_r_tran*density - 2.0_r_tran*val_right
           coeffs(2) = -6.0_r_tran*density + 6.0_r_tran*val_right
           coeffs(3) =  3.0_r_tran*density - 3.0_r_tran*val_right
         end if
      end if
    end if

  end subroutine parabola_mono

  !----------------------------------------------------------------------------------
  !> @brief computes the 3 coefficients of a cell parabola rho(x) which satifies
  !!        (i)   rho(0) = density_left,
  !!        (ii)  rho(1) = density_right,
  !!        (iii) integral(rho(x),x=0...1) = density
  !! @param[out]  coeffs          The cell parabola coefficients
  !! @param[in]   density_left    Parabola left value at x=0
  !! @param[in]   density         Density of the cell
  !! @param[in]   density_right   Parabola right value at x=1
  !----------------------------------------------------------------------------------
  subroutine calculate_parabola_coeffs(density_left,density_right,density_cell,coeffs)

    implicit none

    real(kind=r_tran),    intent(in)  :: density_left
    real(kind=r_tran),    intent(in)  :: density_right
    real(kind=r_tran),    intent(in)  :: density_cell
    real(kind=r_tran),    intent(out) :: coeffs(1:3)

    ! Calculate coefficients
    coeffs(1) = density_left
    coeffs(2) = -4.0_r_tran*density_left - 2.0_r_tran*density_right + 6.0_r_tran*density_cell
    coeffs(3) =  3.0_r_tran*density_left + 3.0_r_tran*density_right - 6.0_r_tran*density_cell

  end subroutine calculate_parabola_coeffs

end module subgrid_horizontal_support_mod
