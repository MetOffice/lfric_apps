! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_cape_mod

implicit none

contains

! Subroutine to add contribution from the current convective level-step
! to the CAPE (and other vertical integrals).
! The integrals use a piecewise constant method, assuming a discontinuity
! half-way between neighbouring grid-levels
subroutine calc_cape( n_points, sublevs, i_next, l_within_bl,                  &
                      index_ij, ij_first, ij_last,                             &
                      fields_2d )

use comorph_constants_mod, only: real_cvprec, zero, one, half, gravity,        &
                                 l_calc_cape, l_calc_mfw_cape
use fields_2d_mod, only: n_fields_2d, i_cape,                                  &
                         i_mfw_cape, i_int_m_dz, i_m_ref
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_mean_buoy, j_massflux_d, j_env_tv

implicit none

! Number of points
integer, intent(in) :: n_points

! Super-array storing parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Previous model-level interface
! b) Next model-level interface
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: sublevs                              &
                     ( n_points, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(in) :: i_next(n_points)

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Indices used to reference a collapsed horizontal coordinate from the parcel
integer, intent(in) :: index_ij(n_points)

! Collapsed horizontal indices of the first and last point available
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Super-array containing CAPE and other 2D fields to be incremented
real(kind=real_cvprec), intent(in out) :: fields_2d                            &
                                          ( ij_first:ij_last, n_fields_2d )

! Height interval used in vertical integrals.
! NOTE: this goes negative when integrating downwards (for downdrafts),
! changing the sign of the resulting integral.  So e.g. a negatively buoyant
! downdraft will have positive CAPE.
real(kind=real_cvprec) :: dz

! Increments to vertical integrals
real(kind=real_cvprec) :: b_dz(n_points)
real(kind=real_cvprec) :: mb_dz(n_points)
real(kind=real_cvprec) :: m_dz(n_points)
real(kind=real_cvprec) :: mm_dz(n_points)

! Buoyancy differences, used for vertical interpolation of mass-flux
real(kind=real_cvprec) :: dbuoy_lev, dbuoy_next
! Interpolation weight
real(kind=real_cvprec) :: interp

! Loop counters
integer :: ic, i_lev, ij


! See if any points are above the BL-top (nothing to do otherwise)
if ( .not. all( l_within_bl ) ) then

  ! Initially calculate integrals assuming no sub-level steps...

  if ( l_calc_cape ) then

    ! Calculate CAPE contribution, g/Tv Tv' dz
    ! (only at points where mass-flux is nonzero)
    do ic = 1, n_points
      b_dz(ic) = zero
      if ( sublevs(ic,j_massflux_d,i_prev) > zero ) then
        b_dz(ic) = b_dz(ic) + sublevs(ic,j_mean_buoy,i_prev)                   &
                            / sublevs(ic,j_env_tv,i_prev)
      end if
      if ( sublevs(ic,j_massflux_d,i_prev+1) > zero ) then
        b_dz(ic) = b_dz(ic) + sublevs(ic,j_mean_buoy,i_prev+1)                 &
                            / sublevs(ic,j_env_tv,i_prev+1)
      end if
    end do
    do ic = 1, n_points
      dz = sublevs(ic,j_height,i_prev+1) - sublevs(ic,j_height,i_prev)
      b_dz(ic) = b_dz(ic) * gravity * half * dz
    end do

  end if  ! ( l_calc_cape )

  if ( l_calc_mfw_cape ) then

    do ic = 1, n_points
      dz = sublevs(ic,j_height,i_prev+1) - sublevs(ic,j_height,i_prev)

      ! Calculate contribution g/Tv Tv' M dz
      mb_dz(ic) = gravity * half * dz                                          &
       * ( sublevs(ic,j_massflux_d,i_prev) * sublevs(ic,j_mean_buoy,i_prev)    &
                                           / sublevs(ic,j_env_tv,i_prev)       &
         + sublevs(ic,j_massflux_d,i_prev+1) * sublevs(ic,j_mean_buoy,i_prev+1)&
                                             / sublevs(ic,j_env_tv,i_prev+1) )

      ! Contribution M dz
      m_dz(ic) = half * dz * ( sublevs(ic,j_massflux_d,i_prev)                 &
                             + sublevs(ic,j_massflux_d,i_prev+1) )

      ! Contribution M^2 dz
      mm_dz(ic) = half * dz                                                    &
       *( sublevs(ic,j_massflux_d,i_prev) * sublevs(ic,j_massflux_d,i_prev)    &
        + sublevs(ic,j_massflux_d,i_prev+1)*sublevs(ic,j_massflux_d,i_prev+1) )

    end do

  end if

  ! Now overwrite with corrected contributions where there are sub-level
  ! steps to account for...
  do ic = 1, n_points
    if ( i_next(ic) > i_prev + 1 .and. ( .not. l_within_bl(ic) ) ) then
      ! If there are any sub-levels between prev and next...
      ! (also don't bother with these expensive calculations if
      !  within the boundary-layer, where the answers won't be used).

      ! TEMPORARY CODE RETAINED TO PRESERVE KGO;
      ! (sub-level mass-fluxes should be calculated more accurately in set_det)

      ! Find mass-flux
      ! at each sub-level step in the buoyancy super-array...

      ! Linearly interpolate mass-flux in buoyancy space
      ! (assuming wiggles in buoyancy will cause wiggles in mass-flux;
      !  want to downweight the mass-flux weighting of the saturation
      !  height if the minimum buoyancy occurred there, as that will have
      !  caused detrainment seen as a reduced mass-flux at next).
      dbuoy_next = sublevs(ic,j_mean_buoy,i_next(ic))                          &
                 - sublevs(ic,j_mean_buoy,i_prev)
      do i_lev = i_prev+1, i_next(ic)-1
        dbuoy_lev = sublevs(ic,j_mean_buoy,i_lev)                              &
                  - sublevs(ic,j_mean_buoy,i_prev)
        if ( dbuoy_lev > zero ) then
          if ( dbuoy_lev < dbuoy_next ) then
            interp = dbuoy_lev / dbuoy_next
          else
            interp = one
          end if
        else if (dbuoy_lev < zero ) then
          if ( dbuoy_lev > dbuoy_next ) then
            interp = dbuoy_lev / dbuoy_next
          else
            interp = one
          end if
        else
          interp = zero
        end if
        sublevs(ic,j_massflux_d,i_lev)                                         &
          = (one-interp) * sublevs(ic,j_massflux_d,i_prev)                     &
          +      interp  * sublevs(ic,j_massflux_d,i_next(ic))
      end do

      if ( l_calc_cape ) then

        ! Calculate CAPE contribution, g/Tv Tv' dz
        ! (only at points where mass-flux is nonzero)
        b_dz(ic) = zero
        do i_lev = i_prev, i_next(ic)
          if ( sublevs(ic,j_massflux_d,i_lev) > zero ) then
            dz = zero
            if ( i_lev > i_prev ) then
              dz = dz + half * ( sublevs(ic,j_height,i_lev)                    &
                               - sublevs(ic,j_height,i_lev-1) )
            end if
            if ( i_lev < i_next(ic) ) then
              dz = dz + half * ( sublevs(ic,j_height,i_lev+1)                  &
                               - sublevs(ic,j_height,i_lev) )
            end if
            b_dz(ic) = b_dz(ic) + ( sublevs(ic,j_mean_buoy,i_lev)              &
                                  / sublevs(ic,j_env_tv,i_lev) )               &
                                * gravity * dz
          end if
        end do

      end if  ! ( l_calc_cape )

      if ( l_calc_mfw_cape ) then

        mb_dz(ic) = zero
        m_dz(ic) = zero
        mm_dz(ic) = zero
        do i_lev = i_prev, i_next(ic)

          dz = zero
          if ( i_lev > i_prev ) then
            dz = dz + half * ( sublevs(ic,j_height,i_lev)                      &
                             - sublevs(ic,j_height,i_lev-1) )
          end if
          if ( i_lev < i_next(ic) ) then
            dz = dz + half * ( sublevs(ic,j_height,i_lev+1)                    &
                             - sublevs(ic,j_height,i_lev) )
          end if

          ! Calculate contribution g/Tv Tv' M dz
          mb_dz(ic) = mb_dz(ic) + gravity * dz                                 &
             * sublevs(ic,j_massflux_d,i_lev) * sublevs(ic,j_mean_buoy,i_lev)  &
                                              / sublevs(ic,j_env_tv,i_lev)

          ! Contribution M dz
          m_dz(ic) = m_dz(ic) + dz * sublevs(ic,j_massflux_d,i_lev)

          ! Contribution M^2 dz
          mm_dz(ic) = mm_dz(ic) + dz                                           &
           * sublevs(ic,j_massflux_d,i_lev) * sublevs(ic,j_massflux_d,i_lev)

        end do  ! i_lev = i_prev, i_next(ic)

      end if  ! ( l_calc_mfw_cape )

    end if  ! ( i_next(ic) > i_prev + 1 )
  end do  ! ic = 1, n_points

  ! Add final contributions onto the main arrays only where NOT
  ! within the boundary-layer...

  if ( l_calc_cape ) then
    do ic = 1, n_points
      if ( .not. l_within_bl(ic) ) then
        ij = index_ij(ic)
        fields_2d(ij,i_cape) = fields_2d(ij,i_cape) + b_dz(ic)
      end if
    end do
  end if

  if ( l_calc_mfw_cape ) then
    do ic = 1, n_points
      if ( .not. l_within_bl(ic) ) then
        ij = index_ij(ic)
        fields_2d(ij,i_mfw_cape) = fields_2d(ij,i_mfw_cape) + mb_dz(ic)
        fields_2d(ij,i_int_m_dz) = fields_2d(ij,i_int_m_dz) + abs(m_dz(ic))
        fields_2d(ij,i_m_ref)    = fields_2d(ij,i_m_ref)    + abs(mm_dz(ic))
        ! Note: for the buoyancy integrals, we deliberately have opposite
        ! sign for downdrafts by having negative dz (so that negatively-buoyant
        ! downdrafts have positive CAPE).
        ! However, for the mass-flux integrals used for normalising the
        ! mass-flux-weighted CAPE, we wish the integrals to always be positive,
        ! so using ABS above to change sign back to positive for downdrafts.
      end if
    end do
  end if

end if  ! ( .NOT. ALL( l_within_bl ) )


return
end subroutine calc_cape

end module calc_cape_mod
