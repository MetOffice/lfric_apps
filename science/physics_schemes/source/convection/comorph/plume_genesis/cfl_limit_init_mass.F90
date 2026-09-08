! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cfl_limit_init_mass_mod

implicit none

contains

! Subroutine to impose a CFL limit on the mass-flux initiating
! from each model-level / sub-grid region
subroutine cfl_limit_init_mass( n_points, nc, index_ic, n_conv_types,          &
                                max_ent_frac, l_within_bl,                     &
                                layer_mass_k, frac_r_k ,                       &
                                init_mass_t )

use comorph_constants_mod, only: real_cvprec, zero, one, comorph_timestep,     &
                                 i_cfl_local, i_cfl_local_nobl

implicit none

! Total number of points in super-arrays
integer, intent(in) :: n_points

! Points where initiation mass-sources occur in current region
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Number of convection types
integer, intent(in) :: n_conv_types

! Max fraction of mass in region which can initiate (CFL limit)
real(kind=real_cvprec), intent(in) :: max_ent_frac

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Total dry-mass on level k
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Fraction occupied by current region
real(kind=real_cvprec), intent(in) :: frac_r_k(n_points)

! Initiating mass-flux from current region, for each convection type
real(kind=real_cvprec), intent(in out) :: init_mass_t(n_points,n_conv_types)

! Scaling factor needed to reduce initiating mass below limit
real(kind=real_cvprec) :: factor(nc)

! Sum of initiating mass over convection types
real(kind=real_cvprec) :: init_mass_sum

! Loop counters
integer :: ic, ic2, i_type


if ( n_conv_types == 1 ) then
  ! Calculation is simpler when there is only 1 convection type
  ! TEMPORARY CODE TO PRESERVE KGO:
  ! the multi-type implementation is perfectly correct with only 1 type,
  ! but changes answers at bit-level compared to what we had before so
  ! keeping the original code in the case where there is only 1 type.
  if ( i_cfl_local==i_cfl_local_nobl ) then
    ! Only apply vertically local limit above the BL-top
    do ic2 = 1, nc
      ic = index_ic(ic2)
      if ( .not. l_within_bl(ic) ) then
        init_mass_t(ic2,1) = min( init_mass_t(ic2,1),                          &
                                  max_ent_frac * layer_mass_k(ic)              &
                                  * frac_r_k(ic) / comorph_timestep )
      end if
    end do
  else
    ! Always impose vertically local limit
    do ic2 = 1, nc
      ic = index_ic(ic2)
      init_mass_t(ic2,1) = min( init_mass_t(ic2,1),                            &
                                max_ent_frac * layer_mass_k(ic)                &
                                * frac_r_k(ic) / comorph_timestep )
    end do
  end if

else
  ! Multiple convection types; need to apply CFL limit to the sum over
  ! types and then scale down each type proportionally...

  ! Calculate scaling factor needed to reduce init_mass to the CFL limit
  ! (defaults to one when mass is already below limit)
  do ic2 = 1, nc
    ic = index_ic(ic2)
    init_mass_sum = zero
    do i_type = 1, n_conv_types
      init_mass_sum = init_mass_sum + init_mass_t(ic2,i_type)
    end do
    factor(ic2) = min( max_ent_frac * layer_mass_k(ic) * frac_r_k(ic)          &
                      / ( comorph_timestep * init_mass_sum ), one )
  end do

  if ( i_cfl_local==i_cfl_local_nobl ) then
    ! Options to only apply local limit when above the boundary-layer top
    ! (so reset factor to 1 when within the BL)
    do ic2 = 1, nc
      ic = index_ic(ic2)
      if ( l_within_bl(ic) )  factor(ic2) = one
    end do
  end if

  ! Apply the CFL limit scaling
  do i_type = 1, n_conv_types
    do ic2 = 1, nc
      init_mass_t(ic2,i_type) = init_mass_t(ic2,i_type) * factor(ic2)
    end do
  end do

end if


return
end subroutine cfl_limit_init_mass

end module cfl_limit_init_mass_mod
