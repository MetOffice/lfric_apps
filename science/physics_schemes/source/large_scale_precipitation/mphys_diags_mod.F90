! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module mphys_diags_mod

! Description:
! Holds diagnostic switches and values required by the large-scale
! precipitation scheme that have been integrated over the particle size
! distribution
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use um_types, only: real_lsprec

implicit none

save

! Switches for aggregate fraction and points used
logical :: l_aggfr_diag = .false. ! Aggregate Fraction

logical :: l_point_diag = .false. ! Points used

! Microphysical process rate diagnostic logical switches
logical :: l_psdep_diag = .false.
                     ! Deposition of vapour to snow agg.
logical :: l_psaut_diag = .false.
                     ! Autoconversion of aggregates from cry
logical :: l_psacw_diag = .false.
                     ! Accretion of liq. water by snow agg.
logical :: l_psacr_diag = .false.
                    ! Collection of rain by snow aggregates
logical :: l_psaci_diag = .false.
                     ! Collection of ice crystals by agg.
logical :: l_psmlt_diag = .false.
                     ! Melting of snow aggregates
logical :: l_psmltevp_diag = .false.
                     ! Evaporation of melting aggregates
logical :: l_praut_diag = .false.
                     ! Autoconversion of cloud drops to rain
logical :: l_pracw_diag = .false.
                     ! Accretion of liq. water by rain
logical :: l_prevp_diag = .false.
                     ! Evaporation of rain
logical :: l_pgaut_diag = .false.
                     ! Autoconversion of graupel from agg.
logical :: l_pgacw_diag = .false.
                     ! Accretion of liq. water by graupel
logical :: l_pgacs_diag = .false.
                     ! Collection of snow agg. by graupel
logical :: l_pgmlt_diag = .false.
                     ! Melting of graupel
logical :: l_pifrw_diag = .false.
                     ! Homogeneous freezing nucleation
logical :: l_piprm_diag = .false.
                     ! Heterogeneous (primary) nucleation
logical :: l_pidep_diag = .false.
                     ! Deposition of vapour to ice crystals
logical :: l_piacw_diag = .false.
                     ! Accretion of liq. water by ice cry.
logical :: l_piacr_diag = .false.
                     ! Collection of rain by ice crystals
logical :: l_pimlt_diag = .false.
                     ! Melting of ice crystals
logical :: l_pimltevp_diag = .false.
                     ! Evaporation of melting ice crystals
logical :: l_pifall_diag = .false.
                     ! Sedimentation of ice crystals
logical :: l_psfall_diag = .false.
                     ! Sedimentation of aggregates
logical :: l_prfall_diag = .false.
                     ! Sedimentation of rain
logical :: l_pgfall_diag = .false.
                     ! Sedimentation of graupel
logical :: l_plset_diag = .false.
                     ! Droplet settling of liquid water
logical :: l_plevpset_diag = .false.
                     ! Evaporated settled droplets
logical :: l_pifrr_diag = .false.
                     ! Homogeneous freezing of rain
logical :: l_piprr_diag = .false.
                     ! Heterogeneous freezing of rain
logical :: l_vtbranch_diag = .false.
                     ! Flag for vt-D branch choice
logical :: l_vm_cry_diag = .false.
                     ! Fallspeed with crystal parameters
logical :: l_vm_agg_diag = .false.
                     ! Fallspeed with aggregate parameters
logical :: l_vm_used_diag = .false.
                     ! Fallspeed used


logical :: l_refl_tot  = .false.
                     ! Total reflectivity diagnostic
logical :: l_refl_gr   = .false.
                     ! Graupel reflectivity diagnostic
logical :: l_refl_ice  = .false.
                     ! Ice Aggregate reflectivity diagnostic
logical :: l_refl_ice2 = .false.
                     ! Ice Crystal reflectivity diagnostic
logical :: l_refl_rain = .false.
                     ! Rain reflectivity diagnostic
logical :: l_refl_qcl  = .false.
                     ! Liquid cloud reflectivity diagnostic

logical :: l_refl_1km  = .false.
                     ! Total reflectivity at 1km AGL
logical :: l_refl_surf = .false.
                     ! Total reflectivity at surface (lowest model level)
logical :: l_refl_max  = .false.
                     ! Maximum reflectivity in the column

logical :: l_cth_dbZ  = .false.
                     ! Reflectivity echo-top height in metres AGL

! Switches for seeder feeder scheme
logical :: l_sfwater_diag = .false.
                     ! Subgrid orographic water mixing ratio
logical :: l_sfrain_diag = .false.
                     ! Subgrid orographic rain production
logical :: l_sfsnow_diag = .false.
                     ! Subgrid orographic snow production

real(kind=real_lsprec), allocatable ::                                         &
  psdep(:,:,:),  psaut(:,:,:),   psacw(:,:,:),  psacr(:,:,:), psaci(:,:,:),    &
  psmlt(:,:,:),  psmltevp(:,:,:), praut(:,:,:),  pracw(:,:,:), prevp(:,:,:),   &
  frac_agg(:,:,:), pgaut(:,:,:),  pgacw(:,:,:),  pgacs(:,:,:), pgmlt(:,:,:),   &
  pifrw(:,:,:),    pifrr(:,:,:),  piprm(:,:,:),  pidep(:,:,:), piacw(:,:,:),   &
  piacr(:,:,:), pimlt(:,:,:), pimltevp(:,:,:), pifall(:,:,:), psfall(:,:,:),   &
  prfall(:,:,:), pgfall(:,:,:), plset(:,:,:), plevpset(:,:,:), piprr(:,:,:),   &
  vm_cry(:,:,:), vm_agg(:,:,:), vtbranch_flag(:,:,:), vm_used(:,:,:)


real(kind=real_lsprec), allocatable :: dbz_tot(:,:,:) ! Total reflectivity (dBZ)
real(kind=real_lsprec), allocatable :: dbz_g(:,:,:)
                                    ! Graupel reflectivity (dBZ)
real(kind=real_lsprec), allocatable :: dbz_i(:,:,:)
                                    ! Ice Agg. reflectivity (dBZ)
real(kind=real_lsprec), allocatable :: dbz_i2(:,:,:)
                                    ! Ice Cry. reflectivity (dBZ)
real(kind=real_lsprec), allocatable :: dbz_l(:,:,:)
                                    ! Cloud liquid reflectivity (dBZ)
real(kind=real_lsprec), allocatable :: dbz_r(:,:,:)   ! Rain reflectivity (dBZ)

real(kind=real_lsprec), allocatable ::                                         &
  sfwater(:,:,:), sfrain(:,:,:), sfsnow(:,:,:)

logical, allocatable :: mphys_pts(:,:,:) ! Flag for whether grid box is used

! Logical flags for warnings when user outputs similar graupel and snowfall
! diagnostics on the same timestep.
! It would be rather annoying if these were to output a warning message
! on every single timestep, so these flags control that the warning is only
! produced once, if required.
! Their values are initialised to .true. and if the STASH flags are set
! a warning message is output and the flag is then set to .false.
logical :: l_warn_gr_am  = .true. ! Flag for warning when graupel and
                                  ! snow amounts output together
logical :: l_warn_gr_ra  = .true. ! Flag for warning when graupel and
                                  ! snow rates output together
logical :: l_warn_gr_ra3 = .true. ! Flag for warning when graupel and
                                  ! snow rates output together

end module mphys_diags_mod
