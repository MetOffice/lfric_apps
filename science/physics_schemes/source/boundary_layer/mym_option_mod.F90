! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
!  Module mym_option_mod----------------------------------------------

!  Purpose: To define options and symbols in the Mellor-Yamada model

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_option_mod

use  missing_data_mod, only: rmdi, imdi
use control_max_sizes, only: max_bl_levels
use um_types, only: real_umphys

implicit none

!=======================================================================
! TKE namelist options - ordered as in meta-data, with namelist option
! and possible values
!=======================================================================

! 01 Switch to choose the TKE schemes
integer :: bdy_tke = imdi ! suggested mymodel3
integer, parameter ::                                                          &
! the first order scheme
     deardorff = 1,                                                            &
! the improved Mellor-Yamada level 2.5 model
     mymodel25 = 2,                                                            &
! the improved Mellor-Yamada level 3 model
     mymodel3  = 3

! 02 Switch to use code implemented for 3D TKE
logical :: l_3dtke = .false.

!
! 03 Maximum level to predict the prognostic variables in the TKE scheme.
integer :: tke_levels = imdi

! 04 If true and the parameter "tke_levels" is less than bl_levels, diffusion
! coefficients between tke_levels + 1 and bl_levels are evaluated with the
! stability function (i.e. local scheme)
logical :: l_local_above_tkelvs = .false.

! 05 If true, the prognostic variables in the TKE schemes are initialized
! assuming balance between production and dissipation. If missing values
! are set to the prognostic variables through the reconfiguration, the
! initialization will be automatically conducted.
logical :: l_my_initialize = .false. ! suggested true

! 06 In the case that l_my_initialize == .true. or missing values are set to
! the prognostic variables through the reconfiguration, they are initialized
! to zeros if it is true.
logical :: l_my_ini_zero = .false.

! 07 The lower limit for dbdz in the initialization to avoid to diagnose huge
! initial values.
real(kind=real_umphys) :: my_ini_dbdz_min = rmdi ! suggested 1.0e-5

! 08 A switch to turn on advection of the prognostic variables in the TKE
! scheme (E_TRB, TSQ, QSQ, COV).
logical :: l_adv_turb_field = .false. ! suggested true

! 09 If true, buoyancy parameters appearing in TKE production by buoyancy are
! evaluated with predicted variants assuming that fluctuation of heat and
! moisture can be described by the bi-normal probability distribution
! function. If false, buoyancy parameters calculated in the large scale
! clouds scheme (e.g. Smith, PC2), are employed to evaluate the buoyancy flux.
logical :: l_my_condense = .false. ! suggested true

! 10 If true, non-gradient buoyancy flux based on Lock and Mailhot(2006)
! is added.
logical :: l_shcu_buoy = .false. ! suggested true

! 11 Maximum level to evaluate the non-gradient buoyancy flux
integer :: shcu_levels = imdi ! suggested -1 to set to TKE levels

! 12 The maximum limit for the non-gradient buoyancy flux
real(kind=real_umphys) :: wb_ng_max = rmdi  ! suggested 0.05

! 13 Switch related to production terms at the lowest levels
integer :: my_lowest_pd_surf = imdi ! suggested bh1991
integer, parameter ::                                                          &
! No use of surface fluxes
     no_pd_surf = 0,                                                           &
! use the gradient function by Businger
     businger = 1,                                                             &
! use the gradient function by Beljaars and Holtslag
     bh1991 = 2

! 14 If true, production terms of covariance brought by the counter gradient
! terms are adjusted so that stability in integration can be secured.
logical :: l_my_prod_adj = .false. ! suggested true

! 15 if Z_TQ > MY_z_limit_elb, elb is limited less than vertical grid spacing.
real(kind=real_umphys) :: my_z_limit_elb = rmdi ! suggested 1.0e10 to not use

! 16 If true, the maximum values of the prognostic variables are printed.
logical :: l_print_max_tke = .false.

! 17 A proportional coef CM below the top of mixed layer K = CM * sqrt(E) * L
real(kind=real_umphys) :: tke_cm_mx = rmdi ! suggested 0.1

! 18 A proportional coef CM above the top of mixed layer K = CM * sqrt(E) * L
real(kind=real_umphys) :: tke_cm_fa = rmdi ! suggested 0.1

! 19 Switch to choose mixing length in the first order model
integer :: tke_dlen = imdi ! suggested my_length
integer, parameter ::                                                          &
! use the same mixing length in the MY model
     my_length = 1,                                                            &
! original mixing length suggested by Deardorff(1980)
     ddf_length = 2,                                                           &
! with correction based on Sun and Chang (1986)
     non_local_like_length = 3

!=======================================================================
! TKE options not in a namelist, ordered by type
!=======================================================================

! Choice of the higher order SL advection schemes for the prognostic
! variables in the TKE scheme.
integer, parameter :: high_order_scheme_adv_turb = 1

! Choice of the monotone SL advection schemes for the prognostic variables
! in the TKE scheme.
integer, parameter :: monotone_scheme_adv_turb = 1

! if L_MY_EXTRA_LEVEL == .true., the extra level is assigned at
! Z_TQ(:,:,1) * MY_Z_EXTRA_FACT above the surface.
! not a parameter as set in mym_initialise
real(kind=real_umphys) :: my_z_extra_fact = 0.5

! Factor in production term adjustment related to diffusion.
! A smaller factor makes the adjustment activate more often,
! but too strong adjustment might adversely affect the accuracy of
! forecasts. With the value 0.225, sufficient computational stability
! is secured in the UKV.
real(kind=real_umphys), parameter :: my_prod_adj_fact(1:max_bl_levels) = 0.225

! Use the correction to the mixing length by Blackadar (valid only in the
! first order model)
logical, parameter :: l_tke_dlen_blackadar = .false.

! Use the improved closure constants in the MY model by Nakanishi and Niino.
! If False, use the original value from Mellor-Yamada(1982)
logical, parameter :: l_my3_improved_closure = .true.

! A switch for higher order schemes of the SL advection of the prognostic
! variables in the TKE scheme.
logical, parameter :: l_high_adv_turb = .true.

! A switch for monotone schemes of the SL advection of the prognostic
! variables in the TKE scheme.
logical, parameter :: l_mono_adv_turb = .true.

! A switch for the conservative SL advection of the prognostic variables
! in the TKE scheme.
logical, parameter :: l_conserv_adv_turb = .false.

! If true, the extra level below the lowest level in the atmosphere is
! generated for the prognostic variables in the TKE scheme.
! not a parameter as set in mym_initialise
logical :: l_my_extra_level = .false.

! If true, the production terms for the covariaces are evaluated with surface
! fluxes. It is valid only when MY_lowest_pd_surf > 0.
logical, parameter :: l_my_lowest_pd_surf_tqc = .false.

end module mym_option_mod
