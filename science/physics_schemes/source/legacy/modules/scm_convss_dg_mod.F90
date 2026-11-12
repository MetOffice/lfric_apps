! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module containing a type definition for SCM diagnostics from convection,
!  and subroutines to initialise, output and deallocate the diagnostic arrays.

module scm_convss_dg_mod


implicit none

!
! Description:
!   Defines a derived type structure containing all the SCM diagnostics that
!   need to be passed up from inside the convection scheme (e.g. diagnostics
!   based on variables which are only available inside deep_conv_mod).
!   SCMoutput can't be called from inside the convection scheme, because
!   it is only called on some timesteps at some points, using compressed
!   arrays.  SCMoutput assumes the data passed to it are of dimension
!   (row_length, rows).  Therefore, the type defined here is used to declare
!   a structure, containing numerous diagnostics, with the full dimensions
!   in ni_conv_ctl, which is passed down into compressed dimensions in
!   glueconv, deep_conv_mod, etc.
!   This module also contains subroutines to allocate and initialise the
!   diagnostics to zero, and output them via scmoutput.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: scm
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

character(len=*), parameter, private :: ModuleName =                           &
                                        'SCM_CONVSS_DG_MOD'



!------------------------------------------------------------------------
! Type definition for SCM convection diagnostics
!------------------------------------------------------------------------
! NOTE: If adding more diagnostics to this list, they also need to be
! added to the equivalent type definition in the stub routine
! src/scm/stub/scm_convss_dg_mod_stub.F90.
! Otherwise, where the new components are referenced in the UM code,
! the full UM will fail at compile time due to the stated fields not
! actually being components of this type.
type :: scm_convss_dg_type

  ! Note that making all the fields below allocatable is unnecessary,
  ! given that the first routine that uses this type already knows that
  ! all the arrays should have the same size (model_levels).
  ! This structure is itself declared as an allocatable array; no
  ! memory is used until that is allocated, so there is no advantage
  ! to using allocatable fields in here.
  ! This would be an ideal time to use the fortran 2003 parameterised
  ! derived type feature (i.e. it allows all the arrays to be dimensioned
  ! with model_levels from the outset, by making model_levels an argument
  ! to this type).  We then wouldn't need the allocate / deallocate
  ! subroutines below.  However, UMDP3 states that f2003 parameterised
  ! derived types are not allowed, so we'll just have to do it this way
  ! for now.

  ! Environment profile seen by the convective main ascent
  real, allocatable :: env_theta(:)
  real, allocatable :: env_q(:)
  real, allocatable :: env_qcl(:)
  real, allocatable :: env_qcf(:)
  real, allocatable :: env_thetav(:)

  ! Parcel properties from the main ascent
  real, allocatable :: par_theta(:)
  real, allocatable :: par_q(:)
  real, allocatable :: par_qcl(:)
  real, allocatable :: par_qcf(:)
  real, allocatable :: par_thetav(:)

  ! Parcel buoyancy excess
  real, allocatable :: par_thetav_excess(:)

  ! Sub-step mass-flux profile before and after closure
  real, allocatable :: up_flx_guess(:)
  real, allocatable :: up_flx(:)

  ! Fractional entrainment and mixing detrainment rates in layer
  real, allocatable :: ekp14(:)
  real, allocatable :: ekp34(:)
  real, allocatable :: amdetk(:)

  ! Diagnostics for adaptive detrainment
  real, allocatable :: deltak(:)
  real, allocatable :: rbuoy_star(:)
  real, allocatable :: xsbmin(:)
  real, allocatable :: thrk(:)
  real, allocatable :: qrk(:)
  real, allocatable :: thvrk_excess(:)

  ! Status of deep, shallow and mid-level convection
  ! (0 = not diagnosed, 1 = failed ascent, 2 = zero closure, 3 = real convection
  integer :: status_deep
  integer :: status_shallow
  integer :: status_mid

  ! Surface precip from each convection type
  real :: precip_deep
  real :: precip_shallow
  real :: precip_mid

end type scm_convss_dg_type


contains




end module scm_convss_dg_mod
