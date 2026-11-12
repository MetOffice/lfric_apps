! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash

module cstash_mod

use version_mod, only: ndiagpm, nproftp, ntimep, nprofdp, NTimSerP,            &
                       nlevp, npslevp,npslistp,nprofup,nprofpp

use profilename_length_mod, only: profilename_length,packagename_length

implicit none

! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.

! Global scalars:
integer ::   ndiag   ! No. of diagnostics
integer ::   ntprof  ! No. of time profiles
integer ::   nseries ! No. of stash time series
integer ::   ndprof  ! No. of domain profiles
integer ::   nuprof  ! No. of useage profiles
integer ::   npprof  ! No. of excluded package profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
integer ::   modl_b(ndiagpm)  ! Internal model no.
integer ::   isec_b(ndiagpm)  ! Section
integer ::   item_b(ndiagpm)  ! Item
integer ::   itim_b(ndiagpm)  ! Time profile number
integer ::   idom_b(ndiagpm)  ! Domain profile number
integer ::   iuse_b(ndiagpm)  ! Useage profile number

!   Time profile information:

character(len=profilename_length) :: timpro(nproftp)   ! Name of time profile
integer ::   ityp_t(nproftp)         ! Type of profile
integer ::   intv_t(nproftp)         ! Time Interval
integer ::   unt1_t(nproftp)         ! Units for time interval
integer ::   isam_t(nproftp)         ! Sampling period
integer ::   unt2_t(nproftp)         ! Units for sampling period
integer ::   iopt_t(nproftp)         ! Output option
integer ::   istr_t(nproftp)         ! Output Start time
integer ::   iend_t(nproftp)         ! Output End time
integer ::   isdt_t(6, nproftp)      ! Output Start date
integer ::   iedt_t(6, nproftp)      ! Output End date
integer ::   ifre_t(nproftp)         ! Output frequency
integer ::   ioff_t(nproftp)         ! Offset for sampling
integer ::   unt3_t(nproftp)         ! Units for output times
integer ::   itim_t(nproftp)         ! No. of times in times table
integer ::   iser_t(ntimep ,nproftp) ! Times table (with units)
integer ::   modl_t(nproftp)         ! Indicates internal model
                                     !  for each times table
logical ::   lts0_t(nproftp)         ! Logical for PP output at timestep zero

!   Domain profile information:
character(len=profilename_length) :: dompro  (nprofdp) ! Name of domain profile
integer ::  iopl_d  (nprofdp)           ! Levels option
integer ::  levb_d  (nprofdp)           ! Bottom level
integer ::  levt_d  (nprofdp)           ! Top level
integer ::  iopa_d  (nprofdp)           ! Area option
integer ::  inth_d  (nprofdp)           ! North boundary
integer ::  isth_d  (nprofdp)           ! South boundary
integer ::  iest_d  (nprofdp)           ! East boundary
integer ::  iwst_d  (nprofdp)           ! West boundary
integer ::  imsk_d  (nprofdp)           ! Mask type
integer ::  imn_d   (nprofdp)           ! Meaning option
integer ::  iwt_d   (nprofdp)           ! Weighting option
logical ::  ts_d    (nprofdp)           ! Time series profile
integer ::  ig_ts
integer ::  i1_ts
integer ::  i51_ts
integer ::  blim_ts (NTimSerP)
integer ::  tlim_ts (NTimSerP)
real ::     blimr_ts(NTimSerP)
real ::     tlimr_ts(NTimSerP)
integer ::  nlim_ts (NTimSerP)
integer ::  slim_ts (NTimSerP)
integer ::  elim_ts (NTimSerP)
integer ::  wlim_ts (NTimSerP)
integer ::  ilev_d  (nprofdp)           ! Output levels code
integer ::  levlst_d(nlevp   ,nprofdp ) ! Levels list
real ::    rlevlst_d(nlevp   ,nprofdp ) ! Levels list
integer ::  plt_d   (nprofdp)
integer ::  pllen_d (nprofdp)
integer ::  plpos_d (nprofdp)
integer ::  pslist_d(npslevp ,npslistp)
integer ::  npslists

! Useage information:

character(len=profilename_length) :: usepro(nprofup) ! Name of useage profile
integer ::  locn_u(nprofup)   ! Storage location of profile
integer ::  iunt_u(nprofup)   ! Unit no.
logical ::  lnetcdf_u(nprofup)! Is output format netCDF?

!Package information
character(len=packagename_length) :: pckpro(nprofpp)

! Information from ppxref file:

integer ::   model_st       ! Internal model number
integer ::   ispace         ! Space code
integer ::   itima          ! Time availability code
integer ::   igp            ! Grid of data code
integer ::   ilev           ! Level type code
integer ::   ibot           ! First level code
integer ::   itop           ! Last level code
integer ::   iflag          ! Level compression flag
integer ::   iopn(6)        ! Sectional option code
integer ::   vmsk           ! Integer equiv of bin vers mask
integer ::   ipseudo        ! Pseudo dimension type
integer ::   ipfirst        ! First pseudo dim code
integer ::   iplast         ! Last pseudo dim code
integer ::   ptr_prog       ! Section zero point back
integer ::   halo_type      ! Type of halo the field has

! Available time units

integer, parameter ::  stsh_timesteps = 1
integer, parameter ::  stsh_hours     = 2
integer, parameter ::  stsh_days      = 3
integer, parameter ::  stsh_dumps     = 4
integer, parameter ::  stsh_minutes   = 5
integer, parameter ::  stsh_seconds   = 6
integer, parameter ::  stsh_realmonth = 7

end module cstash_mod
