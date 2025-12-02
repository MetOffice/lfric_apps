! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: STASH related variables/arrays for describing output requests
!          and space management.
!
! Programming Standard : Unified Model Documentation paper number 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH
!
!
module stash_array_mod

use submodel_mod, only: n_internal_model
use version_mod,  only: nelemp
use parkind1,     only: jprb, jpim
use yomhook,      only: lhook, dr_hook
use ereport_mod,  only: ereport
use umPrintMgr,   only: umMessage, newline

implicit none

save

integer, parameter :: len_stlist   = nelemp ! No of items in stlist

! No of items per timeseries recd
integer, parameter :: time_series_rec_len = 9

integer :: nsects                ! Max no of diagnostic sections
integer :: n_req_items           ! Max item number in any section
integer :: nitems                ! No of distinct items requested
integer :: n_ppxrecs             ! No of PP_XREF records this run
integer :: totitems              ! Total no of processing requests
integer :: nsttims               ! Max no of STASHtimes in a table
integer :: nsttabl               ! No of STASHtimes tables
integer :: num_stash_levels      ! Max no of levels in a levelslist
integer :: num_level_lists       ! No of levels lists
integer :: num_stash_pseudo      ! Max no of pseudo-levs in a list
integer :: num_pseudo_lists      ! No of pseudo-level lists
integer :: nstash_series_block   ! No of blocks of timeseries recds
integer :: nstash_series_records ! Total no of timeseries records


! This file is needed to get ppxref_codelen to dimension PP_XREF
! sizes in STASH used for defining local array dimensions at a
! lower level.
integer :: max_stash_levs  ! Max no of output levels for any diag
integer :: pp_len2_lookup  ! Max no of lookups needed in stwork


! stashflag (.true. for processing this timestep). sf(0,is) .false.
! if no flags on for section is.
logical, allocatable :: sf (:,:)

! Whether a calculation is needed for sf above
logical, allocatable :: sf_calc (:,:)

! STASH list index
integer, allocatable :: stindex (:,:,:,:)

! List of STASH output requests
integer, allocatable :: stlist (:,:)

! Address of item from generating plug compatible routine (often workspace)
integer, allocatable :: si (:,:,:)

! Address of the end of item to allow slices to be passed to copydiag
integer, allocatable :: si_last (:,:,:)

! STASH times tables
integer, allocatable :: sttabl (:,:)

! Length of STASH workspace required in each section
integer, allocatable :: stash_maxlen        (:,:)
integer, allocatable :: ppindex             (:,:)
integer, allocatable :: stash_levels        (:,:)
integer, allocatable :: stash_pseudo_levels (:,:)
integer, allocatable :: stash_series        (:,:)
integer, allocatable :: stash_series_index  (:,:)

character(len=*), parameter, private :: ModuleName = 'STASH_ARRAY_MOD'

private :: n_internal_model

contains

subroutine allocate_stash_arrays()

! Purpose: Allocate STASH arrays based on current values of
!          the dimensioning variables!

implicit none

integer                       :: icode
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName = 'ALLOCATE_STASH_ARRAYS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

if (allocated( sf                  )) deallocate(sf)
if (allocated( sf_calc             )) deallocate(sf_calc)
if (allocated( stindex             )) deallocate(stindex)
if (allocated( stlist              )) deallocate(stlist)
if (allocated( si                  )) deallocate(si)
if (allocated( si_last             )) deallocate(si_last)
if (allocated( sttabl              )) deallocate(sttabl)
if (allocated( stash_maxlen        )) deallocate(stash_maxlen)
if (allocated( ppindex             )) deallocate(ppindex)
if (allocated( stash_levels        )) deallocate(stash_levels)
if (allocated( stash_pseudo_levels )) deallocate(stash_pseudo_levels)
if (allocated( stash_series        )) deallocate(stash_series)
if (allocated( stash_series_index  )) deallocate(stash_series_index)


allocate( sf (0:nitems, 0:nsects)                                              &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf]'//                    newline//&
              trim(umMessage) )


allocate( sf_calc (0:nitems, 0:nsects)                                         &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf_calc]'//               newline//&
              trim(umMessage) )


allocate( stindex (2, nitems, 0:nsects, n_internal_model)                      &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf_calc]'//               newline//&
              trim(umMessage) )


allocate( stlist (len_stlist, totitems)                                        &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stlist]'//                newline//&
              trim(umMessage) )


allocate( si (nitems, 0:nsects, n_internal_model)                              &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [si]'//                    newline//&
              trim(umMessage) )

allocate( si_last (nitems, 0:nsects, n_internal_model)                         &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [si_last]'//               newline//&
              trim(umMessage) )


allocate( sttabl (nsttims, nsttabl)                                            &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sttabl]'//                newline//&
              trim(umMessage) )


allocate( stash_maxlen (0:nsects, n_internal_model)                            &
         , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_maxlen]'//          newline//&
              trim(umMessage) )


allocate( ppindex (nitems, n_internal_model)                                   &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [ppindex]'//               newline//&
              trim(umMessage) )


allocate( stash_levels (num_stash_levels+1, num_level_lists)                   &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_levels]'//          newline//&
              trim(umMessage) )


allocate( stash_pseudo_levels (num_stash_pseudo+1, num_pseudo_lists)           &
        , stat=icode, errmsg=umMessage)
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_pseudo_levels]'//   newline//&
              trim(umMessage) )


allocate( stash_series (time_series_rec_len, nstash_series_records)            &
        , stat=icode, errmsg=umMessage)
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_series]'//          newline//&
              trim(umMessage) )


allocate( stash_series_index (2, nstash_series_block)                          &
        , stat=icode, errmsg=umMessage )
if (icode /= 0)                                                                &
call ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_series_index]'//    newline//&
              trim(umMessage) )


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine allocate_stash_arrays

end module stash_array_mod
