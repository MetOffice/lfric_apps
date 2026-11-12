! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Defines the derived-type, SCMop_type, necessary for declaring
!  SCMop - the structure that carries all the SCM diagnostic
!  information.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: scm

module scmoptype_defn

! Description:
! The NetCDF library routines were written to work in 32 bit and
! so we may have to call 32 bit routines from the 64 bit
! UM. Using the Intel compiler on the Linux platform this can
! easily be done by explicitly declaring the integers variables
! that will be inputs to or outputs from the NetCDF routines to
! be 32 bit, and the UM routines that call them
! (DUMP_STREAMS_INIT, DUMP_STREAMS and DUMP_STREAMS_END) can
! still be compiled in 64 bits.

! The NEC compiler however, will ignore any explicit requests for
! 32 bit integers in a routine that is being compiled at 64 bits
! - everything will be promoted to 64 bits. Thus on the NEC the UM
! routines that call the NetCDF library routines have to be compiled
! at 32 bits and have the inputs that they are recieving from the
! rest of the 64 bit UM explicitly declared as 64 bit. The upshot is
! that, in order to support both these methods, we will explictly
! declare the precision of every variable being passed to or
! from DUMP_STREAMS_INIT, DUMP_STREAMS, DUMP_STREAMS_END or any
! NetCDF library routine. This includes the contents of the
! SCMop structure, defined here.

use um_types, only: integer_64, real_64, logical_64, integer_32
implicit none

! Short-hands for integer_64, real_64 and logical_64
integer, parameter :: i64 = integer_64
integer, parameter :: r64 = real_64
integer, parameter :: l64 = logical_64

! Precision of integers passed in to & out of NetCDF routines
integer, parameter :: incdf = integer_32

! Some limits
integer, parameter :: maxndomprof = 100

! The lengths of certain character variables
integer, parameter ::                                                          &
  lfilename = 100                                                              &
, lsname    = 30                                                               &
, llname_old= 50      & ! Long-name field width used to write text file output,
                        ! retained at 50 to ensure backwards compatibility
                        ! with old programs that read SCM text-format files.
, llname    = 100     & ! Long-name field width used internally, and for
                        ! writing the long-name in netcdf output. This has been
                        ! increased to prevent truncation of some long-names.
, lunits    = 15

! If the filename of a stream is set to this the name will be
! ignored and the stream will be opened with the default name
! given the unit number.
character(len=lfilename), parameter :: default_stream='<default>'

! The maximum no. of streams must be such that
! Stream(maxnstreams+1) gives an integer which is not
! too high for the machine precision
integer, parameter :: maxnstreams  = 30
integer, parameter :: inot_written = maxnstreams+1
integer, parameter :: listlength   = 5000

!-----------------------------------------------------------------------
! A diagnostic's dump array has its own type so that we can
! declare an array of arrays and allocate memory to each
! separately
!-----------------------------------------------------------------------
type :: allocatable_array
  real(r64), pointer :: dump(:)
end type allocatable_array

!-----------------------------------------------------------------------
! A stream has its own type to clearly separate its entries in
! SCMop from those corresponding to individual diagnostics
!-----------------------------------------------------------------------
type :: astream

  ! The unit the stream will write to
  integer(i64) :: op_unit

  ! The name of the file it will create (if not set to default)
  character (len=lfilename) :: filename
       ! The dumping period of the diagnostics sent to this stream
       ! (i.e. if a diagnostic is an average it is the number of
       ! timesteps it is averaged over, if it is an accumulation it
       ! is the number of timesteps it is accumulated over, if it is
       ! a maximum it is the number of timesteps it is the maximum
       ! over, etc.)

  integer(i64) :: dump_step
       ! Flags whether diagnostics sent to this stream by the
       ! respective input to routine SCMoutput will actually be
       ! sent to this stream.

  integer(i64) :: heed_hardwired
       ! Flags whether diagnostics sent to this stream by namelist
       ! request will actually be sent to this stream.

  integer(i64) :: heed_acceptlist
       ! Flags whether diagnostics prevented from going to this
       ! stream by namelist request will actually be prevented.

  integer(i64) :: heed_rejectlist

  character (len=lsname), pointer :: accept_list(:)
  character (len=lsname), pointer :: reject_list(:)

  integer(i64) :: switch ! If zero, stream is not active.

  integer(i64) :: output_format    ! Determines format of output file.
                                   ! 0 = format intended for subsequent reading
                                   !     by PV-wave routine scmread2.pro
                                   !     (used by scmoutput.pro)
                                   ! 1 = format geared to easy perusal by eye
                                   ! 2 = format suitable for FSSI database
                                   ! 3 = new PV-wave format designed to replace
                                   !     format 0. Can be read by same routines.

  integer(i64) :: n_output  ! The number of diagnostics that
                            ! will be output to this stream
end type astream


!-----------------------------------------------------------------------
! Define the derived type for SCMop. SCMop carries all
! necessary diagnostic information from the top(ish) level
! down to wherever any diagnostic is actually calculated,
! and then back up to the top for output.
!-----------------------------------------------------------------------
type :: SCMop_type

       ! Flags whether diagnostic system is "switched on"
  logical(l64) :: on

       ! first_pass will be true during all calls to SCMoutput in
       ! the first timestep (a formative stage for the list of
       ! diagnostics), and false thereafter (when the creation of
       ! new diagnostics will not be allowed).
  logical(l64) :: first_pass

       ! Pointers to daycount and stepcount, and knowledge of
       ! full_daysteps (all declared and defined in scm_main) are
       ! required so SCMop can tell the time.
  integer(i64), pointer :: daycount,stepcount

  integer(i64) :: full_daysteps

       ! Knowledge of ntrad1 and ntrad (the first timestep on which
       ! radiation is called and the number of timesteps between calls
       ! thereafter) is required by SCMop to make sense of diagnostics
       ! only calculated on radiation timesteps.
  integer(i64) :: ntrad1,ntrad

       ! An encoded integer representing which output streams are
       ! open (i.e. which streams will be output to file).
  integer(i64) :: openstreams

       ! We want a certain number of streams, we will allocate
       ! exactly how many at runtime.
  type(astream), pointer :: strm(:)

       ! maxnentries will be the size of all the arrays associated
       ! with the diagnostic entries (sname, etc.) once allocated.
       ! It can be increased at runtime in routine expand_SCMop.
  integer(i64) :: maxnentries

       ! nentries, n_output and nSCMoutput will be the total number
       ! of diagnostic entries in SCMop, the number of those entries
       ! being output to any stream (with no multiple counting of
       ! entries resulting from the same call to SCMoutput), and
       ! the number of calls made to SCMoutput so far this timestep
       ! respectively.
  integer(i64) :: nentries,n_output,nSCMoutput

       ! Will hold the total number of expected/observed sub-steps
       ! and the current sub-step number
  integer(i64) :: num_substeps,substep_number

       ! The diagnostic entries, 1:nentries have been set by newdiag
       ! via a call to SCMoutput. if A NEW ARRAY is ADDED HERE, THERE
       ! MUST BE AN associated SECTION OF CODE in ROUTINE EXPAND_SCMOP.
  character(len=lsname), pointer :: sname(:) ! Short name
  character(len=llname), pointer :: lname(:) ! Long name
  character(len=lunits), pointer :: units(:) ! Units

  integer(i64), pointer ::                                                     &
    domprof(:)       &! Domain profile
  , timprof(:)       &! Time profile
  , streams(:)       &! List of streams to write to
  , dump_step(:)     &! Dumping period
  , nadd2dump(:)     &! Number of calls to add2dump this period
  , ncols(:)         &! No. of columns
  , nrows(:)         &! No. of rows
  , nlevs(:)         &! No. of levels
  , nelements(:)     &! elements in total (given the domain profile)
  , sname_id(:)      &! An integer unique to each sname
  , wd(:)            &! Index of another entry upon which this entry depends
  , lastencounter(:) &! Last timestep this entry was seen by SCMoutput
  , substep(:)        ! Substep this entry pertains to

  logical(l64), pointer ::                                                     &
    only_radsteps(:)  ! Only defined on radiation timesteps?

  integer(incdf), pointer ::                                                   &
    netcdf_id(:)      ! A NetCDF integer id

  type(allocatable_array), pointer ::                                          &
    diag(:)           ! Dump array

  ! The domain profiles, set by define_domprof.
  character (len=15) ::                                                        &
    d_name(maxndomprof)

  integer(i64) ::                                                              &
    d_rowa1 (maxndomprof)                                                      &
  , d_rowa2 (maxndomprof)                                                      &
  , d_rowb1 (maxndomprof)                                                      &
  , d_rowb2 (maxndomprof)                                                      &
  , d_lev1  (maxndomprof)                                                      &
  , d_lev2  (maxndomprof)

       ! This array will be used to memorise the order in which the
       ! calls to SCMoutput take place, which saves time in SCMoutput
       ! by avoiding the translation between the inputs and the
       ! corresponding diagnostic entries in this structure
  integer(i64), pointer :: diag_mem(:,:)

end type SCMop_type

end module scmoptype_defn

