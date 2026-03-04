! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
! Program does a standalone single-timestep test-run of the
! entire convection scheme, using idealised initial profiles
program test_comorph

use comorph_constants_mod, only: nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 k_top_init, n_tracers
use standalone_test_mod, only: standalone_test

implicit none

! Set array sizes etc
nx_full = 10
ny_full = 10
k_bot_conv = 1
k_top_conv = 50
k_top_init = 40
n_tracers = 2

! Call standalone test interface routine to CoMorph
call standalone_test()

end program test_comorph
#endif
