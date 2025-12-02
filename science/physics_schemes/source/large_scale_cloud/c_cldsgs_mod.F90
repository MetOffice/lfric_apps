! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
module c_cldsgs_mod

use um_types, only: real_umphys

implicit none


! Minimum ice content to perform calculations
real(kind=real_umphys),parameter:: qcfmin=1.0e-8

! end C_CLDSGS

end module c_cldsgs_mod
