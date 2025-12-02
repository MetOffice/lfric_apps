!----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief UKCA Photolysis parameter module.
!!        Currently only contains hardwired photolysis rates (factors).

module ukca_photol_param_mod

  use constants_mod, only : r_def, i_um

  ! UKCA API module
  use ukca_api_mod,  only : ukca_photol_varname_len

  implicit none

  public

  ! Number of species in photolysis - has to match UKCA definitions
  integer(kind=i_um), parameter :: jppj = 60

  ! Set up photolytically active species and scaling factor for each species
  !  w.r.t the 'standard' (i.e. highest) rate found in an example FastJX run.

  ! Names of species - should match ratj reactions in ukca_chem_master
  !
  ! NOTE: Some reaction names are repeated and it is not trivial to separate
  ! out the rates for each. Hence specifying the same rate for each copy
  ! assuming they represent similar molecules.

  character(len=ukca_photol_varname_len), parameter, dimension(jppj)::         &
              !   1           2(R)        3           4           5
   jp_names = ['jaceta   ','jaceta   ','jacetb   ','jaceto   ','jbrcl    ',    &
              !   6           7           8           9           10
               'jbrna    ','jbrnb    ','jbro     ','jcfc2    ','jcfcl3   ',    &
              !   11          12          13          14          15
               'jch4     ', 'jcl2o2   ','jclna    ','jclnb    ','jco2     ',   &
              !   16          17          18          19          20
               'jcos     ','jcs2     ','jetcho   ','jh2o     ','jh2o2    ',    &
              !   21          22          23          24          25
               'jh2so4   ','jhacet   ','jhchoa   ','jhchob   ','jhcl     ',    &
              !   26          27          28          29          30
               'jhobr    ','jhocl    ','jhono    ','jhono2   ','jiprn    ',    &
              !   31          32          33(R)       34          35
               'jmacr    ','jmacro   ','jmacro   ','jmebr    ','jmeco3h  ',    &
              !   36          37          38(R)       39(R)       40(R)
               'jmena    ','jmhp     ','jmhp     ','jmhp     ','jmhp     ',    &
              !   41(R)       42(R)       43          44          45
               'jmhp     ','jmhp     ','jmkal    ','jn2o     ','jn2o5    ',    &
              !   46          47          48          49          50
               'jno      ','jno2     ','jno3a    ','jno3b    ','jo2      ',    &
              !   51          52          53          54          55
               'jo2b     ','jo3a     ','jo3b     ','joclo    ','jpan     ',    &
              !   56(R)       57(R)       58          59          60
               'jpan     ','jpan     ','jpna33   ','jpna67   ','jso3     ']

  real(kind=r_def), parameter, dimension(jppj):: jrate_fac = [                  &
   !     1              2(R)           3              4              5
   2.246e-05_r_def,2.246e-05_r_def,1.738e-07_r_def,2.485e-05_r_def,3.031e-02_r_def, &
   !     6              7              8              9              10
   1.000e+00_r_def,1.762e-01_r_def,1.192e-01_r_def,3.638e-07_r_def,3.431e-06_r_def, &
   !     11             12             13             14             15
   0.000e+00_r_def,4.369e-03_r_def,1.108e-04_r_def,5.469e-05_r_def,1.677e-11_r_def, &
   !     16             17             18             19             20
   5.838e-06_r_def,5.038e-04_r_def,1.515e-04_r_def,5.362e-10_r_def,2.269e-05_r_def, &
   !     21             22             23             24             25
   0.000e+00_r_def,2.515e-05_r_def,1.015e-04_r_def,1.538e-04_r_def,1.331e-07_r_def, &
   !     26             27             28             29             30
   6.700e-03_r_def,8.692e-04_r_def,5.792e-03_r_def,3.038e-05_r_def,9.231e-05_r_def, &
   !     31             32             33(R)          34             35
   1.285e-05_r_def,8.385e-06_r_def,8.385e-06_r_def,1.246e-05_r_def,7.023e-06_r_def, &
   !     36             37             38(R)          39(R)          40(R)
   4.862e-06_r_def,1.677e-05_r_def,1.677e-05_r_def,1.677e-05_r_def,1.677e-05_r_def, &
   !     41(R)          42(R)          43             44             45
   1.677e-05_r_def,1.677e-05_r_def,2.569e-03_r_def,1.746e-07_r_def,1.731e-04_r_def, &
   !     46             47             48             49             50
   1.654e-07_r_def,2.715e-02_r_def,6.823e-02_r_def,5.069e-01_r_def,1.438e-10_r_def, &
   !     51             52             53             54             55
   0.000e+00_r_def,6.400e-04_r_def,1.277e-03_r_def,2.492e-01_r_def,3.915e-05_r_def, &
   !     56(R)          57(R)          58             59             60
   3.915e-05_r_def,3.915e-05_r_def,1.431e-05_r_def,7.846e-05_r_def,0.000e+00_r_def ]

end module ukca_photol_param_mod
