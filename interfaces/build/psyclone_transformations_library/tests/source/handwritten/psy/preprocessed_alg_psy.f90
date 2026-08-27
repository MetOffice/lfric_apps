  MODULE preprocessed_alg_mod_psy
    USE constants_mod, ONLY: r_def, i_def
    USE field_mod, ONLY: field_type, field_proxy_type
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE invoke_0(one)
      implicit none
      TYPE(field_type), intent(in) :: one
      INTEGER :: df
      TYPE(field_proxy_type) :: one_proxy
      !
      ! Initialise field and/or operator proxies
      !
      one_proxy = one%get_proxy()
      !
      ! Call kernels and communication routines
      !
      DO df=1,one_proxy%vspace%get_last_dof_annexed()
        one_proxy%data(df) = 3.0_r_def
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL one_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_0
  END MODULE preprocessed_alg_mod_psy
