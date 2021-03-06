!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines useful for iterative matrix calculations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE iterate_matrix
   USE arnoldi_api,                     ONLY: arnoldi_data_type,&
                                              arnoldi_ev,&
                                              arnoldi_extremal,&
                                              deallocate_arnoldi_data,&
                                              get_selected_ritz_val,&
                                              setup_arnoldi_data
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_filter, &
        dbcsr_frobenius_norm, dbcsr_gershgorin_norm, dbcsr_get_info, dbcsr_get_matrix_type, &
        dbcsr_get_occupation, dbcsr_multiply, dbcsr_norm, dbcsr_norm_maxabsnorm, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_trace, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: ifac
   USE mathlib,                         ONLY: abnormal_value
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'iterate_matrix'

   INTERFACE purify_mcweeny
      MODULE PROCEDURE purify_mcweeny_orth, purify_mcweeny_nonorth
   END INTERFACE

   PUBLIC :: invert_Hotelling, matrix_sign_Newton_Schulz, matrix_sqrt_Newton_Schulz, &
             purify_mcweeny

CONTAINS

! **************************************************************************************************
!> \brief invert a symmetric positive definite matrix by Hotelling's method
!>        explicit symmetrization makes this code not suitable for other matrix types
!>        Currently a bit messy with the options, to to be cleaned soon
!> \param matrix_inverse ...
!> \param matrix ...
!> \param threshold convergence threshold nased on the max abs
!> \param use_inv_as_guess logical whether input can be used as guess for inverse
!> \param norm_convergence convergence threshold for the 2-norm, useful for approximate solutions
!> \param filter_eps filter_eps for matrix multiplications, if not passed nothing is filteres
!> \param silent ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2011.10 guess option added [Rustam Z Khaliullin]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE invert_Hotelling(matrix_inverse, matrix, threshold, use_inv_as_guess, &
                               norm_convergence, filter_eps, silent)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix_inverse, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_inv_as_guess
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: norm_convergence, filter_eps
      LOGICAL, INTENT(IN), OPTIONAL                      :: silent

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_Hotelling', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, nrow, unit_nr
      INTEGER(KIND=int_8)                                :: flop1, flop2
      LOGICAL                                            :: use_inv_guess
      REAL(KIND=dp)                                      :: convergence, frob_matrix, &
                                                            gershgorin_norm, max_ev, &
                                                            maxnorm_matrix, min_ev, occ_matrix, &
                                                            t1, t2
      TYPE(arnoldi_data_type)                            :: my_arnoldi
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(1)                   :: mymat
      TYPE(dbcsr_type), TARGET                           :: tmp1, tmp2

! turn this off for the time being

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF
      IF (PRESENT(silent)) THEN
         IF (silent) unit_nr = -1
      END IF

      convergence = threshold
      IF (PRESENT(norm_convergence)) convergence = norm_convergence

      use_inv_guess = .FALSE.
      IF (PRESENT(use_inv_as_guess)) use_inv_guess = use_inv_as_guess
      IF (.NOT. use_inv_guess) THEN
         ! initialize matrix to unity and use arnoldi to scale it into the convergence range
         gershgorin_norm = dbcsr_gershgorin_norm(matrix)
         frob_matrix = dbcsr_frobenius_norm(matrix)
         CALL dbcsr_set(matrix_inverse, 0.0_dp)
         CALL dbcsr_add_on_diag(matrix_inverse, 1.0_dp)
         ! everything commutes, therefor our all products will be symmetric
         CALL dbcsr_create(tmp1, template=matrix_inverse)
      ELSE
         ! It is unlikely that our guess will commute with the matrix, therefore the first product will
         ! be non symmetric
         CALL dbcsr_create(tmp1, template=matrix_inverse, matrix_type=dbcsr_type_no_symmetry)
      ENDIF

      CALL dbcsr_get_info(matrix, nfullrows_total=nrow)
      CALL dbcsr_create(tmp2, template=matrix_inverse)

      IF (unit_nr > 0) WRITE (unit_nr, *)

      ! scale the approximate inverse to be within the convergence radius
      t1 = m_walltime()
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, &
                          0.0_dp, tmp1, flop=flop1, filter_eps=filter_eps)

      mymat(1)%matrix => tmp1
      CALL setup_arnoldi_data(my_arnoldi, mymat, max_iter=30, threshold=1.0E-3_dp, selection_crit=1, &
                              nval_request=2, nrestarts=2, generalized_ev=.FALSE., iram=.TRUE.)
      CALL arnoldi_ev(mymat, my_arnoldi)
      max_eV = REAL(get_selected_ritz_val(my_arnoldi, 2), dp)
      min_eV = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
      CALL deallocate_arnoldi_data(my_arnoldi)

      occ_matrix = dbcsr_get_occupation(matrix_inverse)
      ! 2.0 would be the correct scaling howver, we should make sure here, that we are in the convergence radius
      CALL dbcsr_scale(tmp1, 1.9_dp/(min_ev+max_ev))
      CALL dbcsr_scale(matrix_inverse, 1.9_dp/(min_ev+max_ev))
      min_ev = min_ev*1.9_dp/(min_ev+max_ev)

      DO i = 1, 100

         ! tmp1 = S^-1 S

         ! for the convergence check
         !frob_matrix_base=dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)

         CALL dbcsr_norm(tmp1, &
                         dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix)

         CALL dbcsr_add_on_diag(tmp1, +1.0_dp)

         ! tmp2 = S^-1 S S^-1
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_inverse, 0.0_dp, tmp2, &
                             flop=flop2, filter_eps=filter_eps)
         ! S^-1_{n+1} = 2 S^-1 - S^-1 S S^-1
         CALL dbcsr_add(matrix_inverse, tmp2, 2.0_dp, -1.0_dp)

         CALL dbcsr_filter(matrix_inverse, threshold)
         t2 = m_walltime()
         occ_matrix = dbcsr_get_occupation(matrix_inverse)

         ! use the scalar form of the algorithm to trace the EV
         min_ev = min_ev*(2.0_dp-min_ev)
         IF (PRESENT(norm_convergence)) maxnorm_matrix = ABS(min_eV-1.0_dp)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "Hotelling iter", i, occ_matrix, &
               maxnorm_matrix, t2-t1, &
               (flop1+flop2)/(1.0E6_dp*MAX(0.001_dp, t2-t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (maxnorm_matrix < convergence) EXIT

         ! scale the matrix for improved convergence
         min_ev = min_ev*2.0_dp/(min_ev+1.0_dp)
         CALL dbcsr_scale(matrix_inverse, 2.0_dp/(min_ev+1.0_dp))

         t1 = m_walltime()
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, &
                             0.0_dp, tmp1, flop=flop1, filter_eps=filter_eps)

      ENDDO

      ! try to symmetrize the output matrix
      IF (dbcsr_get_matrix_type(matrix_inverse) == dbcsr_type_no_symmetry) THEN
         CALL dbcsr_transposed(tmp2, matrix_inverse)
         CALL dbcsr_add(matrix_inverse, tmp2, 0.5_dp, 0.5_dp)
      END IF

      IF (unit_nr > 0) THEN
!           WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final Hotelling ",i,occ_matrix,&
!              !frob_matrix/frob_matrix_base
!              maxnorm_matrix
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)

      CALL timestop(handle)

   END SUBROUTINE invert_Hotelling

! **************************************************************************************************
!> \brief compute the sign a matrix using Newton-Schulz iterations
!> \param matrix_sign ...
!> \param matrix ...
!> \param threshold ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE matrix_sign_Newton_Schulz(matrix_sign, matrix, threshold)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sign, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_Newton_Schulz', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: count, handle, i, unit_nr
      INTEGER(KIND=int_8)                                :: flop1, flop2
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base, &
                                                            gersh_matrix, occ_matrix, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tmp1, tmp2

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(tmp1, template=matrix_sign)

      CALL dbcsr_create(tmp2, template=matrix_sign)

      CALL dbcsr_copy(matrix_sign, matrix)
      CALL dbcsr_filter(matrix_sign, threshold)

      ! scale the matrix to get into the convergence range
      frob_matrix = dbcsr_frobenius_norm(matrix_sign)
      gersh_matrix = dbcsr_gershgorin_norm(matrix_sign)
      CALL dbcsr_scale(matrix_sign, 1/MIN(frob_matrix, gersh_matrix))

      IF (unit_nr > 0) WRITE (unit_nr, *)

      count = 0
      DO i = 1, 100

         t1 = m_walltime()
         ! tmp1 = X * X
         CALL dbcsr_multiply("N", "N", -1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1, &
                             filter_eps=threshold, flop=flop1)

         ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, +1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)

         ! update the above to 3*I-X*X
         CALL dbcsr_add_on_diag(tmp1, +2.0_dp)
         occ_matrix = dbcsr_get_occupation(matrix_sign)

         ! tmp2 = 0.5 * X * (3*I-X*X)
         CALL dbcsr_multiply("N", "N", 0.5_dp, matrix_sign, tmp1, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flop2)

         ! done iterating
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sign, tmp2)
         t2 = m_walltime()

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sign iter ", i, occ_matrix, &
               frob_matrix/frob_matrix_base, t2-t1, &
               (flop1+flop2)/(1.0E6_dp*MAX(0.001_dp, t2-t1))
            CALL m_flush(unit_nr)
         ENDIF

         ! frob_matrix/frob_matrix_base < SQRT(threshold)
         IF (frob_matrix*frob_matrix < (threshold*frob_matrix_base*frob_matrix_base)) EXIT

      ENDDO

      ! this check is not really needed
      CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1, &
                          filter_eps=threshold)
      frob_matrix_base = dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(tmp1)
      occ_matrix = dbcsr_get_occupation(matrix_sign)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sign iter", i, occ_matrix, &
            frob_matrix/frob_matrix_base
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)

      CALL timestop(handle)

   END SUBROUTINE matrix_sign_Newton_Schulz

! **************************************************************************************************
!> \brief compute the sqrt of a matrix via the sign function and the corresponding Newton-Schulz iterations
!>        the order of the algorithm should be 2..5, 3 or 5 is recommended
!> \param matrix_sqrt ...
!> \param matrix_sqrt_inv ...
!> \param matrix ...
!> \param threshold ...
!> \param order ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param converged ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE matrix_sqrt_Newton_Schulz(matrix_sqrt, matrix_sqrt_inv, matrix, threshold, order, &
                                        eps_lanczos, max_iter_lanczos, converged)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_sqrt, matrix_sqrt_inv, matrix
      REAL(KIND=dp), INTENT(IN)                          :: threshold
      INTEGER, INTENT(IN)                                :: order
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos
      LOGICAL, OPTIONAL                                  :: converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sqrt_Newton_Schulz', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, unit_nr
      INTEGER(KIND=int_8)                                :: flop1, flop2, flop3, flop4, flop5
      LOGICAL                                            :: arnoldi_converged
      REAL(KIND=dp)                                      :: a, b, c, conv, d, frob_matrix, &
                                                            frob_matrix_base, max_ev, min_ev, oa, &
                                                            ob, oc, occ_matrix, od, scaling, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: tmp1, tmp2, tmp3

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (PRESENT(converged)) converged = .FALSE.

      ! for stability symmetry can not be assumed
      CALL dbcsr_create(tmp1, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(tmp2, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      IF (order .GE. 4) THEN
         CALL dbcsr_create(tmp3, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      ENDIF

      CALL dbcsr_set(matrix_sqrt_inv, 0.0_dp)
      CALL dbcsr_add_on_diag(matrix_sqrt_inv, 1.0_dp)
      CALL dbcsr_filter(matrix_sqrt_inv, threshold)
      CALL dbcsr_copy(matrix_sqrt, matrix)

      ! scale the matrix to get into the convergence range
      CALL arnoldi_extremal(matrix_sqrt, max_ev, min_ev, threshold=eps_lanczos, &
                            max_iter=max_iter_lanczos, converged=arnoldi_converged)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T6,A,1X,L1,A,E12.3)') "Lanczos converged: ", arnoldi_converged, " threshold:", eps_lanczos
         WRITE (unit_nr, '(T6,A,1X,E12.3,E12.3)') "Est. extremal eigenvalues:", max_ev, min_ev
         WRITE (unit_nr, '(T6,A,1X,E12.3)') "Est. condition number :", max_ev/MAX(min_ev, EPSILON(min_ev))
      ENDIF
      ! conservatively assume we get a relatively large error (100*threshold_lanczos) in the estimates
      ! and adjust the scaling to be on the safe side
      scaling = 2/(max_ev+min_ev+100*eps_lanczos)

      CALL dbcsr_scale(matrix_sqrt, scaling)
      CALL dbcsr_filter(matrix_sqrt, threshold)

      DO i = 1, 100

         t1 = m_walltime()

         ! tmp1 = Zk * Yk - I
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1, &
                             filter_eps=threshold, flop=flop1)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)

         ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
         frob_matrix = dbcsr_frobenius_norm(tmp1)

         flop4 = 0; flop5 = 0
         SELECT CASE (order)
         CASE (2)
            ! update the above to 0.5*(3*I-Zk*Yk)
            CALL dbcsr_add_on_diag(tmp1, -2.0_dp)
            CALL dbcsr_scale(tmp1, -0.5_dp)
         CASE (3)
            ! tmp2 = tmp1 ** 2
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            ! tmp1 = 1/16 * (16*I-8*tmp1+6*tmp1**2-5*tmp1**3)
            CALL dbcsr_add(tmp1, tmp2, -4.0_dp, 3.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 8.0_dp)
            CALL dbcsr_scale(tmp1, 0.125_dp)
         CASE (4) ! as expensive as case(5), so little need to use it
            ! tmp2 = tmp1 ** 2
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            ! tmp3 = tmp2 * tmp1
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp1, 0.0_dp, tmp3, &
                                filter_eps=threshold, flop=flop5)
            CALL dbcsr_scale(tmp1, -8.0_dp)
            CALL dbcsr_add_on_diag(tmp1, 16.0_dp)
            CALL dbcsr_add(tmp1, tmp2, 1.0_dp, 6.0_dp)
            CALL dbcsr_add(tmp1, tmp3, 1.0_dp, -5.0_dp)
            CALL dbcsr_scale(tmp1, 1/16.0_dp)
         CASE (5)
            ! Knuth's reformulation to evaluate the polynomial of 4th degree in 2 multiplications
            ! p = y4+A*y3+B*y2+C*y+D
            ! z := y * (y+a); P := (z+y+b) * (z+c) + d.
            ! a=(A-1)/2 ; b=B*(a+1)-C-a*(a+1)*(a+1)
            ! c=B-b-a*(a+1)
            ! d=D-bc
            oa = -40.0_dp/35.0_dp
            ob = 48.0_dp/35.0_dp
            oc = -64.0_dp/35.0_dp
            od = 128.0_dp/35.0_dp
            a = (oa-1)/2
            b = ob*(a+1)-oc-a*(a+1)**2
            c = ob-b-a*(a+1)
            d = od-b*c
            ! tmp2 = tmp1 ** 2 + a * tmp1
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2, &
                                filter_eps=threshold, flop=flop4)
            CALL dbcsr_add(tmp2, tmp1, 1.0_dp, a)
            ! tmp3 = tmp2 + tmp1 + b
            CALL dbcsr_copy(tmp3, tmp2)
            CALL dbcsr_add(tmp3, tmp1, 1.0_dp, 1.0_dp)
            CALL dbcsr_add_on_diag(tmp3, b)
            ! tmp2 = tmp2 + c
            CALL dbcsr_add_on_diag(tmp2, c)
            ! tmp1 = tmp2 * tmp3
            CALL dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp1, &
                                filter_eps=threshold, flop=flop5)
            ! tmp1 = tmp1 + d
            CALL dbcsr_add_on_diag(tmp1, d)
            ! final scale
            CALL dbcsr_scale(tmp1, 35.0_dp/128.0_dp)
         CASE DEFAULT
            CPABORT("")
         END SELECT

         ! tmp2 = Yk * tmp1 = Y(k+1)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt, tmp1, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flop2)
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sqrt, tmp2)

         ! tmp2 = tmp1 * Zk = Z(k+1)
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_sqrt_inv, 0.0_dp, tmp2, &
                             filter_eps=threshold, flop=flop3)
         ! CALL dbcsr_filter(tmp2,threshold)
         CALL dbcsr_copy(matrix_sqrt_inv, tmp2)

         occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)

         ! done iterating
         t2 = m_walltime()

         conv = frob_matrix/frob_matrix_base

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sqrt iter ", i, occ_matrix, &
               conv, t2-t1, &
               (flop1+flop2+flop3+flop4+flop5)/(1.0E6_dp*MAX(0.001_dp, t2-t1))
            CALL m_flush(unit_nr)
         ENDIF

         IF (abnormal_value(conv)) &
            CPABORT("conv is an abnormal value (NaN/Inf).")

         ! conv < SQRT(threshold)
         IF ((conv*conv) < threshold) THEN
            IF (PRESENT(converged)) converged = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! symmetrize the matrices as this is not guaranteed by the algorithm
      CALL dbcsr_transposed(tmp1, matrix_sqrt_inv)
      CALL dbcsr_add(matrix_sqrt_inv, tmp1, 0.5_dp, 0.5_dp)
      CALL dbcsr_transposed(tmp1, matrix_sqrt)
      CALL dbcsr_add(matrix_sqrt, tmp1, 0.5_dp, 0.5_dp)

      ! this check is not really needed
      CALL dbcsr_multiply("N", "N", +1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1, &
                          filter_eps=threshold)
      frob_matrix_base = dbcsr_frobenius_norm(tmp1)
      CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(tmp1)
      occ_matrix = dbcsr_get_occupation(matrix_sqrt_inv)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sqrt iter ", i, occ_matrix, &
            frob_matrix/frob_matrix_base
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      ! scale to proper end results
      CALL dbcsr_scale(matrix_sqrt, 1/SQRT(scaling))
      CALL dbcsr_scale(matrix_sqrt_inv, SQRT(scaling))

      CALL dbcsr_release(tmp1)
      CALL dbcsr_release(tmp2)
      IF (order .GE. 4) THEN
         CALL dbcsr_release(tmp3)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE matrix_sqrt_Newton_Schulz

! **************************************************************************************************
!> \brief ...
!> \param matrix_exp ...
!> \param matrix ...
!> \param omega ...
!> \param alpha ...
!> \param threshold ...
! **************************************************************************************************
   SUBROUTINE matrix_exponential(matrix_exp, matrix, omega, alpha, threshold)
      ! compute matrix_exp=omega*exp(alpha*matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_exp, matrix
      REAL(KIND=dp), INTENT(IN)                          :: omega, alpha, threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_exponential', &
         routineP = moduleN//':'//routineN
      REAL(dp), PARAMETER                                :: one = 1.0_dp, toll = 1.E-17_dp, &
                                                            zero = 0.0_dp

      INTEGER                                            :: handle, i, k, unit_nr
      REAL(dp)                                           :: factorial, norm_C, norm_D, norm_scalar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: B, B_square, C, D, D_product

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! Calculate the norm of the matrix alpha*matrix, and scale it until it is less than 1.0
      norm_scalar = ABS(alpha)*dbcsr_frobenius_norm(matrix)

      ! k=scaling parameter
      k = 1
      DO
         IF ((norm_scalar/2.0_dp**k) <= one) EXIT
         k = k+1
      END DO

      ! copy and scale the input matrix in matrix C and in matrix D
      CALL dbcsr_create(C, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(C, matrix)
      CALL dbcsr_scale(C, alpha_scalar=alpha/2.0_dp**k)

      CALL dbcsr_create(D, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(D, C)

      !   write(*,*)
      !   write(*,*)
      !   CALL dbcsr_print(D, nodata=.FALSE., matlab_format=.TRUE., variable_name="D", unit_nr=6)

      ! set the B matrix as B=Identity+D
      CALL dbcsr_create(B, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(B, D)
      CALL dbcsr_add_on_diag(B, alpha_scalar=one)

      !   CALL dbcsr_print(B, nodata=.FALSE., matlab_format=.TRUE., variable_name="B", unit_nr=6)

      ! Calculate the norm of C and moltiply by toll to be used as a threshold
      norm_C = toll*dbcsr_frobenius_norm(matrix)

      ! iteration for the trucated taylor series expansion
      CALL dbcsr_create(D_product, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      i = 1
      DO
         i = i+1
         ! compute D_product=D*C
         CALL dbcsr_multiply("N", "N", one, D, C, &
                             zero, D_product, filter_eps=threshold)

         ! copy D_product in D
         CALL dbcsr_copy(D, D_product)

         ! calculate B=B+D_product/fat(i)
         factorial = ifac(i)
         CALL dbcsr_add(B, D_product, one, factorial)

         ! check for convergence using the norm of D (copy of the matrix D_product) and C
         norm_D = factorial*dbcsr_frobenius_norm(D)
         IF (norm_D < norm_C) EXIT
      END DO

      ! start the k iteration for the squaring of the matrix
      CALL dbcsr_create(B_square, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      DO i = 1, k
         !compute B_square=B*B
         CALL dbcsr_multiply("N", "N", one, B, B, &
                             zero, B_square, filter_eps=threshold)
         ! copy Bsquare in B to iterate
         CALL dbcsr_copy(B, B_square)
      END DO

      ! copy B_square in matrix_exp and
      CALL dbcsr_copy(matrix_exp, B_square)

      ! scale matrix_exp by omega, matrix_exp=omega*B_square
      CALL dbcsr_scale(matrix_exp, alpha_scalar=omega)
      ! write(*,*) alpha,omega

      CALL timestop(handle)

   END SUBROUTINE matrix_exponential

! **************************************************************************************************
!> \brief McWeeny purification of a matrix in the orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE purify_mcweeny_orth(matrix_p, threshold, max_steps)
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_orth', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_pp, matrix_tmp

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_pp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_trace(matrix_p(1), trace)

      DO ispin = 1, SIZE(matrix_p)
         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_p(ispin), &
                                0.0_dp, matrix_pp, filter_eps=threshold)

            ! test convergence
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_add(matrix_tmp, matrix_p(ispin), 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_tmp) ! tmp = PP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_pp, matrix_p(ispin), &
                                3.0_dp, matrix_tmp, filter_eps=threshold)
            CALL dbcsr_copy(matrix_p(ispin), matrix_tmp) ! tmp = 3PP - 2PPP

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO
      END DO

      CALL dbcsr_release(matrix_pp)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)
   END SUBROUTINE purify_mcweeny_orth

! **************************************************************************************************
!> \brief McWeeny purification of a matrix in the non-orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param matrix_s Overlap-Matrix
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE purify_mcweeny_nonorth(matrix_p, matrix_s, threshold, max_steps)
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p
      TYPE(dbcsr_type)                                   :: matrix_s
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_nonorth', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_ps, matrix_psp, matrix_test

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_ps, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_psp, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_test, template=matrix_p(1), matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, SIZE(matrix_p)
         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_s, &
                                0.0_dp, matrix_ps, filter_eps=threshold)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p(ispin), &
                                0.0_dp, matrix_psp, filter_eps=threshold)
            IF (i == 1) CALL dbcsr_trace(matrix_ps, trace)

            ! test convergence
            CALL dbcsr_copy(matrix_test, matrix_psp)
            CALL dbcsr_add(matrix_test, matrix_p(ispin), 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_test) ! test = PSP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,2f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_p(ispin), matrix_psp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_ps, matrix_psp, &
                                3.0_dp, matrix_p(ispin), filter_eps=threshold)

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO
      END DO

      CALL dbcsr_release(matrix_ps)
      CALL dbcsr_release(matrix_psp)
      CALL dbcsr_release(matrix_test)
      CALL timestop(handle)
   END SUBROUTINE purify_mcweeny_nonorth

END MODULE iterate_matrix
