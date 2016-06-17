!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt4_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_copy,&
                                              cp_dbcsr_deallocate_matrix,&
                                              cp_dbcsr_init_p,&
                                              cp_dbcsr_p_type,&
                                              cp_dbcsr_set,&
                                              cp_dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_schur_product,&
                                              cp_fm_symm,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_invert, cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type, cp_fm_set_element
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: cholesky_dbcsr,&
                                              cholesky_inverse,&
                                              cholesky_off,&
                                              cholesky_restore
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: tddfpt
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt4_methods'

! **************************************************************************************************
!> \brief Single-determinant guess vector descriptor
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_state_descriptor
      !> index of occupied MO to deoccupy
      INTEGER :: imo_occ
      !> index of virtual MO to occupy
      INTEGER :: imo_virt
   END TYPE tddfpt_state_descriptor

! **************************************************************************************************
!> \brief Time-Dependent DFT guess vectors
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_guess_vectors
      !> collection of single-determinant guess vectors (imo_occ -> imo_virt)
      !> ordered by their transition energies
      TYPE(tddfpt_state_descriptor), ALLOCATABLE, DIMENSION(:) :: descrs
      !> transition energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals
   END TYPE tddfpt_guess_vectors

! **************************************************************************************************
!> \brief Time-Dependent DFT action matrix
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
!
!      |----------|-------------|
!      | . d                    |
!      |g  . i      lower_      |
!      |  u  . a    triangle    |
!      |    e  . g    ^T        |
!      |      s  .|-------------|
!      |        s |.            |
!      | lower_   |% .   UNDEF  |
!      | triangle |  e .        |
!      |          |    v .      |
!      |          |      a .    |
!      |          | UNDEF  l .  |
!      |          |          s .|
!      |----------|-------------|
!      <- nvects ->
!
!      <----- nocc * nvirt ----->
! **************************************************************************************************
   TYPE tddfpt_action_matrix
      !> rectangular matrix of dimensions (nocc*nvirt, nvects) which stores
      !> the lower triangular part of the action matrix. All diagonal and
      !> upper-triangular elements are set to zero
      TYPE(cp_fm_type), POINTER                :: lower_triangle
      !> diagonal part of the action matrix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: diag
   END TYPE tddfpt_action_matrix

CONTAINS

   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: ispin = 1
INTEGER, PARAMETER :: max_krylov_vects_param = 2000

      CHARACTER(len=20)                                  :: navail_str, nstates_str
      INTEGER :: handle, nspins_save, log_unit, istate, nstates, nstates_total, nstates_converged, &
                 iter, niters, nvects_new, nvects_exists, nao, nmo_occ, energy_unit, max_krylov_vects
      LOGICAL                                            :: do_hfx
      REAL(kind=dp)                                      :: C_hf, conv, t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_occ, evals_virt, evals_last, evals_prev
      TYPE(cp_fm_type), POINTER                          :: mos_occ, mos_virt, rho_ij_fm
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, rho_ij_ao
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, tddfpt_section, &
                                                            xc_section


      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      TYPE(tddfpt_guess_vectors)                         :: guess_vectors
      TYPE(tddfpt_action_matrix)                         :: action_matrix
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: krylov_vectors

      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: evects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(qs_rho_type), POINTER                         :: rho_ij_struct
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:,:)    :: K_ij_munu
!      TYPE(tddfpt_K_ij_munu_storage)                     :: K_ij_munu
INTEGER :: npoisson, ix, iy

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env, dft_control, input, matrix_ks, matrix_s, mos, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, para_env=para_env, pw_env=pw_env)
      tddfpt_control => dft_control%tddfpt2_control

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requestes for "// &
                      TRIM(nstates_str)//" excited states: nothing to do.")
         CALL timestop(handle)
         RETURN
      END IF

      tddfpt_section => section_vals_get_subs_vals(input, "DFT%TDDFPT2")
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")

      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
         do_hfx = (C_hf /= 0.0_dp)
      END IF

      IF (dft_control%nspins > 1) &
         CPABORT("TDDFPT is currently implemented for restricted molecular orbitals only")

      nspins_save = dft_control%nspins
      dft_control%nspins = 1

      ! generate all virtual molecular orbitals by diagonalising the Kohn-Sham matrix one more time
      !
      ! at this point the scf_env has been partially deallocated;
      ! so do not rely on it and assume scf_env%cholesky_method == cholesky_restore
      NULLIFY (mos_occ, mos_virt)
      CALL tddfpt_lumos(mos_occ=mos_occ, evals_occ=evals_occ, &
                        mos_virt=mos_virt, evals_virt=evals_virt, &
                        mo_set=mos(ispin)%mo_set, matrix_ks=matrix_ks(ispin)%matrix, &
                        matrix_s=matrix_s(1)%matrix, cholesky_method=cholesky_restore, &
                        blacs_env=blacs_env, para_env=para_env)

      ! check we have enough molecular orbitals to contruct the requested number of excited states
      IF (INT(tddfpt_control%nstates, kind=int_8) > SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL integer_to_string(SIZE(evals_occ)*SIZE(evals_virt), navail_str)
         CALL cp_warn(__LOCATION__, TRIM(nstates_str)//" excited states were requested, however only "// &
                      TRIM(navail_str)//" excited states can be constructed using the atomic basis set given.")
         tddfpt_control%nstates = SIZE(evals_occ)*SIZE(evals_virt)
      END IF

      ! adjust the number of trial vectors
      IF (tddfpt_control%added_states < 0 .OR. &
          INT(tddfpt_control%added_states, kind=int_8)+INT(tddfpt_control%nstates, kind=int_8) > &
          SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)) THEN

         tddfpt_control%added_states = INT(SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)- &
                                           INT(tddfpt_control%nstates, kind=int_8))
      END IF
      nstates = tddfpt_control%nstates

      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_guess(guess_vectors=guess_vectors, &
                        nstates=nstates, nstates_active=nstates+tddfpt_control%added_states, &
                        evals_occ=evals_occ, evals_virt=evals_virt, &
                        is_add_degenerate=tddfpt_control%add_degenerate, &
                        degenerate_eps=tddfpt_control%degenerate_eps, &
                        log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_section, "PRINT%GUESS_VECTORS")

      nstates_total = SIZE(guess_vectors%evals)

     ! ++ initialise xc_rho_set and xc_deriv_set (needed for kernel in adiabatic approximation)
      NULLIFY (pw_env)
      CALL get_qs_env(qs_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL cp_fm_get_info(mos_occ, nrow_global=nao, ncol_global=nmo_occ)
      NULLIFY (fm_struct, rho_ij_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ij_fm, fm_struct)

      NULLIFY (rho_ij_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ij_ao, 1)
      CALL cp_dbcsr_init_p(rho_ij_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ij_ao(1)%matrix, matrix_s(1)%matrix)

      NULLIFY (rho_ij_struct, rho_ij_r)
      CALL qs_rho_create(rho_ij_struct)
      CALL qs_rho_set(rho_ij_struct, rho_ao=rho_ij_ao)
      CALL qs_rho_rebuild(rho_ij_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ij_struct, rho_r=rho_ij_r)

      ! alpha == 2.0 due to spin-restricted density
      CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 2.0_dp, mos_occ, mos_occ, 0.0_dp, rho_ij_fm)
      CALL copy_fm_to_dbcsr(rho_ij_fm, rho_ij_ao(1)%matrix, keep_sparsity=.TRUE.)
      CALL qs_rho_update_rho(rho_ij_struct, qs_env)

      NULLIFY (xc_rho_set, xc_deriv_set)

      CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r, auxbas_pw_pool, xc_section=xc_section)

      CALL qs_rho_release(rho_ij_struct)
      CALL cp_fm_release(rho_ij_fm)

      CALL cp_fm_struct_release(fm_struct)

      ! ++ Krylov vectors
      ALLOCATE(evals_last(nstates))
      ALLOCATE(evals_prev(nstates))
      evals_prev(1:nstates) = guess_vectors%evals(1:nstates)

max_krylov_vects = max_krylov_vects_param
IF (max_krylov_vects < nstates) max_krylov_vects = nstates
IF (max_krylov_vects > nstates_total) max_krylov_vects = nstates_total

      CALL cp_fm_get_info(mos_occ, matrix_struct=fm_struct)
      ALLOCATE(evects(max_krylov_vects))
      DO istate = 1,max_krylov_vects
         NULLIFY(evects(istate)%matrix)
         CALL cp_fm_create(evects(istate)%matrix, fm_struct)
      END DO
      NULLIFY(fm_struct)

      ALLOCATE(krylov_vectors(nstates_total, max_krylov_vects+nstates))
      krylov_vectors(:,:) = 0.0_dp
      ! guess first Krylov's vectors
      DO istate = 1,nstates
         krylov_vectors(istate, istate) = 1.0_dp
      END DO

      ! hfx storage
      IF (do_hfx) THEN
         ALLOCATE(K_ij_munu(nmo_occ, nmo_occ))
         DO iter = 1,nmo_occ
            DO istate = 1,nmo_occ
               NULLIFY (K_ij_munu(istate, iter)%matrix)
            END DO
         END DO
      END IF

      ! start the wave-function optimization
      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%ITERATION_INFO", extension=".tddfptLog")
      energy_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%DETAILED_ENERGY", extension=".tddfptLog")

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') 'TDDFT WAVEFUNCTION OPTIMIZATION'
         WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
         WRITE (log_unit, '(1X,50("-"))')
      END IF

      CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

      NULLIFY(action_matrix%lower_triangle)
      nvects_exists = 0
      t1 = m_walltime()
      niters = tddfpt_control%niters

      DO iter = 1,niters
         CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter)

         nvects_new = tddfpt_am_expand(action_matrix, nstates_total=nstates_total, nvects_new=nstates, &
                                    nvects_max=max_krylov_vects, blacs_env=blacs_env, para_env=para_env)
         nvects_exists = nvects_exists + nvects_new

         IF (nvects_new == 0) THEN
            ! maximum number of Krylov vectors has been reached, restart the Davidson eigensolver
            CALL tddfpt_davidson_restart(krylov_vectors=krylov_vectors(:,1:nvects_exists), &
                                         nstates=nstates, action_matrix=action_matrix, &
                                         diag_guess=guess_vectors%evals, &
                                         blacs_env=blacs_env, para_env=para_env)
            CALL tddfpt_am_release(action_matrix)
            nvects_new = tddfpt_am_expand(action_matrix, nstates_total=nstates_total, nvects_new=nstates, &
                                          nvects_max=max_krylov_vects, blacs_env=blacs_env, para_env=para_env)
            nvects_exists = nvects_new
            IF (log_unit > 0) &
               WRITE(log_unit, '(1X,10("-"),1X,A,1X,11("-"))') "Restart Davidson iterations"
         END IF

         CALL tddfpt_assemble_evects(evects=evects(nvects_exists-nvects_new+1:nvects_exists), &
                                     krylov_vectors=krylov_vectors(:,nvects_exists-nvects_new+1:nvects_exists), &
                                     guess_vectors=guess_vectors, mos_virt=mos_virt)
CALL print_densities(log_unit, nvects_exists-nvects_new+1, evects(nvects_exists-nvects_new+1:nvects_exists), &
     mos_occ, qs_env, matrix_s(1)%matrix)

         CALL tddfpt_am_compute(action_matrix=action_matrix, nvects_new=nvects_new, &
                                evects=evects(1:nvects_exists), guess_vectors=guess_vectors, &
                                mos_occ=mos_occ, evals_occ=evals_occ, mos_virt=mos_virt, &
                                C_hf=C_hf, K_ij_munu=K_ij_munu, &
                                qs_env=qs_env, matrix_ks=matrix_ks(1)%matrix, matrix_s=matrix_s(1)%matrix, &
                                xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, xc_section=xc_section)


         CALL tddfpt_davidson_iteration(residual_vectors=krylov_vectors(:,nvects_exists+1:nvects_exists+nstates), &
                                        evals=evals_last, action_matrix=action_matrix, &
                                        diag_guess=guess_vectors%evals, &
                                        krylov_vectors=krylov_vectors(:,1:nvects_exists), &
                                        blacs_env=blacs_env, para_env=para_env)

         t2 = m_walltime()

         IF (energy_unit > 0) THEN
            WRITE (energy_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
            DO istate = 1, nstates
               WRITE (energy_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
                  evals_last(istate)*evolt, (evals_last(istate)-evals_prev(istate))*evolt
            END DO
            WRITE (energy_unit, *)
            CALL m_flush(energy_unit)

 IF (do_hfx) THEN
 npoisson = 0
 DO ix = 1,nmo_occ
    DO iy = 1,nmo_occ
       IF (ASSOCIATED(K_ij_munu(iy, ix)%matrix)) npoisson = npoisson + 1
    END DO
 END DO
 WRITE(energy_unit, '(1X,A,T70,I8,/)') "Total number of Poisson equations solved to compute (i j || mu nu)", npoisson
 END IF
         END IF

         conv = MAXVAL(ABS(evals_last(1:nstates)-evals_prev(1:nstates)))

         IF (log_unit > 0) THEN
            nstates_converged = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate)-evals_prev(istate)) <= tddfpt_control%conv) &
                  nstates_converged = nstates_converged+1
            END DO

            WRITE (log_unit, '(1X,I8,T12,F7.1,T24,ES11.4,T42,I8)') iter, t2-t1, conv, nstates_converged
            CALL m_flush(log_unit)
         END IF

         t1 = t2
         evals_prev(1:nstates) = evals_last(1:nstates)

         IF (conv <= tddfpt_control%conv) THEN
            CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter)
            EXIT
         END IF
      END DO

      CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

      IF (log_unit > 0) THEN
         IF (iter <= niters) THEN
            CALL integer_to_string(iter, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFT run converged in "//TRIM(nstates_str)//" iteration(s) ***"
         ELSE
            CALL integer_to_string(iter-1, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ***"
         END IF
      END IF

      CALL cp_print_key_finished_output(energy_unit, logger, tddfpt_section, "PRINT%DETAILED_ENERGY")
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_section, "PRINT%ITERATION_INFO")


      CALL tddfpt_davidson_restart(krylov_vectors=krylov_vectors(:,1:nvects_exists), &
                                   nstates=nstates, action_matrix=action_matrix, &
                                   diag_guess=guess_vectors%evals, &
                                   blacs_env=blacs_env, para_env=para_env)
      CALL tddfpt_assemble_evects(evects=evects(1:nstates), &
                                  krylov_vectors=krylov_vectors(:,1:nstates), &
                                  guess_vectors=guess_vectors, mos_virt=mos_virt)

      log_unit = cp_logger_get_default_io_unit()
      CALL tddfpt_print_summary(log_unit=log_unit, nstates=nstates, &
                                evects=evects(1:nstates), evals=evals_last, nmo_occ_total=nmo_occ, &
                                mos_occ=mos_occ, evals_occ=evals_occ, mos_virt=mos_virt, evals_virt=evals_virt, &
                                matrix_s=matrix_s, min_amplitude=tddfpt_control%min_excitation_amplitude)

      IF (ALLOCATED(K_ij_munu)) THEN
         DO iter = 1,nmo_occ
            DO istate = 1,nmo_occ
               IF (ASSOCIATED(K_ij_munu(istate, iter)%matrix)) &
                  CALL cp_fm_release(K_ij_munu(istate, iter)%matrix)
            END DO
         END DO
         DEALLOCATE(K_ij_munu)
      END IF

      DEALLOCATE(krylov_vectors)
      DO istate = SIZE(evects),1,-1
          CALL cp_fm_release(evects(istate)%matrix)
      END DO
      DEALLOCATE(evects)

      DEALLOCATE(evals_last)
      CALL tddfpt_am_release(action_matrix)

      IF (ALLOCATED(guess_vectors%descrs)) &
         DEALLOCATE(guess_vectors%descrs)
      IF (ALLOCATED(guess_vectors%evals)) &
         DEALLOCATE(guess_vectors%evals)

      CALL xc_dset_release(xc_deriv_set)
      CALL xc_rho_set_release(xc_rho_set)

      CALL cp_fm_release(mos_virt)
      CALL cp_fm_release(mos_occ)
      DEALLOCATE (evals_virt)
      DEALLOCATE (evals_occ)

      dft_control%nspins = nspins_save
      CALL timestop(handle)
   END SUBROUTINE tddfpt

   SUBROUTINE tddfpt_am_compute(action_matrix, nvects_new, evects, guess_vectors, mos_occ, evals_occ, mos_virt, &
                               C_hf, K_ij_munu, qs_env, matrix_ks, matrix_s, xc_rho_set, xc_deriv_set, xc_section)
      TYPE(tddfpt_action_matrix), INTENT(inout)          :: action_matrix
      INTEGER, INTENT(in) :: nvects_new
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
!      REAL(kind=dp), DIMENSION(:,:), INTENT(in)          :: krylov_vectors
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), INTENT(in)                          :: C_hf
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: K_ij_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section

      INTEGER :: nvects_exists, nstates_total, nvects_old, nvects_remain
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: block_nn, block_gn, block_on

      nvects_exists = SIZE(evects)
      nstates_total = SIZE(guess_vectors%evals)

      CPASSERT(SIZE(action_matrix%diag) == nvects_exists)
      nvects_old    = nvects_exists - nvects_new
      nvects_remain = nstates_total - nvects_exists

      CPASSERT(nvects_new > 0)
      CPASSERT(nvects_old >= 0)
      CPASSERT(nvects_remain >= 0)

      ALLOCATE(block_nn(nvects_new, nvects_new))
      block_nn(:,:) = 0.0_dp

      IF (nvects_old > 0) THEN
         ALLOCATE(block_on(nvects_old, nvects_new))
         block_on(:,:) = 0.0_dp
      END IF

      IF (nvects_remain > 0) THEN
         ALLOCATE(block_gn(nvects_remain, nvects_new))
         block_gn(:,:) = 0.0_dp
      END IF

      CALL tddfpt_init_by_energy_diff(block_nn=block_nn, &
                                      evects=evects(nvects_exists-nvects_new+1:nvects_exists), evals_occ=evals_occ, &
                                      matrix_ks=matrix_ks, matrix_s=matrix_s)

      ! add contributions from coulomb terms
      ! alpha == 2 due to singlet state
      CALL tddfpt_apply_coulomb(block_on=block_on, block_nn=block_nn, block_gn=block_gn, alpha=2.0_dp, &
                                evects=evects, mos_occ=mos_occ, guess_vectors=guess_vectors, &
                                mos_virt=mos_virt, qs_env=qs_env, matrix_s=matrix_s)

      ! add contributions from the adiabatic TDDFT kernel
      ! alpha == 2 due to relation K_alpha,alpha == K_alpha,beta for singlet states
      CALL tddfpt_apply_kernel(block_on=block_on, block_nn=block_nn, block_gn=block_gn, alpha=2.0_dp, &
                               evects=evects, mos_occ=mos_occ, guess_vectors=guess_vectors, &
                               mos_virt=mos_virt, xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, &
                               qs_env=qs_env, matrix_s=matrix_s, xc_section=xc_section)

      IF (ALLOCATED(K_ij_munu)) &
         CALL tddfpt_apply_hfx(block_on=block_on, block_nn=block_nn, K_ij_munu=K_ij_munu, alpha=C_hf, &
                            evects=evects, mos_occ=mos_occ, qs_env=qs_env, matrix_s=matrix_s)

      CALL tddfpt_am_update(action_matrix=action_matrix, block_on=block_on, block_nn=block_nn, block_gn=block_gn)

      IF (ALLOCATED(block_gn)) &
         DEALLOCATE(block_gn)
      IF (ALLOCATED(block_on)) &
         DEALLOCATE(block_on)
      DEALLOCATE(block_nn)
   END SUBROUTINE tddfpt_am_compute

   SUBROUTINE tddfpt_davidson_iteration(residual_vectors, evals, action_matrix, diag_guess, &
                                        krylov_vectors, blacs_env, para_env)
      REAL(kind=dp), DIMENSION(:,:), INTENT(out)         :: residual_vectors
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(tddfpt_action_matrix), INTENT(in)             :: action_matrix
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: diag_guess
      REAL(kind=dp), DIMENSION(:,:), INTENT(in)          :: krylov_vectors
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      REAL(kind=dp), PARAMETER :: threshold = 16.0_dp * EPSILON(1.0_dp), eref_scale = 0.99_dp

      INTEGER :: nstates, nstates_total, nvects, nvects_residual, ivect, jvect
      REAL(kind=dp)                                      :: eref, lambda
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: evals_Atilde
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: krylov_vectors_fm, ritz_vectors_fm, &
                                                            A_v_fm, vT_A_v_fm, evects_Atilde_fm

      nstates = SIZE(evals)
      nstates_total = SIZE(krylov_vectors, 1)
      nvects = SIZE(krylov_vectors, 2)

      CPASSERT(nstates <= nvects .AND. nvects <= nstates_total)
      CPASSERT(ALLOCATED(action_matrix%diag))
      CPASSERT(SIZE(action_matrix%diag) == nvects)
      CPASSERT(SIZE(diag_guess) == nstates_total)

      NULLIFY(fm_struct, krylov_vectors_fm, A_v_fm, vT_A_v_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nvects, &
                                  para_env=para_env, context=blacs_env)
      CALL cp_fm_create(krylov_vectors_fm, fm_struct)
      CALL cp_fm_create(A_v_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nvects, ncol_global=nvects, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(vT_A_v_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_set_submatrix(krylov_vectors_fm, krylov_vectors)
      CALL tddfpt_am_multiply(action_matrix=action_matrix, diag_guess=diag_guess, matrix_b=krylov_vectors_fm, matrix_c=A_v_fm)

      CALL cp_gemm('T', 'N', nvects, nvects, nstates_total, 1.0_dp, krylov_vectors_fm, A_v_fm, 0.0_dp, vT_A_v_fm)

      ALLOCATE(evals_Atilde(nvects))
      CALL choose_eigv_solver(vT_A_v_fm, evects_Atilde_fm, evals_Atilde)

      evals(1:nstates) = evals_Atilde(1:nstates)
      DEALLOCATE(evals_Atilde)

      CALL cp_fm_release(vT_A_v_fm)
      CALL cp_fm_release(A_v_fm)

      NULLIFY(ritz_vectors_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nstates, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(ritz_vectors_fm, fm_struct)
      CALL cp_fm_create(A_v_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_gemm('N', 'N', nstates_total, nstates, nvects, 1.0_dp, krylov_vectors_fm, &
                   evects_Atilde_fm, 0.0_dp, ritz_vectors_fm)
      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(krylov_vectors_fm)

      ! residual = A * U - evals * U
      CALL tddfpt_am_multiply(action_matrix=action_matrix, diag_guess=diag_guess, matrix_b=ritz_vectors_fm, matrix_c=A_v_fm)
      CALL cp_fm_column_scale(ritz_vectors_fm, evals(1:nstates))
      CALL cp_fm_scale_and_add(1.0_dp, A_v_fm, -1.0_dp, ritz_vectors_fm)

      nvects_residual = SIZE(residual_vectors, 2)
      IF (nvects_residual > nstates) &
         nvects_residual = nstates
      CALL cp_fm_get_submatrix(A_v_fm, residual_vectors(:,1:nvects_residual), n_cols=nvects_residual)

      CALL cp_fm_release(A_v_fm)
      CALL cp_fm_release(ritz_vectors_fm)

      DO jvect = 1,nvects_residual
         lambda = evals(jvect)
         DO ivect = 1,nvects
            eref = action_matrix%diag(ivect) - lambda
            IF (ABS(eref) < threshold) &
               eref = action_matrix%diag(ivect) - eref_scale * lambda

            residual_vectors(ivect, jvect) = residual_vectors(ivect, jvect) / eref
         END DO

         DO ivect = nvects+1,nstates_total
            eref = diag_guess(ivect) - lambda
            IF (ABS(eref) < threshold) &
               eref = diag_guess(ivect) - eref_scale * lambda

            residual_vectors(ivect, jvect) = residual_vectors(ivect, jvect) / eref
         END DO
      END DO

      CALL tddfpt_orthonormalize_residual_vectors(residual_vectors, krylov_vectors)
   END SUBROUTINE tddfpt_davidson_iteration

   SUBROUTINE tddfpt_davidson_restart(krylov_vectors, nstates, action_matrix, diag_guess, blacs_env, para_env)
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout)       :: krylov_vectors
      INTEGER, INTENT(in)                                :: nstates
      TYPE(tddfpt_action_matrix), INTENT(in)             :: action_matrix
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: diag_guess
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER :: nstates_total, nvects
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: evals_Atilde
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: krylov_vectors_fm, ritz_vectors_fm, &
                                                            A_v_fm, vT_A_v_fm, evects_Atilde_fm

      nstates_total = SIZE(krylov_vectors, 1)
      nvects = SIZE(krylov_vectors, 2)

      CPASSERT(nstates <= nvects .AND. nvects <= nstates_total)
      CPASSERT(ALLOCATED(action_matrix%diag))
      CPASSERT(SIZE(action_matrix%diag) == nvects)
      CPASSERT(SIZE(diag_guess) == nstates_total)

      NULLIFY(fm_struct, krylov_vectors_fm, A_v_fm, vT_A_v_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nvects, &
                                  para_env=para_env, context=blacs_env)
      CALL cp_fm_create(krylov_vectors_fm, fm_struct)
      CALL cp_fm_create(A_v_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nvects, ncol_global=nvects, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(vT_A_v_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_set_submatrix(krylov_vectors_fm, krylov_vectors)
      CALL tddfpt_am_multiply(action_matrix=action_matrix, diag_guess=diag_guess, matrix_b=krylov_vectors_fm, matrix_c=A_v_fm)

      CALL cp_gemm('T', 'N', nvects, nvects, nstates_total, 1.0_dp, krylov_vectors_fm, A_v_fm, 0.0_dp, vT_A_v_fm)

      ALLOCATE(evals_Atilde(nvects))
      CALL choose_eigv_solver(vT_A_v_fm, evects_Atilde_fm, evals_Atilde)
      DEALLOCATE(evals_Atilde)

      CALL cp_fm_release(vT_A_v_fm)
      CALL cp_fm_release(A_v_fm)

      NULLIFY(ritz_vectors_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nstates, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(ritz_vectors_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_gemm('N', 'N', nstates_total, nstates, nvects, 1.0_dp, krylov_vectors_fm, &
                   evects_Atilde_fm, 0.0_dp, ritz_vectors_fm)
      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(krylov_vectors_fm)

!      ALLOCATE(ritz_vectors(nstates_total, nstates))
      CALL cp_fm_get_submatrix(ritz_vectors_fm, krylov_vectors(:,1:nstates))
      CALL cp_fm_release(ritz_vectors_fm)

      CALL tddfpt_orthonormalize_residual_vectors(krylov_vectors(:,1:nstates), krylov_vectors(:,1:0))
   END SUBROUTINE tddfpt_davidson_restart

   SUBROUTINE tddfpt_orthonormalize_residual_vectors(residual_vectors, krylov_vectors)
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: residual_vectors
      REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: krylov_vectors

      INTEGER :: ikv, irv, nkv, nrv, nstates_total
      REAL(kind=dp) :: norm

      REAL(kind=dp), EXTERNAL :: ddot, dnrm2

      nstates_total = SIZE(krylov_vectors, 1)
      nkv = SIZE(krylov_vectors, 2)
      nrv = SIZE(residual_vectors, 2)

      DO irv = 1, nrv
         DO ikv = 1,nkv
!            norm = DOT_PRODUCT(residual_vectors(:,irv), krylov_vectors(:,ikv))
            norm = ddot(nstates_total, residual_vectors(1:nstates_total,irv), 1, krylov_vectors(1:nstates_total,ikv), 1)
            residual_vectors(:,irv) = residual_vectors(:,irv) - norm * krylov_vectors(:,ikv)
         END DO

         DO ikv = 1,irv-1
!            norm = DOT_PRODUCT(residual_vectors(:,irv), residual_vectors(:,ikv))
            norm = ddot(nstates_total, residual_vectors(1:nstates_total,irv), 1, residual_vectors(1:nstates_total,ikv), 1)
            residual_vectors(:,irv) = residual_vectors(:,irv) - norm * residual_vectors(:,ikv)
         END DO

!         norm = 1.0_dp / SQRT(DOT_PRODUCT(residual_vectors(:,irv), residual_vectors(:,irv)))
         norm = 1.0_dp / dnrm2(nstates_total, residual_vectors(1:nstates_total, irv), 1)
         residual_vectors(:,irv) = norm * residual_vectors(:,irv)
      END DO
   END SUBROUTINE tddfpt_orthonormalize_residual_vectors

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param mos_occ         occupied molecular orbitals (allocated on exit)
!> \param evals_occ       occupied orbital energies (allocated on exit)
!> \param mos_virt        virtual molecular orbitals (allocated on exit)
!> \param evals_virt      virtual orbital energies (allocated on exit)
!> \param mo_set          ground state molecular orbitals for a given spin
!> \param matrix_ks       Kohn-Sham matrix for a given spin
!> \param matrix_s        overlap matrix
!> \param cholesky_method Cholesky method to compute the inverse overlap matrix
!> \param blacs_env       BLACS environment for newly allocated distributed matrices
!> \param para_env        parallel environment for newly allocated distributed matrices
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_lumos(mos_occ, evals_occ, mos_virt, evals_virt, &
                           mo_set, matrix_ks, matrix_s, cholesky_method, blacs_env, para_env)
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: evals_occ
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: evals_virt
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_lumos', routineP = moduleN//':'//routineN

      INTEGER                                            :: cholesky_method_inout, handle, homo, &
                                                            iocc, nao, nelectrons, nmo_occ_actual, &
                                                            nmo_scf, nmo_virt_actual
      REAL(kind=dp)                                      :: maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      TYPE(cp_fm_pool_type), POINTER                     :: ao_ao_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      CPASSERT((.NOT. ASSOCIATED(mos_occ)) .AND. (.NOT. ASSOCIATED(mos_virt)))
      CPASSERT((.NOT. ALLOCATED(evals_occ)) .AND. (.NOT. ALLOCATED(evals_virt)))

      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=homo, maxocc=maxocc, nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_occ_actual = homo
      nmo_virt_actual = nao-homo

      IF (nmo_virt_actual <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'Unable to generate virtual molecular orbitals. Please use a larger atomic basis set.')

      ! structures of distributed dense matrices
      NULLIFY (ao_ao_fm_struct, ao_mo_occ_fm_struct, ao_mo_virt_fm_struct)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ_actual, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt_actual, &
                               para_env=para_env, context=blacs_env)

      ! storage space
      CALL cp_fm_create(mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(mos_virt, ao_mo_virt_fm_struct)
      ALLOCATE (evals_occ(nmo_occ_actual))
      ALLOCATE (evals_virt(nmo_virt_actual))

      ! set of molecular orbitals
      NULLIFY (ao_ao_fm_pool, mos_extended)
      CALL fm_pool_create(ao_ao_fm_pool, ao_ao_fm_struct)
      CALL allocate_mo_set(mos_extended, nao, nao, nelectrons, &
                           REAL(nelectrons, dp), maxocc, flexible_electron_count=0.0_dp)
      CALL init_mo_set(mos_extended, fm_pool=ao_ao_fm_pool, name="mos-extended")
      CALL fm_pool_release(ao_ao_fm_pool)

      CALL get_mo_set(mos_extended, mo_coeff=mo_coeff_extended, &
                      eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)

      ! use the explicit loop in order to avoid temporary arrays.
      !
      ! The assignment statement : mo_occ_extended(1:nmo_scf) = mo_occ_scf(1:nmo_scf)
      ! implies temporary arrays as a compiler does not know in advance that the pointers
      ! on both sides of the statement point to non-overlapped memory regions
      DO iocc = 1, nmo_scf
         mo_occ_extended(iocc) = mo_occ_scf(iocc)
      END DO
      mo_occ_extended(nmo_scf+1:) = 0.0_dp

      ! ++ create temporary matrices
      NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
      CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
      CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
      CALL cp_fm_create(work_fm, ao_ao_fm_struct)

      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)
      CALL cp_fm_struct_release(ao_ao_fm_struct)

      ! some stuff from the subroutine general_eigenproblem()
      CALL copy_dbcsr_to_fm(matrix_s, ortho_fm)
      CALL copy_dbcsr_to_fm(matrix_ks, matrix_ks_fm)

      IF (cholesky_method == cholesky_dbcsr) THEN
         CPABORT('CHOLESKY DBCSR_INVERSE is not implemented in TDDFT.')
      ELSE IF (cholesky_method == cholesky_off) THEN
         CPABORT('CHOLESKY OFF is not implemented in TDDFT.')
      ELSE
         CALL cp_fm_cholesky_decompose(ortho_fm)
         IF (cholesky_method == cholesky_inverse) THEN
            CALL cp_fm_triangular_invert(ortho_fm)
         END IF

         ! need this temporary variable, as the subroutine eigensolver() is going to update it.
         cholesky_method_inout = cholesky_method
         CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=mos_extended, ortho=ortho_fm, &
                          work=work_fm, cholesky_method=cholesky_method_inout, &
                          do_level_shift=.FALSE., level_shift=0.0_dp, matrix_u_fm=null(), use_jacobi=.FALSE.)
      END IF

      ! -- clean up needless matrices
      CALL cp_fm_release(ortho_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(matrix_ks_fm)

      ! return the requested occupied and virtual molecular orbitals and corresponding orbital energies
      CALL cp_fm_to_fm_submat(msource=mo_coeff_extended, mtarget=mos_occ, nrow=nao, ncol=nmo_occ_actual, &
                              s_firstrow=1, s_firstcol=homo-nmo_occ_actual+1, t_firstrow=1, t_firstcol=1)
      evals_occ(1:nmo_occ_actual) = mo_evals_extended(homo-nmo_occ_actual+1:homo)

      CALL cp_fm_to_fm_submat(msource=mo_coeff_extended, mtarget=mos_virt, nrow=nao, ncol=nmo_virt_actual, &
                              s_firstrow=1, s_firstcol=homo+1, t_firstrow=1, t_firstcol=1)
      evals_virt(1:nmo_virt_actual) = mo_evals_extended(homo+1:homo+nmo_virt_actual)

      CALL deallocate_mo_set(mos_extended)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_lumos

   SUBROUTINE tddfpt_guess(guess_vectors, nstates, nstates_active, evals_occ, evals_virt, &
                           is_add_degenerate, degenerate_eps, log_unit)
      TYPE(tddfpt_guess_vectors), INTENT(inout)          :: guess_vectors
      INTEGER, INTENT(in)                                :: nstates
      INTEGER, INTENT(in)                                :: nstates_active
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_virt
      LOGICAL, INTENT(in)                                :: is_add_degenerate
      REAL(kind=dp), INTENT(in)                          :: degenerate_eps
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_guess', routineP = moduleN//':'//routineN

      INTEGER :: handle, nmo_occ_avail, nmo_virt_avail, nstates_total, nmo_occ_selected, nmo_virt_selected, &
                 nstates_selected, imo_occ, imo_virt, istate, jstate
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ALLOCATED(guess_vectors%descrs))
      CPASSERT(.NOT. ALLOCATED(guess_vectors%evals))

      nmo_occ_avail  = SIZE(evals_occ)
      nmo_virt_avail = SIZE(evals_virt)
      nstates_total = nstates_active

      DO nmo_occ_selected = MIN(nstates_total, nmo_occ_avail), nmo_occ_avail-1
         IF (evals_occ(nmo_occ_avail-nmo_occ_selected+1)-evals_occ(nmo_occ_avail-nmo_occ_selected) > degenerate_eps) EXIT
      END DO

      DO nmo_virt_selected = MIN(nstates_total, nmo_virt_avail), nmo_virt_avail-1
         IF (evals_virt(nmo_virt_selected+1)-evals_virt(nmo_virt_selected) > degenerate_eps) EXIT
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = nmo_occ_selected*nmo_virt_selected

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO imo_occ = 1, nmo_occ_selected
         e_occ = evals_occ(nmo_occ_avail-imo_occ+1)
         DO imo_virt = 1, nmo_virt_selected
            istate = istate+1
            e_virt_minus_occ(istate) = evals_virt(imo_virt)-e_occ
         END DO
      END DO

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,18X,A)') "TDDFT initial guess"
         WRITE (log_unit, '(1X,52("-"))')
         WRITE (log_unit, '(5X,A)') "State     Occupied  ->  Virtual      Excitation"
         WRITE (log_unit, '(5X,A)') "number     orbital      orbital      energy (eV)"
         WRITE (log_unit, '(1X,52("-"))')

         DO istate = 1, nstates
            imo_occ = (inds(istate)-1)/nmo_virt_selected+1
            imo_virt = MOD(inds(istate)-1, nmo_virt_selected)+1
               WRITE (log_unit, '(1X,I8,4X,I8,5X,I8,5X,F14.5)') istate, nmo_occ_avail-imo_occ+1, &
               imo_virt+nmo_occ_avail, e_virt_minus_occ(istate)*evolt
         END DO
      END IF

      ! detect degenerate states outside the active space
      DO istate = nstates_total+1, nstates_selected
         IF (e_virt_minus_occ(istate)-e_virt_minus_occ(istate-1) > degenerate_eps) EXIT
      END DO

      IF (istate > nstates_total+1) THEN
         istate = istate - 1

         IF (is_add_degenerate) THEN
            nstates_total = istate
         ELSE
            IF (log_unit > 0) THEN
               DO jstate = nstates_total, 2, -1
                  IF (e_virt_minus_occ(jstate)-e_virt_minus_occ(jstate-1) > degenerate_eps) EXIT
               END DO

               IF (jstate <= nstates) &
                  jstate = nstates + 1

               IF (jstate <= nstates_total) THEN
                  WRITE (log_unit, '(1X,15("-"),1X,A,1X,16("-"))') "Other active states"

                  DO jstate = jstate, nstates_total
                     imo_occ = (inds(jstate)-1)/nmo_virt_selected
                     imo_virt = MOD(inds(jstate)-1, nmo_virt_selected)+1

                     WRITE (log_unit, '(1X,I8,4X,I8,5X,I8,5X,F14.5)') jstate, nmo_occ_avail-imo_occ, &
                     imo_virt+nmo_occ_avail, e_virt_minus_occ(jstate)*evolt
                  END DO
               END IF

               WRITE (log_unit, '(1X,17("-"),1X,A,1X,18("-"))') "Inactive states"

               DO jstate = nstates_total+1,istate
                  imo_occ = (inds(jstate)-1)/nmo_virt_selected
                  imo_virt = MOD(inds(jstate)-1, nmo_virt_selected)+1

                  WRITE (log_unit, '(1X,I8,4X,I8,5X,I8,5X,F14.5)') jstate, nmo_occ_avail-imo_occ, &
                  imo_virt+nmo_occ_avail, e_virt_minus_occ(jstate)*evolt
               END DO
            END IF

            CALL cp_warn(__LOCATION__, "It is advised to include all proposed"// &
                         " degenerate excited states into TDDFT active space")
         END IF
      END IF

      ALLOCATE (guess_vectors%evals(nstates_total))
      guess_vectors%evals(1:nstates_total) = e_virt_minus_occ(1:nstates_total)
      DEALLOCATE (e_virt_minus_occ)

      IF (log_unit > 0) &
         WRITE (log_unit, '(/,1X,A,T46,I8)') 'Number of active states:', nstates_total

      ALLOCATE (guess_vectors%descrs(nstates_total))
      DO istate = 1,nstates_total
         guess_vectors%descrs(istate)%imo_occ  = nmo_occ_avail-(inds(istate)-1)/nmo_virt_selected
         guess_vectors%descrs(istate)%imo_virt = MOD(inds(istate)-1, nmo_virt_selected)+1
      END DO

      DEALLOCATE(inds)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_guess

   SUBROUTINE tddfpt_assemble_evects(evects, krylov_vectors, guess_vectors, mos_virt)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      REAL(kind=dp), DIMENSION(:,:), INTENT(in)          :: krylov_vectors
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_fm_type), POINTER                          :: mos_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_assemble_evects', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, ivect, jvect, nvects_krylov, nstates, nao, nmo_occ, nmo_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_p_type), DIMENSION(:), ALLOCATABLE :: mos_virt_combined, mos_virt_columns

      CALL timeset(routineN, handle)

      nvects_krylov = SIZE(evects)
      nstates = SIZE(guess_vectors%descrs)

!      CPASSERT(nvects_krylov > 0 .AND. nstates > 0)
      CPASSERT(SIZE(krylov_vectors, 1) == nstates)
      CPASSERT(SIZE(krylov_vectors, 2) == nvects_krylov)

!      NULLIFY (fm_struct, evect_descr)
      CALL cp_fm_get_info(matrix=evects(1)%matrix, ncol_global=nmo_occ, context=blacs_env, para_env=para_env)
      CALL cp_fm_get_info(matrix=mos_virt, nrow_global=nao, ncol_global=nmo_virt)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=1, context=blacs_env, para_env=para_env)

      ALLOCATE (mos_virt_columns(nmo_virt))
      DO ivect = 1,nmo_virt
         NULLIFY(mos_virt_columns(ivect)%matrix)
         CALL cp_fm_create(mos_virt_columns(ivect)%matrix, fm_struct)
         CALL cp_fm_to_fm(mos_virt, mos_virt_columns(ivect)%matrix, ncol=1, source_start=ivect, target_start=1)
      END DO

      ALLOCATE (mos_virt_combined(nmo_occ))
      DO ivect = 1,nmo_occ
         NULLIFY(mos_virt_combined(ivect)%matrix)
         CALL cp_fm_create(mos_virt_combined(ivect)%matrix, fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      DO ivect = 1,nvects_krylov
         DO jvect = 1,nmo_occ
            CALL cp_fm_set_all(mos_virt_combined(jvect)%matrix, 0.0_dp)
         END DO

         DO jvect = 1, nstates
            CALL cp_fm_scale_and_add(1.0_dp, mos_virt_combined(guess_vectors%descrs(jvect)%imo_occ)%matrix, &
                                     krylov_vectors(jvect, ivect), mos_virt_columns(guess_vectors%descrs(jvect)%imo_virt)%matrix)
         END DO

         DO jvect = 1,nmo_occ
            CALL cp_fm_to_fm(mos_virt_combined(jvect)%matrix, evects(ivect)%matrix, ncol=1, source_start=1, target_start=jvect)
         END DO
      END DO

      DO ivect = SIZE(mos_virt_combined),1,-1
         CALL cp_fm_release(mos_virt_combined(ivect)%matrix)
      END DO
      DEALLOCATE (mos_virt_combined)

      DO ivect = SIZE(mos_virt_columns),1,-1
         CALL cp_fm_release(mos_virt_columns(ivect)%matrix)
      END DO
      DEALLOCATE (mos_virt_columns)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_assemble_evects

! **************************************************************************************************
!> \brief Initialise the action_matrix using the following expression:
!>        trace { C_{1,i}^T * ( KS * C_{1,j} - S * C_{1,j} * diag(evals_occ) ) } .
!> \param a_kv      the next KK block of the action_matrix (initialised on exit)
!> \param evects    TDDFPT trial vectors
!> \param evals_occ occupied orbital energies
!> \param matrix_ks Kohn-Sham matrix for a given spin
!> \param matrix_s  overlap matrix
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_by_energy_diff(block_nn, evects, evals_occ, matrix_ks, matrix_s)
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: block_nn
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_by_energy_diff', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, istate, jstate, nmo, nstates
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ks_psi1_fm, s_psi1_fm

      CALL timeset(routineN, handle)

      nstates = SIZE(evects)
      CPASSERT(SIZE(block_nn, 1) == nstates)
      CPASSERT(SIZE(block_nn, 2) == nstates)

      NULLIFY (fm_struct, s_psi1_fm)
      CALL cp_fm_get_info(matrix=evects(1)%matrix, ncol_global=nmo, matrix_struct=fm_struct)
      CALL cp_fm_create(ks_psi1_fm, fm_struct)
      CALL cp_fm_create(s_psi1_fm, fm_struct)

      DO jstate = 1, nstates
         CALL cp_dbcsr_sm_fm_multiply(matrix_ks, evects(jstate)%matrix, ks_psi1_fm, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(jstate)%matrix, s_psi1_fm, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_fm_column_scale(s_psi1_fm, evals_occ)
         !  KS * C_1 - S * C_1 * occupied_orbital_energies
         CALL cp_fm_scale_and_add(1.0_dp, ks_psi1_fm, -1.0_dp, s_psi1_fm)

         DO istate = 1, jstate
            CALL cp_fm_trace(evects(istate)%matrix, ks_psi1_fm, block_nn(istate, jstate))
         END DO
      END DO

      ! on exit : a_kv(i, j) = \sum_{occ, virt} (C_virt * epsilon_virt - C_occ * epsilon_occ)
      CALL cp_fm_release(s_psi1_fm)
      CALL cp_fm_release(ks_psi1_fm)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_by_energy_diff

   SUBROUTINE tddfpt_apply_coulomb(block_on, block_nn, block_gn, alpha, evects, mos_occ, guess_vectors, mos_virt, qs_env, matrix_s)
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: block_on
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout)       :: block_nn
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: block_gn
!      TYPE(tddfpt_action_matrix), INTENT(inout)          :: action_matrix
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikv, nao, nstates_total, &
             nvects_exists, nkv_new, ivect, jvect, jkv, nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: J_ia_jb
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: J_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: J_ia_jnu, rho_ia_fm
      TYPE(cp_fm_p_type), DIMENSION(:), ALLOCATABLE :: mos_occ_columns, mos_virt_columns
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_gspace, v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      nkv_new = SIZE(block_nn, 1)
      nvects_exists = nkv_new
      IF (ALLOCATED(block_on)) THEN
         nvects_exists = nvects_exists + SIZE(block_on, 1)
      END IF
      nstates_total = nvects_exists
      IF (ALLOCATED(block_gn)) THEN
         nstates_total = nstates_total + SIZE(block_gn, 1)
      END IF
      CPASSERT(SIZE(evects) == nvects_exists)
      CPASSERT(SIZE(guess_vectors%evals) == nstates_total)

      CALL cp_fm_get_info(mos_occ,  nrow_global=nao, ncol_global=nmo_occ)
      CALL cp_fm_get_info(mos_virt, ncol_global=nmo_virt)

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_g=rho_ia_g)

      NULLIFY (J_ia_munu%matrix)
      CALL cp_dbcsr_init_p(J_ia_munu%matrix)
      ! TO DO: is there a better way to initialise a DBCSR matrix?
      ! The subroutine cp_dbcsr_create() alone does not work here,
      ! as block distribution of the matrix J_ia_munu should be identical with the matrix_s
      ! (the subroutine integrate_v_rspace() will fail otherwise) .
      CALL cp_dbcsr_copy(J_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (J_ia_jnu)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(J_ia_jnu, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      IF (ALLOCATED(block_gn)) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=1, para_env=para_env, context=blacs_env)

         ALLOCATE (mos_occ_columns(nmo_occ))
         DO ivect = 1,nmo_occ
            NULLIFY(mos_occ_columns(ivect)%matrix)
            CALL cp_fm_create(mos_occ_columns(ivect)%matrix, fm_struct)
         END DO

         ALLOCATE (mos_virt_columns(nmo_virt))
         DO ivect = 1,nmo_virt
            NULLIFY(mos_virt_columns(ivect)%matrix)
            CALL cp_fm_create(mos_virt_columns(ivect)%matrix, fm_struct)
            CALL cp_fm_to_fm(mos_virt, mos_virt_columns(ivect)%matrix, ncol=1, source_start=ivect, target_start=1)
         END DO

         CALL cp_fm_struct_release(fm_struct)
      END IF

      NULLIFY (v_gspace%pw, v_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw, use_data=REALDATA3D, in_space=REALSPACE)

      DO ikv = 1,nkv_new
         ivect = nvects_exists - nkv_new + ikv

         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, mos_occ, evects(ivect)%matrix, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, evects(ivect)%matrix, mos_occ, 1.0_dp, rho_ia_fm)

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         ! here J_ia_jb is used as a temporary variable to store needless pair energy
         CALL pw_poisson_solve(poisson_env, rho_ia_g(1)%pw, J_ia_jb, v_gspace%pw)
         CALL pw_transfer(v_gspace%pw, v_rspace%pw)
         CALL pw_scale(v_rspace%pw, v_rspace%pw%pw_grid%dvol)

         CALL cp_dbcsr_set(J_ia_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=v_rspace, hmat=J_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

         CALL cp_dbcsr_sm_fm_multiply(J_ia_munu%matrix, mos_occ, J_ia_jnu, ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)

         IF (ALLOCATED(block_on)) THEN
            DO jvect = 1,nvects_exists-nkv_new
               CALL cp_fm_trace(evects(jvect)%matrix, J_ia_jnu, J_ia_jb)
               block_on(jvect, ikv) = block_on(jvect, ikv) + alpha*J_ia_jb
            END DO
         END IF

         DO jkv = 1,ikv
            jvect = nvects_exists - nkv_new + jkv

            CALL cp_fm_trace(evects(jvect)%matrix, J_ia_jnu, J_ia_jb)
            block_nn(jkv, ikv) = block_nn(jkv, ikv) + alpha*J_ia_jb
         END DO

         IF (ALLOCATED(block_gn)) THEN
            DO jvect = 1,nmo_occ
               CALL cp_fm_to_fm(J_ia_jnu, mos_occ_columns(jvect)%matrix, ncol=1, source_start=jvect, target_start=1)
            END DO

            DO jvect = nvects_exists+1,nstates_total
               jkv = jvect - nvects_exists

               CALL cp_fm_trace(mos_virt_columns(guess_vectors%descrs(jvect)%imo_virt)%matrix, &
                    mos_occ_columns(guess_vectors%descrs(jvect)%imo_occ)%matrix, J_ia_jb)
               block_gn(jkv, ikv) = block_gn(jkv, ikv) + alpha*J_ia_jb
            END DO
         END IF
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)

      IF (ALLOCATED(mos_virt_columns)) THEN
         DO ivect = SIZE(mos_virt_columns),1,-1
            CALL cp_fm_release(mos_virt_columns(ivect)%matrix)
         END DO
         DEALLOCATE (mos_virt_columns)
      END IF

      IF (ALLOCATED(mos_occ_columns)) THEN
         DO ivect = SIZE(mos_occ_columns),1,-1
            CALL cp_fm_release(mos_occ_columns(ivect)%matrix)
         END DO
         DEALLOCATE (mos_occ_columns)
      END IF

      CALL cp_fm_release(J_ia_jnu)
      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(J_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_coulomb

   SUBROUTINE tddfpt_apply_kernel(block_on, block_nn, block_gn, alpha, evects, mos_occ, guess_vectors, mos_virt, &
                                  xc_rho_set, xc_deriv_set, qs_env, matrix_s, xc_section)
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: block_on
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout)       :: block_nn
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: block_gn
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikv, nao, nstates_total, &
             xc_deriv_method_id, xc_rho_smooth_id, &
             nvects_exists, nkv_new, ivect, jvect, jkv, nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: F_ia_jb
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: F_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: F_ia_jnu, rho_ia_fm
      TYPE(cp_fm_p_type), DIMENSION(:), ALLOCATABLE      :: mos_occ_columns, mos_virt_columns
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_r, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set


      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nkv_new = SIZE(block_nn, 1)
      nvects_exists = nkv_new
      IF (ALLOCATED(block_on)) THEN
         nvects_exists = nvects_exists + SIZE(block_on, 1)
      END IF
      nstates_total = nvects_exists
      IF (ALLOCATED(block_gn)) THEN
         nstates_total = nstates_total + SIZE(block_gn, 1)
      END IF
      CPASSERT(SIZE(evects) == nvects_exists)
      CPASSERT(SIZE(guess_vectors%evals) == nstates_total)

      CALL cp_fm_get_info(mos_occ,  nrow_global=nao, ncol_global=nmo_occ)
      CALL cp_fm_get_info(mos_virt, ncol_global=nmo_virt)

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct, rho_ia_r, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho_ia_r, rho_g=rho_ia_g)

      NULLIFY (F_ia_munu%matrix)
      CALL cp_dbcsr_init_p(F_ia_munu%matrix)
      CALL cp_dbcsr_copy(F_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (fm_struct, F_ia_jnu)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(F_ia_jnu, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (v_xc(1))
      NULLIFY (v_xc(1)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(1)%pw, use_data=REALDATA3D, in_space=REALSPACE)

      NULLIFY (xc_rho1_set)
      CALL xc_rho_set_create(xc_rho1_set, rho_ia_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      xc_deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      xc_rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      ! get arguments needed by the given functional
      needs = xc_functionals_get_needs(functionals=xc_fun_section, lsd=.FALSE., add_basic_components=.TRUE.)

      IF (ALLOCATED(block_gn)) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=1, para_env=para_env, context=blacs_env)

         ALLOCATE (mos_occ_columns(nmo_occ))
         DO ivect = 1,nmo_occ
            NULLIFY(mos_occ_columns(ivect)%matrix)
            CALL cp_fm_create(mos_occ_columns(ivect)%matrix, fm_struct)
         END DO

         ALLOCATE (mos_virt_columns(nmo_virt))
         DO ivect = 1,nmo_virt
            NULLIFY(mos_virt_columns(ivect)%matrix)
            CALL cp_fm_create(mos_virt_columns(ivect)%matrix, fm_struct)
            CALL cp_fm_to_fm(mos_virt, mos_virt_columns(ivect)%matrix, ncol=1, source_start=ivect, target_start=1)
         END DO

         CALL cp_fm_struct_release(fm_struct)
      END IF

      DO ikv = 1,nkv_new
         ivect = nvects_exists - nkv_new + ikv

         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, mos_occ, evects(ivect)%matrix, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, evects(ivect)%matrix, mos_occ, 1.0_dp, rho_ia_fm)

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         CALL xc_rho_set_update(rho_set=xc_rho1_set, rho_r=rho_ia_r, rho_g=rho_ia_g, &
                                tau=null(), needs=needs, xc_deriv_method_id=xc_deriv_method_id, &
                                xc_rho_smooth_id=xc_rho_smooth_id, pw_pool=auxbas_pw_pool)

         CALL pw_zero(v_xc(1)%pw)
         CALL xc_calc_2nd_deriv(v_xc=v_xc, deriv_set=xc_deriv_set, rho_set=xc_rho_set, &
                                rho1_set=xc_rho1_set, pw_pool=auxbas_pw_pool, &
                                xc_section=xc_section, gapw=.FALSE., tddfpt_fac=0.0_dp)

         CALL pw_scale(v_xc(1)%pw, v_xc(1)%pw%pw_grid%dvol)

         CALL cp_dbcsr_set(F_ia_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=v_xc(1), hmat=F_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

         CALL cp_dbcsr_sm_fm_multiply(F_ia_munu%matrix, mos_occ, F_ia_jnu, ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)

         IF (ALLOCATED(block_on)) THEN
            DO jvect = 1,nvects_exists-nkv_new
               CALL cp_fm_trace(evects(jvect)%matrix, F_ia_jnu, F_ia_jb)
               block_on(jvect, ikv) = block_on(jvect, ikv) + alpha*F_ia_jb
            END DO
         END IF

         DO jkv = 1,ikv
            jvect = nvects_exists - nkv_new + jkv

            CALL cp_fm_trace(evects(jvect)%matrix, F_ia_jnu, F_ia_jb)
            block_nn(jkv, ikv) = block_nn(jkv, ikv) + alpha*F_ia_jb
         END DO

         IF (ALLOCATED(block_gn)) THEN
            DO jvect = 1,nmo_occ
               CALL cp_fm_to_fm(F_ia_jnu, mos_occ_columns(jvect)%matrix, ncol=1, source_start=jvect, target_start=1)
            END DO

            DO jvect = nvects_exists+1,nstates_total
               jkv = jvect - nvects_exists

               CALL cp_fm_trace(mos_virt_columns(guess_vectors%descrs(jvect)%imo_virt)%matrix, &
                    mos_occ_columns(guess_vectors%descrs(jvect)%imo_occ)%matrix, F_ia_jb)
               block_gn(jkv, ikv) = block_gn(jkv, ikv) + alpha*F_ia_jb
            END DO
         END IF
      END DO

      IF (ALLOCATED(mos_virt_columns)) THEN
         DO ivect = SIZE(mos_virt_columns),1,-1
            CALL cp_fm_release(mos_virt_columns(ivect)%matrix)
         END DO
         DEALLOCATE (mos_virt_columns)
      END IF

      IF (ALLOCATED(mos_occ_columns)) THEN
         DO ivect = SIZE(mos_occ_columns),1,-1
            CALL cp_fm_release(mos_occ_columns(ivect)%matrix)
         END DO
         DEALLOCATE (mos_occ_columns)
      END IF

      CALL xc_rho_set_release(xc_rho1_set)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(1)%pw)
      DEALLOCATE (v_xc)

      CALL cp_fm_release(F_ia_jnu)
      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(F_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_kernel

   SUBROUTINE tddfpt_apply_hfx(block_on, block_nn, K_ij_munu, alpha, evects, mos_occ, qs_env, matrix_s)
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: block_on
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout)       :: block_nn
      TYPE(cp_fm_p_type), DIMENSION(:,:), INTENT(inout)  :: K_ij_munu
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_hfx', &
         routineP = moduleN//':'//routineN

REAL(kind=dp), PARAMETER :: threshold = 1.0e-8_dp

      INTEGER                                            :: handle, nkv_new, nvects_exists, ikv, jkv, &
               ivect, jvect, nao, nmo_occ!, iocc, jocc
!      REAL(kind=dp)                                      :: K_ij_ab
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: weights_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_a, mo_coeff_b, S_ivect, weights_occ_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env

      TYPE(pw_p_type)                                    :: psi0_i, psi0_j, rho_g, v_gspace
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_p_type)                              :: hfx_munu
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
!INTEGER :: ix

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, blacs_env, cell, dft_control, para_env, particle_set, pw_env, qs_kind_set)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, blacs_env=blacs_env, &
                      cell=cell, dft_control=dft_control, para_env=para_env, &
                      particle_set=particle_set, pw_env=pw_env, qs_kind_set=qs_kind_set)

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      NULLIFY (fm_struct, S_ivect)
      CALL cp_fm_get_info(matrix=mos_occ, nrow_global=nao, ncol_global=nmo_occ, matrix_struct=fm_struct)

      NULLIFY (psi0_i%pw, psi0_j%pw, v_gspace%pw, rho_g%pw, hfx_munu%matrix, work_fm)

      nkv_new = SIZE(block_nn, 1)
      nvects_exists = nkv_new
      IF (ALLOCATED(block_on)) THEN
         nvects_exists = nvects_exists + SIZE(block_on, 1)
      END IF
      CPASSERT(SIZE(evects) == nvects_exists)

      CALL cp_fm_create(S_ivect, fm_struct)

      NULLIFY (fm_struct, mo_coeff_a, mo_coeff_b, weights_occ_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=1, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(mo_coeff_a, fm_struct)
      CALL cp_fm_create(mo_coeff_b, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ, ncol_global=nmo_occ, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(weights_occ_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      ALLOCATE(weights_occ(nmo_occ, nmo_occ))

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)

      DO ikv = 1,nkv_new
         ivect = nvects_exists - nkv_new + ikv

         CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ivect)%matrix, S_ivect, ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)

         DO jkv = 1,ikv
            jvect = nvects_exists - nkv_new + jkv
            CALL cp_gemm('T', 'N', nmo_occ, nmo_occ, nao, 1.0_dp, evects(jvect)%matrix, S_ivect, 0.0_dp, weights_occ_fm)
            CALL cp_fm_get_submatrix(weights_occ_fm, weights_occ)

            block_nn(jkv, ikv) = block_nn(jkv, ikv) - alpha * tddfpt_compute_K_ij_ab(evects(ivect)%matrix, evects(jvect)%matrix)
         END DO

!         DO jvect = 1,nvects_exists-nkv_new
!            CALL cp_gemm('T', 'N', nmo_occ, nmo_occ, nao, 1.0_dp, evects(jvect)%matrix, S_ivect, 0.0_dp, weights_occ_fm)
!            CALL cp_fm_get_submatrix(weights_occ_fm, weights_occ)
!
!            block_on(jvect, ikv) = block_on(jvect, ikv) - &
!                                   alpha * tddfpt_compute_K_ij_ab(evects(ivect)%matrix, evects(jvect)%matrix)
!         END DO
      END DO

      IF (ASSOCIATED(work_fm)) THEN
         CALL cp_fm_release(work_fm)
      END IF
      IF (ASSOCIATED(hfx_munu%matrix)) &
         CALL cp_dbcsr_deallocate_matrix(hfx_munu%matrix)

      IF (ASSOCIATED(rho_g%pw)) &
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)

      IF (ASSOCIATED(v_gspace%pw)) &
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)

      IF (ASSOCIATED(psi0_j%pw)) &
         CALL pw_pool_give_back_pw(auxbas_pw_pool, psi0_j%pw)

      IF (ASSOCIATED(psi0_i%pw)) &
         CALL pw_pool_give_back_pw(auxbas_pw_pool, psi0_i%pw)

      CALL cp_fm_struct_release(fm_struct)

      DEALLOCATE(weights_occ)
      CALL cp_fm_release(weights_occ_fm)
      CALL cp_fm_release(mo_coeff_b)
      CALL cp_fm_release(mo_coeff_a)
      CALL cp_fm_release(S_ivect)

      CALL timestop(handle)

      CONTAINS

      FUNCTION tddfpt_compute_K_ij_ab(evect_i, evect_j) RESULT(K_ij_ab)
         TYPE(cp_fm_type), POINTER                          :: evect_i, evect_j
         REAL(kind=dp) :: K_ij_ab

         INTEGER :: iocc, jocc
         REAL(kind=dp) :: pair_energy

         K_ij_ab = 0.0_dp
         DO iocc = 1,nmo_occ
            DO jocc = 1,iocc
               IF (ABS(weights_occ(jocc, iocc)) >= threshold .AND. .NOT. ASSOCIATED(K_ij_munu(jocc, iocc)%matrix)) THEN
                  ! compute missed element ((i j || mu nu))

                  IF (.NOT. ASSOCIATED(psi0_i%pw)) &
                     CALL pw_pool_create_pw(auxbas_pw_pool, psi0_i%pw, use_data=REALDATA3D, in_space=REALSPACE)
                  IF (.NOT. ASSOCIATED(psi0_j%pw)) &
                     CALL pw_pool_create_pw(auxbas_pw_pool, psi0_j%pw, use_data=REALDATA3D, in_space=REALSPACE)
                  IF (.NOT. ASSOCIATED(v_gspace%pw)) &
                     CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  IF (.NOT. ASSOCIATED(rho_g%pw)) &
                     CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  IF (.NOT. ASSOCIATED(hfx_munu%matrix)) THEN
                     CALL cp_dbcsr_init_p(hfx_munu%matrix)
                     CALL cp_dbcsr_copy(hfx_munu%matrix, matrix_s)
                  END IF
                  IF (.NOT. ASSOCIATED(work_fm)) THEN
                     CALL cp_fm_create(work_fm, fm_struct)
                  END IF
                  CALL calculate_wavefunction(mos_occ, iocc, psi0_i, rho_g, atomic_kind_set, &
                                qs_kind_set, cell, dft_control, particle_set, pw_env) ! TO DO: move 1 level top-ward

                  CALL calculate_wavefunction(mos_occ, jocc, psi0_j, rho_g, atomic_kind_set, &
                                   qs_kind_set, cell, dft_control, particle_set, pw_env)
                  psi0_j%pw%cr3d = psi0_j%pw%cr3d*psi0_i%pw%cr3d

                  CALL pw_transfer(psi0_j%pw, rho_g%pw)
                  CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, v_gspace%pw)
                  CALL pw_transfer(v_gspace%pw, psi0_j%pw)
                  CALL pw_scale(psi0_j%pw, psi0_j%pw%pw_grid%dvol)

                  CALL cp_dbcsr_set(hfx_munu%matrix, 0.0_dp)
                  CALL integrate_v_rspace(v_rspace=psi0_j, hmat=hfx_munu, &
                       qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

                  CALL cp_fm_create(K_ij_munu(jocc, iocc)%matrix, fm_struct)
                  CALL copy_dbcsr_to_fm(hfx_munu%matrix, K_ij_munu(jocc, iocc)%matrix)
                  CALL cp_fm_upper_to_full(K_ij_munu(jocc, iocc)%matrix, work_fm)
               END IF

               IF (ASSOCIATED(K_ij_munu(jocc, iocc)%matrix)) THEN
                  ! count existing element ((i j || mu nu))
                  CALL cp_fm_to_fm(evect_i, mo_coeff_b, ncol=1, source_start=iocc, target_start=1)
                  CALL cp_gemm('N', 'N', nao, 1, nao, 1.0_dp, K_ij_munu(jocc, iocc)%matrix, &
                               mo_coeff_b, 0.0_dp, mo_coeff_a)
                  CALL cp_fm_to_fm(evect_j, mo_coeff_b, ncol=1, source_start=jocc, target_start=1)
                  CALL cp_fm_trace(mo_coeff_b, mo_coeff_a, pair_energy)

                  K_ij_ab = K_ij_ab + pair_energy
               END IF

            END DO
         END DO

      END FUNCTION tddfpt_compute_K_ij_ab
   END SUBROUTINE tddfpt_apply_hfx

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and analysis.
!> \param log_unit output unit
!> \param nstates  number of excited states to print
!> \param evects   TDDFPT trial vectors
!> \param evals    TDDFPT eigenvalues
!> \param nmo_occ_total total number of occupied molecular orbitals
!> \param mos_occ occupied molecular orbitals
!> \param evals_occ  occupied orbital energies
!> \param mos_virt virtual molecular orbitals
!> \param evals_virt  virtual orbital energies
!> \param matrix_s overlap matrix and their firts derivatives over x, y, and z directions
!> \param min_amplitude the smallest excitation amplitude to print
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       To compute transition dipole integrals we use the fact that the ground state wave-function
!>       and all response functions (which are in fact a linear combination of virtual molecular
!>       orbitals) are eigenfunctions of the Hamiltonian operator. Using the identity:
!>       \f[\vec{r}\hat{H} - \hat{H}\vec{r} = [\vec{r},\hat{H}] = [\vec{r},-1/2 \nabla^2] = \nabla\f]
!>       we can derive a relationship between transition dipole integrals and momentum integrals:
!>       \f[<\psi_i|\nabla|\psi_a> = <\psi_i|\vec{r}|\hat{H}\psi_a> - <\hat{H}\psi_i|\vec{r}|\psi_a>=
!>       (\epsilon_a - \epsilon_i) <\psi_i|\vec{r}|\psi_a> .\f]
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_print_summary(log_unit, nstates, evects, evals, nmo_occ_total, &
                                   mos_occ, evals_occ, mos_virt, evals_virt, matrix_s, min_amplitude)
      INTEGER, INTENT(in)                                :: log_unit, nstates
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: nmo_occ_total
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_virt
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_summary', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nderivs = 3

      INTEGER                                            :: handle, ideriv, iexc, imo_occ, imo_virt, &
                                                            istate, nao, nmo_occ, nmo_virt, &
                                                            nmo_virt_occ
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp)                                      :: oscillator_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_abs_1d
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_2d
      REAL(kind=dp), DIMENSION(nderivs)                  :: trans_dipole
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(nderivs)             :: dS_mos_occ
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ediff_inv, ediff_inv_weights, &
                                                            S_mos_virt, scaled_evect, weights_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env

      EXTERNAL :: dcopy

      CALL timeset(routineN, handle)

      ! ensure we have the overlap matrix itself and its first derivatives over x, y and z
      CPASSERT(SIZE(matrix_s) >= nderivs+1)

      NULLIFY (fm_struct, blacs_env, para_env, S_mos_virt)
      CALL cp_fm_get_info(matrix=mos_virt, nrow_global=nao, ncol_global=nmo_virt, &
                          matrix_struct=fm_struct, context=blacs_env, para_env=para_env)

      CALL cp_fm_create(S_mos_virt, fm_struct)

      NULLIFY (fm_struct)
      CALL cp_fm_get_info(matrix=evects(1)%matrix, ncol_global=nmo_occ, matrix_struct=fm_struct)
      DO ideriv = 1, nderivs
         NULLIFY (dS_mos_occ(ideriv)%matrix)
         CALL cp_fm_create(dS_mos_occ(ideriv)%matrix, fm_struct)
      END DO

      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mos_virt, S_mos_virt, ncol=nmo_virt, alpha=1.0_dp, beta=0.0_dp)
      DO ideriv = 1, nderivs
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(ideriv+1)%matrix, mos_occ, dS_mos_occ(ideriv)%matrix, &
                                      ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)
      END DO

      NULLIFY (fm_struct, weights_fm, ediff_inv, ediff_inv_weights)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt, ncol_global=nmo_occ, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(weights_fm, fm_struct)
      CALL cp_fm_create(ediff_inv, fm_struct)
      CALL cp_fm_create(ediff_inv_weights, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(matrix=evects(1)%matrix, matrix_struct=fm_struct)
      NULLIFY (scaled_evect)
      CALL cp_fm_create(scaled_evect, fm_struct)

      ALLOCATE (weights_2d(nmo_virt, nmo_occ))
      DO imo_occ = 1, nmo_occ
         DO imo_virt = 1, nmo_virt
            weights_2d(imo_virt, imo_occ) = 1.0_dp/(evals_virt(imo_virt)-evals_occ(imo_occ))
         END DO
      END DO
      CALL cp_fm_set_submatrix(ediff_inv, weights_2d)

      ! summary information
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') "TDDFT singlet states"

         WRITE (log_unit, '(3X,A,T14,A,T31,A,T63,A)') "State", "Excitation", "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(3X,A,T14,A,T33,A,T43,A,T53,A,T61,A)') "number", "energy (eV)", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(1X,74("-"))')
      END IF

      DO istate = 1, nstates
         CALL cp_gemm('T', 'N', nmo_virt, nmo_occ, nao, 1.0_dp, S_mos_virt, evects(istate)%matrix, 0.0_dp, weights_fm)

         CALL cp_fm_schur_product(weights_fm, ediff_inv, ediff_inv_weights)
         CALL cp_gemm('N', 'N', nao, nmo_occ, nmo_virt, 1.0_dp, mos_virt, ediff_inv_weights, 0.0_dp, scaled_evect)

         DO ideriv = 1, nderivs
            CALL cp_fm_trace(dS_mos_occ(ideriv)%matrix, scaled_evect, trans_dipole(ideriv))
         END DO
         oscillator_strength = 2.0_dp/3.0_dp*evals(istate)*DOT_PRODUCT(trans_dipole, trans_dipole)

         IF (log_unit > 0) &
            WRITE (log_unit, '(1X,I8,T14,F11.5,T28,3(1X,F9.4),T63,F10.5)') &
            istate, evals(istate)*evolt, trans_dipole(1:nderivs), oscillator_strength
      END DO

      CALL cp_fm_release(scaled_evect)
      CALL cp_fm_release(ediff_inv_weights)
      CALL cp_fm_release(ediff_inv)

      ! excitation analysis
      IF (min_amplitude < 1.0_dp) THEN
         nmo_virt_occ = nmo_virt*nmo_occ
         ALLOCATE (weights_abs_1d(nmo_virt_occ))
         ALLOCATE (inds(nmo_virt_occ))

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') "Excitation analysis"

            WRITE (log_unit, '(3X,A,T15,A,T28,A,T41,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(3X,A,T16,A,T28,A,T41,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,49("-"))')
         END IF

         DO istate = 1, nstates
            CALL cp_gemm('T', 'N', nmo_virt, nmo_occ, nao, 1.0_dp, S_mos_virt, evects(istate)%matrix, 0.0_dp, weights_fm)
            CALL cp_fm_get_submatrix(weights_fm, weights_2d)
            CALL dcopy(nmo_virt_occ, weights_2d, 1, weights_abs_1d, 1)
            weights_abs_1d = ABS(weights_abs_1d)
            CALL sort(weights_abs_1d, nmo_virt_occ, inds)

            IF (log_unit > 0) &
               WRITE (log_unit, '(1X,I8)') istate

            DO iexc = nmo_virt_occ, 1, -1
               IF (weights_abs_1d(iexc) < min_amplitude) EXIT

               imo_occ = (inds(iexc)-1)/nmo_virt+1
               imo_virt = MOD(inds(iexc)-1, nmo_virt)+1

               IF (log_unit > 0) &
                  WRITE (log_unit, '(T14,I8,T26,I8,T42,F9.6)') nmo_occ_total-nmo_occ+imo_occ, &
                  nmo_occ_total+imo_virt, weights_2d(imo_virt, imo_occ)
            END DO
         END DO

         DEALLOCATE (inds)
         DEALLOCATE (weights_abs_1d)
      END IF

      DEALLOCATE (weights_2d)
      CALL cp_fm_release(weights_fm)
      DO ideriv = 1, nderivs
         CALL cp_fm_release(dS_mos_occ(ideriv)%matrix)
      END DO
      CALL cp_fm_release(S_mos_virt)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary


! action matrix related stuff
   SUBROUTINE tddfpt_am_release(action_matrix)
      TYPE(tddfpt_action_matrix), INTENT(inout)          :: action_matrix

      IF (ASSOCIATED(action_matrix%lower_triangle)) &
         CALL cp_fm_release(action_matrix%lower_triangle)

      IF (ALLOCATED(action_matrix%diag)) &
         DEALLOCATE (action_matrix%diag)
   END SUBROUTINE tddfpt_am_release

   FUNCTION tddfpt_am_expand(action_matrix, nstates_total, nvects_new, nvects_max, &
                             blacs_env, para_env) RESULT(added_vectors)
      TYPE(tddfpt_action_matrix), INTENT(inout)          :: action_matrix
      INTEGER, INTENT(in)                                :: nstates_total
      INTEGER, INTENT(in)                                :: nvects_new
      INTEGER, INTENT(in)                                :: nvects_max
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER :: added_vectors

      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: diag_moved
      TYPE(cp_fm_type), POINTER                          :: tril_moved
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      INTEGER :: nvects_exists, nvects_max_local

      IF (ALLOCATED(action_matrix%diag)) THEN
         nvects_exists = SIZE(action_matrix%diag)
      ELSE
         nvects_exists = 0
      END IF

      CPASSERT(nvects_exists == 0 .OR. ASSOCIATED(action_matrix%lower_triangle))

      IF (nvects_max > 0) THEN
         nvects_max_local = nvects_max
      ELSE
         nvects_max_local = nstates_total
      END IF

      added_vectors = nvects_max_local - nvects_exists
      IF (added_vectors > nvects_new) THEN
         added_vectors = nvects_new
      ELSE IF (added_vectors < 0) THEN
         added_vectors = 0
      END IF

      IF (added_vectors > 0) THEN
         IF (nvects_exists > 0) THEN
            ! -- %lower_triangle
            tril_moved => action_matrix%lower_triangle
!            CALL cp_fm_get_info(tril_moved, context=blacs_env, para_env=para_env)
            NULLIFY(action_matrix%lower_triangle, fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nvects_exists+added_vectors, &
                                     context=blacs_env, para_env=para_env)
            CALL cp_fm_create(action_matrix%lower_triangle, fm_struct)
            CALL cp_fm_struct_release(fm_struct)

            CALL cp_fm_set_all(action_matrix%lower_triangle, 0.0_dp)
            CALL cp_fm_to_fm(msource=tril_moved, mtarget=action_matrix%lower_triangle, ncol=nvects_exists)
            CALL cp_fm_release(tril_moved)

            ! -- %diag
            CALL MOVE_ALLOC(action_matrix%diag, diag_moved)
            ALLOCATE(action_matrix%diag(nvects_exists + added_vectors))
            action_matrix%diag(1:nvects_exists) = diag_moved(1:nvects_exists)
            action_matrix%diag(nvects_exists+1:nvects_exists+added_vectors) = 0.0_dp
            DEALLOCATE(diag_moved)
         ELSE
            ! -- %lower_triangle
            NULLIFY(fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=added_vectors, &
                                     context=blacs_env, para_env=para_env)
            CALL cp_fm_create(action_matrix%lower_triangle, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_fm_set_all(action_matrix%lower_triangle, 0.0_dp)

            ! -- %diag
            ALLOCATE(action_matrix%diag(added_vectors))
            action_matrix%diag(:) = 0.0_dp
         END IF

         nvects_exists = nvects_exists + added_vectors
      END IF
   END FUNCTION tddfpt_am_expand

   !  [NO] = new_krylov_vector x old_krylov_vectors
   !  [NN] = new_krylov_vector x new_krylov_vector (upper triangular matrix)
   !  [GN] = guess_vector x new_krylov_vector
   !  1     first_vector
   !  |-----|-----|-----| <- 1
   !  |     | ON= |     |
   !  |     |NO^T |     |
   !  |=====|=====|-----| <- first_vector
   !  | NO  | .NN |     |
   !  |     |   . |     |
   !  |=====|=====|-----|
   !  |     |     |     |
   !  |     |     |     |
   !  |     | GN  |     |
   !  |     |     |     |
   !  |-----|-----|-----|
   SUBROUTINE tddfpt_am_update(action_matrix, block_on, block_nn, block_gn)
      TYPE(tddfpt_action_matrix), INTENT(inout)          :: action_matrix
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in) :: block_on
      REAL(kind=dp), DIMENSION(:,:), INTENT(inout)       :: block_nn
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in) :: block_gn

      INTEGER :: nstates_total, first_vector, nvects_new, ivect, jvect


      nvects_new = SIZE(block_nn, 1)

      first_vector = 1
      IF (ALLOCATED(block_on)) THEN
         first_vector = first_vector + SIZE(block_on, 1)
         CPASSERT(SIZE(block_on, 2) == nvects_new)
      END IF

      nstates_total = first_vector + nvects_new - 1
      CPASSERT(ALLOCATED(action_matrix%diag))
      ! SIZE(action_matrix%diag) == current number of Krylov's vectors
      CPASSERT(SIZE(action_matrix%diag) == nstates_total)

      IF (ALLOCATED(block_gn)) THEN
         nstates_total = nstates_total + SIZE(block_gn)
         CPASSERT(SIZE(block_gn, 2) == nvects_new)
      END IF

      ! update corresponding diagonal elements
      DO ivect = 1,nvects_new
         action_matrix%diag(first_vector+ivect-1) = block_nn(ivect, ivect)
      END DO

      ! zeroise diagonal elements of the block_nn as well its lower triangle
      ! (the lower triangle should already be zero, but do it anyway)
      DO jvect = 1,nvects_new
         DO ivect = jvect,nvects_new
            block_nn(ivect, jvect) = 0.0_dp
         END DO
      END DO

      CALL cp_fm_set_submatrix(action_matrix%lower_triangle, block_nn, &
                               start_row=first_vector, start_col=first_vector, &
                               n_rows=nvects_new, n_cols=nvects_new, transpose=.TRUE.)

      IF (ALLOCATED(block_on)) &
         CALL cp_fm_set_submatrix(action_matrix%lower_triangle, block_on, &
                                  start_row=first_vector, start_col=1, &
                                  n_rows=nvects_new, n_cols=first_vector-1, transpose=.TRUE.)

      IF (ALLOCATED(block_gn)) &
         CALL cp_fm_set_submatrix(action_matrix%lower_triangle, block_gn, &
                                  start_row=first_vector+nvects_new, start_col=first_vector, &
                                  n_rows=SIZE(block_gn, 1), n_cols=nvects_new, transpose=.FALSE.)
   END SUBROUTINE tddfpt_am_update

  ! compute matrix product: C = A * B, where A is the compressed action matrix,
  ! B and C are Krylov's vectors
   SUBROUTINE tddfpt_am_multiply(action_matrix, diag_guess, matrix_b, matrix_c)
      TYPE(tddfpt_action_matrix), INTENT(in)             :: action_matrix
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: diag_guess
      TYPE(cp_fm_type), POINTER                          :: matrix_b, matrix_c

      INTEGER :: nvects_exists, nstates_total, offset, nstates
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: weight_fm
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: a_diag
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY(fm_struct, blacs_env, para_env)
      CALL cp_fm_get_info(matrix_b, ncol_global=nstates, matrix_struct=fm_struct, para_env=para_env, context=blacs_env)

      nstates_total = SIZE(diag_guess)
      IF (ALLOCATED(action_matrix%diag)) THEN
         nvects_exists = SIZE(action_matrix%diag)
      ELSE
         nvects_exists = 0
      END IF

      ALLOCATE(a_diag(nstates_total, 1))
      IF (nvects_exists > 0) &
         a_diag(1:nvects_exists,1) = action_matrix%diag(1:nvects_exists)
      a_diag(nvects_exists+1:nstates_total,1) = diag_guess(nvects_exists+1:nstates_total)

      ! diag(action_matrix) * matrix_b : scale rows of the matrix_b by a_diag(:,1) .
      !
      ! In the absence of the subroutine cp_fm_row_scale()
      ! set all the rows of the matrix weight_fm(i,:) = a_diag(i,1)
      ! and compute the element-wise matrix product.
      NULLIFY(weight_fm)
      CALL cp_fm_create(weight_fm, fm_struct)
      DO offset = 1,nstates
         CALL cp_fm_set_submatrix(weight_fm, a_diag, start_col=offset, n_cols=1)
      END DO

      CALL cp_fm_schur_product(weight_fm, matrix_b, matrix_c)

      CALL cp_fm_release(weight_fm)
      DEALLOCATE(a_diag)

      ! non-diagonal elements
      IF (nvects_exists > 0) THEN
         CALL cp_gemm('N', 'N', nstates_total, nstates, nvects_exists, 1.0_dp, &
                      action_matrix%lower_triangle, matrix_b, 1.0_dp, matrix_c)
         CALL cp_gemm('T', 'N', nvects_exists, nstates, nstates_total, 1.0_dp, &
                      action_matrix%lower_triangle, matrix_b, 1.0_dp, matrix_c)
      END IF
   END SUBROUTINE tddfpt_am_multiply

   SUBROUTINE print_densities(log_unit, first_vector, evects, mos_occ, qs_env, matrix_s)
   USE kahan_sum,                       ONLY: accurate_sum
   USE pw_methods,                      ONLY: pw_integrate_function
      INTEGER, INTENT(in) :: log_unit, first_vector
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evects
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      INTEGER :: ivect, nvects, nao, nmo_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_r, rho_ia_g
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm

      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      REAL(KIND=dp) :: total_rho_rspace, total_rho_gspace

      nvects = SIZE(evects)
      CALL cp_fm_get_info(mos_occ, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env, para_env=para_env)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct, rho_ia_r, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)
      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho_ia_r, rho_g=rho_ia_g, tot_rho_r=tot_rho_r)

      IF (log_unit > 0) &
         WRITE(log_unit, '(/,1X,A,T25,A,T50,A)') 'Vector no.', 'total_rho1_rspace', 'total_rho1_gspace'
      DO ivect = 1,nvects
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, evects(ivect)%matrix, mos_occ, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ, 0.5_dp, mos_occ, evects(ivect)%matrix, 1.0_dp, rho_ia_fm)

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         total_rho_gspace = pw_integrate_function(rho_ia_g(1)%pw, isign=-1)
         total_rho_rspace = accurate_sum(tot_rho_r)

         IF (log_unit > 0) THEN
            WRITE(log_unit, '(1X,I8,T25,ES20.10,T50,ES20.10)') ivect+first_vector-1, total_rho_rspace, total_rho_gspace
         END IF
      END DO

      CALL qs_rho_release(rho_ia_struct)
      CALL cp_fm_release(rho_ia_fm)

   END SUBROUTINE print_densities
END MODULE qs_tddfpt4_methods
