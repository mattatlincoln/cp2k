!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09-JUL-2002, TCH, development started
! *****************************************************************************
MODULE qs_tddfpt_utils
  USE cp_control_types,                ONLY: dft_control_type,&
                                             tddfpt_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                             tddfpt_env_deallocate,&
                                             tddfpt_env_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_utils'
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

! *****************************************************************************
  TYPE simple_solution_sorter
     INTEGER                               :: orbit
     INTEGER                               :: lumo
     REAL(KIND=DP)                        :: value
     TYPE(simple_solution_sorter), POINTER :: next
  END TYPE simple_solution_sorter

  PRIVATE

  ! METHODS
  PUBLIC :: tddfpt_cleanup,&
            tddfpt_init,&
            co_initial_guess,&
            find_contributions,&
            normalize,&
            reorthogonalize

CONTAINS

! *****************************************************************************
!> \brief Initialize some necessary structures for a tddfpt calculation.
!> \param p_env perturbation environment to be initialized
!> \param t_env tddfpt environment to be initialized
!> \param qs_env Quickstep environment with the results of a
!>                   ground state calcualtion
! *****************************************************************************
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_env_type), INTENT(out)       :: t_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init', &
      routineP = moduleN//':'//routineN

!------------------!
! create the p_env !
!------------------!

    CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.)
    CALL p_env_psi0_changed(p_env, qs_env) ! update the m_epsilon matrix

    !------------------!
    ! create the t_env !
    !------------------!
    CALL tddfpt_env_allocate(t_env, p_env, qs_env)
    CALL tddfpt_env_init(t_env, qs_env)

  END SUBROUTINE tddfpt_init

! *****************************************************************************
!> \brief Initialize t_env with meaningfull values.
!> \param t_env ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE tddfpt_env_init(t_env, qs_env)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_spins, spin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(matrix_s, dft_control)

    CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control)
    n_spins = dft_control%nspins
    IF (dft_control%tddfpt_control%invert_S) THEN
       DO spin=1, n_spins
          CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, t_env%invS(spin)%matrix)
          CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix)
          CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix)
       END DO
    END IF

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************
!> \brief ...
!> \param t_env ...
!> \param p_env ...
! *****************************************************************************
  SUBROUTINE tddfpt_cleanup(t_env, p_env)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_cleanup', &
      routineP = moduleN//':'//routineN

    CALL tddfpt_env_deallocate(t_env)
    CALL p_env_release(p_env)

  END SUBROUTINE tddfpt_cleanup

! *****************************************************************************
!> \brief ...
!> \param matrices ...
!> \param energies ...
!> \param n_v ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE co_initial_guess(matrices, energies, n_v, qs_env)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: matrices
    REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: energies
    INTEGER, INTENT(IN)                      :: n_v
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'co_initial_guess', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_cols, n_lumos, n_orbits, &
                                                n_rows, n_spins, oo, spin, vo
    REAL(KIND=DP)                            :: evd
    REAL(KIND=DP), ALLOCATABLE, &
      DIMENSION(:, :)                        :: guess, lumos
    REAL(KIND=DP), DIMENSION(:), POINTER     :: orbital_eigenvalues
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(simple_solution_sorter), POINTER    :: sorter_iterator, &
                                                sorter_pointer, sorter_start
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

! number of vectors to initialize

    NULLIFY(dft_control)

    CALL get_qs_env(qs_env, dft_control=dft_control)
    tddfpt_control => dft_control%tddfpt_control
    n_spins = dft_control%nspins
    energies = 0.0_dp

    IF (.NOT.ASSOCIATED(tddfpt_control%lumos)) THEN
       CPABORT("LUMOS missing")
    END IF

    DO spin=1, n_spins

       n_cols = matrices(1,spin)%matrix%matrix_struct%ncol_global
       n_rows = matrices(1,spin)%matrix%matrix_struct%nrow_global

       DO i=1, n_v
          CALL cp_fm_set_all(matrices(i,spin)%matrix, 0.0_dp)
       END DO

       CALL get_mo_set(qs_env%mos(spin)%mo_set, eigenvalues=orbital_eigenvalues)

       n_lumos = tddfpt_control%lumos(spin)%matrix%matrix_struct%ncol_global

       n_orbits = SIZE(orbital_eigenvalues)

       !-----------------------------------------!
       ! create a SORTED list of initial guesses !
       !-----------------------------------------!
       ! first element
       evd = tddfpt_control%lumos_eigenvalues(1,spin)-orbital_eigenvalues(n_orbits)
       ALLOCATE(sorter_start)
       sorter_start%orbit=n_orbits
       sorter_start%lumo=1
       sorter_start%value=evd
       NULLIFY(sorter_start%next)
       ! rest of the elements
       DO oo=n_orbits, 1, -1
          DO vo=1, n_lumos

             IF (oo==n_orbits .AND. vo==1) CYCLE ! already in list

             evd = tddfpt_control%lumos_eigenvalues(vo,spin)-orbital_eigenvalues(oo)

             sorter_iterator => sorter_start
             NULLIFY(sorter_pointer)
             DO WHILE (ASSOCIATED(sorter_iterator%next))
                IF (sorter_iterator%next%value > evd) THEN
                   sorter_pointer => sorter_iterator%next
                   EXIT
                END IF
                sorter_iterator => sorter_iterator%next
             END DO

             ALLOCATE(sorter_iterator%next)
             sorter_iterator%next%orbit = oo
             sorter_iterator%next%lumo  = vo
             sorter_iterator%next%value = evd
             sorter_iterator%next%next  => sorter_pointer

          END DO
       END DO

       ALLOCATE(lumos(n_rows, n_lumos), guess(n_rows, n_orbits))
       CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                start_col=1, n_cols=n_lumos)

       !-------------------!
       ! fill the matrices !
       !-------------------!
       sorter_iterator=>sorter_start
       DO i=1, MIN(n_v,n_orbits*n_lumos)
          CALL dcopy(n_rows*n_orbits, 0.0_dp, 0, guess, 1)
          CALL dcopy(n_rows, lumos(:,sorter_iterator%lumo), 1, &
                             guess(:,sorter_iterator%orbit),1)
          CALL cp_fm_set_submatrix(matrices(i,spin)%matrix, &
                                   guess(:,1:n_orbits))
          energies(i) = energies(i) + sorter_iterator%value/REAL(n_spins,dp)
          sorter_iterator=>sorter_iterator%next
       END DO
       IF (n_v > n_orbits*n_lumos) THEN
          DO i=n_orbits*n_lumos+1, n_v
             CALL cp_fm_init_random(matrices(i,spin)%matrix, n_orbits)
             energies(i)=1.0E38_dp
          END DO
       END IF

       !--------------!
       ! some cleanup !
       !--------------!
       DEALLOCATE(lumos, guess)
       sorter_iterator=>sorter_start
       DO WHILE(ASSOCIATED(sorter_iterator))
          sorter_pointer=>sorter_iterator
          sorter_iterator=>sorter_iterator%next
          DEALLOCATE(sorter_pointer)
       END DO

    END DO

  END SUBROUTINE co_initial_guess

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param t_env ...
! *****************************************************************************
  SUBROUTINE find_contributions(qs_env, t_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(tddfpt_env_type), INTENT(IN)        :: t_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_contributions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, n_ev, n_spins, occ, &
                                                output_unit, spin, virt, n_lumos, nmo,& 
                                                nao, component
    INTEGER, DIMENSION(2)                    :: nhomos, nlumos, nrows
    REAL(KIND=dp)                            :: contribution, &
                                                summed_contributions, xint, yint, zint,&
                                                osx,osy,osz
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: homo_coeff_col, lumo_coeff_col
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: S_lumos,  matrix_elements
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(tddfpt_control_type)                :: t_control
    TYPE(cp_fm_struct_type), POINTER         :: me_struct, template
    TYPE(cp_fm_type), POINTER                :: mo_coeff, buf, test

    NULLIFY(S_lumos, matrix_s, dft_control, mo_coeff, buf, me_struct, matrix_elements, test)
    output_unit= cp_logger_get_default_io_unit()
    CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control)

    IF (output_unit>0) WRITE (output_unit,*)
    IF (output_unit>0) WRITE (output_unit,*)

    t_control = dft_control%tddfpt_control
    n_ev      = t_control%n_ev
    n_spins   = dft_control%nspins

    ALLOCATE(matrix_elements(3*n_spins))

    DO spin=1, n_spins

       CALL get_mo_set(qs_env%mos(spin)%mo_set, nao=nao, nmo=nmo, mo_coeff=mo_coeff)
       n_lumos = t_control%lumos(spin)%matrix%matrix_struct%ncol_global

       IF (output_unit>0) THEN
          WRITE (output_unit,*) "Calculating momentum operator matrix elements for spin", spin
          WRITE (output_unit,*) "Found", nmo,"occupied and ", n_lumos,"unoccupied orbitals"
       END IF

       template => t_control%lumos(spin)%matrix%matrix_struct
       CALL cp_fm_create(buf,template,"buffer")
       CALL cp_fm_struct_create(me_struct,para_env=template%para_env,context=template%context, &
         nrow_global=nmo,ncol_global=n_lumos)

       DO component=1,3
          CALL cp_fm_create(matrix_elements(component+(spin-1)*3)%matrix,me_struct,"matrix_elements")
       END DO

       DO component=1,3
          IF (output_unit>0) THEN
             WRITE (output_unit,*) "Calculating momentum operator matrix elements for spin", spin, &
               "component", component
          END IF
          CALL cp_dbcsr_sm_fm_multiply(matrix_s(component+1)%matrix, t_control%lumos(spin)%matrix ,buf,n_lumos)
          CALL cp_gemm("T","N",nmo,n_lumos,nao,1.0_dp,mo_coeff,buf,0.0_dp,matrix_elements(component+(spin-1)*3)%matrix)
       END DO

       CALL cp_fm_struct_release(me_struct)
       CALL cp_fm_release(buf)
       NULLIFY(mo_coeff)

    END DO


    ALLOCATE(S_lumos(n_spins))

    DO spin=1, n_spins
       nrows(spin) = t_control%lumos(spin)%matrix%matrix_struct%nrow_global
       nhomos(spin) = t_env%evecs(1,spin)%matrix%matrix_struct%ncol_global
       nlumos(spin) = t_control%lumos(spin)%matrix%matrix_struct%ncol_global
       CALL cp_fm_create(S_lumos(spin)%matrix, t_control%lumos(spin)%matrix%matrix_struct, &
                         "S times lumos")
       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, t_control%lumos(spin)%matrix, &
                              S_lumos(spin)%matrix, nlumos(spin), 1.0_dp, 0.0_dp)
    END DO

    ALLOCATE(homo_coeff_col(MAXVAL(nrows(1:n_spins)),1), &
             lumo_coeff_col(MAXVAL(nrows(1:n_spins)),1))
    DO i=1, n_ev
       IF (output_unit>0) THEN
          WRITE (output_unit,'(A,I3,5X,F15.6)') "  excited state : ", i, t_env%evals(i)*evolt
          WRITE (output_unit,*)
       END IF
       summed_contributions = 0.0_dp
       osx = 0.0
       osy = 0.0
       osz = 0.0
       DO spin=1, n_spins
          IF (n_spins == 2) THEN
             IF (spin == 1) THEN
                IF (output_unit>0) WRITE (output_unit,*) 'alpha:'
             ELSE
                IF (output_unit>0) WRITE (output_unit,*) 'beta:'
             END IF
             template => t_control%lumos(spin)%matrix%matrix_struct
          END IF
          searchloop: DO occ=nhomos(spin), 1, -1
             CALL cp_fm_get_submatrix(t_env%evecs(i,spin)%matrix, homo_coeff_col, &
                                      1, occ, nrows(spin), 1)
             DO virt=1, nlumos(spin)
                CALL cp_fm_get_submatrix(S_lumos(spin)%matrix, lumo_coeff_col, &
                                         1, virt, nrows(spin), 1)
                contribution = 0.0_dp
                DO j=1, nrows(spin)
                   contribution = contribution + homo_coeff_col(j,1)*lumo_coeff_col(j,1)
                END DO
                summed_contributions = summed_contributions + (contribution)**2
                CALL cp_fm_get_element(matrix_elements(1+(spin-1)*3)%matrix,occ,virt, xint)
                CALL cp_fm_get_element(matrix_elements(2+(spin-1)*3)%matrix,occ,virt, yint)
                CALL cp_fm_get_element(matrix_elements(3+(spin-1)*3)%matrix,occ,virt, zint)
                osx = osx + ((xint)**2)*(contribution)**2 
                osy = osy + ((yint)**2)*(contribution)**2 
                osz = osz + ((zint)**2)*(contribution)**2 
                IF (ABS(contribution) > 5.0e-2_dp) THEN
                   IF (output_unit>0) THEN
                      WRITE (output_unit,'(14X,I5,A,I5,10X,F8.3,5X,F8.3)') &
                        occ, " ->", nhomos(spin)+virt, ABS(contribution), summed_contributions

                      WRITE (output_unit,'(14X,A,5X,F8.3,5X,F8.3,5X,F8.3)') &
                           "oscillator strength", xint**2, yint**2, zint**2
                   END IF
                END IF
                IF (ABS(summed_contributions - 1.0_dp) < 1.0e-3_dp) THEN 
                   CYCLE searchloop
                END IF
             END DO
          END DO searchloop
          WRITE (output_unit,'(14X,A,5X,F8.3,5X,F8.3,5X,F8.3)') &
             "total oscillator strength", osx, osy,osz
       END DO
       IF (output_unit>0) WRITE (output_unit,*)
    END DO

    !
    ! punch a checksum for the regs
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E14.6)') ' TDDFPT : CheckSum  =',SQRT(SUM(t_env%evals**2))
    ENDIF

    DO spin=1, n_spins
       CALL cp_fm_release(S_lumos(spin)%matrix)
    END DO
    DEALLOCATE(homo_coeff_col, lumo_coeff_col)

    DEALLOCATE(S_lumos)

  END SUBROUTINE find_contributions

! *****************************************************************************
!> \brief ...
!> \param X ...
!> \param tmp_vec ...
!> \param metric ...
! *****************************************************************************
  SUBROUTINE normalize(X, tmp_vec, metric)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: x, tmp_vec
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: metric

    INTEGER                                  :: n_spins, spin
    REAL(KIND=dp)                            :: norm, tmp

    n_spins = SIZE(x)
    norm = 0.0_dp

    DO spin=1, n_spins
       tmp = 0.0_dp
       CALL cp_dbcsr_sm_fm_multiply(metric(1)%matrix, X(spin)%matrix, &
            tmp_vec(spin)%matrix, &
            X(spin)%matrix%matrix_struct%ncol_global, &
            1.0_dp, 0.0_dp)
       CALL cp_fm_trace(X(spin)%matrix, tmp_vec(spin)%matrix, tmp)
       norm = norm + tmp
    END DO

    norm = SQRT(norm)
    DO spin=1, n_spins
       CALL cp_fm_scale((1.0_dp/norm), X(spin)%matrix)
    END DO

  END SUBROUTINE normalize

  !---------------------------------------!
  ! x must not be changed in this routine !
  ! tmp_vec may be changed                !
  !---------------------------------------!
! *****************************************************************************
!> \brief ...
!> \param X ...
!> \param V_set ...
!> \param SV_set ...
!> \param work ...
!> \param n ...
! *****************************************************************************
  SUBROUTINE reorthogonalize(X, V_set, SV_set, work, n)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: X
    TYPE(cp_fm_p_type), DIMENSION(:, :)      :: V_set, SV_set
    TYPE(cp_fm_p_type), DIMENSION(:)         :: work
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reorthogonalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, n_spins, spin
    REAL(DP)                                 :: dot_product, tmp

    CALL timeset(routineN,handle)

    IF (n>0) THEN

       n_spins = SIZE(X)
       DO spin=1, n_spins
          CALL cp_fm_to_fm(X(spin)%matrix, work(spin)%matrix)
       END DO

       DO i=1, n
          dot_product = 0.0_dp
          DO spin=1, n_spins
             CALL cp_fm_trace(SV_set(i,spin)%matrix, work(spin)%matrix, tmp)
             dot_product = dot_product + tmp
          END DO
          DO spin=1, n_spins
             CALL cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                      -1.0_dp*dot_product, V_set(i,spin)%matrix)
          END DO
       END DO

    END IF

    CALL timestop(handle)

  END SUBROUTINE reorthogonalize

! *****************************************************************************

END MODULE qs_tddfpt_utils
