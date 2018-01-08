!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data needed to print additional quantities based on force_env
! **************************************************************************************************
MODULE force_print_env_types
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_print_env_types'

   PUBLIC :: force_print_env_type

   PUBLIC :: force_print_env_create, &
             force_print_env_retain, &
             force_print_env_release

! **************************************************************************************************
!> \brief data needed to print additional quantities.
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param atom_list_fixed  list of selected atoms
!> \param forces_fixed     sum of the forces on selected atoms
!> \author SKC
! **************************************************************************************************
   TYPE force_print_env_type
      INTEGER :: ref_count
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_list_fixed
      REAL(kind=dp), DIMENSION(3)                      :: forces_fixed
   END TYPE force_print_env_type

CONTAINS

! **************************************************************************************************
!> \brief create a new force_print env
!> \param force_print_env   force_print environment to create
!> \par History
!>      01.2018 created [SKC]
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE force_print_env_create(force_print_env)
      TYPE(force_print_env_type), POINTER                :: force_print_env

      CHARACTER(len=*), PARAMETER :: routineN = 'force_print_env_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(force_print_env))
      ALLOCATE (force_print_env)

      force_print_env%ref_count = force_print_env%ref_count+1
   END SUBROUTINE force_print_env_create

! **************************************************************************************************
!> \brief retains the given force_print env
!> \param force_print_env   force_print environment to retain
!> \par History
!>      01.2018 created [SKC]
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE force_print_env_retain(force_print_env)
      TYPE(force_print_env_type), POINTER                :: force_print_env

      CHARACTER(len=*), PARAMETER :: routineN = 'force_print_env_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(force_print_env))
      CPASSERT(force_print_env%ref_count > 0)
      force_print_env%ref_count = force_print_env%ref_count+1
   END SUBROUTINE force_print_env_retain

! **************************************************************************************************
!> \brief releases the given force_print env
!> \param force_print_env    force_print environment to release
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE force_print_env_release(force_print_env)
      TYPE(force_print_env_type), POINTER                :: force_print_env

      CHARACTER(len=*), PARAMETER :: routineN = 'force_print_env_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(force_print_env)) THEN
         CPASSERT(force_print_env%ref_count > 0)
         force_print_env%ref_count = force_print_env%ref_count-1
         IF (force_print_env%ref_count == 0) THEN
            IF (ALLOCATED(force_print_env%atom_list_fixed)) &
               DEALLOCATE (force_print_env%atom_list_fixed)

            DEALLOCATE (force_print_env)
         END IF
      END IF
      NULLIFY (force_print_env)
   END SUBROUTINE force_print_env_release
END MODULE force_print_env_types

