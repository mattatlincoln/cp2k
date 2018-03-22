!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Part of the NEGF environment that needs to be stored within qs_environment.
!>        It is used to update the Hartree potential and the Kohn-Sham matrix within qs_scf().
! **************************************************************************************************
MODULE negf_qs_env_types
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE negf_alloc_types,                ONLY: negf_allocatable_rmatrix
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED,&
                                              pw_grid_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_qs_env_types'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

! **************************************************************************************************
!> \brief ...
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_matrix_block_map_type
      !> column and row atomic indices
      INTEGER                                            :: iatom_col, iatom_row
      !> cell replica
      INTEGER                                            :: image
      !> absolute value is the index of the block within the block data storage;
      !>  0  all matrix elements of the given block are equal to zero,
      !> >0  block is stored simply as block_storage(block_index),
      !> <0  block is stored as TRANSPOSE(block_storage(-block_index))
      INTEGER                                            :: block_index
   END TYPE negf_matrix_block_map_type

! **************************************************************************************************
!> \brief Part of the NEGF environment stored within the QuickStep environment.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_qs_env_contact_type
      REAL(kind=dp), DIMENSION(3)                        :: direction_vector, origin_scattering
      REAL(kind=dp)                                      :: v_external
      !> atom_row, atom_col, image, and the index of the actual data block in 'block_storage'
      !> all indices are identical for both spins
      TYPE(negf_matrix_block_map_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ks_block_map
      !> collection of DBSCR Kohn-Sham matrix blocks ([index, spin])
      TYPE(negf_allocatable_rmatrix), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: ks_block_storage
   END TYPE negf_qs_env_contact_type

! **************************************************************************************************
!> \brief data to update the Hartree potential and the Kohn-Sham matrix within qs_scf().
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_qs_env_type
      LOGICAL                                            :: apply_bias
      LOGICAL                                            :: update_ks_matrix
      TYPE(negf_qs_env_contact_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: contacts
      !> reference Hartree potential (with no bias) !!! no longer needed !!!
      !TYPE(pw_type), POINTER                             :: v_hartree_ref
      !> pw_neighbour_blocks_mepos(bound, axis), bound = 1 (lower), 2 (upper); axis = 1 (x), 2 (y), 3 (z)
      INTEGER, DIMENSION(2, 3)                            :: pw_neighbour_blocks_mepos
      REAL(kind=dp)                                      :: conv_jacobi
   END TYPE negf_qs_env_type

   PUBLIC :: negf_qs_env_type, negf_qs_env_contact_type, negf_matrix_block_map_type
   PUBLIC :: negf_qs_env_create, negf_qs_env_neighbour_mepos, negf_qs_env_release

CONTAINS

! **************************************************************************************************
!> \brief Allocate control options for Non-equilibrium Green's Function calculation.
!> \param negf_qs_env an object to create
!> \par History
!>    * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_qs_env_create(negf_qs_env)
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_qs_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(negf_qs_env))
      CALL timeset(routineN, handle)

      ALLOCATE (negf_qs_env)
      negf_qs_env%apply_bias = .FALSE.
      negf_qs_env%update_ks_matrix = .FALSE.
      !NULLIFY (negf_qs_env%v_hartree_ref)

      CALL timestop(handle)
   END SUBROUTINE negf_qs_env_create

! **************************************************************************************************
!> \brief ...
!> \param negf_qs_env ...
!> \param pw_pool ...
! **************************************************************************************************
   SUBROUTINE negf_qs_env_neighbour_mepos(negf_qs_env, pw_pool)
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_qs_env_neighbour_mepos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group_size, handle, igroup
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: bounds_local
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CALL timeset(routineN, handle)
      pw_grid => pw_pool%pw_grid

      negf_qs_env%pw_neighbour_blocks_mepos(:, :) = -1

      IF (pw_grid%para%mode == PW_MODE_DISTRIBUTED) THEN
         group_size = pw_grid%para%group_size
         ALLOCATE (bounds_local(2, 3, 0:group_size))

         bounds_local(:, :, :) = 0
         bounds_local(1:2, 1:3, pw_grid%para%my_pos) = pw_grid%bounds_local
         CALL mp_sum(bounds_local, pw_grid%para%group)

         DO igroup = 0, group_size-1
            IF (bounds_local(1, 1, igroup) == pw_grid%bounds_local(1, 1) .AND. &
                bounds_local(2, 1, igroup) == pw_grid%bounds_local(2, 1)) THEN

               IF (bounds_local(1, 2, igroup) == pw_grid%bounds_local(1, 2) .AND. &
                   bounds_local(2, 2, igroup) == pw_grid%bounds_local(2, 2)) THEN

                  ! same xy plane
                  IF (bounds_local(2, 3, igroup) == pw_grid%bounds_local(1, 3)-1) &
                     negf_qs_env%pw_neighbour_blocks_mepos(1, 3) = igroup ! -z

                  IF (bounds_local(1, 3, igroup) == pw_grid%bounds_local(2, 3)+1) &
                     negf_qs_env%pw_neighbour_blocks_mepos(2, 3) = igroup ! +z
               END IF

               IF (bounds_local(1, 3, igroup) == pw_grid%bounds_local(1, 3) .AND. &
                   bounds_local(2, 3, igroup) == pw_grid%bounds_local(2, 3)) THEN

                  ! same xz plane
                  IF (bounds_local(2, 2, igroup) == pw_grid%bounds_local(1, 2)-1) &
                     negf_qs_env%pw_neighbour_blocks_mepos(1, 2) = igroup ! -y

                  IF (bounds_local(1, 2, igroup) == pw_grid%bounds_local(2, 2)+1) &
                     negf_qs_env%pw_neighbour_blocks_mepos(2, 2) = igroup ! +y
               END IF
            ELSE IF (bounds_local(2, 1, igroup) == pw_grid%bounds_local(1, 1)-1) THEN
               ! -x if yz plane is the same

               IF (bounds_local(1, 2, igroup) == pw_grid%bounds_local(1, 2) .AND. &
                   bounds_local(2, 2, igroup) == pw_grid%bounds_local(2, 2) .AND. &
                   bounds_local(1, 3, igroup) == pw_grid%bounds_local(1, 3) .AND. &
                   bounds_local(2, 3, igroup) == pw_grid%bounds_local(2, 3)) &
                  negf_qs_env%pw_neighbour_blocks_mepos(1, 1) = igroup ! -x

            ELSE IF (bounds_local(1, 1, igroup) == pw_grid%bounds_local(2, 1)+1) THEN
               ! +x if yz plane is the same

               IF (bounds_local(1, 2, igroup) == pw_grid%bounds_local(1, 2) .AND. &
                   bounds_local(2, 2, igroup) == pw_grid%bounds_local(2, 2) .AND. &
                   bounds_local(1, 3, igroup) == pw_grid%bounds_local(1, 3) .AND. &
                   bounds_local(2, 3, igroup) == pw_grid%bounds_local(2, 3)) &
                  negf_qs_env%pw_neighbour_blocks_mepos(2, 1) = igroup ! +x
            END IF
         END DO

         DEALLOCATE (bounds_local)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_qs_env_neighbour_mepos

! **************************************************************************************************
!> \brief Release memory allocated for NEGF control options.
!> \param negf_qs_env an object to release
!> \par History
!>    * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_qs_env_release(negf_qs_env)
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_qs_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock, icontact, ispin

!      TYPE(pw_env_type), POINTER                         :: pw_env
!      TYPE(pw_pool_type), POINTER                        :: pw_pool

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(negf_qs_env)) THEN
         IF (ALLOCATED(negf_qs_env%contacts)) THEN
            DO icontact = SIZE(negf_qs_env%contacts), 1, -1
               IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_map)) &
                  DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_map)

               IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_storage)) THEN
                  DO ispin = SIZE(negf_qs_env%contacts(icontact)%ks_block_storage, 2), 1, -1
                     DO iblock = SIZE(negf_qs_env%contacts(icontact)%ks_block_storage, 1), 1, -1
                        IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_storage(iblock, ispin)%matrix)) &
                           DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_storage(iblock, ispin)%matrix)
                     END DO
                  END DO

                  DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_storage)
               END IF
            END DO

            DEALLOCATE (negf_qs_env%contacts)
         END IF

         !IF (ASSOCIATED(negf_qs_env%v_hartree_ref)) THEN
         !   CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
         !   CALL pw_pool_give_back_pw(pw_pool, negf_qs_env%v_hartree_ref)
         !END IF

         DEALLOCATE (negf_qs_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_qs_env_release
END MODULE negf_qs_env_types
