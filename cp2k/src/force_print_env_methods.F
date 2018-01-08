!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface for force_print_env
! **************************************************************************************************
MODULE force_print_env_methods
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE force_print_env_types,           ONLY: force_print_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: get_molecule,&
                                              molecule_type
   USE particle_types,                  ONLY: particle_type
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_print_env_methods'

   PUBLIC :: read_force_print_env, update_print_sum_forces

CONTAINS

! **************************************************************************************************
!> \brief Initialise conponents of the given force_print environment based on input parameters.
!> \param force_env force environment
! **************************************************************************************************
   SUBROUTINE read_force_print_env(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'read_force_print_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(force_print_env_type), POINTER                :: force_print_env
      TYPE(section_vals_type), POINTER                   :: root_section, section

      CALL timeset(routineN, handle)

      NULLIFY (force_print_env, root_section, subsys)
      CALL force_env_get(force_env, print_env=force_print_env, root_section=root_section, subsys=subsys)

      section => section_vals_get_subs_vals(root_section, "MOTION%PRINT%SUM_FORCES")

      CALL parse_atom_list(force_print_env%atom_list_fixed, section, subsys)

      force_print_env%forces_fixed(:) = 0.0_dp

      CALL timestop(handle)
   END SUBROUTINE read_force_print_env

! **************************************************************************************************
!> \brief Keep forces acting on selected atoms.
!> \param force_env  force environment
! **************************************************************************************************
   SUBROUTINE update_print_sum_forces(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'update_print_sum_forces', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, natoms
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(force_print_env_type), POINTER                :: force_print_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (force_print_env, subsys)
      CALL force_env_get(force_env, print_env=force_print_env, subsys=subsys)

      IF (ALLOCATED(force_print_env%atom_list_fixed)) THEN
         natoms = SIZE(force_print_env%atom_list_fixed)

         CALL cp_subsys_get(subsys, particle_set=particle_set)

         force_print_env%forces_fixed(:) = 0.0_dp
         DO iatom = 1, natoms
            force_print_env%forces_fixed(:) = force_print_env%forces_fixed(:)+ &
                                              particle_set(force_print_env%atom_list_fixed(iatom))%f(:)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE update_print_sum_forces

! **************************************************************************************************
!> \brief Parse lists of atoms
!> \param atom_list       list of atoms to create
!> \param section         input section that contains LIST and MOLNAME keywords
!> \param subsys          subsystem environment
!> \par History
!>      01.2018 created [SKC]
! **************************************************************************************************
   SUBROUTINE parse_atom_list(atom_list, section, subsys)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: atom_list
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(cp_subsys_type), POINTER                      :: subsys

      CHARACTER(len=*), PARAMETER :: routineN = 'parse_atom_list', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: index_str, natoms_str
      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER                           :: cptr
      INTEGER :: first_atom, handle, iatom, ikind, imol, iname, irep, last_atom, natoms_current, &
         natoms_max, natoms_total, nkinds, nmols, nnames, nrep_list, nrep_molname
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(:), POINTER                     :: iptr
      LOGICAL                                            :: is_list, is_molname
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(molecule_type), POINTER                       :: molecule
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL cp_subsys_get(subsys, particle_set=particle_set, &
                         molecule_set=molecule_set, &
                         molecule_kind_set=molecule_kind_set)
      natoms_max = SIZE(particle_set)
      nkinds = SIZE(molecule_kind_set)

      CALL section_vals_val_get(section, "LIST", explicit=is_list, n_rep_val=nrep_list)
      CALL section_vals_val_get(section, "MOLNAME", explicit=is_molname, n_rep_val=nrep_molname)

      ! compute the number of atoms, and check the validity of atomic indices
      natoms_total = 0
      IF (is_list .AND. nrep_list > 0) THEN
         DO irep = 1, nrep_list
            CALL section_vals_val_get(section, "LIST", i_rep_val=irep, i_vals=iptr)

            natoms_current = SIZE(iptr)
            DO iatom = 1, natoms_current
               IF (iptr(iatom) > natoms_max) THEN
                  CALL integer_to_string(iptr(iatom), index_str)
                  CALL integer_to_string(natoms_max, natoms_str)
                  CALL cp_abort(__LOCATION__, &
                                "Atomic index "//TRIM(index_str)// &
                                " exceeds the maximum number of atoms ("//TRIM(natoms_str)//").")
               END IF
            END DO

            natoms_total = natoms_total+natoms_current
         END DO
      END IF

      IF (is_molname .AND. nrep_molname > 0) THEN
         DO irep = 1, nrep_molname
            CALL section_vals_val_get(section, "MOLNAME", i_rep_val=irep, c_vals=cptr)
            nnames = SIZE(cptr)

            DO iname = 1, nnames
               DO ikind = 1, nkinds
                  IF (molecule_kind_set(ikind)%name .EQ. cptr(iname)) EXIT
               END DO

               IF (ikind <= nkinds) THEN
                  molecule_kind => molecule_kind_set(ikind)
                  CALL get_molecule_kind(molecule_kind, nmolecule=nmols, molecule_list=iptr)

                  DO imol = 1, nmols
                     molecule => molecule_set(iptr(imol))
                     CALL get_molecule(molecule, first_atom=first_atom, last_atom=last_atom)
                     natoms_current = last_atom-first_atom+1
                     natoms_total = natoms_total+natoms_current
                  END DO
               ELSE
                  CALL cp_abort(__LOCATION__, &
                                "A molecule with the name '"//TRIM(cptr(iname))// &
                                " has not been defined. Note that names are case sensitive.")
               END IF
            END DO
         END DO
      END IF

      ! create a list of atomic indices
      IF (natoms_total > 0) THEN
         ALLOCATE (atom_list(natoms_total))

         natoms_total = 0

         IF (is_list .AND. nrep_list > 0) THEN
            DO irep = 1, nrep_list
               CALL section_vals_val_get(section, "LIST", i_rep_val=irep, i_vals=iptr)

               natoms_current = SIZE(iptr)
               atom_list(natoms_total+1:natoms_total+natoms_current) = iptr(1:natoms_current)
               natoms_total = natoms_total+natoms_current
            END DO
         END IF

         IF (is_molname .AND. nrep_molname > 0) THEN
            DO irep = 1, nrep_molname
               CALL section_vals_val_get(section, "MOLNAME", i_rep_val=irep, c_vals=cptr)
               nnames = SIZE(cptr)

               DO iname = 1, nnames
                  DO ikind = 1, nkinds
                     IF (molecule_kind_set(ikind)%name .EQ. cptr(iname)) EXIT
                  END DO

                  IF (ikind <= nkinds) THEN
                     molecule_kind => molecule_kind_set(ikind)
                     CALL get_molecule_kind(molecule_kind, nmolecule=nmols, molecule_list=iptr)

                     DO imol = 1, nmols
                        molecule => molecule_set(iptr(imol))
                        CALL get_molecule(molecule, first_atom=first_atom, last_atom=last_atom)

                        DO natoms_current = first_atom, last_atom
                           natoms_total = natoms_total+1
                           atom_list(natoms_total) = natoms_current
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF

         ! remove duplicated atoms
         ALLOCATE (inds(natoms_total))
         CALL sort(atom_list, natoms_total, inds)
         DEALLOCATE (inds)

         natoms_current = 1
         DO iatom = natoms_current+1, natoms_total
            IF (atom_list(iatom) /= atom_list(natoms_current)) THEN
               natoms_current = natoms_current+1
               atom_list(natoms_current) = atom_list(iatom)
            END IF
         END DO

         IF (natoms_current < natoms_total) THEN
            CALL MOVE_ALLOC(atom_list, inds)

            ALLOCATE (atom_list(natoms_current))
            atom_list(1:natoms_current) = inds(1:natoms_current)
            DEALLOCATE (inds)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE parse_atom_list

END MODULE force_print_env_methods

