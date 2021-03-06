!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/orbital_pointers [1.0] *
!!
!!   NAME
!!     orbital_pointers
!!
!!   FUNCTION
!!      Provides Cartesian and spherical orbital pointers and indices
!!
!!   AUTHOR
!!     Matthias Krack (07.06.2000)
!!
!!   MODIFICATION HISTORY
!!     - reallocate eliminated (17.07.2002,MK)
!!     - restructured and cleaned (20.05.2004,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE orbital_pointers

! *****************************************************************************

! co    : Cartesian orbital pointer for a orbital shell.
! coset : Cartesian orbital pointer for a set of orbitals.
! nco   : Number of Cartesian orbitals for the angular momentum quantum
!         number l.
! ncoset: Number of Cartesian orbitals up to the angular momentum quantum
!         number l.
! nso   : Number of spherical orbitals for the angular momentum quantum
!         number l.
! nsoset: Number of spherical orbitals up to the angular momentum quantum
!         number l.

! *****************************************************************************

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'orbital_pointers'

  INTEGER, SAVE :: current_maxl = -1

  INTEGER, DIMENSION(:), ALLOCATABLE     :: nco,ncoset,nso,nsoset
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: indco,indso
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: so,soset
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: co,coset

! *** Public subroutines ***

  PUBLIC :: deallocate_orbital_pointers,&
            init_orbital_pointers

! *** Public variables ***

  PUBLIC :: co,&
            coset,&
            current_maxl,&
            indco,&
            indso,&
            nco,&
            ncoset,&
            nso,&
            nsoset,&
            so,&
            soset

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE create_orbital_pointers(maxl)

!   Purpose: Allocate and initialize the orbital pointers.

!   History: - Creation (20.05.2004,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: maxl

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_orbital_pointers'

    INTEGER                                  :: iso, istat, l, lx, ly, lz, m

!   ---------------------------------------------------------------------------

    IF (current_maxl > -1) THEN
    END IF

    IF (maxl < 0) THEN
    END IF

!   *** Number of Cartesian orbitals for each l ***

    ALLOCATE (nco(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    nco(-1) = 0

    DO l=0,maxl
      nco(l) = (l + 1)*(l + 2)/2
    END DO

!   *** Number of Cartesian orbitals up to l ***

    ALLOCATE (ncoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    ncoset(-1) = 0

    DO l=0,maxl
      ncoset(l) = ncoset(l-1) + nco(l)
    END DO

!   *** Build the Cartesian orbital pointer and the shell orbital pointer ***

    ALLOCATE (co(0:maxl,0:maxl,0:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    co(:,:,:) = 0

    ALLOCATE (coset(-1:maxl,-1:maxl,-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    coset(:,:,:) = 0

    coset(-1,:,:) = 1
    coset(:,-1,:) = 1
    coset(:,:,-1) = 1

    DO lx=0,maxl
      DO ly=0,maxl
        DO lz=0,maxl
          l = lx + ly + lz
          IF (l > maxl) CYCLE
          co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
          coset(lx,ly,lz) = ncoset(l-1) + co(lx,ly,lz)
        END DO
      END DO
    END DO

    ALLOCATE (indco(3,ncoset(maxl)),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    indco(:,:) = 0

    DO l=0,maxl
      DO lx=0,l
        DO ly=0,l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
        END DO
      END DO
    END DO

!   *** Number of spherical orbitals for each l ***

    ALLOCATE (nso(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    nso(-1) = 0

    DO l=0,maxl
      nso(l) = 2*l + 1
    END DO

!   *** Number of spherical orbitals up to l ***

    ALLOCATE (nsoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    nsoset(-1) = 0

    DO l=0,maxl
      nsoset(l) = nsoset(l-1) + nso(l)
    END DO

    ALLOCATE (indso(2,nsoset(maxl)),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    indso(:,:) = 0

    iso = 0
    DO l=0,maxl
      DO m=-l,l
        iso = iso + 1
        indso(1:2,iso) = (/l,m/)
      END DO
    END DO

    ALLOCATE (so(0:maxl,-maxl:maxl),soset(0:maxl,-maxl:maxl),STAT=istat)
    IF (istat /= 0) THEN
    END IF

    soset(:,:) = 0
    DO l = 0,maxl
      DO m = -l,l
        so(l,m) = nso(l) - (l-m)
        soset(l,m) = nsoset(l-1) + nso(l) - (l-m)
      END DO
    END DO

!   *** Save initialization status ***

    current_maxl = maxl

  END SUBROUTINE create_orbital_pointers

! *****************************************************************************

  SUBROUTINE deallocate_orbital_pointers()

!   Purpose: Deallocate the orbital pointers.

!   History: - Creation (20.05.2005,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_orbital_pointers'

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (current_maxl > -1) THEN

      DEALLOCATE (co,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (coset,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (indco,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (indso,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (nco,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (ncoset,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (nso,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (nsoset,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (so,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      DEALLOCATE (soset,STAT=istat)
      IF (istat /= 0) THEN
      END IF

      current_maxl = -1

    END IF

  END SUBROUTINE deallocate_orbital_pointers

! *****************************************************************************

  SUBROUTINE init_orbital_pointers(maxl)

!   Purpose: Initialize or update the orbital pointers.

!   History: - Creation (07.06.2000,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: maxl

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_orbital_pointers'

!   ---------------------------------------------------------------------------

    IF (maxl < 0) THEN
    END IF

!   *** Check, if the current initialization is sufficient ***

    IF (maxl > current_maxl) THEN
      CALL deallocate_orbital_pointers()
      CALL create_orbital_pointers(maxl)
    END IF

  END SUBROUTINE init_orbital_pointers

! *****************************************************************************

END MODULE orbital_pointers
