! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This module contains all subroutines used in the main program                  !
!                            to search for overlapping configurations.                            !
!                                                                                                 !
! Version number: 1.4.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                         May 15th, 2024                                          !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE OverlapCheck

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! Make subroutines public
PUBLIC :: ListOverlapCheck, FullListOverlapCheck, ListOverlapCheckInitialConfiguration

CONTAINS

! *********************************************************************************************** !
! This subroutine checks if a random displacement (translation or rotation) of a fixed particle i !
!                           causes any overlaps with other particles j                            !
! *********************************************************************************************** !
SUBROUTINE ParticleOverlapCheck( iComponent, iParticle, iQuaternion, iOrientation, iPosition, ContactDistance, CurrentBoxLength, &
&                                CurrentBoxLengthInverse, Overlap )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap        = .FALSE.
ParallelSPC    = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

! Anisomorphic and isomorphic components
DO jComponent = 1, nComponents
  ! Stop if a thread founds an overlap (should never happen)
  IF( SharedOverlap ) STOP "ERROR: Mishap in the overlap check algorithm. Overlap detected but did not return control."
  !###############################################################################################!
  !$OMP PARALLEL DO DEFAULT( Shared ) &                                                           !
  !$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, VectorDistance ) &              !
  !$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, SquaredCutoffDistance ) &               !
  !$OMP PRIVATE( ContactDistance, OverlapSPC, ParallelSPC ) &                                     !
  !$OMP REDUCTION( .OR. : PrivateOverlap )                                                        !
  !###############################################################################################!
  DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
    ! Cycle if a thread founds an overlap
    IF( SharedOverlap ) CYCLE
    ! Cycle if particle index j = i (prevents self counting) 
    IF( jParticle == iParticle ) CYCLE
    ! Position of particle j
    jPosition(1) = pPositionMC(1,jParticle)
    jPosition(2) = pPositionMC(2,jParticle)
    jPosition(3) = pPositionMC(3,jParticle)
    ! Orientation of particle j
    jOrientation(1) = pOrientationMC(1,jParticle)
    jOrientation(2) = pOrientationMC(2,jParticle)
    jOrientation(3) = pOrientationMC(3,jParticle)
    ! Quaternion of particle j
    jQuaternion(0) = pQuaternionMC(0,jParticle)
    jQuaternion(1) = pQuaternionMC(1,jParticle)
    jQuaternion(2) = pQuaternionMC(2,jParticle)
    jQuaternion(3) = pQuaternionMC(3,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(1) = jPosition(1) - iPosition(1)
    VectorDistance(2) = jPosition(2) - iPosition(2)
    VectorDistance(3) = jPosition(3) - iPosition(3)
    ! Minimum image convention
    CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Cutoff distance (squared)
    SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
    SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
    ! Preliminary test (circumscribing spheres)
    IF( SquaredDistance <= SquaredCutoffDistance ) THEN
      ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
      IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
        CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, ContactDistance, &
        &                     PrivateOverlap )
        ! Overlap found
        IF( PrivateOverlap ) SharedOverlap = .TRUE.
      ! Overlap test for spherocylinders (Vega-Lago method)
      ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
        &                     ContactDistance, ParallelSPC, PrivateOverlap )
        ! Overlap found
        IF( PrivateOverlap ) SharedOverlap = .TRUE.
      ! Overlap test for cylinders and/or spheres
      ELSE IF( GeometryType(3) ) THEN ! Cylinders
        IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
          ! Initialization
          OverlapSPC = .FALSE.
          ! Preliminary test (circumscribing spherocylinders)
          CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
          &                     ContactDistance, ParallelSPC, OverlapSPC )
          ! Overlap criterion
          IF( OverlapSPC ) THEN
            ! Apply periodic boundary conditions on the position of particle j
            jPosition(1) = iPosition(1) + VectorDistance(1)
            jPosition(2) = iPosition(2) + VectorDistance(2)
            jPosition(3) = iPosition(3) + VectorDistance(3)
            ! Overlap test for cylinders (modified algorithm of Lopes et al.)
            CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
            &                     iComponent, jComponent, ParallelSPC, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          END IF
        ! Overlap test for cylinders and spheres
        ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
          ! Apply periodic boundary conditions on the position of particle j
          jPosition(1) = iPosition(1) + VectorDistance(1)
          jPosition(2) = iPosition(2) + VectorDistance(2)
          jPosition(3) = iPosition(3) + VectorDistance(3)
          CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, PrivateOverlap )
          ! Overlap found
          IF( PrivateOverlap ) SharedOverlap = .TRUE.
        ! Overlap test for cylinders and spheres
        ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
          ! Apply periodic boundary conditions on the position of particle j
          jPosition(1) = iPosition(1) + VectorDistance(1)
          jPosition(2) = iPosition(2) + VectorDistance(2)
          jPosition(3) = iPosition(3) + VectorDistance(3)
          CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, PrivateOverlap )
          ! Overlap found
          IF( PrivateOverlap ) SharedOverlap = .TRUE.
        ! Overlap test for spheres
        ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
          PrivateOverlap = .TRUE.
          ! Overlap found
          IF( PrivateOverlap ) SharedOverlap = .TRUE.
        END IF
      END IF
    END IF
  END DO
  !###############################################################################################!
  !$OMP END PARALLEL DO                                                                           !
  !###############################################################################################!
  ! Returns control to the calling program unit if an overlap is detected
  IF( SharedOverlap ) THEN
    Overlap = .TRUE.
    RETURN
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ParticleOverlapCheck

! *********************************************************************************************** !
!   This subroutine uses lists to check if a random displacement (translation or rotation) of a   !
!                   fixed particle i causes any overlaps with other particles j                   !
! *********************************************************************************************** !
SUBROUTINE ListOverlapCheck( iComponent, iParticle, iQuaternion, iOrientation, iPosition, ContactDistance, CurrentBoxLength, &
&                            CurrentBoxLengthInverse, Overlap, CellHalfLogical )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: Neighbours, CellIndex
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                          :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 )                          :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 )                          :: jList                  ! Counters (neighbour list)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndex             ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList         ! List of neighbour j particles

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapHER      ! Detects overlap between two ellipsoids of revolution: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC      ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL      ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC     ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
jList       = 0
Overlap     = .FALSE.
OverlapHER  = .FALSE.
OverlapCYL  = .FALSE.
ParallelSPC = .FALSE.

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndex = CellIndex( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogical .AND. ANY( iCellIndex(:) /= pCellIndex(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndex(:), "] and [", pCellIndex(:,iParticle), "]"
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! List of neighbours
jNeighbourList = Neighbours( iParticle, iCellIndex, CellHalfLogical )

! Loop until no more entries in the neighbour list
DO
  ! Next entry
  jList = jList + 1
  ! Neighbour particle index
  jParticle = jNeighbourList(jList)
  ! List exhausted
  IF( jParticle == 0 ) EXIT
  ! Skip self (should never happen)
  IF( iParticle == jParticle ) CYCLE
  ! Neighbour component index
  jComponent = pComponents(jParticle)
  ! Position of particle j
  jPosition(1) = pPositionMC(1,jParticle)
  jPosition(2) = pPositionMC(2,jParticle)
  jPosition(3) = pPositionMC(3,jParticle)
  ! Orientation of particle j
  jOrientation(1) = pOrientationMC(1,jParticle)
  jOrientation(2) = pOrientationMC(2,jParticle)
  jOrientation(3) = pOrientationMC(3,jParticle)
  ! Quaternion of particle j
  jQuaternion(0) = pQuaternionMC(0,jParticle)
  jQuaternion(1) = pQuaternionMC(1,jParticle)
  jQuaternion(2) = pQuaternionMC(2,jParticle)
  jQuaternion(3) = pQuaternionMC(3,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(1) = jPosition(1) - iPosition(1)
  VectorDistance(2) = jPosition(2) - iPosition(2)
  VectorDistance(3) = jPosition(3) - iPosition(3)
  ! Minimum image convention
  CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Cutoff distance (squared)
  SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
  SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
  ! Preliminary test (circumscribing spheres)
  IF( SquaredDistance <= SquaredCutoffDistance ) THEN
    ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
    IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
      CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, ContactDistance, &
      &                     OverlapHER )
      ! Overlap criterion
      IF( OverlapHER ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    ! Overlap test for spherocylinders (Vega-Lago method)
    ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
      &                     ContactDistance, ParallelSPC, OverlapSPC )
      ! Overlap criterion
      IF( OverlapSPC ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    ! Overlap test for cylinders and/or spheres
    ELSE IF( GeometryType(3) ) THEN ! Cylinders
      IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Initialization
        OverlapSPC = .FALSE.
        ! Preliminary test (circumscribing spherocylinders)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
        &                     ContactDistance, ParallelSPC, OverlapSPC )
        ! Overlap criterion
        IF( OverlapSPC ) THEN
          ! Apply periodic boundary conditions on the position of particle j
          jPosition(1) = iPosition(1) + VectorDistance(1)
          jPosition(2) = iPosition(2) + VectorDistance(2)
          jPosition(3) = iPosition(3) + VectorDistance(3)
          ! Overlap test for cylinders (modified algorithm of Lopes et al.)
          CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
          &                     iComponent, jComponent, ParallelSPC, OverlapCYL )
          ! Overlap criterion
          IF( OverlapCYL ) THEN
            ! Returns control to the calling program unit if an overlap is detected
            Overlap = .TRUE.
            RETURN
          END IF
        END IF
      ! Overlap test for cylinders and spheres
      ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, OverlapCYL )
        IF( OverlapCYL ) THEN
          ! Returns control to the calling program unit if an overlap is detected
          Overlap = .TRUE.
          RETURN
        END IF
      ! Overlap test for cylinders and spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, OverlapCYL )
        IF( OverlapCYL ) THEN
          ! Returns control to the calling program unit if an overlap is detected
          Overlap = .TRUE.
          RETURN
        END IF
      ! Overlap test for spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    END IF
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ListOverlapCheck

! *********************************************************************************************** !
!          This subroutine checks if a random volume scaling (isotropic or anisotropic)           !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullOverlapCheck( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, Overlap )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap        = .FALSE.
OverlapSPC     = .FALSE.
ParallelSPC    = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

! Anisomorphic and isomorphic components
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Stop if a thread founds an overlap (should never happen)
    IF( SharedOverlap ) STOP "ERROR: Mishap in the overlap check algorithm. Overlap detected but did not return control."
    ! Cycle if component index j < i (prevents double counting) 
    IF( jComponent < iComponent ) CYCLE
    !#############################################################################################!
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                           !
    !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &     !
    !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &           !
    !$OMP PRIVATE( SquaredDistance, SquaredCutoffDistance, ContactDistance ) &                    !
    !$OMP PRIVATE( OverlapSPC, ParallelSPC ) &                                                    !
    !$OMP REDUCTION( .OR. : PrivateOverlap )                                                      !
    !#############################################################################################!
    ! First loop represents all particles with indexes i of component Ci
    DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Second loop represents all particles with indexes j of component Cj
      DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
        ! Cycle if a thread founds an overlap
        IF( SharedOverlap ) CYCLE
        ! Cycle if particle index j <= i when components are isomorphic (prevents double counting) 
        IF( jComponent == iComponent .AND. jParticle <= iParticle ) CYCLE
        ! Position of particle i
        iPosition(1) = pPositionMC(1,iParticle)
        iPosition(2) = pPositionMC(2,iParticle)
        iPosition(3) = pPositionMC(3,iParticle)
        ! Position of particle j
        jPosition(1) = pPositionMC(1,jParticle)
        jPosition(2) = pPositionMC(2,jParticle)
        jPosition(3) = pPositionMC(3,jParticle)
        ! Orientation of particle i
        iOrientation(1) = pOrientationMC(1,iParticle)
        iOrientation(2) = pOrientationMC(2,iParticle)
        iOrientation(3) = pOrientationMC(3,iParticle)
        ! Orientation of particle j
        jOrientation(1) = pOrientationMC(1,jParticle)
        jOrientation(2) = pOrientationMC(2,jParticle)
        jOrientation(3) = pOrientationMC(3,jParticle)
        ! Quaternion of particle i
        iQuaternion(0) = pQuaternionMC(0,iParticle)
        iQuaternion(1) = pQuaternionMC(1,iParticle)
        iQuaternion(2) = pQuaternionMC(2,iParticle)
        iQuaternion(3) = pQuaternionMC(3,iParticle)
        ! Quaternion of particle j
        jQuaternion(0) = pQuaternionMC(0,jParticle)
        jQuaternion(1) = pQuaternionMC(1,jParticle)
        jQuaternion(2) = pQuaternionMC(2,jParticle)
        jQuaternion(3) = pQuaternionMC(3,jParticle)
        ! Vector distance between particles i and j
        VectorDistance(1) = jPosition(1) - iPosition(1)
        VectorDistance(2) = jPosition(2) - iPosition(2)
        VectorDistance(3) = jPosition(3) - iPosition(3)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (squared)
        SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
        SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
        ! Preliminary test (circumscribing spheres)
        IF( SquaredDistance <= SquaredCutoffDistance ) THEN
          ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
          IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
            CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          ! Overlap test for spherocylinders (Vega-Lago method)
          ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, ParallelSPC, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          ! Overlap test for cylinders and/or spheres
          ELSE IF( GeometryType(3) ) THEN ! Cylinders
            IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Initialization
              OverlapSPC = .FALSE.
              ! Preliminary test (circumscribing spherocylinders)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
              &                     ContactDistance, ParallelSPC, OverlapSPC )
              ! Overlap criterion
              IF( OverlapSPC ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                jPosition(1) = iPosition(1) + VectorDistance(1)
                jPosition(2) = iPosition(2) + VectorDistance(2)
                jPosition(3) = iPosition(3) + VectorDistance(3)
                ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
                &                     iComponent, jComponent, ParallelSPC, PrivateOverlap )
                ! Overlap found
                IF( PrivateOverlap ) SharedOverlap = .TRUE.
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            ! Overlap test for cylinders and spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            ! Overlap test for spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              PrivateOverlap = .TRUE.
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            END IF
          END IF
        END IF
      END DO
    END DO
    !#############################################################################################!
    !$OMP END PARALLEL DO                                                                         !
    !#############################################################################################!
    ! Returns control to the calling program unit if an overlap is detected
    IF( SharedOverlap ) THEN
      Overlap = .TRUE.
      RETURN
    END IF
  END DO
END DO

RETURN ! No overlaps detected

END SUBROUTINE FullOverlapCheck

! *********************************************************************************************** !
!          This subroutine checks if a random volume scaling (isotropic or anisotropic)           !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullOverlapCheckFixPFraction( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, Overlap )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap        = .FALSE.
OverlapSPC     = .FALSE.
ParallelSPC    = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

! Anisomorphic components
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Stop if a thread founds an overlap (should never happen)
    IF( SharedOverlap ) STOP "ERROR: Mishap in the overlap check algorithm. Overlap detected but did not return control."
    ! Cycle if component index j < i (prevents double counting) 
    IF( jComponent < iComponent ) CYCLE
    !#############################################################################################!
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                           !
    !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &     !
    !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &           !
    !$OMP PRIVATE( SquaredDistance, SquaredCutoffDistance, ContactDistance ) &                    !
    !$OMP PRIVATE( OverlapSPC, ParallelSPC ) &                                                    !
    !$OMP REDUCTION( .OR. : PrivateOverlap )                                                      !
    !#############################################################################################!
    ! First loop represents all particles with indexes i of component Ci
    DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Second loop represents all particles with indexes j of component Cj
      DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
        ! Cycle if a thread founds an overlap
        IF( SharedOverlap ) CYCLE
        ! Cycle if particle index j <= i when components are isomorphic (prevents double counting) 
        IF( jComponent == iComponent .AND. jParticle <= iParticle ) CYCLE
        ! Position of particle i
        iPosition(1) = pPosition(1,iParticle)
        iPosition(2) = pPosition(2,iParticle)
        iPosition(3) = pPosition(3,iParticle)
        ! Position of particle j
        jPosition(1) = pPosition(1,jParticle)
        jPosition(2) = pPosition(2,jParticle)
        jPosition(3) = pPosition(3,jParticle)
        ! Orientation of particle i
        iOrientation(1) = pOrientation(1,iParticle)
        iOrientation(2) = pOrientation(2,iParticle)
        iOrientation(3) = pOrientation(3,iParticle)
        ! Orientation of particle j
        jOrientation(1) = pOrientation(1,jParticle)
        jOrientation(2) = pOrientation(2,jParticle)
        jOrientation(3) = pOrientation(3,jParticle)
        ! Quaternion of particle i
        iQuaternion(0) = pQuaternion(0,iParticle)
        iQuaternion(1) = pQuaternion(1,iParticle)
        iQuaternion(2) = pQuaternion(2,iParticle)
        iQuaternion(3) = pQuaternion(3,iParticle)
        ! Quaternion of particle j
        jQuaternion(0) = pQuaternion(0,jParticle)
        jQuaternion(1) = pQuaternion(1,jParticle)
        jQuaternion(2) = pQuaternion(2,jParticle)
        jQuaternion(3) = pQuaternion(3,jParticle)
        ! Vector distance between particles i and j
        VectorDistance(1) = jPosition(1) - iPosition(1)
        VectorDistance(2) = jPosition(2) - iPosition(2)
        VectorDistance(3) = jPosition(3) - iPosition(3)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (squared)
        SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
        SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
        ! Preliminary test (circumscribing spheres)
        IF( SquaredDistance <= SquaredCutoffDistance ) THEN
          ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
          IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
            CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          ! Overlap test for spherocylinders (Vega-Lago method)
          ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, ParallelSPC, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          ! Overlap test for cylinders and/or spheres
          ELSE IF( GeometryType(3) ) THEN ! Cylinders
            IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Initialization
              OverlapSPC = .FALSE.
              ! Preliminary test (circumscribing spherocylinders)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
              &                     ContactDistance, ParallelSPC, OverlapSPC )
              ! Overlap criterion
              IF( OverlapSPC ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                jPosition(1) = iPosition(1) + VectorDistance(1)
                jPosition(2) = iPosition(2) + VectorDistance(2)
                jPosition(3) = iPosition(3) + VectorDistance(3)
                ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
                &                     iComponent, jComponent, ParallelSPC, PrivateOverlap )
                ! Overlap found
                IF( PrivateOverlap ) SharedOverlap = .TRUE.
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            ! Overlap test for cylinders and spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            ! Overlap test for spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              PrivateOverlap = .TRUE.
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            END IF
          END IF
        END IF
      END DO
    END DO
    !#############################################################################################!
    !$OMP END PARALLEL DO                                                                         !
    !#############################################################################################!
    ! Returns control to the calling program unit if an overlap is detected
    IF( SharedOverlap ) THEN
      Overlap = .TRUE.
      RETURN
    END IF
  END DO
END DO

RETURN ! No overlaps detected

END SUBROUTINE FullOverlapCheckFixPFraction

! *********************************************************************************************** !
!    This subroutine uses lists to check if a random volume scaling (isotropic or anisotropic)    !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullListOverlapCheck( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, Overlap, HalfNeighbours )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: MakeList
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )           :: iParticle      ! Counter (particle)
INTEGER( Kind= Int64 )           :: iComponent     ! Counter (component)
INTEGER( Kind= Int64 ), OPTIONAL :: HalfNeighbours ! Checks whether a cell and its 26 surrounding cells are searched (PRESENT) or just its 13 neighbour cells (NOT PRESENT)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: ContactDistance         ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition               ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation            ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion             ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap   ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
Overlap         = .FALSE.
SharedOverlap   = .FALSE.
PrivateOverlap  = .FALSE.
CellHalfLogical = .FALSE.

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphereDiameter / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphereDiameter / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphereDiameter / CurrentBoxLength(9)
CALL MakeList( BoxCutoff, pPositionMC, CurrentBoxLengthInverse )

! Check whether the number of cells along any direction (x, y, or z) is less than 3
IF( .NOT. CellListControl ) RETURN

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( iParticle, iComponent, iPosition, iOrientation, iQuaternion, ContactDistance ) &   !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                          !
!#################################################################################################!
! Loop over all particles (consider half of neighbours)
DO iParticle = 1, nParticles
  ! Cycle if a thread founds an overlap
  IF( SharedOverlap ) CYCLE
  ! Component index of particle i
  iComponent = pComponents(iParticle)
  ! Position of particle i
  iPosition = pPositionMC(:,iParticle)
  ! Quaternion of particle i
  iQuaternion = pQuaternionMC(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientationMC(:,iParticle)
  ! Overlap check between particle i and its neighbours
  CALL ListOverlapCheck( iComponent, iParticle, iQuaternion, iOrientation, iPosition, ContactDistance, CurrentBoxLength, &
  &                      CurrentBoxLengthInverse, PrivateOverlap, CellHalfLogical )
  ! Overlap found
  IF( PrivateOverlap ) SharedOverlap = .TRUE.
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN ! No overlaps detected

END SUBROUTINE FullListOverlapCheck

! *********************************************************************************************** !
!                    This subroutine checks there is any overlapping particles                    !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE OverlapCheckInitialConfiguration(  )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle                 ! Counters (particle)
INTEGER( Kind= Int64 ) :: iParticleMessage, jParticleMessage   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent               ! Counters (component)
INTEGER( Kind= Int64 ) :: iComponentMessage, jComponentMessage ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap        = .FALSE.
OverlapSPC     = .FALSE.
ParallelSPC    = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

! Anisomorphic and isomorphic components
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Stop if a thread founds an overlap (should never happen)
    IF( SharedOverlap ) STOP "ERROR: Mishap in the overlap check algorithm. Overlap detected but did not return control."
    ! Cycle if component index j < i (prevents double counting) 
    IF( jComponent < iComponent ) CYCLE
    !#############################################################################################!
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                           !
    !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &     !
    !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &           !
    !$OMP PRIVATE( SquaredDistance, SquaredCutoffDistance, ContactDistance ) &                    !
    !$OMP PRIVATE( OverlapSPC, ParallelSPC ) &                                                    !
    !$OMP REDUCTION( .OR. : PrivateOverlap )                                                      !
    !#############################################################################################!
    ! First loop represents all particles with indexes i of component Ci
    DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Second loop represents all particles with indexes j of component Cj
      DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
        ! Cycle if a thread founds an overlap
        IF( SharedOverlap ) CYCLE
        ! Cycle if particle index j <= i when components are isomorphic (prevents double counting) 
        IF( jComponent == iComponent .AND. jParticle <= iParticle ) CYCLE
        ! Position of particle i
        iPosition(1) = pPosition(1,iParticle)
        iPosition(2) = pPosition(2,iParticle)
        iPosition(3) = pPosition(3,iParticle)
        ! Position of particle j
        jPosition(1) = pPosition(1,jParticle)
        jPosition(2) = pPosition(2,jParticle)
        jPosition(3) = pPosition(3,jParticle)
        ! Orientation of particle i
        iOrientation(1) = pOrientation(1,iParticle)
        iOrientation(2) = pOrientation(2,iParticle)
        iOrientation(3) = pOrientation(3,iParticle)
        ! Orientation of particle j
        jOrientation(1) = pOrientation(1,jParticle)
        jOrientation(2) = pOrientation(2,jParticle)
        jOrientation(3) = pOrientation(3,jParticle)
        ! Quaternion of particle i
        iQuaternion(0) = pQuaternion(0,iParticle)
        iQuaternion(1) = pQuaternion(1,iParticle)
        iQuaternion(2) = pQuaternion(2,iParticle)
        iQuaternion(3) = pQuaternion(3,iParticle)
        ! Quaternion of particle j
        jQuaternion(0) = pQuaternion(0,jParticle)
        jQuaternion(1) = pQuaternion(1,jParticle)
        jQuaternion(2) = pQuaternion(2,jParticle)
        jQuaternion(3) = pQuaternion(3,jParticle)
        ! Vector distance between particles i and j
        VectorDistance(1) = jPosition(1) - iPosition(1)
        VectorDistance(2) = jPosition(2) - iPosition(2)
        VectorDistance(3) = jPosition(3) - iPosition(3)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT(ScalingDistanceUnitBox)
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (squared)
        SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
        SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
        ! Preliminary test (circumscribing spheres)
        IF( SquaredDistance <= SquaredCutoffDistance ) THEN
          ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
          IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
            CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) THEN
              SharedOverlap = .TRUE.
              iParticleMessage  = iParticle
              jParticleMessage  = jParticle
              iComponentMessage = iComponent
              jComponentMessage = jComponent
            END IF
          ! Overlap test for spherocylinders (Vega-Lago method)
          ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
            &                     ContactDistance, ParallelSPC, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) THEN
              SharedOverlap = .TRUE.
              iParticleMessage  = iParticle
              jParticleMessage  = jParticle
              iComponentMessage = iComponent
              jComponentMessage = jComponent
            END IF
          ! Overlap test for cylinders and/or spheres
          ELSE IF( GeometryType(3) ) THEN ! Cylinders
            IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Initialization
              OverlapSPC = .FALSE.
              ! Preliminary test (circumscribing spherocylinders)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
              &                     ContactDistance, ParallelSPC, OverlapSPC )
              ! Overlap criterion
              IF( OverlapSPC ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                jPosition(1) = iPosition(1) + VectorDistance(1)
                jPosition(2) = iPosition(2) + VectorDistance(2)
                jPosition(3) = iPosition(3) + VectorDistance(3)
                ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
                &                     iComponent, jComponent, ParallelSPC, PrivateOverlap )
                ! Overlap found
                IF( PrivateOverlap ) THEN
                  SharedOverlap = .TRUE.
                  iParticleMessage  = iParticle
                  jParticleMessage  = jParticle
                  iComponentMessage = iComponent
                  jComponentMessage = jComponent
                END IF
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) THEN
                SharedOverlap = .TRUE.
                iParticleMessage  = iParticle
                jParticleMessage  = jParticle
                iComponentMessage = iComponent
                jComponentMessage = jComponent
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              jPosition(1) = iPosition(1) + VectorDistance(1)
              jPosition(2) = iPosition(2) + VectorDistance(2)
              jPosition(3) = iPosition(3) + VectorDistance(3)
              CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) THEN
                SharedOverlap = .TRUE.
                iParticleMessage  = iParticle
                jParticleMessage  = jParticle
                iComponentMessage = iComponent
                jComponentMessage = jComponent
              END IF
            ! Overlap test for spheres
            ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
              PrivateOverlap = .TRUE.
              ! Overlap found
              IF( PrivateOverlap ) THEN
                SharedOverlap = .TRUE.
                iParticleMessage  = iParticle
                jParticleMessage  = jParticle
                iComponentMessage = iComponent
                jComponentMessage = jComponent
              END IF
            END IF
          END IF
        END IF
      END DO
    END DO
    !#############################################################################################!
    !$OMP END PARALLEL DO                                                                         !
    !#############################################################################################!
    ! Overlap detected in the initial configuration
    IF( SharedOverlap ) THEN
      IF( nComponents > 1 ) THEN
        WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", iParticleMessage, &
        &                   " of component ", iComponentMessage, " and ", jParticleMessage, " of component ", jComponentMessage, &
        &                   ". Exiting..."
      ELSE IF( nComponents == 1 ) THEN
        WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", iParticleMessage, " and ", &
        &                   jParticleMessage, ". Exiting..."
      END IF
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
  END DO
END DO

! No overlaps
RETURN

END SUBROUTINE OverlapCheckInitialConfiguration

! *********************************************************************************************** !
!         This subroutine uses lists to check whether there are any overlapping particles         !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE FullListOverlapCheckInitialConfiguration( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, HalfNeighbours )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: MakeList
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )           :: iParticle, jParticle                 ! Counter (particle)
INTEGER( Kind= Int64 )           :: iParticleMessage, jParticleMessage   ! Counters (particle)
INTEGER( Kind= Int64 )           :: iComponent, jComponent               ! Counter (component)
INTEGER( Kind= Int64 )           :: iComponentMessage, jComponentMessage ! Counters (component)
INTEGER( Kind= Int64 ), OPTIONAL :: HalfNeighbours                       ! Checks whether a cell and its 26 surrounding cells are searched (PRESENT) or just its 13 neighbour cells (NOT PRESENT)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: ContactDistance         ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition               ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation            ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion             ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap   ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
Overlap         = .FALSE.
SharedOverlap   = .FALSE.
PrivateOverlap  = .FALSE.
CellHalfLogical = .FALSE.

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphereDiameter / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphereDiameter / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphereDiameter / CurrentBoxLength(9)
CALL MakeList( BoxCutoff, pPosition, CurrentBoxLengthInverse )

! Check whether the number of cells along any direction (x, y, or z) is less than 3
IF( .NOT. CellListControl ) RETURN

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( iParticle, jParticle, iComponent, jComponent, iPosition, iQuaternion ) &           !
!$OMP PRIVATE( iOrientation, ContactDistance ) &                                                  !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                          !
!#################################################################################################!
! Loop over all particles (consider half of neighbours)
DO iParticle = 1, nParticles
  ! Cycle if a thread founds an overlap
  IF( SharedOverlap ) CYCLE
  ! Component index of particle i
  iComponent = pComponents(iParticle)
  ! Position of particle i
  iPosition = pPosition(:,iParticle)
  ! Quaternion of particle i
  iQuaternion = pQuaternion(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientation(:,iParticle)
  ! Overlap check between particle i and its neighbours
  CALL ListOverlapCheckInitialConfiguration( iComponent, iParticle, iQuaternion, iOrientation, iPosition, ContactDistance, &
  &                                          CurrentBoxLength, CurrentBoxLengthInverse, PrivateOverlap, jParticle, jComponent, &
  &                                          CellHalfLogical )
  ! Overlap found
  IF( PrivateOverlap ) THEN
    SharedOverlap = .TRUE.
    iParticleMessage  = iParticle
    jParticleMessage  = jParticle
    iComponentMessage = iComponent
    jComponentMessage = jComponent
  END IF
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Overlap detected in the initial configuration
IF( SharedOverlap ) THEN
  IF( nComponents > 1 ) THEN
    WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", iParticleMessage, &
    &                   " of component ", iComponentMessage, " and ", jParticleMessage, " of component ", jComponentMessage, &
    &                   ". Exiting..."
  ELSE IF( nComponents == 1 ) THEN
    WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", iParticleMessage, " and ", &
    &                   jParticleMessage, ". Exiting..."
  END IF
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! No overlaps
RETURN

END SUBROUTINE FullListOverlapCheckInitialConfiguration

! *********************************************************************************************** !
!         This subroutine uses lists to check whether there are any overlapping particles         !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE ListOverlapCheckInitialConfiguration( iComponent, iParticle, iQuaternion, iOrientation, iPosition, ContactDistance, &
&                                                CurrentBoxLength, CurrentBoxLengthInverse, Overlap, jParticle, jComponent, &
&                                                CellHalfLogical )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: Neighbours, CellIndex
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                          :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 )                          :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 )                          :: jList                  ! Counters (neighbour list)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndex             ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList         ! List of neighbour j particles

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                   :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: SquaredCutoffDistance      ! Cutoff distance (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapHER      ! Detects overlap between two ellipsoids of revolution: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC      ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL      ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC     ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
jList       = 0
Overlap     = .FALSE.
OverlapHER  = .FALSE.
OverlapCYL  = .FALSE.
ParallelSPC = .FALSE.

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndex = CellIndex( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogical .AND. ANY( iCellIndex(:) /= pCellIndex(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndex(:), "] and [", pCellIndex(:,iParticle), "]"
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! List of neighbours
jNeighbourList = Neighbours( iParticle, iCellIndex, CellHalfLogical )

! Loop until no more entries in the neighbour list
DO
  ! Next entry
  jList = jList + 1
  ! Neighbour particle index
  jParticle = jNeighbourList(jList)
  ! List exhausted
  IF( jParticle == 0 ) EXIT
  ! Skip self (should never happen)
  IF( iParticle == jParticle ) CYCLE
  ! Neighbour component index
  jComponent = pComponents(jParticle)
  ! Position of particle j
  jPosition(1) = pPosition(1,jParticle)
  jPosition(2) = pPosition(2,jParticle)
  jPosition(3) = pPosition(3,jParticle)
  ! Orientation of particle j
  jOrientation(1) = pOrientation(1,jParticle)
  jOrientation(2) = pOrientation(2,jParticle)
  jOrientation(3) = pOrientation(3,jParticle)
  ! Quaternion of particle j
  jQuaternion(0) = pQuaternion(0,jParticle)
  jQuaternion(1) = pQuaternion(1,jParticle)
  jQuaternion(2) = pQuaternion(2,jParticle)
  jQuaternion(3) = pQuaternion(3,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(1) = jPosition(1) - iPosition(1)
  VectorDistance(2) = jPosition(2) - iPosition(2)
  VectorDistance(3) = jPosition(3) - iPosition(3)
  ! Minimum image convention
  CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Cutoff distance (squared)
  SquaredCutoffDistance = 0.5D0 * ( cCircumscribingSphereDiameter(iComponent) + cCircumscribingSphereDiameter(jComponent) )
  SquaredCutoffDistance = SquaredCutoffDistance * SquaredCutoffDistance
  ! Preliminary test (circumscribing spheres)
  IF( SquaredDistance <= SquaredCutoffDistance ) THEN
    ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
    IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
      CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, ContactDistance, &
      &                     OverlapHER )
      ! Overlap criterion
      IF( OverlapHER ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    ! Overlap test for spherocylinders (Vega-Lago method)
    ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
      &                     ContactDistance, ParallelSPC, OverlapSPC )
      ! Overlap criterion
      IF( OverlapSPC ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    ! Overlap test for cylinders and/or spheres
    ELSE IF( GeometryType(3) ) THEN ! Cylinders
      IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Initialization
        OverlapSPC = .FALSE.
        ! Preliminary test (circumscribing spherocylinders)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, &
        &                     ContactDistance, ParallelSPC, OverlapSPC )
        ! Overlap criterion
        IF( OverlapSPC ) THEN
          ! Apply periodic boundary conditions on the position of particle j
          jPosition(1) = iPosition(1) + VectorDistance(1)
          jPosition(2) = iPosition(2) + VectorDistance(2)
          jPosition(3) = iPosition(3) + VectorDistance(3)
          ! Overlap test for cylinders (modified algorithm of Lopes et al.)
          CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
          &                     iComponent, jComponent, ParallelSPC, OverlapCYL )
          ! Overlap criterion
          IF( OverlapCYL ) THEN
            ! Returns control to the calling program unit if an overlap is detected
            Overlap = .TRUE.
            RETURN
          END IF
        END IF
      ! Overlap test for cylinders and spheres
      ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, OverlapCYL )
        IF( OverlapCYL ) THEN
          ! Returns control to the calling program unit if an overlap is detected
          Overlap = .TRUE.
          RETURN
        END IF
      ! Overlap test for cylinders and spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, OverlapCYL )
        IF( OverlapCYL ) THEN
          ! Returns control to the calling program unit if an overlap is detected
          Overlap = .TRUE.
          RETURN
        END IF
      ! Overlap test for spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        ! Returns control to the calling program unit if an overlap is detected
        Overlap = .TRUE.
        RETURN
      END IF
    END IF
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ListOverlapCheckInitialConfiguration

END MODULE