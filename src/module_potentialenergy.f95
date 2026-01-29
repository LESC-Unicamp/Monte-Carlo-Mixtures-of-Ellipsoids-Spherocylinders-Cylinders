! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                         to compute the potential energy of the system.                          !
!                                                                                                 !
! Version number: 2.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       January 28th, 2026                                        !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE ComputePotential

! Uses two modules: global variables and force fields
USE GlobalVar
USE ForceFields

IMPLICIT NONE

! Make subroutines public
PUBLIC :: ListComputeTotalPotentialEnergy, ListComputeParticlePotentialEnergy, ListComputeTotalPotentialEnergyInitialConfig, &
&         ListComputeParticlePotentialEnergyInitialConfig

CONTAINS

! *********************************************************************************************** !
!   This subroutine computes the total potential energy for the initial molecular configuration   !
! *********************************************************************************************** !
SUBROUTINE ComputeTotalPotentialEnergyInitialConfig(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: iParticle, jParticle   ! Counters
INTEGER( KIND= INT64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox     ! Position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( KIND= REAL64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy
REAL( KIND= REAL64 ), DIMENSION( nRange ) :: PairPotentialEnergyShared  ! Pair potential energy (shared)

! Initialization
TotalPotentialEnergy = 0.D0
PairPotentialEnergyShared = 0.D0

! Anisomorphic components
DO iComponent = 1, nComponents - 1
  DO jComponent = iComponent + 1, nComponents
    !#############################################################################################!
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                           !
    !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &     !
    !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &           !
    !$OMP PRIVATE( SquaredDistance, PairPotentialEnergy ) &                                       !
    !$OMP REDUCTION( + : PairPotentialEnergyShared )                                              !
    !#############################################################################################!
    ! First loop represents all particles with indexes i of component Ci
    DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Second loop represents all particles with indexes j of component Cj
      DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
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
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Compute pair potential
        IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
          CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
        ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
          CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
          &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
        END IF
        ! Increment total potential energy
        PairPotentialEnergyShared = PairPotentialEnergyShared + PairPotentialEnergy
      END DO
    END DO
    !#############################################################################################!
    !$OMP END PARALLEL DO                                                                         !
    !#############################################################################################!
    ! Update total potential energy
    TotalPotentialEnergy = TotalPotentialEnergy + PairPotentialEnergyShared
    PairPotentialEnergyShared = 0.D0
  END DO
END DO

! Isomorphic components
DO iComponent = 1, nComponents
  jComponent = iComponent
  !###############################################################################################!
  !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                             !
  !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &       !
  !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &             !
  !$OMP PRIVATE( SquaredDistance, PairPotentialEnergy ) &                                         !
  !$OMP REDUCTION( + : PairPotentialEnergyShared )                                                !
  !###############################################################################################!
  ! First loop represents a particle with an index i of component Ci
  DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
    ! Second loop represents a particle with an index j of component Cj = Ci
    DO jParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Cycle if index j <= i (prevents double counting)
      IF( jParticle <= iParticle ) CYCLE
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
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Compute pair potential
      IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
        CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
      ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
        CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
        &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
      END IF
      ! Increment total potential energy
      PairPotentialEnergyShared = PairPotentialEnergyShared + PairPotentialEnergy
    END DO
  END DO
  !###############################################################################################!
  !$OMP END PARALLEL DO                                                                           !
  !###############################################################################################!
  ! Update total potential energy
  TotalPotentialEnergy = TotalPotentialEnergy + PairPotentialEnergyShared
  PairPotentialEnergyShared = 0.D0
END DO

RETURN

END SUBROUTINE ComputeTotalPotentialEnergyInitialConfig

! *********************************************************************************************** !
!            This subroutine uses lists to compute the total potential energy for the             !
!                                 initial molecular configuration                                 !
! *********************************************************************************************** !
SUBROUTINE ListComputeTotalPotentialEnergyInitialConfig( CurrentBoxLength, CurrentBoxLengthInverse, HalfNeighbours )

! Uses one module: global variables
USE LinkedLists, ONLY: MakeListPotential

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
REAL( Kind= Real64 ), DIMENSION( 3 )      :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition               ! Position of particle i
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation            ! Orientation of particle i
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion             ! Quaternion of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy        ! Potential energy of particle i
REAL( KIND= REAL64 ), DIMENSION( nRange ) :: iPotentialEnergyShared  ! Potential energy of particle i (shared)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphericalWell / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphericalWell / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphericalWell / CurrentBoxLength(9)
CALL MakeListPotential( BoxCutoff, pPosition, CurrentBoxLengthInverse )

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

! Initialization
TotalPotentialEnergy = 0.D0
iPotentialEnergyShared = 0.D0

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( iParticle, iComponent, iPosition, iOrientation, iQuaternion, iPotentialEnergy ) &  !
!$OMP REDUCTION( + : iPotentialEnergyShared )                                                     !
!#################################################################################################!
! Loop over all particles (consider neighbours)
DO iParticle = 1, nParticles
  ! Component index of particle i
  iComponent = pComponents(iParticle)
  ! Position of particle i
  iPosition = pPosition(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientation(:,iParticle)
  ! Quaternion of particle i
  iQuaternion = pQuaternion(:,iParticle)
  ! Potential between particle i and its neighbours
  CALL ListComputeParticlePotentialEnergyInitialConfig( iComponent, iParticle, iQuaternion, iOrientation, iPosition, &
  &                                                     iPotentialEnergy, CurrentBoxLength, CurrentBoxLengthInverse, &
  &                                                     CellHalfLogical )
  ! Increment total potential energy
  iPotentialEnergyShared = iPotentialEnergyShared + iPotentialEnergy
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Total potential energy
TotalPotentialEnergy = iPotentialEnergyShared
iPotentialEnergyShared = 0.D0

RETURN

END SUBROUTINE ListComputeTotalPotentialEnergyInitialConfig

! *********************************************************************************************** !
!        This subroutine uses lists to compute the potential energy of a random particle i        !
! *********************************************************************************************** !
SUBROUTINE ListComputeParticlePotentialEnergyInitialConfig( iComponent, iParticle, iQuaternion, iOrientation, iPosition, &
&                                                           iPotentialEnergy, CurrentBoxLength, CurrentBoxLengthInverse, &
&                                                           CellHalfLogicalPotential )

! Uses one module: linked lists
USE LinkedLists, ONLY: NeighboursPotential, CellIndexPotential

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                          :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 )                          :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 )                          :: jList                  ! Counters (neighbour list)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndexPotential    ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList         ! List of neighbour j particles

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse    ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy           ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: CellHalfLogicalPotential ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
iPotentialEnergy = 0.D0

! Initialization
jList = 0

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndexPotential = CellIndexPotential( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogicalPotential .AND. ANY( iCellIndexPotential(:) /= pCellIndexPotential(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndexPotential(:), "] and [", &
  &                                                pCellIndexPotential(:,iParticle), "]"
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! List of neighbours
jNeighbourList = NeighboursPotential( iParticle, iCellIndexPotential, CellHalfLogicalPotential )

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
  ! Compute pair potential
  IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
    CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
    &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO

RETURN

END SUBROUTINE ListComputeParticlePotentialEnergyInitialConfig

! *********************************************************************************************** !
!           This subroutine computes the total potential energy for the NPT simulation            !
! *********************************************************************************************** !
SUBROUTINE ComputeTotalPotentialEnergy( CurrentBoxLength, CurrentBoxLengthInverse, SystemPotentialEnergy )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: iParticle, jParticle   ! Counters
INTEGER( KIND= INT64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                      :: SystemPotentialEnergy      ! System potential energy
REAL( KIND= REAL64 )                      :: PairPotentialEnergyShared  ! Pair potential energy (shared)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox     ! Position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( KIND= REAL64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy

! Initialization
SystemPotentialEnergy = 0.D0
PairPotentialEnergyShared = 0.D0

! Anisomorphic components
DO iComponent = 1, nComponents - 1
  DO jComponent = iComponent + 1, nComponents
    !#############################################################################################!
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                           !
    !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &     !
    !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &           !
    !$OMP PRIVATE( SquaredDistance, PairPotentialEnergy ) &                                       !
    !$OMP REDUCTION( + : PairPotentialEnergyShared )                                              !
    !#############################################################################################!
    ! First loop represents all particles with indexes i of component Ci
    DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Second loop represents all particles with indexes j of component Cj
      DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
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
        ! Compute pair potential
        IF( FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
          CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
        ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
          CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
          &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
        END IF
        ! Increment total potential energy
        PairPotentialEnergyShared = PairPotentialEnergyShared + PairPotentialEnergy(1)
      END DO
    END DO
    !#############################################################################################!
    !$OMP END PARALLEL DO                                                                         !
    !#############################################################################################!
    ! Update total potential energy
    SystemPotentialEnergy = SystemPotentialEnergy + PairPotentialEnergyShared
    PairPotentialEnergyShared = 0.D0
  END DO
END DO

! Isomorphic components
DO iComponent = 1, nComponents
  jComponent = iComponent
  !###############################################################################################!
  !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                             !
  !$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &       !
  !$OMP PRIVATE( iQuaternion, jQuaternion, VectorDistance, ScalingDistanceUnitBox ) &             !
  !$OMP PRIVATE( SquaredDistance, PairPotentialEnergy ) &                                         !
  !$OMP REDUCTION( + : PairPotentialEnergyShared )                                                !
  !###############################################################################################!
  ! First loop represents a particle with an index i of component Ci
  DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
    ! Second loop represents a particle with an index j of component Cj = Ci
    DO jParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) )
      ! Cycle if index j <= i (prevents double counting)
      IF( jParticle <= iParticle ) CYCLE
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
      ! Compute pair potential
      IF( FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
        CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
      ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
        CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
        &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
      END IF
      ! Increment total potential energy
      PairPotentialEnergyShared = PairPotentialEnergyShared + PairPotentialEnergy(1)
    END DO
  END DO
  !###############################################################################################!
  !$OMP END PARALLEL DO                                                                           !
  !###############################################################################################!
  ! Update total potential energy
  SystemPotentialEnergy = SystemPotentialEnergy + PairPotentialEnergyShared
  PairPotentialEnergyShared = 0.D0
END DO

RETURN

END SUBROUTINE ComputeTotalPotentialEnergy

! *********************************************************************************************** !
!              This subroutine computes the potential energy of a random particle i               !
! *********************************************************************************************** !
SUBROUTINE ComputeParticlePotentialEnergy( iComponent, iParticle, iQuaternion, iOrientation, iPosition, iPotentialEnergy, &
&                                          CurrentBoxLength, CurrentBoxLengthInverse )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse    ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy           ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergyShared     ! Potential energy of particle i (shared)
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy

! Initialization
iPotentialEnergy = 0.D0
iPotentialEnergyShared = 0.D0

! Anisomorphic components I
DO jComponent = 1, iComponent - 1
  !###############################################################################################!
  !$OMP PARALLEL DO DEFAULT( Shared ) &                                                           !
  !$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, VectorDistance ) &              !
  !$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, PairPotentialEnergy ) &                 !
  !$OMP REDUCTION( + : iPotentialEnergyShared )                                                   !
  !###############################################################################################!
  ! Unique loop takes only particles whose component indexes are less than Ci
  DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
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
    ! Compute pair potential
    IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
      CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
    ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
      CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
      &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
    END IF
    ! Increment potential energy of particle i
    iPotentialEnergyShared = iPotentialEnergyShared + PairPotentialEnergy
  END DO
  !###############################################################################################!
  !$OMP END PARALLEL DO                                                                           !
  !###############################################################################################!
  ! Update potential energy of particle i
  iPotentialEnergy = iPotentialEnergy + iPotentialEnergyShared
  iPotentialEnergyShared = 0.D0
END DO

! Anisomorphic components II
DO jComponent = iComponent + 1, nComponents
  !###############################################################################################!
  !$OMP PARALLEL DO DEFAULT( Shared ) &                                                           !
  !$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, VectorDistance ) &              !
  !$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, PairPotentialEnergy ) &                 !
  !$OMP REDUCTION( + : iPotentialEnergyShared )                                                   !
  !###############################################################################################!
  ! Unique loop takes only particles whose component indexes are greater than Ci
  DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
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
    ! Compute pair potential
    IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
      CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
    ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
      CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
      &                                    SquaredDistance, iComponent, jComponent, &
      &                                    PairPotentialEnergy )
    END IF
    ! Increment potential energy of particle i
    iPotentialEnergyShared = iPotentialEnergyShared + PairPotentialEnergy
  END DO
  !###############################################################################################!
  !$OMP END PARALLEL DO                                                                           !
  !###############################################################################################!
  ! Update potential energy of particle i
  iPotentialEnergy = iPotentialEnergy + iPotentialEnergyShared
  iPotentialEnergyShared = 0.D0
END DO

! Isomorphic components
jComponent = iComponent

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, VectorDistance ) &                !
!$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, PairPotentialEnergy ) &                   !
!$OMP REDUCTION( + : iPotentialEnergyShared )                                                     !
!#################################################################################################!
! First loop takes only particles whose j-indexes are below the i-index of the particles of the component Ci
DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, iParticle - 1
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
  ! Compute pair potential
  IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
    CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
    &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy of particle i
  iPotentialEnergyShared = iPotentialEnergyShared + PairPotentialEnergy
END DO
!#################################################################################################!
!$OMP END PARALLEL DO
!#################################################################################################!

! Update potential energy of particle i
iPotentialEnergy = iPotentialEnergy + iPotentialEnergyShared
iPotentialEnergyShared = 0.D0

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, VectorDistance ) &                !
!$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, PairPotentialEnergy ) &                   !
!$OMP REDUCTION( + : iPotentialEnergyShared )                                                     !
!#################################################################################################!
! Second loop takes only particles whose j-indexes are above the i-index of the particles of the component Ci
DO jParticle = iParticle + 1, SUM( cParticles(0:jComponent) )
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
  ! Compute pair potential
  IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
    CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
    &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy of particle i
  iPotentialEnergyShared = iPotentialEnergyShared + PairPotentialEnergy
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Update potential energy of particle i
iPotentialEnergy = iPotentialEnergy + iPotentialEnergyShared
iPotentialEnergyShared = 0.D0

RETURN

END SUBROUTINE ComputeParticlePotentialEnergy

! *********************************************************************************************** !
!            This subroutine uses lists to compute the total potential energy for the             !
!                                         NPT simulation                                          !
! *********************************************************************************************** !
SUBROUTINE ListComputeTotalPotentialEnergy( CurrentBoxLength, CurrentBoxLengthInverse, SystemPotentialEnergy, HalfNeighbours )

! Uses one module: global variables
USE LinkedLists, ONLY: MakeListPotential

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
REAL( KIND= REAL64 )                      :: SystemPotentialEnergy   ! System potential energy
REAL( Kind= Real64 )                      :: iPotentialEnergyShared  ! Potential energy of particle i (shared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition               ! Position of particle i
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation            ! Orientation of particle i
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion             ! Quaternion of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy        ! Potential energy of particle i

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphericalWell / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphericalWell / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphericalWell / CurrentBoxLength(9)
CALL MakeListPotential( BoxCutoff, pPositionMC, CurrentBoxLengthInverse )

! Check whether the number of cells along any direction (x, y, or z) is less than 3
IF( .NOT. CellListControlPotential ) RETURN

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

! Initialization
SystemPotentialEnergy = 0.D0
iPotentialEnergyShared = 0.D0

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( iParticle, iComponent, iPosition, iOrientation, iQuaternion, iPotentialEnergy ) &  !
!$OMP REDUCTION( + : iPotentialEnergyShared )                                                     !
!#################################################################################################!
! Loop over all particles (consider neighbours)
DO iParticle = 1, nParticles
  ! Component index of particle i
  iComponent = pComponents(iParticle)
  ! Position of particle i
  iPosition = pPositionMC(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientationMC(:,iParticle)
  ! Quaternion of particle i
  iQuaternion = pQuaternionMC(:,iParticle)
  ! Potential between particle i and its neighbours
  CALL ListComputeParticlePotentialEnergy( iComponent, iParticle, iQuaternion, iOrientation, iPosition, iPotentialEnergy, &
  &                                        CurrentBoxLength, CurrentBoxLengthInverse, CellHalfLogical )
  ! Increment total potential energy
  iPotentialEnergyShared = iPotentialEnergyShared + iPotentialEnergy(1)
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Update total potential energy
SystemPotentialEnergy = iPotentialEnergyShared
iPotentialEnergyShared = 0.D0

RETURN

END SUBROUTINE ListComputeTotalPotentialEnergy

! *********************************************************************************************** !
!        This subroutine uses lists to compute the potential energy of a random particle i        !
! *********************************************************************************************** !
SUBROUTINE ListComputeParticlePotentialEnergy( iComponent, iParticle, iQuaternion, iOrientation, iPosition, iPotentialEnergy, &
&                                              CurrentBoxLength, CurrentBoxLengthInverse, CellHalfLogicalPotential )

! Uses one module: linked lists
USE LinkedLists, ONLY: NeighboursPotential, CellIndexPotential

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                          :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 )                          :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 )                          :: jList                  ! Counters (neighbour list)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndexPotential    ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList         ! List of neighbour j particles

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse    ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy           ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: CellHalfLogicalPotential ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
iPotentialEnergy = 0.D0

! Initialization
jList = 0

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndexPotential = CellIndexPotential( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogicalPotential .AND. ANY( iCellIndexPotential(:) /= pCellIndexPotential(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndexPotential(:), "] and [", &
  &                                                pCellIndexPotential(:,iParticle), "]"
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! List of neighbours
jNeighbourList = NeighboursPotential( iParticle, iCellIndexPotential, CellHalfLogicalPotential )

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
  ! Compute pair potential
  IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
    CALL AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
    &                                    SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO

RETURN

END SUBROUTINE ListComputeParticlePotentialEnergy

END MODULE ComputePotential