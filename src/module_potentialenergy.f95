! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                         to compute the potential energy of the system.                          !
!                                                                                                 !
! Version number: 1.2.1                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       November 22nd, 2023                                       !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE ComputePotential

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! Make subroutines public
PUBLIC :: ListComputeTotalPotentialEnergyInitialConfiguration, ListComputeParticlePotentialEnergyInitialConfiguration

CONTAINS

! *********************************************************************************************** !
!   This subroutine computes the total potential energy for the initial molecular configuration   !
! *********************************************************************************************** !
SUBROUTINE ComputeTotalPotentialEnergy(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: iParticle, jParticle   ! Counters
INTEGER( KIND= INT64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                       :: SquaredDistance        ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 ), DIMENSION( 3 )       :: iPosition, jPosition   ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )       :: VectorDistance         ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )       :: ScalingDistanceUnitBox ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( nRange ) :: PairPotentialEnergy     ! Pair potential energy

! Initialization
TotalPotentialEnergy = 0.D0

! Anisomorphic components
DO iComponent = 1, nComponents - 1
  DO jComponent = iComponent + 1, nComponents
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
        IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
          CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
        END IF
        ! Increment total potential energy
        TotalPotentialEnergy = TotalPotentialEnergy + PairPotentialEnergy
      END DO
    END DO
  END DO
END DO
! Isomorphic components
DO iComponent = 1, nComponents
  jComponent = iComponent
  ! First loop represents a particle with an index i of component Ci
  DO iParticle = SUM( cParticles(0:(iComponent-1)) ) + 1, SUM( cParticles(0:iComponent) ) - 1
    ! Second loop represents all other particles with indexes j > i of component Cj = Ci
    DO jParticle = iParticle + 1, SUM( cParticles(0:jComponent) )
      ! Position of particle i
      iPosition(1) = pPosition(1,iParticle)
      iPosition(2) = pPosition(2,iParticle)
      iPosition(3) = pPosition(3,iParticle)
      ! Position of particle j
      jPosition(1) = pPosition(1,jParticle)
      jPosition(2) = pPosition(2,jParticle)
      jPosition(3) = pPosition(3,jParticle)
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
      IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
        CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
      END IF
      ! Increment total potential energy
      TotalPotentialEnergy = TotalPotentialEnergy + PairPotentialEnergy
    END DO
  END DO
END DO

RETURN

END SUBROUTINE ComputeTotalPotentialEnergy

! *********************************************************************************************** !
!            This subroutine uses lists to compute the total potential energy for the             !
!                                 initial molecular configuration                                 !
! *********************************************************************************************** !
SUBROUTINE ListComputeTotalPotentialEnergyInitialConfiguration( CurrentBoxLength, CurrentBoxLengthInverse, HalfNeighbours )

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
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition               ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy        ! Potential energy of particle i

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

! Loop over all particles (consider neighbours)
DO iParticle = 1, nParticles
  ! Component index of particle i
  iComponent = pComponents(iParticle)
  ! Position of particle i
  iPosition = pPosition(:,iParticle)
  ! Potential between particle i and its neighbours
  CALL ListComputeParticlePotentialEnergyInitialConfiguration( iComponent, iParticle, iPosition, iPotentialEnergy, &
  &                                                            CurrentBoxLength, CurrentBoxLengthInverse, CellHalfLogical )
  ! Increment total potential energy
  TotalPotentialEnergy = TotalPotentialEnergy + iPotentialEnergy
END DO

RETURN

END SUBROUTINE ListComputeTotalPotentialEnergyInitialConfiguration

! *********************************************************************************************** !
!              This subroutine computes the potential energy of a random particle i               !
! *********************************************************************************************** !
SUBROUTINE ComputeParticlePotentialEnergy( iComponent, iParticle, iPosition, iPotentialEnergy, CurrentBoxLength, &
&                                          CurrentBoxLengthInverse )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance         ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition    ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance          ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox  ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy        ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy     ! Pair potential energy

! Initialization
iPotentialEnergy = 0.D0

! Anisomorphic components I
DO jComponent = 1, iComponent - 1
  ! Unique loop takes only particles whose component indexes are less than Ci
  DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
    ! Position of particle j
    jPosition(1) = pPositionMC(1,jParticle)
    jPosition(2) = pPositionMC(2,jParticle)
    jPosition(3) = pPositionMC(3,jParticle)
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
    IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
      CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
    END IF
    ! Increment potential energy of particle i
    iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
  END DO
END DO

! Anisomorphic components II
DO jComponent = iComponent + 1, nComponents
  ! Unique loop takes only particles whose component indexes are greater than Ci
  DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, SUM( cParticles(0:jComponent) )
    ! Position of particle j
    jPosition(1) = pPositionMC(1,jParticle)
    jPosition(2) = pPositionMC(2,jParticle)
    jPosition(3) = pPositionMC(3,jParticle)
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
    IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
      CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
    END IF
    ! Increment potential energy of particle i
    iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
  END DO
END DO

! Isomorphic components
jComponent = iComponent
! First loop takes only particles whose j-indexes are below the i-index of the particles of the component Ci
DO jParticle = SUM( cParticles(0:(jComponent-1)) ) + 1, iParticle - 1
  ! Position of particle j
  jPosition(1) = pPositionMC(1,jParticle)
  jPosition(2) = pPositionMC(2,jParticle)
  jPosition(3) = pPositionMC(3,jParticle)
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
  IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy of particle i
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO
! Second loop takes only particles whose j-indexes are above the i-index of the particles of the component Ci
DO jParticle = iParticle + 1, SUM( cParticles(0:jComponent) )
  ! Position of particle j
  jPosition(1) = pPositionMC(1,jParticle)
  jPosition(2) = pPositionMC(2,jParticle)
  jPosition(3) = pPositionMC(3,jParticle)
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
  IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy of particle i
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO

RETURN

END SUBROUTINE ComputeParticlePotentialEnergy

! *********************************************************************************************** !
!        This subroutine uses lists to compute the potential energy of a random particle i        !
! *********************************************************************************************** !
SUBROUTINE ListComputeParticlePotentialEnergyInitialConfiguration( iComponent, iParticle, iPosition, iPotentialEnergy, &
&                                                                  CurrentBoxLength, CurrentBoxLengthInverse, &
&                                                                  CellHalfLogicalPotential )

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
REAL( Kind= Real64 )                      :: SquaredDistance         ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition    ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance          ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox  ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy         ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy      ! Pair potential energy

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
  IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO

RETURN

END SUBROUTINE ListComputeParticlePotentialEnergyInitialConfiguration

! *********************************************************************************************** !
!        This subroutine uses lists to compute the potential energy of a random particle i        !
! *********************************************************************************************** !
SUBROUTINE ListComputeParticlePotentialEnergy( iComponent, iParticle, iPosition, iPotentialEnergy, CurrentBoxLength, &
&                                              CurrentBoxLengthInverse, CellHalfLogicalPotential )

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
REAL( Kind= Real64 )                      :: SquaredDistance         ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition    ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance          ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: ScalingDistanceUnitBox  ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )      :: CurrentBoxLengthInverse ! Inverse of box length
REAL( Kind= Real64 ), DIMENSION( nRange ) :: iPotentialEnergy        ! Potential energy of particle i
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy     ! Pair potential energy

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
  IF( PotentialTypeLogical(2) ) THEN ! Spherical square-well potential
    CALL SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )
  END IF
  ! Increment potential energy
  iPotentialEnergy = iPotentialEnergy + PairPotentialEnergy
END DO

RETURN

END SUBROUTINE ListComputeParticlePotentialEnergy

! *********************************************************************************************** !
!            This subroutine computes the pair potential between particles i and j by             !
!             applying a discrete square-well potential to compute the pair potential             !
! *********************************************************************************************** !
SUBROUTINE SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 ) :: rRange                 ! Counter (range)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance       ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                      :: SquaredPotentialRange ! Effective range of attraction (squared)
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy   ! Pair potential energy

! Compute pair potential for every attractive range
DO rRange = 1, nRange
  ! Effective range of attraction (squared)
  SquaredPotentialRange = 0.5D0 * ( cPotentialRange(iComponent,rRange) + cPotentialRange(jComponent,rRange) )
  SquaredPotentialRange = SquaredPotentialRange * SquaredPotentialRange
  ! Overlap in the region of attraction
  IF( SquaredDistance <= SquaredPotentialRange ) THEN
    ! Pair potential (reduced units)
    PairPotentialEnergy(rRange) = - 1.D0
  ELSE
    PairPotentialEnergy(rRange) = 0.D0
  END IF
END DO

RETURN

END SUBROUTINE SquareWellPotential

END MODULE ComputePotential