! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!   This module create and manipulates linked cell lists of particles according to the textbook   !
!                                   Allen and Tildesley (2017).                                   !
!            Check https://github.com/Allen-Tildesley/examples/ for more information.             !
!                                                                                                 !
! Version number: 1.2.1                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       October 31st, 2023                                        !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE LinkedLists

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! Make subroutine public
PUBLIC :: MakeList

CONTAINS

! *********************************************************************************************** !
!                         This subroutine allocates the cell list arrays                          !
! *********************************************************************************************** !
SUBROUTINE InitializeList(  )

IMPLICIT NONE

IF( .NOT. ALLOCATED( pCellList ) ) ALLOCATE( pCellList(nParticles) )
IF( .NOT. ALLOCATED( pCellIndex ) ) ALLOCATE( pCellIndex(3,nParticles) )
IF( .NOT. ALLOCATED( pCellHead ) ) THEN
  ALLOCATE( pCellHead(0:(pCells(1)-1),0:(pCells(2)-1),0:(pCells(3)-1)) )
  ! Initialize head of the list (0 means the list is exhausted)
  pCellHead = 0
END IF

RETURN

END SUBROUTINE InitializeList

! *********************************************************************************************** !
!                        This subroutine deallocates the cell list arrays                         !
! *********************************************************************************************** !
SUBROUTINE FinalizeList(  )

IMPLICIT NONE

IF( ALLOCATED( pCellList ) ) DEALLOCATE( pCellList )
IF( ALLOCATED( pCellIndex ) ) DEALLOCATE( pCellIndex )
IF( ALLOCATED( pCellHead ) ) DEALLOCATE( pCellHead )

RETURN

END SUBROUTINE FinalizeList

! *********************************************************************************************** !
!                               This subroutine creates a cell list                               !
! *********************************************************************************************** !
SUBROUTINE MakeList( BoxCutoff, Position, bLengthInverse, FixPFraction )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )           :: pParticle    ! Counter (particles)
INTEGER( Kind= Int64 ), OPTIONAL :: FixPFraction ! Optional control variable for the algorithm that fixes the packing fraction

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 )             :: ScalingDistanceUnitBox ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )             :: BoxCutoff              ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 9 )             :: bLengthInverse         ! Length of the box (inverse)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles ) :: Position               ! Position of particles

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FixPFractionLogical ! Control variable for the algorithm that fixes the packing fraction

! Initialization
CellListControl = .TRUE.

! Number of cells in each dimension
pCells = FLOOR( 1.D0 / BoxCutoff )

! Disable cell list
IF( ANY( pCells < 3 ) ) THEN
  CALL FinalizeList(  )
  CellListControl = .FALSE.
  RETURN
END IF

! Initialize array
CALL InitializeList(  )

! Optional control
FixPFractionLogical = PRESENT( FixPFraction )

! Check whether the new number of cells is different from the previous number of cells (NPT simulation)
IF( pCells(1) /= SIZE( pCellHead, Dim= 1 ) .OR. pCells(2) /= SIZE( pCellHead, Dim= 2 ) .OR. &
&   pCells(3) /= SIZE( pCellHead, Dim= 3 ) .OR. FixPFractionLogical ) THEN
  ! Recreate list
  DEALLOCATE( pCellHead )
  ALLOCATE( pCellHead(0:(pCells(1)-1),0:(pCells(2)-1),0:(pCells(3)-1)) )
  ! Initialize head of the list (0 means the list is exhausted)
  pCellHead = 0
END IF

! No need to recreate lists
IF( ANY( pCellHead /= 0 ) ) RETURN

! Allocate particles to cells
DO pParticle = 1, nParticles
  ! Spatial transformation
  CALL MatrixVectorMultiplication( bLengthInverse, Position(:,pParticle), ScalingDistanceUnitBox )
  ! Minimum image convention
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  ! Index function allocating particle i to cell
  pCellIndex(:,pParticle) = CellIndex( ScalingDistanceUnitBox )
  ! Add particle i to list
  CALL CreateParticleList( pParticle, pCellIndex(:,pParticle) )
END DO

RETURN

END SUBROUTINE MakeList

! *********************************************************************************************** !
!               This function creates an index that allocates a particle to a cell                !
! *********************************************************************************************** !
FUNCTION CellIndex( iPosition ) RESULT( iCellIndex )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: iCellIndex ! 3D cell index

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iPosition ! Scaled position of particles

! Debug condition
IF( ANY( DABS( iPosition ) > 0.5D0 ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0.10,', '),G0.10,G0)" ) "Particle not in main box: [", iPosition, "]. Exiting..."
  CALL Sleep( 1 )
  CALL Exit(  )
END IF

! Cell index formula
iCellIndex(:) = FLOOR( ( iPosition(:) + 0.5D0 ) * DBLE( pCells(:) ) )

! Guard against small chance of roundoff error
WHERE( iCellIndex(:) < 0 ) iCellIndex(:) = 0
WHERE( iCellIndex(:) > ( pCells(:) - 1 ) ) iCellIndex(:) = pCells(:) - 1

RETURN

END FUNCTION CellIndex

! *********************************************************************************************** !
!                         This subroutine creates a particle i in a cell                          !
! *********************************************************************************************** !
SUBROUTINE CreateParticleList( pParticle, iCellIndex )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: pParticle  ! Counter (particles)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: iCellIndex ! 3D index of cell in which particle i lies

! Transfer old head to list (the 'head' is always the largest index in the list)
pCellList(pParticle) = pCellHead(iCellIndex(1),iCellIndex(2),iCellIndex(3))

! Make particle i the new head for this list
pCellHead(iCellIndex(1),iCellIndex(2),iCellIndex(3)) = pParticle 

! Store 3D index in array (necessary when moving particles among cells)
pCellIndex(:,pParticle) = iCellIndex(:)

RETURN

END SUBROUTINE CreateParticleList

! *********************************************************************************************** !
!                         This subroutine destroys a particle i in a cell                         !
! *********************************************************************************************** !
SUBROUTINE DestroyParticleList( pParticle, iCellIndex )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: pParticle       ! Counter (particles)
INTEGER( Kind= Int64 )                 :: pCurrent, pNext ! Counter (particles)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: iCellIndex      ! 3D index of cell in which particle i lies

! Locate head of list corresponding to cell
pCurrent = pCellHead(iCellIndex(1),iCellIndex(2),iCellIndex(3))

! Particle i is the head atom in that cell
IF( pCurrent == pParticle ) THEN
  ! Simply point head at next atom, we're done (previous head of list)
  pCellHead(iCellIndex(1),iCellIndex(2),iCellIndex(3)) = pCellList(pParticle) ! Simply point head at next atom, we're done (old head da lista)
! Particle i lies further down in the list
ELSE
  ! Loop traversing a link-list
  FindParticleInList: DO
    ! Look ahead to the next entry
    pNext = pCellList(pCurrent)
    ! Particle has been found
    IF ( pNext == pParticle ) THEN
      ! Link over particle
      pCellList(pCurrent) = pCellList(pParticle)
      EXIT FindParticleInList
    ! Debug check
    ELSE IF ( pNext == 0 ) THEN ! This should never happen
      WRITE( *, "(3G0,2(G0,', '),2G0)" ) "Could not find particle #", pParticle, " in its cell: [", iCellIndex, "]"
      CALL Sleep( 1 )
      CALL Exit(  )
    ! Move on to the next
    ELSE
      ! Keep index for next iteration
      pCurrent = pNext
    END IF
  END DO FindParticleInList
END IF

RETURN

END SUBROUTINE DestroyParticleList

! *********************************************************************************************** !
!            This subroutine moves a particle i from the current cell to another cell             !
! *********************************************************************************************** !
SUBROUTINE MoveParticleList( pParticle, iCellIndex )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: pParticle  ! Counter (particles)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: iCellIndex ! 3D index of cell in which particle i lies

! This is only necessary when a particle is transferred, leaving its cell towards another cell
IF( ALL( iCellIndex(:) == pCellIndex(:,pParticle) ) ) RETURN

! Remove particle i from old cell
CALL DestroyParticleList( pParticle, pCellIndex(:,pParticle) )

! Add particle i to a new cell
CALL CreateParticleList( pParticle, iCellIndex(:) )

! *********************************************************************************************** !
! OBSERVATIONS                                                                                    !
! *********************************************************************************************** !
!  Rotations do not move particles between cells. Volume changes also do not move particles but   !
!  can change the number of cells. However, in this case, another list is created.                !
! *********************************************************************************************** !

RETURN

END SUBROUTINE MoveParticleList

! *********************************************************************************************** !
!              This subroutine checks whether translating a particle i will move it               !
!                              from its current cell to another cell                              !
! *********************************************************************************************** !
SUBROUTINE ParticleTranslationNVT( pParticle, ScalingDistanceUnitBox )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: pParticle  ! Counter (particles)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: iCellIndex ! 3D index of cell in which particle i lies

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: ScalingDistanceUnitBox

! New cell index
iCellIndex(:) = CellIndex( ScalingDistanceUnitBox )

! Move particle from a cell to another cell, only if necessary
CALL MoveParticleList( pParticle, iCellIndex(:) )

RETURN

END SUBROUTINE ParticleTranslationNVT

! *********************************************************************************************** !
!            This subroutine checks whether changing the volume of the simulation box             !
!                   changes the number of cells in the x-, y-, or z-directions                    !
! *********************************************************************************************** !
SUBROUTINE BoxCheckNPT( Position, bLengthOld, bLengthNew, bLengthInverseNew )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: pCellsOld ! Number of cells (old)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: pCellsNew ! Number of cells (new)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 )             :: BoxCutoffOld      ! Old box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )             :: BoxCutoffNew      ! New box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 9 )             :: bLengthNew        ! New length of the box
REAL( Kind= Real64 ), DIMENSION( 9 )             :: bLengthOld        ! Old length of the box
REAL( Kind= Real64 ), DIMENSION( 9 )             :: bLengthInverseNew ! New length of the box (inverse)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles ) :: Position          ! Position of particles

! Box cutoff
BoxCutoffOld(1) = cLargestSphereDiameter / bLengthOld(1)
BoxCutoffOld(2) = cLargestSphereDiameter / bLengthOld(5)
BoxCutoffOld(3) = cLargestSphereDiameter / bLengthOld(9)

! Old number of cells
pCellsOld = FLOOR( 1.D0 / BoxCutoffOld )

! Box cutoff
BoxCutoffNew(1) = cLargestSphereDiameter / bLengthNew(1)
BoxCutoffNew(2) = cLargestSphereDiameter / bLengthNew(5)
BoxCutoffNew(3) = cLargestSphereDiameter / bLengthNew(9)

! New number of cells
pCellsNew = FLOOR( 1.D0 / BoxCutoffNew )

! Compare old number of cells with new number of cells (remake cell list)
IF( ANY( pCellsNew /= pCellsOld ) ) CALL MakeList( BoxCutoffNew, Position, bLengthInverseNew )

RETURN

END SUBROUTINE BoxCheckNPT

! *********************************************************************************************** !
!               This routine uses the link-list cell structure to fill out an array               !
!                   with possible neighbours of a particle, padding with zeroes                   !
! *********************************************************************************************** !
FUNCTION Neighbours( iParticle, iCellIndex, CellHalfLogical ) RESULT( jNeighbourList )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                          :: iParticle            ! Particle whose neighbours are required
INTEGER( Kind= Int64 )                          :: jParticle            ! Neighbour particles
INTEGER( Kind= Int64 )                          :: sCell, cCell, nCells ! Counter (cells)
INTEGER( Kind= Int64 )                          :: jCounter             ! Counter (particle)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndex           ! Cell of particle of interest
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: jCellIndex           ! Cell of neighbour particles
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList       ! Resulting list of indices

! *********************************************************************************************** !
! INTEGER VARIABLES (PARAMETERS)                                                                  !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER                                                   :: nNeighbourCells = 13  ! Number of neighbouring cells (half)
INTEGER( Kind= Int64 ), DIMENSION( 3, -nNeighbourCells:nNeighbourCells ), PARAMETER :: cBasis = RESHAPE( [ & ! Reference position (basis)
&                       [-1,-1,-1], [0,-1,-1], [1,-1,-1], [-1, 1,-1], [0, 1,-1], [1, 1,-1], &
&                       [-1, 0,-1], [1, 0,-1], [0, 0,-1], [0,-1, 0], [1,-1, 0], [-1,-1, 0], &
&                       [-1, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [-1, 1, 0], [0, 1, 0], &
&                       [0, 0, 1], [-1, 0, 1], [1, 0, 1], [-1,-1, 1], [0,-1, 1], [1,-1, 1], &
&                       [-1, 1, 1], [0, 1, 1], [1, 1, 1] ], [ 3, INT( 2 * nNeighbourCells + 1 ) ] )

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Check half neighbour cells and j down-list from i in current cell
IF( CellHalfLogical ) THEN
  sCell  = 0
  nCells = nNeighbourCells
  ! Debug check
  IF( ANY( iCellIndex(:) /= pCellIndex(:,iParticle) ) ) THEN
    WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndex(:), "] and [", pCellIndex(:,iParticle), "]"
    CALL Sleep( 1 )
    CALL Exit(  )
  END IF
! Check every particle other than i in all cells
ELSE
  sCell  = - nNeighbourCells
  nCells =   nNeighbourCells
END IF

! Initialize with zero values everywhere
jNeighbourList = 0

! Next position in list to be filled
jCounter = 0

! Begin loop over neighbouring cells
DO cCell = sCell, nCells
  ! Neighbour cell index
  jCellIndex(:) = iCellIndex(:) + cBasis(:,cCell)
  ! Periodic boundary correction
  jCellIndex(:) = MODULO( jCellIndex(:), pCells )
  ! Check down-list from i in central cell
  IF ( cCell == 0 .AND. CellHalfLogical ) THEN
    jParticle = pCellList(iParticle)
  ! Check all cells, including the central cell
  ELSE
    jParticle = pCellHead(jCellIndex(1),jCellIndex(2),jCellIndex(3))
  END IF
  ! Begin loop over j particles in list
  ListLoop: DO
    ! Exhausted list
    IF ( jParticle == 0 ) EXIT ListLoop
    ! Skip self
    IF ( jParticle /= iParticle ) THEN
      ! Increment count of j particles
      jCounter = jCounter + 1
      ! Debug check
      IF ( jCounter >= nParticles ) THEN
        WRITE( *, "(5G0)" ) "Error in the list of neighbours: current counter is ", jCounter, " but the number of particles is ", &
        &                   nParticles, "."
        CALL Sleep( 1 )
        CALL Exit(  )
      END IF
      ! Store new neighbour particle
      jNeighbourList(jCounter) = jParticle
    END IF
    ! Next particle in current cell
    jParticle = pCellList(jParticle) 
  END DO ListLoop
END DO

RETURN

END FUNCTION Neighbours

END MODULE LinkedLists