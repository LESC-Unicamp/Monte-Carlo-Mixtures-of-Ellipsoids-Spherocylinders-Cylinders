! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!         This code contains the algorithms of the force fields used in the main program          !
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
MODULE ForceFields

! Uses one module: global variables
USE GlobalVar
USE OverlapCheckAlgorithms, ONLY: OverlapCheckEOR, OverlapCheckSPC, OverlapCheckCYL, OverlapCheckCYLSPH

IMPLICIT NONE

CONTAINS

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

! *********************************************************************************************** !
!            This subroutine computes the pair potential between particles i and j by             !
!       applying a discrete anisotropic square-well potential to compute the pair potential       !
! *********************************************************************************************** !
SUBROUTINE AnisotropicSquareWellPotential( iPosition, iOrientation, jOrientation, iQuaternion, jQuaternion, VectorDistance, &
&                                          SquaredDistance, iComponent, jComponent, PairPotentialEnergy )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 ) :: rRange                 ! Counter (range)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                      :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                      :: SquaredPotentialRange      ! Effective range of attraction (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )      :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )      :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )    :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy        ! Pair potential energy

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap     ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Compute pair potential for every attractive range
DO rRange = 1, nRange
  ! Effective range of attraction (squared)
  SquaredPotentialRange = 0.5D0 * ( cCircumscribingPotentialRange(iComponent,rRange) + &
  &                                 cCircumscribingPotentialRange(jComponent,rRange) )
  SquaredPotentialRange = SquaredPotentialRange * SquaredPotentialRange
  ! Possible overlap in the region of attraction
  IF( SquaredDistance <= SquaredPotentialRange ) THEN
    IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
      CALL OverlapCheckEOR( iQuaternion, jQuaternion, VectorDistance, SquaredDistance, iComponent, jComponent, ContactDistance, &
      &                     Overlap, .TRUE., rRange )
    ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, iComponent, jComponent, ContactDistance, &
      &                     ParallelSPC, Overlap, .TRUE., rRange )
    ELSE IF( GeometryType(3) ) THEN ! Cylinders
      IF( .NOT. SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        ! Overlap test for cylinders (modified algorithm of Lopes et al.)
        CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
        &                     iComponent, jComponent, .FALSE., Overlap, .TRUE., rRange )
      ! Overlap test for cylinders and spheres
      ELSE IF( .NOT. SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( iComponent, jComponent, iQuaternion, iPosition, jPosition, Overlap, .TRUE., rRange )
      ! Overlap test for cylinders and spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. .NOT. SphericalComponentLogical(jComponent) ) THEN
        ! Apply periodic boundary conditions on the position of particle j
        jPosition(1) = iPosition(1) + VectorDistance(1)
        jPosition(2) = iPosition(2) + VectorDistance(2)
        jPosition(3) = iPosition(3) + VectorDistance(3)
        CALL OverlapCheckCYLSPH( jComponent, iComponent, jQuaternion, jPosition, iPosition, Overlap, .TRUE., rRange )
      ! Overlap test for spheres
      ELSE IF( SphericalComponentLogical(iComponent) .AND. SphericalComponentLogical(jComponent) ) THEN
        Overlap = .TRUE.
      END IF
    END IF
    ! Attraction occurs only when anisotropic shells overlap
    IF( Overlap ) THEN
      ! Pair potential (reduced units)
      PairPotentialEnergy(rRange) = - 1.D0
    ELSE
      PairPotentialEnergy(rRange) = 0.D0
    END IF
  ELSE
    PairPotentialEnergy(rRange) = 0.D0
  END IF
END DO

RETURN

END SUBROUTINE AnisotropicSquareWellPotential

END MODULE ForceFields
