! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                  to eliminate any undesirable rotations of the simulation box.                  !
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
! Recommended Reference:  V. Del Tatto, P. Raiteri, M. Bernetti, G. Bussi                         !
!                         Journal of Applied Sciences, 12(3), 1139 (2022)                         !
!                             DOI: https://doi.org/10.3390/app12031139                            !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!   This subroutine undoes any rotation of the box vectors after performing a lattice reduction   !
!                                        of the box matrix                                        !
! *********************************************************************************************** !
SUBROUTINE UndoBoxRotation( BoxLengthMC, BoxLengthInverseMC, BoxVolumeMC )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: Particle

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                               :: BoxVolumeBoxRotation        ! Box volume (after undoing box rotation)
REAL( Kind= Real64 )                               :: ThetaAngle                  ! Angle between a box vector and a coordinate axis
REAL( Kind= Real64 )                               :: RotationAxisMagnitude       ! Magnitude of rotation axis
REAL( Kind= Real64 )                               :: BoxVolumeMC                 ! Box volume (NPT Simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )               :: BoxLengthMC                 ! Box length (NPT Simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )               :: BoxLengthInverseMC          ! Inverse of box length (NPT Simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )               :: OldBoxLength                ! Box length (NPT Simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )               :: OldBoxLengthInverse         ! Inverse of box length (NPT Simulation)
REAL( Kind= Real64 ), DIMENSION( 3 )               :: ProjectionYXY               ! Projection of the y-vector of the box onto the ZY-plane and the unit vector of the y-axis
REAL( Kind= Real64 ), DIMENSION( 3 )               :: RotationAxis                ! Rotation axis
REAL( Kind= Real64 ), DIMENSION( 3 )               :: RotatedVector               ! Rotated vector
REAL( Kind= Real64 ), DIMENSION( 3 )               :: BoxVectorX                  ! Box vector, X
REAL( Kind= Real64 ), DIMENSION( 3 )               :: BoxVectorY                  ! Box vector, Y
REAL( Kind= Real64 ), DIMENSION( 3 )               :: BoxVectorZ                  ! Box vector, Z
REAL( Kind= Real64 ), DIMENSION( 9 )               :: BoxLengthBoxRotation        ! Box length (after undoing box rotation)
REAL( Kind= Real64 ), DIMENSION( 9 )               :: BoxLengthInverseBoxRotation ! Inverse of box length (after undoing box rotation)
REAL( Kind= Real64 ), DIMENSION( 0:3 )             :: RotationQuaternion          ! Rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 3 )               :: ScalingDistanceUnitBox      ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles )   :: pPositionBoxRotation        ! Position of the center of mass (after undoing box rotation)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles )   :: pOrientationBoxRotation     ! Orientation of the center of mass (after undoing box rotation)
REAL( Kind= Real64 ), DIMENSION( 0:3, nParticles ) :: pQuaternionBoxRotation      ! Quaternion of the center of mass (after undoing box rotation)

! Initialization
OldBoxLength = BoxLengthMC
OldBoxLengthInverse = BoxLengthInverseMC

! Box vectors
BoxVectorX = BoxLengthMC(1:3)
BoxVectorY = BoxLengthMC(4:6)
BoxVectorZ = BoxLengthMC(7:9)

! Angle between x-vector and x-axis
ThetaAngle = DACOS( DOT_PRODUCT( BoxVectorX, [1.D0,0.D0,0.D0] ) / DSQRT( DOT_PRODUCT( BoxVectorX, BoxVectorX ) ) )

! Cross product between x-vector and x-axis (rotation axis)
RotationAxis(1) = 0.D0
RotationAxis(2) = BoxVectorX(3)
RotationAxis(3) = - BoxVectorX(2)

! Magnitude of rotation axis
RotationAxisMagnitude = DSQRT( DOT_PRODUCT( RotationAxis, RotationAxis ) )

! Avoid null vectors
IF( DABS( RotationAxisMagnitude - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  RotationAxis = 0.D0
ELSE
  RotationAxis = RotationAxis / RotationAxisMagnitude
END IF

! Rotation quaternion
RotationQuaternion(0) = DCOS( ThetaAngle * 0.5D0 )                   ! Real part
RotationQuaternion(1) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(1) ! Imaginary part (Vector)
RotationQuaternion(2) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(2) ! Imaginary part (Vector)
RotationQuaternion(3) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(3) ! Imaginary part (Vector)

! Make box x-vector parallel to the x-axis of the coordinate system
IF( DABS( RotationAxisMagnitude - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
  ! Rotated x-vector
  CALL ActiveTransformation( BoxVectorX / ( DSQRT( DOT_PRODUCT( BoxVectorX, BoxVectorX ) ) ), RotationQuaternion, RotatedVector )
  ! New x-vector of the simulation box
  BoxLengthBoxRotation(1:3) = DSQRT( DOT_PRODUCT( BoxLengthMC(1:3), BoxLengthMC(1:3) ) ) * RotatedVector
  BoxLengthBoxRotation(2:3) = 0.D0
  ! Rotated y-vector
  CALL ActiveTransformation( BoxVectorY / ( DSQRT( DOT_PRODUCT( BoxVectorY, BoxVectorY ) ) ), RotationQuaternion, RotatedVector )
  ! New y-vector of the simulation box
  BoxLengthBoxRotation(4:6) = DSQRT( DOT_PRODUCT( BoxLengthMC(4:6), BoxLengthMC(4:6) ) ) * RotatedVector
  ! Rotated z-vector
  CALL ActiveTransformation( BoxVectorZ / ( DSQRT( DOT_PRODUCT( BoxVectorZ, BoxVectorZ ) ) ), RotationQuaternion, RotatedVector )
  ! New z-vector of the simulation box
  BoxLengthBoxRotation(7:9) = DSQRT( DOT_PRODUCT( BoxLengthMC(7:9), BoxLengthMC(7:9) ) ) * RotatedVector
  ! Calculate the new reciprocal box basis vectors
  CALL InverseMatrixCofactorVec( BoxLengthBoxRotation, BoxLengthInverseBoxRotation, BoxVolumeBoxRotation )
  ! Rescale positions and orientations of particles accordingly (the spatial distribution of particles remains unchanged)
  DO Particle = 1, nParticles
    ! Transform spatial coordinates using old box dimensions
    CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionMC(:,Particle), ScalingDistanceUnitBox )
    ! New spatial coordinates using new box dimensions
    CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, pPositionBoxRotation(:,Particle) )
    ! Reorient particles
    CALL QuaternionMultiplication( RotationQuaternion, pQuaternionMC(:,Particle), pQuaternionBoxRotation(:,Particle) )
    ! Active transformation (rotation)
    CALL ActiveTransformation( zAxis, pQuaternionBoxRotation(:,Particle), pOrientationBoxRotation(:,Particle) )
  END DO
! Box x-vector already parallel to the x-axis of the coordinate system
ELSE
  ! Retrieve old box properties
  BoxLengthBoxRotation        = BoxLengthMC
  BoxLengthInverseBoxRotation = BoxLengthInverseMC
  BoxVolumeBoxRotation        = BoxVolumeMC
  ! Retrieve old particle properties
  pPositionBoxRotation    = pPositionMC
  pQuaternionBoxRotation  = pQuaternionMC
  pOrientationBoxRotation = pOrientationMC
END IF

! Initialization
OldBoxLength = BoxLengthBoxRotation
OldBoxLengthInverse = BoxLengthInverseBoxRotation

! Box vectors
BoxVectorX = BoxLengthBoxRotation(1:3)
BoxVectorY = BoxLengthBoxRotation(4:6)
BoxVectorZ = BoxLengthBoxRotation(7:9)

! Axis of rotation
RotationAxis = BoxVectorX

! Magnitude of rotation axis
RotationAxisMagnitude = DSQRT( DOT_PRODUCT( RotationAxis, RotationAxis ) )

! Avoid null vectors
IF( DABS( RotationAxisMagnitude - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  RotationAxis = 0.D0
ELSE
  RotationAxis = RotationAxis / RotationAxisMagnitude
END IF

! Projection of the y-vector of the box onto the ZY-plane
ProjectionYXY = BoxVectorY - ( DOT_PRODUCT( BoxVectorY, RotationAxis ) ) * RotationAxis

! Versor of the projection of the y-vector of the box onto the ZY-plane
ProjectionYXY = ProjectionYXY / DSQRT( DOT_PRODUCT( ProjectionYXY, ProjectionYXY ) )

! Angle between the projection of the y-vector of the box and the y-axis of the coordinate system
ThetaAngle = DACOS( DOT_PRODUCT( ProjectionYXY, [0.D0,1.D0,0.D0] ) / DSQRT( DOT_PRODUCT( ProjectionYXY, ProjectionYXY ) ) )

! Direction of rotation
IF( ProjectionYXY(3) < 0.D0 ) THEN
  ! Rotation quaternion (clockwise rotation)
  RotationQuaternion(0) = DCOS( ThetaAngle * 0.5D0 )                   ! Real part
  RotationQuaternion(1) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(1) ! Imaginary part (Vector)
  RotationQuaternion(2) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(2) ! Imaginary part (Vector)
  RotationQuaternion(3) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(3) ! Imaginary part (Vector)
ELSE
  ! Rotation quaternion (counterclockwise rotation)
  RotationQuaternion(0) = DCOS( - ThetaAngle * 0.5D0 )                   ! Real part
  RotationQuaternion(1) = DSIN( - ThetaAngle * 0.5D0 ) * RotationAxis(1) ! Imaginary part (Vector)
  RotationQuaternion(2) = DSIN( - ThetaAngle * 0.5D0 ) * RotationAxis(2) ! Imaginary part (Vector)
  RotationQuaternion(3) = DSIN( - ThetaAngle * 0.5D0 ) * RotationAxis(3) ! Imaginary part (Vector)
END IF

! Make box y-vector parallel to the XY-plane of the coordinate system
IF( DABS( RotationAxisMagnitude - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
  ! Rotated x-vector
  CALL ActiveTransformation( BoxVectorX / ( DSQRT( DOT_PRODUCT( BoxVectorX, BoxVectorX ) ) ), RotationQuaternion, RotatedVector )
  ! New x-vector of the simulation box
  BoxLengthBoxRotation(1:3) = DSQRT( DOT_PRODUCT( BoxLengthBoxRotation(1:3), BoxLengthBoxRotation(1:3) ) ) * RotatedVector
  BoxLengthBoxRotation(2:3) = 0.D0
  ! Rotated y-vector
  CALL ActiveTransformation( BoxVectorY / ( DSQRT( DOT_PRODUCT( BoxVectorY, BoxVectorY ) ) ), RotationQuaternion, RotatedVector )
  ! New y-vector of the simulation box
  BoxLengthBoxRotation(4:6) = DSQRT( DOT_PRODUCT( BoxLengthBoxRotation(4:6), BoxLengthBoxRotation(4:6) ) ) * RotatedVector
  BoxLengthBoxRotation(6)   = 0.D0
  ! Rotated z-vector
  CALL ActiveTransformation( BoxVectorZ / ( DSQRT( DOT_PRODUCT( BoxVectorZ, BoxVectorZ ) ) ), RotationQuaternion, RotatedVector )
  ! New z-vector of the simulation box
  BoxLengthBoxRotation(7:9) = DSQRT( DOT_PRODUCT( BoxLengthBoxRotation(7:9), BoxLengthBoxRotation(7:9) ) ) * RotatedVector
  ! Calculate the new reciprocal box basis vectors
  CALL InverseMatrixCofactorVec( BoxLengthBoxRotation, BoxLengthInverseBoxRotation, BoxVolumeBoxRotation )
  ! Rescale positions and orientations of particles accordingly (the spatial distribution of particles remains unchanged)
  DO Particle = 1, nParticles
    ! Transform spatial coordinates using old box dimensions
    CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionBoxRotation(:,Particle), ScalingDistanceUnitBox )
    ! New spatial coordinates using new box dimensions
    CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, pPositionBoxRotation(:,Particle) )
    ! Reorient particles
    CALL QuaternionMultiplication( RotationQuaternion, pQuaternionBoxRotation(:,Particle), pQuaternionBoxRotation(:,Particle) )
    ! Active transformation (rotation)
    CALL ActiveTransformation( zAxis, pQuaternionBoxRotation(:,Particle), pOrientationBoxRotation(:,Particle) )
  END DO
END IF

! Update box properties
BoxLengthMC        = BoxLengthBoxRotation
BoxLengthInverseMC = BoxLengthInverseBoxRotation
BoxVolumeMC        = BoxVolumeBoxRotation

! Update particle properties
pPositionMC    = pPositionBoxRotation
pQuaternionMC  = pQuaternionBoxRotation
pOrientationMC = pOrientationBoxRotation

RETURN

END SUBROUTINE UndoBoxRotation