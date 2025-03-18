! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!           This module allows the user to choose the initial molecular configuration.            !
!               The only configuration avaiable at the moment is: packed box (PB).                !
!        Molecules will be then allocated in accordance to the selected crystal structure.        !
!     This module also writes out a file containing all particles' positions and quaternions.     !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       February 15th, 2024                                       !
! ############################################################################################### !
! Main References:                 M. P. Allen, D. J. Tildesley                                   !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             doi: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE InitialConfig

! Uses two modules: global variables and overlap check
USE GlobalVar
USE OverlapCheck

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BodyFixedAxis ! Body-fixed axis of rotation (for the initial configuration only)

! *********************************************************************************************** !
!                                      INITIAL CONFIGURATION                                      !
!                                         Packed box = PB                                         !
! *********************************************************************************************** !
CONTAINS

! *********************************************************************************************** !
!          This subroutine allows the user to choose the initial molecular configuration          !
! *********************************************************************************************** !
SUBROUTINE InitialConfigurationSelection(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = PB                                                                                       !
! *********************************************************************************************** !
ConfigurationSelection = .FALSE.

! *********************************************************************************************** !
! MOLECULAR PROPERTIES                                                                            !
! *********************************************************************************************** !
! Initial configuration file
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )
READ( 100, * ) Dummy, ConfigurationInquiry
CALL ToUpper( ConfigurationInquiry, LEN_TRIM( ConfigurationInquiry ), ConfigurationInquiry )
CLOSE( 100 )

! Extended configuration name
IF( ConfigurationInquiry == "PB" ) THEN
  InitialConfiguration = "Packed Box"
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined ", TRIM( ConfigurationInquiry ), " is not an available initial configuration. Exiting... "
  CALL Exit(  )
END IF

! Initial configuration inquiry
IF( ConfigurationInquiry == "PB" ) THEN
  ConfigurationSelection(1) = .TRUE.
END IF

RETURN

END SUBROUTINE InitialConfigurationSelection

! *********************************************************************************************** !
!            This subroutine allocates particles according to the PB configuration                !
! *********************************************************************************************** !
SUBROUTINE PackedBoxConfiguration(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SeedSize ! Seed array size

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: RandomSeed ! Random seed array

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: sSeed                      ! Counter (seed)
INTEGER( Kind= Int64 )                 :: iParticle, pParticle       ! Counters (particles)
INTEGER( Kind= Int64 )                 :: iCell, jCell, kCell, cCell ! Counters (cells)
INTEGER( Kind= Int64 )                 :: Counter                    ! Counter
INTEGER( Kind= Int64 )                 :: cCylinder                  ! Counter (cylinders)
INTEGER( Kind= Int64 )                 :: pImage                     ! Counter (images)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: rParticles                 ! Number of particles along x-, y-, and z-directions

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO)                                                                 !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCycle                      ! Counter of cycles
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation      ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation         ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter    ! Move counter (Rotation)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance                  ! Vega-Lago contact distance
REAL( Kind= Real64 )                    :: pWidth                           ! Width of the molecular geometry
REAL( Kind= Real64 )                    :: RandomFactor                     ! Expansion factor (initial configuration only)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox           ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: rSpacing                         ! Row spacing
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition       ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iOldQuaternion, iNewQuaternion   ! Quaternion (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciOldPosition, ciNewPosition     ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPosition                 ! Rotated position

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: MaxTranslationalDisplacement ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: MaxAngularDisplacement       ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: Ratio                        ! Acceptance ratio

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: ScalingDistanceImageUnitBox          ! Scaled position of images (unit box)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: iimOldOrientation, iimNewOrientation ! Orientation of the images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: iimOldPosition, iimNewPosition       ! Position of center of mass of the images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: iimOldQuaternion, iimNewQuaternion   ! Quaternion of the images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: iimcOldPosition, iimcNewPosition     ! Position of cylinder images (before/after a trial move)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 03 ) :: Ensemble ! Ensemble type

! *********************************************************************************************** !
! LOGICAL VARIABLES (MONTE CARLO)                                                                 !
! *********************************************************************************************** !
LOGICAL :: MovementRotationLogical    ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap ! Detects overlap between two particles : TRUE = overlap detected; FALSE = overlap not detected

! Allocation
ALLOCATE( ScalingDistanceImageUnitBox(3,nImages) )
ALLOCATE( iimOldPosition(3,nImages), iimNewPosition(3,nImages) )
ALLOCATE( iimcOldPosition(3,4,nImages), iimcNewPosition(3,4,nImages) )
ALLOCATE( iimOldOrientation(3,nImages), iimNewOrientation(3,nImages) )
ALLOCATE( iimOldQuaternion(0:3,nImages), iimNewQuaternion(0:3,nImages) )

! Random number generator seed
IF( RNGeneratorLogical(1) ) THEN ! Fortran
  CALL Random_Seed(  )
ELSE IF( RNGeneratorLogical(2) ) THEN ! Bitwise
  SeedSize = 5
  IF( Allocated( SeedValue ) ) DEALLOCATE( SeedValue )
  ALLOCATE( SeedValue(SeedSize) )
  ! Generate random seed
  CALL Random_Seed( Size= SeedSize )
  ALLOCATE( RandomSeed(SeedSize) )
  CALL Random_Seed( Get= RandomSeed )
  ! Attribute seed
  DO sSeed = 1, 5
    SeedValue(sSeed) = RandomSeed(sSeed)
  END DO
END IF

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AxisSelection(1) ) THEN
  BodyFixedAxis = xAxis
ELSE IF( AxisSelection(2) ) THEN
  BodyFixedAxis = yAxis
ELSE IF( AxisSelection(3) ) THEN
  BodyFixedAxis = zAxis
END IF

! Atom ID (required in some visualization and analysis software)
Atom(1) = "R" ! Real
Atom(2) = "I" ! Image

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
pQuaternion(0,:)  = DCOS( QuaternionAngle * 0.5D0 )                    ! Real part
pQuaternion(1,:)  = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:)  = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:)  = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(3) ! Imaginary part (Vector)

! Volume of the simulation box (initial estimative)
BoxVolume  = 1.D0
iBoxVolume = BoxVolume

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 13 )//"CREATE INITIAL CONFIGURATION"//REPEAT( " ", 14 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0.5,G0)" ) "Initial volume estimative: ", iBoxVolume, "Å³"
WRITE( *, "(G0)" ) " "

! Initialization
rParticles(:) = 0

! Molecular width
IF( PetalArrangementLogical(1) ) THEN
  pWidth = 1.5D0 * cDiameter ! Å
ELSE IF( PetalArrangementLogical(2) ) THEN
  pWidth = ( 1.D0 + 0.5D0 * DSQRT( 2.D0 ) ) * cDiameter ! Å
END IF

! Optimization of the volume of the simulation box
BoxVolumeOptimization: DO
  ! Number density
  TotalNumberDensity = DBLE( nParticles ) / BoxVolume ! Å⁻³
  iNumberDensity     = TotalNumberDensity ! Å⁻³
  ! Packing fraction
  PackingFraction  = TotalNumberDensity * pMolecularVolume
  iPackingFraction = PackingFraction
  ! Cubic box (only diagonal elements)
  BoxLength    = 0.D0
  BoxLength(1) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  BoxLength(5) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  BoxLength(9) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  iBoxLength   = BoxLength
  ! Structure Validation
  IF( ( rParticles(1) * rParticles(2) * rParticles(3) ) < nParticles ) THEN
    ! Update maximum number of particles along x, y and z directions
    rParticles(1) = FLOOR( BoxLength(1) / ( pWidth ) )
    rParticles(2) = FLOOR( BoxLength(5) / ( pWidth ) )
    rParticles(3) = FLOOR( BoxLength(9) / ( cLength ) )
    ! Attempt to increase box volume (new estimative)
    BoxVolume = BoxVolume * 1.01D0 ! Å³
    CYCLE BoxVolumeOptimization
  END IF
  ! Length increment
  rSpacing(1) = ( BoxLength(1) - ( ( pWidth ) * rParticles(1) )  ) / rParticles(1) ! Å
  rSpacing(2) = ( BoxLength(5) - ( ( pWidth ) * rParticles(2) )  ) / rParticles(2) ! Å
  rSpacing(3) = ( BoxLength(9) - ( ( cLength ) * rParticles(3) ) ) / rParticles(3) ! Å
  ! Avoid super compact configurations
  DO cCell = 1, 3
    IF( cAspectRatio <= 1.D0 ) THEN
      IF( rSpacing(cCell) < 1.D-1 * cLength ) THEN ! Arbitrary minimum spacing
        BoxVolume = BoxVolume * 1.01D0
        CYCLE BoxVolumeOptimization
      END IF
    ELSE IF( cAspectRatio > 1.D0 ) THEN
      IF( rSpacing(cCell) < 1.D-1 / cLength ) THEN ! Arbitrary minimum spacing
        BoxVolume = BoxVolume * 1.01D0
        CYCLE BoxVolumeOptimization
      END IF
    END IF
  END DO
  EXIT BoxVolumeOptimization
END DO BoxVolumeOptimization

! Make simulation box N times bigger (random initial volume)
IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
RandomFactor = (10.D0 * RandomNumber) + 10.D0
BoxVolume    = RandomFactor * BoxVolume ! Å³
iBoxVolume   = BoxVolume ! Å³
! New number density
TotalNumberDensity = ( DBLE( nParticles ) / BoxVolume ) ! Å⁻³
iNumberDensity = TotalNumberDensity ! Å⁻³
! New packing fraction
PackingFraction  = TotalNumberDensity * pMolecularVolume
iPackingFraction = PackingFraction
! New cubic box (only diagonal elements)
BoxLength    = 0.D0
BoxLength(1) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
BoxLength(5) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
BoxLength(9) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
iBoxLength   = BoxLength
! New length increment
rSpacing(1) = ( BoxLength(1) - ( ( pWidth ) * rParticles(1) )  ) / rParticles(1) ! Å
rSpacing(2) = ( BoxLength(5) - ( ( pWidth ) * rParticles(2) )  ) / rParticles(2) ! Å
rSpacing(3) = ( BoxLength(9) - ( ( cLength ) * rParticles(3) ) ) / rParticles(3) ! Å

! Status
WRITE( *, "(G0,G0.5,G0)" ) "Final volume estimative: ", iBoxVolume, "Å³" ! Å³
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Positioning particles..."
WRITE( *, "(G0)" ) " "

! Placement of centers of mass
Counter = 1
ParticlePositioning: DO iCell = 1, rParticles(1)
  DO jCell = 1, rParticles(2)
    DO kCell = 1, rParticles(3)
      ! Positioning (centers of mass)
      pPosition(1,Counter) = DBLE( iCell - 1 ) * ( pWidth + rSpacing(1)  ) ! Å
      pPosition(2,Counter) = DBLE( jCell - 1 ) * ( pWidth + rSpacing(2)  ) ! Å
      pPosition(3,Counter) = DBLE( kCell - 1 ) * ( cLength + rSpacing(3) ) ! Å
      ! Cylinders
      CALL ActiveTransformation( cReferencePosition(:,1), pQuaternion(:,Counter), cRotatedPosition(:,1) )
      cPosition(:,1,Counter) = pPosition(:,Counter) + cRotatedPosition(:,1) ! First quarter
      CALL ActiveTransformation( cReferencePosition(:,2), pQuaternion(:,Counter), cRotatedPosition(:,2) )
      cPosition(:,2,Counter) = pPosition(:,Counter) + cRotatedPosition(:,2) ! Second quarter
      CALL ActiveTransformation( cReferencePosition(:,3), pQuaternion(:,Counter), cRotatedPosition(:,3) )
      cPosition(:,3,Counter) = pPosition(:,Counter) + cRotatedPosition(:,3) ! Third quarter
      CALL ActiveTransformation( cReferencePosition(:,4), pQuaternion(:,Counter), cRotatedPosition(:,4) )
      cPosition(:,4,Counter) = pPosition(:,Counter) + cRotatedPosition(:,4) ! Fourth quarter
      ! Increment
      Counter = Counter + 1
      ! Condition
      IF( Counter > nParticles ) THEN
        EXIT ParticlePositioning
      END IF
    END DO
  END DO
END DO ParticlePositioning

! Simulation box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO pParticle = 1, nParticles
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,pParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,pParticle) )
  DO cCylinder = 1, 4
    CALL MatrixVectorMultiplication( BoxLengthInverse, cPosition(:,cCylinder,pParticle), ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, cPosition(:,cCylinder,pParticle) )
  END DO
END DO

! Positioning images of particles
DO pParticle = 1, nParticles
  ! Transform spatial coordinates into coordinates of a unit box
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,pParticle), ScalingDistanceUnitBox )
  ! Build image list
  CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
  ! Particle images
  DO pImage = 1, nImages
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceImageUnitBox(:,pImage), imPosition(:,pImage,pParticle) )
  END DO
  DO cCylinder = 1, 4
    ! Transform spatial coordinates into coordinates of a unit box
    CALL MatrixVectorMultiplication( BoxLengthInverse, cPosition(:,cCylinder,pParticle), ScalingDistanceUnitBox )
    ! Build image list
    CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
    ! Cylinder images
    DO pImage = 1, nImages
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceImageUnitBox(:,pImage), &
      &                                imcPosition(:,cCylinder,pImage,pParticle) )
    END DO
  END DO
  ! Quaternions of images
  DO pImage = 1, nImages
    imQuaternion(0,pImage,pParticle) = pQuaternion(0,pParticle)
    imQuaternion(1,pImage,pParticle) = pQuaternion(1,pParticle)
    imQuaternion(2,pImage,pParticle) = pQuaternion(2,pParticle)
    imQuaternion(3,pImage,pParticle) = pQuaternion(3,pParticle)
  END DO
END DO

! Active transformation
DO pParticle = 1, nParticles
  CALL ActiveTransformation( BodyFixedAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
END DO
DO pParticle = 1, nParticles
  DO pImage = 1, nImages
    imOrientation(1,pImage,pParticle) = pOrientation(1,pParticle)
    imOrientation(2,pImage,pParticle) = pOrientation(2,pParticle)
    imOrientation(3,pImage,pParticle) = pOrientation(3,pParticle)
  END DO
END DO

! *********************************************************************************************** !
! Monte Carlo parameters (NVT simulation)                                                         !
! *********************************************************************************************** !
MovementTranslationLogical   = .FALSE.            ! Translational move selector           (initial value)
MovementRotationLogical      = .FALSE.            ! Rotational move selector              (initial value)
MaxTranslationalDisplacement = UserMaxTranslation ! Maximum translational displacement    (initial value)
MaxAngularDisplacement       = UserMaxRotation    ! Maximum rotational displacement       (initial value)
nAcceptanceTranslation       = 0                  ! Translational move acceptance counter (initial value)
nAcceptanceRotation          = 0                  ! Rotational move acceptance counter    (initial value)
nMovementTranslationCounter  = 0                  ! Translational move counter            (initial value)
nMovementRotationCounter     = 0                  ! Rotational move counter               (initial value)
iCycle                       = 0                  ! Number of cycles                      (initial value)
cPositionMC                  = cPosition          ! Position (Cylinders)                  (initial value)
imcPositionMC                = imcPosition        ! Position (Cylinder images)            (initial value)
pQuaternionMC                = pQuaternion        ! Quaternion algebra                    (initial value)
imQuaternionMC               = imQuaternion       ! Quaternion algebra (images)           (initial value)
pPositionMC                  = pPosition          ! Position of particles                 (initial value)
imPositionMC                 = imPosition         ! Position of images                    (initial value)
pOrientationMC               = pOrientation       ! Orientation of particles              (initial value)
imOrientationMC              = imOrientation      ! Orientation of images                 (initial value)
Ensemble                     = "NVT"              ! Ensemble type                         (initial value)

WRITE( *, "(3G0)" ) "Starting an NVT-Monte Carlo simulation with ", InitialConfigurationMaxCycles, &
&                   " cycle(s) to randomly distribute particles..."
WRITE( *, "(G0)" ) " "

! Displace particles as long as the number of cycles is lower than the specified value
DO WHILE( iCycle <= InitialConfigurationMaxCycles )

  ! Progress bar
  CALL ProgressBarMC( iCycle, InitialConfigurationMaxCycles, Ensemble )

  ! Random particle
  IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
  IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
  iParticle = INT( DBLE( nParticles ) * RandomNumber ) + 1

  ! Assignment of previous configuration (Microstate m)
  iOldPosition(:)        = pPositionMC(:,iParticle)       ! Position
  iimOldPosition(:,:)    = imPositionMC(:,:,iParticle)    ! Position (image)
  iOldQuaternion(:)      = pQuaternionMC(:,iParticle)     ! Quaternion
  iimOldQuaternion(:,:)  = imQuaternionMC(:,:,iParticle)  ! Quaternion (image)
  ciOldPosition(:,:)     = cPositionMC(:,:,iParticle)     ! Position (cylinders)
  iimcOldPosition(:,:,:) = imcPositionMC(:,:,:,iParticle) ! Position (cylinder image)
  iOldOrientation(:)     = pOrientationMC(:,iParticle)    ! Orientation
  iimOldOrientation(:,:) = imOrientationMC(:,:,iParticle) ! Orientation (image)

  ! Pseudorandom number generator (uniform distribution)
  IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
  IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

  ! Translation criterion
  IF( RandomNumber < TranslationalProbability ) THEN
    MovementTranslationLogical  = .TRUE.  ! Enable translation
    MovementRotationLogical     = .FALSE. ! Disable rotation
    nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
  ! Rotation criterion
  ELSE IF( RandomNumber >= TranslationalProbability ) THEN
    MovementRotationLogical    = .TRUE.  ! Enable rotation
    MovementTranslationLogical = .FALSE. ! Disable translation
    nMovementRotationCounter   = nMovementRotationCounter + 1 ! Increment move counter
  END IF

  ! Translational movement
  IF( MovementTranslationLogical ) THEN
    ! Random translation along x-axis
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
    iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Random translation along y-axis
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
    iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Random translation along z-axis
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
    iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
    ! Transform spatial coordinates into coordinates of a unit box
    CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
    ! Build image list
    CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
    ! Positioning images of particle i
    DO pImage = 1, nImages
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceImageUnitBox(:,pImage), iimNewPosition(:,pImage) )
    END DO
  ! No translation
  ELSE IF( .NOT. MovementTranslationLogical ) THEN
    iNewPosition   = iOldPosition
    iimNewPosition = iimOldPosition
  END IF

  ! Rotational movement
  IF( MovementRotationLogical ) THEN
    ! Random Composed Unit Quaternion (see 'subroutines' code for more details)
    CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
    ! Active transformation (rotation)
    CALL ActiveTransformation( BodyFixedAxis, iNewQuaternion, iNewOrientation )
    ! Orientation of images of particle i
    DO pImage = 1, nImages
      iimNewQuaternion(:,pImage)  = iNewQuaternion(:)
      iimNewOrientation(:,pImage) = iNewOrientation(:)
    END DO
  ! No rotation
  ELSE IF( .NOT. MovementRotationLogical ) THEN
    iNewQuaternion    = iOldQuaternion
    iNewOrientation   = iOldOrientation
    iimNewQuaternion  = iimOldQuaternion
    iimNewOrientation = iimOldOrientation
  END IF

  ! Position of cylinders (after translation/rotation)
  DO cCylinder = 1, 4
    ! Active transformation (translation)
    CALL ActiveTransformation( cReferencePosition(:,cCylinder), iNewQuaternion, cRotatedPosition(:,cCylinder) )
    ciNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPosition(:,cCylinder)
    ! Transform spatial coordinates into coordinates of a unit box
    CALL MatrixVectorMultiplication( BoxLengthInverse, ciNewPosition(:,cCylinder), ScalingDistanceUnitBox )
    ! Build image list
    CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
    ! Position of cylinder images of particle i
    DO pImage = 1, nImages
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceImageUnitBox(:,pImage), iimcNewPosition(:,cCylinder,pImage) )
    END DO
  END DO

  ! Overlap check
  CALL ParticleOverlapCheck( iParticle, iNewQuaternion, iimNewQuaternion, iNewOrientation, iimNewOrientation, iNewPosition, &
  &                          iimNewPosition, ciNewPosition, iimcNewPosition, ContactDistance, Overlap )

  ! Acceptance criterion
  IF( .NOT. Overlap ) THEN
    ! Update system configuration
    pPositionMC(:,iParticle)       = iNewPosition(:)        ! Update position
    imPositionMC(:,:,iParticle)    = iimNewPosition(:,:)    ! Update position (images)
    cPositionMC(:,:,iParticle)     = ciNewPosition(:,:)     ! Update position (cylinders)
    imcPositionMC(:,:,:,iParticle) = iimcNewPosition(:,:,:) ! Update position (cylinder images)
    pQuaternionMC(:,iParticle)     = iNewQuaternion(:)      ! Update quaternion
    imQuaternionMC(:,:,iParticle)  = iimNewQuaternion(:,:)  ! Update quaternion (images)
    pOrientationMC(:,iParticle)    = iNewOrientation(:)     ! Update orientation
    imOrientationMC(:,:,iParticle) = iimNewOrientation(:,:) ! Update orientation (images)
    ! Update displacement counter
    IF( MovementTranslationLogical ) THEN
      nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
    ELSE IF( MovementRotationLogical ) THEN
      nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
    END IF
  ELSE
    ! Retrieve old configuration
    pPositionMC(:,iParticle)       = iOldPosition(:)        ! Retrieve old position
    imPositionMC(:,:,iParticle)    = iimOldPosition(:,:)    ! Retrieve old position (images)
    cPositionMC(:,:,iParticle)     = ciOldPosition(:,:)     ! Retrieve old position (cylinders)
    imcPositionMC(:,:,:,iParticle) = iimcOldPosition(:,:,:) ! Retrieve old position (cylinder images)
    pQuaternionMC(:,iParticle)     = iOldQuaternion(:)      ! Retrieve old quaternion
    imQuaternionMC(:,:,iParticle)  = iimOldQuaternion(:,:)  ! Retrieve old quaternion (images)
    pOrientationMC(:,iParticle)    = iOldOrientation(:)     ! Retrieve old orientation
    imOrientationMC(:,:,iParticle) = iimOldOrientation(:,:) ! Retrieve old orientation (images)
  END IF

  ! Adjustment of maximum translation
  IF( MOD( iCycle, nAdjustmentMovementFrequency ) == 0 .AND. nMovementTranslationCounter > 0 ) THEN
    ! Acceptance ratio (non-overlapping microstates over sampled microstates)
    Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
    ! Translational adjustment
    IF( Ratio <= AcceptanceRatioTranslation ) THEN
      MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
    ELSE
      MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
    END IF
    ! Avoid large translations
    IF( ( MaxTranslationalDisplacement >= 2.D0 * MAXVAL( BoxLength ) ) ) THEN
      MaxTranslationalDisplacement = MaxTranslationalDisplacement - 1.0D0 * MAXVAL( BoxLength )
    END IF
    ! Reset counter
    nAcceptanceTranslation      = 0
    nMovementTranslationCounter = 0
  END IF

  ! Adjustment of maximum rotation
  IF( MOD( iCycle, nAdjustmentMovementFrequency ) == 0 .AND. nMovementRotationCounter > 0 ) THEN
    ! Acceptance ratio (non-overlapping microstates over sampled microstates)
    Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
    ! Rotational adjustment
    IF( Ratio <= AcceptanceRatioRotation ) THEN
      MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
    ELSE
      MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
    END IF
    ! Avoid 2π-rotations (multiple turns)
    IF( ( MaxAngularDisplacement >= 4.D0 * cPi ) ) THEN
      MaxAngularDisplacement = MaxAngularDisplacement - 2.0D0 * cPi
    END IF
    ! Reset counter
    nAcceptanceRotation      = 0
    nMovementRotationCounter = 0
  END IF

  ! Update number of cycles
  iCycle = iCycle + 1

END DO

! Deallocation
DEALLOCATE( ScalingDistanceImageUnitBox, iimOldPosition, iimNewPosition, iimOldQuaternion, iimNewQuaternion, iimcOldPosition, &
&           iimcNewPosition )
IF( Allocated( RandomSeed ) ) DEALLOCATE( RandomSeed )

WRITE( *, * ) " "
WRITE( *, * ) " "

! *********************************************************************************************** !
! Configuration update                                                                            !
! *********************************************************************************************** !
pPosition     = pPositionMC     ! Update position
imPosition    = imPositionMC    ! Update position (images)
cPosition     = cPositionMC     ! Update position (cylinders)
imcPosition   = imcPositionMC   ! Update position (cylinder images)
pQuaternion   = pQuaternionMC   ! Update quaternion
imQuaternion  = imQuaternionMC  ! Update quaternion (images)
pOrientation  = pOrientationMC  ! Update orientation
imOrientation = imOrientationMC ! Update orientation (images)

RETURN

END SUBROUTINE PackedBoxConfiguration

! *********************************************************************************************** !
!          This subroutine creates a file with all particles' positions and orientations          !
! *********************************************************************************************** !
SUBROUTINE ConfigurationOutput(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle ! Counter (particles)
INTEGER( Kind= Int64 ) :: pImage    ! Counter (images)
INTEGER( Kind= Int64 ) :: cCylinder ! Counter (cylinders)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: TrajectoryPosition ! Position array (cylinders)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Packed-box structure
IF( ConfigurationSelection(1) ) THEN
  OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_initconf_N" &
  &                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//"_pb.xyz" )
  WRITE( 10, "(G0)" ) nParticles * 4 + nParticles * nImages * 4
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', (1.D0 + 2.D0 * DBLE( nLayers) ) * BoxLength(1:9), '" Origin="', -0.5D0 * (1.D0 + 2.D0 &
  &                             * DBLE( nLayers) ) * ( BoxLength(1) + BoxLength(4) + BoxLength(7) ), -0.5D0 * (1.D0 + 2.D0 * &
  &                             DBLE( nLayers) ) * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * (1.D0 + 2.D0 * &
  &                             DBLE( nLayers) ) * ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO pParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of cylinders
      TrajectoryPosition(1) = cPosition(1,cCylinder,pParticle)
      TrajectoryPosition(2) = cPosition(2,cCylinder,pParticle)
      TrajectoryPosition(3) = cPosition(3,cCylinder,pParticle)
      WRITE( 10, * ) Atom(1), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), pQuaternion(1,pParticle), &
      &              pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), 0.5D0 * cDiameter, &
      &              0.5D0 * cDiameter, cLength
      ! Position of cylinder images
      DO pImage = 1, nImages
        TrajectoryPosition(1) = imcPosition(1,cCylinder,pImage,pParticle)
        TrajectoryPosition(2) = imcPosition(2,cCylinder,pImage,pParticle)
        TrajectoryPosition(3) = imcPosition(3,cCylinder,pImage,pParticle)
        WRITE( 10, * ) Atom(2), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), &
        &              imQuaternion(1,pImage,pParticle), imQuaternion(2,pImage,pParticle), imQuaternion(3,pImage,pParticle), &
        &              imQuaternion(0,pImage,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      END DO
    END DO
  END DO
  FLUSH( 10 )
  CLOSE( 10 )
END IF

RETURN

END SUBROUTINE ConfigurationOutput

END MODULE InitialConfig