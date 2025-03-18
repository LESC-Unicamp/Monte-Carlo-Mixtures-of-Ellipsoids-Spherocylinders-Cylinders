! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!  This algorithm employs the FBMC method (J. Chem. Phys. 137, 214101, 2012) to predict crystal   !
!                                  arrangements in solid phases.                                  !
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
! Main References:    J. Graaf, L. Filion, M. Marechal, R. Roij, M. Dijkstra                      !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
!                             ---------------------------------------                             !
!                                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
PROGRAM FloppyBoxMonteCarlo

! Uses five modules: global variables, variable initialization, initial configuration, overlap check, and directory creator
USE GlobalVar
USE VariableInitialization
USE InitialConfig
USE OverlapCheck
USE Folders

! Use intrinsic module for OUTPUT_UNIT
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit

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
INTEGER( Kind= Int64 ) :: sSeed                ! Counter (seed)
INTEGER( Kind= Int64 ) :: iParticle, pParticle ! Counters (particles)
INTEGER( Kind= Int64 ) :: cCylinder            ! Counters (cylinders)
INTEGER( Kind= Int64 ) :: bEdge                ! Counter (box edges)
INTEGER( Kind= Int64 ) :: pImage               ! Counter (images)
INTEGER( Kind= Int64 ) :: rPressure            ! Counter (pressure ramp)
INTEGER( Kind= Int64 ) :: BoxMatrixComponent   ! Box matrix component

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO)                                                                 !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCycle                             ! Counter of cycles
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation             ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation                ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nAcceptanceIsotropicVolumeChange   ! Move acceptance counter: Isotropic volume change
INTEGER( Kind= Int64 ) :: nAcceptanceAnisotropicVolumeChange ! Move acceptance counter: Anisotropic volume change
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter        ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter           ! Move counter (Rotation)
INTEGER( Kind= Int64 ) :: nMovementIsoVolumeChangeCounter    ! Move counter (Isotropic)
INTEGER( Kind= Int64 ) :: nMovementAnisoVolumeChangeCounter  ! Move counter (Anisotropic)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance                  ! Vega-Lago contact distance (variable)
REAL( Kind= Real64 )                    :: pSphereDistance                  ! Cutoff diameter (larger sphere)
REAL( Kind= Real64 )                    :: pSpheroCylinderDistance          ! Cutoff diameter (spherocylinder)
REAL( Kind= Real64 )                    :: OrderParameterMC                 ! Nematic order parameter
REAL( Kind= Real64 )                    :: cSphereDistance                  ! Cutoff diameter (smaller sphere)
REAL( Kind= Real64 )                    :: OldBoxVolume, NewBoxVolume       ! Volume of simulation box (before/after a trial move)
REAL( Kind= Real64 )                    :: EnthalpyChange                   ! Enthalpy criterion (reduced)
REAL( Kind= Real64 )                    :: VolumeScalingFactor              ! Scaling factor (volume)
REAL( Kind= Real64 )                    :: BoxDistortionMC                  ! Box distortion
REAL( Kind= Real64 )                    :: BoxVolumeMC                      ! Volume of simulation box
REAL( Kind= Real64 ), DIMENSION( 3 )    :: TrajectoryPosition               ! Position array (cylinders)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxVectorAngle                   ! Cossine of angle between box vectors
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthMC                      ! Length (x,y,z) of simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverseMC               ! Length (x,y,z) of simulation box (inverse)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxEdgeLength                    ! Length of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxEdgeRatio                     ! Length ratio of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )    :: RescalingFactor                  ! Rescaling distance factor for the cylinders
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox           ! Scaled coordinates (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition       ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: OldBoxLength, NewBoxLength       ! Length of simulation box (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: OldBoxLengthInverse              ! Length of simulation box (before a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: NewBoxLengthInverse              ! Length of simulation box (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciOldPosition, ciNewPosition     ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPosition                 ! Rotated position
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iOldQuaternion, iNewQuaternion   ! Quaternion (before/after a trial move)

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: MaxTranslationalDisplacement   ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: MaxAngularDisplacement         ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: MaxIsoVolumetricDisplacement   ! Maximum displacement [+/-] (Isotropic)
REAL( Kind= Real64 ) :: MaxAnisoVolumetricDisplacement ! Maximum displacement [+/-] (Anisotropic)
REAL( Kind= Real64 ) :: Ratio                          ! Acceptance ratio (Simulation)

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: ScalingDistanceImageUnitBox          ! Position of images (unit box)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: iimOldOrientation, iimNewOrientation ! Orientation of images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: iimOldPosition, iimNewPosition       ! Position of images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: iimOldQuaternion, iimNewQuaternion   ! Quaternion of images (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: pPositionSaveMC                      ! Position array
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: cPositionSaveMC                      ! Position array (cylinders)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: imPositionSaveMC                     ! Position array (images)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: iimcOldPosition, iimcNewPosition     ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( :, :, :, : ), ALLOCATABLE :: imcPositionSaveMC                    ! Position array (cylinder images)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap                 ! Detects overlap between two particles (Vega-Lago) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: LatticeReductionLogical ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: CheckBoxDistortion      ! Detects if a box deformation is valid or not : TRUE = ignore box deformation; FALSE = consider box deformation

! *********************************************************************************************** !
! LOGICAL VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
LOGICAL :: MovementRotationLogical          ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical       ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementIsoVolumeChangeLogical   ! Isotropic volume change selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementAnisoVolumeChangeLogical ! Anisotropic volume change selection : TRUE = movement selected; FALSE = movement not selected

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 02 )  :: VolumeType       ! Expansion/Compression type
CHARACTER( LEN= 03 )  :: Ensemble         ! Ensemble type
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Title
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 11 )//"FLOPPY-BOX MONTE CARLO SIMULATION"//REPEAT( " ", 11 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) REPEAT( " ", 13 )//"Author: Nathan Barros de Souza"
WRITE( *, "(G0)" ) REPEAT( " ", 13 )//"E-mail: n264179@dac.unicamp.br"
WRITE( *, "(G0)" ) " "

! Molecular configuration selection
CALL InitialConfigurationSelection(  )

! X-axis
xAxis(1) = 1.D0
xAxis(2) = 0.D0
xAxis(3) = 0.D0
! Y-axis
yAxis(1) = 0.D0
yAxis(2) = 1.D0
yAxis(3) = 0.D0
! Z-axis
zAxis(1) = 0.D0
zAxis(2) = 0.D0
zAxis(3) = 1.D0

! Surface geometry selection
CALL GeometrySelection(  )

! CPU Clock
CALL DATE_AND_TIME( Values= DateTime )

! Date format (YYYY/MM/DD)
FormatDate = "(I4,2I2.2)"
! Date descriptor
WRITE( DescriptorDate, FormatDate ) DateTime(1), DateTime(2), DateTime(3)
! Time format (HH:MM:SS)
FormatHour = "(3I2.2)"
! Hour descriptor
WRITE( DescriptorHour, FormatHour ) DateTime(5), DateTime(6), DateTime(7)

! Initialization of common variables
CALL CommonVariables(  )

! Output file descriptors (Aspect ratio [1])
FormatFileGeometry = "(F0.5)"
WRITE( DescriptorFileGeometry, FormatFileGeometry ) cAspectRatio

! Output file descriptors (Number of particles [2])
FormatFileParticles = "(I0.3)"
WRITE( DescriptorFileParticles, FormatFileParticles ) nParticles

! Allocation (particles)
ALLOCATE( pPosition(3,nParticles) )
ALLOCATE( pOrientation(3,nParticles) )
ALLOCATE( pQuaternion(0:3,nParticles) )
ALLOCATE( cPosition(3,4,nParticles) )
ALLOCATE( pPositionMC(3,nParticles) )
ALLOCATE( pOrientationMC(3,nParticles) )
ALLOCATE( pQuaternionMC(0:3,nParticles) )
ALLOCATE( cPositionMC(3,4,nParticles) )
ALLOCATE( pPositionSaveMC(3,nParticles) )
ALLOCATE( cPositionSaveMC(3,4,nParticles) )
! Allocation (images)
ALLOCATE( ScalingDistanceImageUnitBox(3,nImages) )
ALLOCATE( imQuaternion(0:3,nImages,nParticles) )
ALLOCATE( imPosition(3,nImages,nParticles) )
ALLOCATE( imOrientation(3,nImages,nParticles) )
ALLOCATE( imcPosition(3,4,nImages,nParticles) )
ALLOCATE( imQuaternionMC(0:3,nImages,nParticles) )
ALLOCATE( imPositionMC(3,nImages,nParticles) )
ALLOCATE( imPositionSaveMC(3,nImages,nParticles) )
ALLOCATE( imOrientationMC(3,nImages,nParticles) )
ALLOCATE( imcPositionMC(3,4,nImages,nParticles) )
ALLOCATE( imcPositionSaveMC(3,4,nImages,nParticles) )
ALLOCATE( iimOldOrientation(3,nImages), iimNewOrientation(3,nImages) )
ALLOCATE( iimOldPosition(3,nImages), iimNewPosition(3,nImages) )
ALLOCATE( iimOldQuaternion(0:3,nImages), iimNewQuaternion(0:3,nImages) )
ALLOCATE( iimcOldPosition(3,4,nImages), iimcNewPosition(3,4,nImages) )

! Initialization of Monte Carlo parameters
CALL MonteCarloVariables(  )

! Initialization of Inquiry/Control variables
CALL VariablesInquest(  )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"FOLDER ORGANIZER"//REPEAT( " ", 20 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Setting up folders..."
WRITE( *, "(G0)" ) " "

! Initial configuration folder
CALL InitialConfigurationFolder(  )

! Create output directories
CALL InquireFolders(  )

! Create date subfolders
CALL DateFolders(  )

! *********************************************************************************************** !
! Cutoff diameter (non-convex body)                                                               !
! *********************************************************************************************** !
!  A cutoff distance equivalent to the diameter of a spherical geometry circumscribing            !
!  the spherocylinder circumscribing the non-convex body.                                         !
! *********************************************************************************************** !
IF( PetalArrangementLogical(1) ) THEN
  pSphereDistance                = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter + cLength )
  pSquaredSphereDistance         = pSphereDistance * pSphereDistance
  pSpheroCylinderDistance        = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter )
  pSquaredSpheroCylinderDistance = pSpheroCylinderDistance * pSpheroCylinderDistance
ELSE IF( PetalArrangementLogical(2) ) THEN
  pSphereDistance                = ( 2.D0 * cDiameter + cLength )
  pSquaredSphereDistance         = pSphereDistance * pSphereDistance
  pSpheroCylinderDistance        = 2.D0 * cDiameter
  pSquaredSpheroCylinderDistance = pSpheroCylinderDistance * pSpheroCylinderDistance
END IF

! Cutoff diameter (cylinders)
cSphereDistance        = ( cDiameter + cLength )
cSquaredSphereDistance = cSphereDistance * cSphereDistance

! Unrotated reference position
IF( PetalArrangementLogical(1) ) THEN ! Arrangement [1]
  ! First quarter
  cReferencePosition(1,1) = 0.25D0 * cDiameter ! Å
  cReferencePosition(2,1) = 0.25D0 * cDiameter ! Å
  cReferencePosition(3,1) = 0.D0               ! Å
  ! Second quarter
  cReferencePosition(1,2) = - 0.25D0 * cDiameter ! Å
  cReferencePosition(2,2) = 0.25D0 * cDiameter   ! Å
  cReferencePosition(3,2) = 0.D0                 ! Å
  ! Third quarter
  cReferencePosition(1,3) = - 0.25D0 * cDiameter ! Å
  cReferencePosition(2,3) = - 0.25D0 * cDiameter ! Å
  cReferencePosition(3,3) = 0.D0                 ! Å
  ! Fourth quarter
  cReferencePosition(1,4) = 0.25D0 * cDiameter   ! Å
  cReferencePosition(2,4) = - 0.25D0 * cDiameter ! Å
  cReferencePosition(3,4) = 0.D0                 ! Å
ELSE IF( PetalArrangementLogical(2) ) THEN ! Arrangement [2]
  ! First quarter
  cReferencePosition(1,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(2,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(3,1) = 0.D0                               ! Å
  ! Second quarter
  cReferencePosition(1,2) = - 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(2,2) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter   ! Å
  cReferencePosition(3,2) = 0.D0                                 ! Å
  ! Third quarter
  cReferencePosition(1,3) = - 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(2,3) = - 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(3,3) = 0.D0                                 ! Å
  ! Fourth quarter
  cReferencePosition(1,4) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter   ! Å
  cReferencePosition(2,4) = - 0.25D0 * DSQRT( 2.D0 ) * cDiameter ! Å
  cReferencePosition(3,4) = 0.D0                                 ! Å
END IF

! *********************************************************************************************** !
! Initial configuration                                                                           !
! *********************************************************************************************** !
! Calls 'config_pb' subroutine if the user chooses a packed box structure
IF( ConfigurationSelection(1) ) THEN
  CALL PackedBoxConfiguration(  )
END IF

! Initial configuration file
CALL ConfigurationOutput(  )

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

! Start simulation timer
CALL CPU_Time( StartTimer )

! Active transformation (rotation)
DO pParticle = 1, nParticles
  CALL ActiveTransformation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
END DO

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 6 )//"OVERLAP CHECK FOR THE INITIAL CONFIGURATION"//REPEAT( " ", 6 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Verifying initial configuration..."
WRITE( *, "(G0)" ) " "

! Search for overlapping particles in initial configuration
CALL OverlapCheckInitialConfiguration(  )

! Status
WRITE( *, "(G0)" ) "No overlaps found in the initial configuration. Resuming..."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"FILE ORGANIZER"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)", Advance= "NO" ) "Creating files..."

! Output file units
CALL FileHandler(  )

! Status
CALL Sleep( 1 )
WRITE( *, "(G0)", Advance= "YES" ) " Done!"
WRITE( *, "(G0)" ) " "

! Initialize box distortion parameter if restoration is not selected
CALL LatticeReduction( BoxLength, BoxDistortionMC, LatticeReductionLogical )

! *********************************************************************************************** !
! Monte Carlo Parameters                                                                          !
! *********************************************************************************************** !
MovementTranslationLogical         = .FALSE.               ! Translational move selector           (initial value)
MovementRotationLogical            = .FALSE.               ! Rotational move selector              (initial value)
MovementIsoVolumeChangeLogical     = .FALSE.               ! Volume move selector                  (initial value)
MovementAnisoVolumeChangeLogical   = .FALSE.               ! Volume move selector                  (initial value)
MaxTranslationalDisplacement       = UserMaxTranslation    ! Maximum translational displacement    (initial value)
MaxAngularDisplacement             = UserMaxRotation       ! Maximum rotational displacement       (initial value)
MaxIsoVolumetricDisplacement       = UserMaxIsotropicVol   ! Maximum volumetric displacement       (initial value)
MaxAnisoVolumetricDisplacement     = UserMaxAnisotropicVol ! Maximum volumetric displacement       (initial value)
nAcceptanceTranslation             = 0                     ! Translational move acceptance counter (initial value)
nAcceptanceRotation                = 0                     ! Rotational move acceptance counter    (initial value)
nAcceptanceIsotropicVolumeChange   = 0                     ! Volumetric move acceptance counter    (initial value)
nAcceptanceAnisotropicVolumeChange = 0                     ! Volumetric move acceptance counter    (initial value)
nMovementTranslationCounter        = 0                     ! Translational move counter            (initial value)
nMovementRotationCounter           = 0                     ! Rotational move counter               (initial value)
nMovementIsoVolumeChangeCounter    = 0                     ! Volume change counter                 (initial value)
nMovementAnisoVolumeChangeCounter  = 0                     ! Volume change counter                 (initial value)
cPositionMC                        = cPosition             ! Position (cylinders)                  (initial value)
imcPositionMC                      = imcPosition           ! Position (cylinders)                  (initial value)
pQuaternionMC                      = pQuaternion           ! Quaternion algebra                    (initial value)
imQuaternionMC                     = imQuaternion          ! Quaternion algebra                    (initial value)
pPositionMC                        = pPosition             ! Position of particles                 (initial value)
imPositionMC                       = imPosition            ! Position of particles                 (initial value)
pOrientationMC                     = pOrientation          ! Orientation of particles              (initial value)
imOrientationMC                    = imOrientation         ! Orientation of particles              (initial value)
BoxLengthMC                        = BoxLength             ! Box length                            (initial value)
BoxLengthInverseMC                 = BoxLengthInverse      ! Box length (inverse)                  (initial value)
BoxVolumeMC                        = BoxVolume             ! Box volume                            (initial value)
Ensemble                           = "NPT"                 ! Ensemble type                         (initial value)

! Initialize system properties
PackingFraction    = ( DBLE( nParticles ) * pMolecularVolume ) / BoxVolume
TotalNumberDensity = ( DBLE( nParticles ) / BoxVolume )

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"FLOPPY-BOX SIMULATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Starting Floppy-Box Monte Carlo simulation..."
WRITE( *, "(G0)" ) " "

DO rPressure = 1, nPressure

  ! Status
  WRITE( *, "(5G0,G0.5)", Advance= "NO" ) "Current Pressure (", rPressure, "/", nPressure, "): ", ReducedPressure(rPressure)
  WRITE( *, "(G0,G0.5)", Advance= "NO" ) " | Packing Fraction: ", PackingFraction
  WRITE( *, "(G0,G0.5)", Advance= "YES" ) " | Number Density: ", TotalNumberDensity

  ! Simulation cycles
  DO iCycle = 1, MaxSimulationCycles

    ! Progress bar
    CALL ProgressBarMC( iCycle, MaxSimulationCycles, Ensemble )

    ! Generates a random number for the NPT-simulation
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

    ! Movement (Translation or Rotation)
    IF( RandomNumber <= MovementProbability ) THEN

      ! Disable volume movement
      MovementIsoVolumeChangeLogical   = .FALSE.
      MovementAnisoVolumeChangeLogical = .FALSE.

      ! Pseudorandom number generator (uniform distribution)
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      iParticle = INT( DBLE( nParticles ) * RandomNumber ) + 1

      ! Assignment of previous configuration (Microstate m)
      iOldPosition(:)        = pPositionMC(:,iParticle)       ! Position
      iimOldPosition(:,:)    = imPositionMC(:,:,iParticle)    ! Position (image)
      iOldQuaternion(:)      = pQuaternionMC(:,iParticle)     ! Quaternion
      iimOldQuaternion(:,:)  = imQuaternionMC(:,:,iParticle)  ! Quaternion (image)
      iOldOrientation(:)     = pOrientationMC(:,iParticle)    ! Orientation
      iimOldOrientation(:,:) = imOrientationMC(:,:,iParticle) ! Orientation (image)
      ciOldPosition(:,:)     = cPositionMC(:,:,iParticle)     ! Position (cylinders)
      iimcOldPosition(:,:,:) = imcPositionMC(:,:,:,iParticle) ! Position (cylindrical image)

      ! Disable translation if only one particle is present
      IF( nParticles /= 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        ! Translation criterion
        IF( RandomNumber <= TranslationalProbability ) THEN
          MovementTranslationLogical  = .TRUE.  ! Enable translation
          MovementRotationLogical     = .FALSE. ! Disable rotation
          nMovementTranslationCounter = nMovementTranslationCounter + 1
        ! Rotation criterion
        ELSE IF( RandomNumber > TranslationalProbability ) THEN
          MovementRotationLogical    = .TRUE.  ! Enable rotation
          MovementTranslationLogical = .FALSE. ! Disable translation
          nMovementRotationCounter   = nMovementRotationCounter + 1
        END IF
      ELSE
        MovementTranslationLogical = .FALSE. ! Disable translation
        MovementRotationLogical    = .TRUE.  ! Enable rotation
        nMovementRotationCounter   = nMovementRotationCounter + 1
      END IF

      ! Translation
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
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, iNewPosition )
        ! Construct images of particle i
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox )
        ! Initialize image list
        CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
        ! Create images
        DO pImage = 1, nImages
          CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceImageUnitBox(:,pImage), iimNewPosition(:,pImage) )
        END DO
      ! No translation
      ELSE IF( .NOT. MovementTranslationLogical ) THEN
        iNewPosition   = iOldPosition
        iimNewPosition = iimOldPosition
      END IF

      ! Rotation
      IF( MovementRotationLogical ) THEN
        ! Random Composed Unit Quaternion
        CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
        ! Active transformation (rotation)
        CALL ActiveTransformation( zAxis, iNewQuaternion, iNewOrientation )
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

      ! Random position of cylinders (after translation or rotation)
      DO cCylinder = 1, 4
        ! Active transformation (translation)
        CALL ActiveTransformation( cReferencePosition(:,cCylinder), iNewQuaternion, cRotatedPosition(:,cCylinder) )
        ciNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPosition(:,cCylinder)
        ! Construc images of cylinders of particle i
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, ciNewPosition(:,cCylinder), ScalingDistanceUnitBox )
        ! Initialize list of images
        CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
        ! Create images
        DO pImage = 1, nImages
          CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceImageUnitBox(:,pImage), &
          &                                iimcNewPosition(:,cCylinder,pImage) )
        END DO
      END DO

      ! Overlap check after displacement of a particle
      CALL ParticleOverlapCheck( iParticle, iNewQuaternion, iimNewQuaternion, iNewOrientation, iimNewOrientation, iNewPosition, &
      &                          iimNewPosition, ciNewPosition, iimcNewPosition, ContactDistance, Overlap )

      ! Acceptance criterion
      IF( .NOT. Overlap ) THEN
        ! Update system configuration
        pPositionMC(:,iParticle)       = iNewPosition(:)        ! Update position
        imPositionMC(:,:,iParticle)    = iimNewPosition(:,:)    ! Update position (image)
        cPositionMC(:,:,iParticle)     = ciNewPosition(:,:)     ! Update position of cylinders
        imcPositionMC(:,:,:,iParticle) = iimcNewPosition(:,:,:) ! Update position of cylinders (image)
        pQuaternionMC(:,iParticle)     = iNewQuaternion(:)      ! Update quaternion
        imQuaternionMC(:,:,iParticle)  = iimNewQuaternion(:,:)  ! Update quaternion (image)
        pOrientationMC(:,iParticle)    = iNewOrientation(:)     ! Update orientation
        imOrientationMC(:,:,iParticle) = iimNewOrientation(:,:) ! Update orientation (image)
        ! Update displacement counter
        IF( MovementTranslationLogical ) THEN
          nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
        ELSE IF( MovementRotationLogical ) THEN
          nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
        END IF
      ! Move rejected
      ELSE
        ! Retrieve old configuration
        pPositionMC(:,iParticle)       = iOldPosition(:)         ! Retrieve old position
        imPositionMC(:,:,iParticle)    = iimOldPosition(:,:)     ! Retrieve old position (images)
        cPositionMC(:,:,iParticle)     = ciOldPosition(:,:)      ! Retrieve old position of cylinders
        imcPositionMC(:,:,:,iParticle) = iimcOldPosition(:,:,:)  ! Retrieve old position of cylinders (images)
        pQuaternionMC(:,iParticle)     = iOldQuaternion(:)       ! Retrieve old quaternion
        imQuaternionMC(:,:,iParticle)  = iimOldQuaternion(:,:)   ! Retrieve old quaternion (images)
        pOrientationMC(:,iParticle)    = iOldOrientation(:)      ! Retrieve old orientation
        imOrientationMC(:,:,iParticle)  = iimOldOrientation(:,:) ! Retrieve old orientation (images)
      END IF

    ! Volume scaling
    ELSE IF( RandomNumber > MovementProbability ) THEN

      ! Disable translation and rotation
      MovementTranslationLogical = .FALSE.
      MovementRotationLogical    = .FALSE.

      ! Assignment of previous configuration (Microstate m)
      OldBoxLength        = BoxLengthMC        ! Box length
      OldBoxLengthInverse = BoxLengthInverseMC ! Box length (inverse)
      OldBoxVolume        = BoxVolumeMC        ! Box volume

      ! Expansion/compression type
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

      ! Isotropic volume change
      IF( RandomNumber < IsoVolumetricProbability .OR. AnisotropicVolumeChangeIndex >= rPressure ) THEN
        ! Random factor
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        VolumeScalingFactor = DLOG( OldBoxVolume ) + (RandomNumber - 0.5D0) * MaxIsoVolumetricDisplacement
        VolumeScalingFactor = DEXP( VolumeScalingFactor )
        VolumeScalingFactor = (VolumeScalingFactor / OldBoxVolume) ** (1.D0 / 3.D0)
        ! Proportional box length
        NewBoxLength = OldBoxLength * VolumeScalingFactor
        VolumeType   = "IS"
        ! Calculate the new reciprocal box basis vectors and the volume of the system
        CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
        ! Movement counter
        nMovementIsoVolumeChangeCounter = nMovementIsoVolumeChangeCounter + 1
        ! Movement type
        MovementIsoVolumeChangeLogical   = .TRUE.
        MovementAnisoVolumeChangeLogical = .FALSE.
      ! Anisotropic volume change
      ELSE IF( RandomNumber >= IsoVolumetricProbability .AND. AnisotropicVolumeChangeIndex < rPressure ) THEN
        ! Random component
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        BoxMatrixComponent = INT( RandomNumber * 6.D0 ) + 1
        IF( BoxMatrixComponent == 1 ) THEN
          BoxMatrixComponent = 1 ! xx component
        ELSE IF( BoxMatrixComponent == 2 ) THEN
          BoxMatrixComponent = 4 ! xy component
        ELSE IF( BoxMatrixComponent == 3 ) THEN
          BoxMatrixComponent = 5 ! yy component
        ELSE IF( BoxMatrixComponent == 4 ) THEN
          BoxMatrixComponent = 7 ! xz component
        ELSE IF( BoxMatrixComponent == 5 ) THEN
          BoxMatrixComponent = 8 ! yz component
        ELSE IF( BoxMatrixComponent == 6 ) THEN
          BoxMatrixComponent = 9 ! zz component
        END IF
        NewBoxLength = OldBoxLength
        ! Random factor
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        NewBoxLength(BoxMatrixComponent) = OldBoxLength(BoxMatrixComponent) + MaxAnisoVolumetricDisplacement * &
        &                                  (2.D0 * RandomNumber - 1.0D0)
        VolumeType = "AN"
        ! Calculate the new reciprocal box basis vectors and the volume of the system
        CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
        ! Movement counter
        nMovementAnisoVolumeChangeCounter = nMovementAnisoVolumeChangeCounter + 1
        ! Movement type
        MovementAnisoVolumeChangeLogical = .TRUE.
        MovementIsoVolumeChangeLogical   = .FALSE.
      END IF

      ! Reset condition of anisotropic volume change
      CheckBoxDistortion = .FALSE.

      ! Condition of anisotropic volume change
      IF( MovementAnisoVolumeChangeLogical ) THEN
        ! Box length
        BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(1:3) ) )
        BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(4:6) ) )
        BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( NewBoxLength(7:9), NewBoxLength(7:9) ) )
        ! Length ratio
        BoxEdgeRatio(1) = BoxEdgeLength(1) / BoxEdgeLength(2)
        BoxEdgeRatio(2) = BoxEdgeLength(1) / BoxEdgeLength(3)
        BoxEdgeRatio(3) = BoxEdgeLength(2) / BoxEdgeLength(3)
        ! Angle between box vectors
        BoxVectorAngle(1) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(4:6) ) / ( BoxEdgeLength(1) * BoxEdgeLength(2) )
        BoxVectorAngle(2) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(7:9) ) / ( BoxEdgeLength(1) * BoxEdgeLength(3) )
        BoxVectorAngle(3) = DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(7:9) ) / ( BoxEdgeLength(2) * BoxEdgeLength(3) )
        ! Avoid big distortions of the simulation box
        DO bEdge = 1, 3
          ! Angle distortion
          IF( BoxVectorAngle(bEdge) < DCOS( (cPi / 2.D0) + BoxVectorMaxAngle ) .OR. &
          &   BoxVectorAngle(bEdge) > DCOS( (cPi / 2.D0) - BoxVectorMaxAngle ) ) THEN
            BoxVolumeMC        = OldBoxVolume
            BoxLengthMC        = OldBoxLength
            BoxLengthInverseMC = OldBoxLengthInverse
            CheckBoxDistortion = .TRUE.
            EXIT
          END IF
          ! Length distortion
          IF( BoxEdgeRatio(bEdge) > BoxEdgeMaxRatio .OR. BoxEdgeRatio(bEdge) < 1.D0 / BoxEdgeMaxRatio ) THEN
            BoxVolumeMC        = OldBoxVolume
            BoxLengthMC        = OldBoxLength
            BoxLengthInverseMC = OldBoxLengthInverse
            CheckBoxDistortion = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! Box not too distorted
      IF( .NOT. CheckBoxDistortion ) THEN

        ! Enthalpy (weighing function)
        IF( MovementIsoVolumeChangeLogical )   EnthalpyChange = ( ReducedPressure(rPressure) * ( NewBoxVolume - OldBoxVolume ) ) &
        &                                                       - ( DBLE( nParticles + 1 ) * DLOG( NewBoxVolume / OldBoxVolume ) )
        IF( MovementAnisoVolumeChangeLogical ) EnthalpyChange = ( ReducedPressure(rPressure) * ( NewBoxVolume - OldBoxVolume ) ) &
        &                                                       - ( DBLE( nParticles ) * DLOG( NewBoxVolume / OldBoxVolume ) )

        ! Random number
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

        ! Enthalpy condition
        IF( DEXP( - EnthalpyChange ) >= RandomNumber ) THEN

          ! System configuration update
          pPositionSaveMC   = pPositionMC   ! Old configuration
          imPositionSaveMC  = imPositionMC  ! Old configuration
          cPositionSaveMC   = cPositionMC   ! Old configuration
          imcPositionSaveMC = imcPositionMC ! Old configuration

          ! Isotropic volume change
          IF( MovementIsoVolumeChangeLogical ) THEN
            ! Rescale positions of all particles and images accordingly
            DO pParticle = 1, nParticles
              ! Rescale positions of all particles accordingly
              pPositionMC(:,pParticle) = pPositionMC(:,pParticle) * VolumeScalingFactor
              ! Transform spatial coordinates into scaled coordinates (unit box)
              CALL MatrixVectorMultiplication( NewBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
              ! Build image list
              CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
              ! Rescale position of all images accordingly
              DO pImage = 1, nImages
                CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceImageUnitBox(:,pImage), &
                &                                imPositionMC(:,pImage,pParticle) )
              END DO
              ! Relative molecular displacement between microstates
              RescalingFactor(:) = pPositionMC(:,pParticle) - pPositionSaveMC(:,pParticle) ! For cylinders
              ! Rescale position of all cylinders accordingly
              DO cCylinder = 1, 4
                cPositionMC(:,cCylinder,pParticle) = cPositionMC(:,cCylinder,pParticle) + RescalingFactor(:)
                ! Transform spatial coordinates into scaled coordinates (unit box)
                CALL MatrixVectorMultiplication( NewBoxLengthInverse, cPositionMC(:,cCylinder,pParticle), ScalingDistanceUnitBox )
                ! Build image list
                CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
                ! Rescale position of all cylindrical images accordingly
                DO pImage = 1, nImages
                  CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceImageUnitBox(:,pImage), &
                  &                                imcPositionMC(:,cCylinder,pImage,pParticle) )
                END DO
              END DO
            END DO
          ! Anisotropic volume change
          ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
            ! Rescale positions of all particles and images accordingly
            DO pParticle = 1, nParticles
              ! Scaling coordinates using the old box length
              CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
              ! New real coordinates using the new box length
              CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceUnitBox, pPositionMC(:,pParticle) )
              ! Transform spatial coordinates into scaled coordinates (unit box)
              CALL MatrixVectorMultiplication( NewBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
              ! Build image list
              CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
              ! Rescale position of all images accordingly
              DO pImage = 1, nImages
                CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceImageUnitBox(:,pImage), &
                &                                imPositionMC(:,pImage,pParticle) )
              END DO
              ! Relative molecular displacement between microstates
              RescalingFactor(:) = pPositionMC(:,pParticle) - pPositionSaveMC(:,pParticle) ! For cylinders
              ! Rescale position of all cylinders accordingly
              DO cCylinder = 1, 4
                cPositionMC(:,cCylinder,pParticle) = cPositionMC(:,cCylinder,pParticle) + RescalingFactor(:)
                ! Transform spatial coordinates into scaled coordinates (unit box)
                CALL MatrixVectorMultiplication( NewBoxLengthInverse, cPositionMC(:,cCylinder,pParticle), ScalingDistanceUnitBox )
                ! Build image list
                CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
                ! Rescale position of all cylindrical images accordingly
                DO pImage = 1, nImages
                  CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceImageUnitBox(:,pImage), &
                  &                                imcPositionMC(:,cCylinder,pImage,pParticle) )
                END DO
              END DO
            END DO
          END IF

          ! Overlap check after expansion/compression of the simulation box
          CALL FullOverlapCheck( ContactDistance, Overlap )

          ! Acceptance criterion
          IF( .NOT. Overlap ) THEN
            ! System configuration update
            BoxVolumeMC        = NewBoxVolume        ! Update volume
            BoxLengthMC        = NewBoxLength        ! Update length
            BoxLengthInverseMC = NewBoxLengthInverse ! Update length (inverse)
            ! Displacement counter update
            IF( MovementIsoVolumeChangeLogical ) THEN
              nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Volumetric move counter
            ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
              nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Volumetric move counter
            END IF
            ! Update packing fraction and reduced number density
            PackingFraction    = ( DBLE( nParticles ) * pMolecularVolume ) / BoxVolumeMC
            TotalNumberDensity = ( DBLE( nParticles ) / BoxVolumeMC )
            ! Re-initialization
            CheckBoxDistortion = .FALSE.
            ! Lattice reduction
            LatticeReductionLogical = .FALSE.
            CALL LatticeReduction( BoxLengthMC, BoxDistortionMC, LatticeReductionLogical )
            IF( LatticeReductionLogical ) THEN
              ! Calculate the new reciprocal box basis vectors
              CALL InverseMatrixCofactorVec( BoxLengthMC, BoxLengthInverseMC, BoxVolumeMC )
              DO pParticle = 1, nParticles
                ! Minimum image convention (the spatial distribution of particles remains unchanged)
                CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pParticle), ScalingDistanceUnitBox ) ! Spatial transformation
                ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
                CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, pPositionMC(:,pParticle) ) ! Spatial transformation
                ! Transform spatial coordinates into scaled coordinates (unit box)
                CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
                ! Build image list of particles
                CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
                ! Position of images of a particle
                DO pImage = 1, nImages
                  CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceImageUnitBox(:,pImage), &
                  &                                imPositionMC(:,pImage,pParticle) )
                END DO
                ! Cylinders
                DO cCylinder = 1, 4
                  ! Active transformation (translation)
                  CALL ActiveTransformation( cReferencePosition(:,cCylinder), pQuaternionMC(:,pParticle), &
                  &                          cRotatedPosition(:,cCylinder) )
                  ! Position of cylinders of a particle
                  cPositionMC(:,cCylinder,pParticle) = pPositionMC(:,pParticle) + cRotatedPosition(:,cCylinder)
                  ! Transform spatial coordinates into scaled coordinates (unit box)
                  CALL MatrixVectorMultiplication( BoxLengthInverseMC, cPositionMC(:,cCylinder,pParticle), ScalingDistanceUnitBox )
                  ! Build image list of cylinders
                  CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
                  ! Position of images of a cylinder
                  DO pImage = 1, nImages
                    CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceImageUnitBox(:,pImage), &
                    &                                imcPositionMC(:,cCylinder,pImage,pParticle) )
                  END DO
                END DO
              END DO
              ! Check orientation of the box (eliminate box rotations)
              IF( DABS( BoxLengthMC(2) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. DABS( BoxLengthMC(3) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. &
              &   DABS( BoxLengthMC(6) - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
                ! Undo box rotation
                CALL UndoBoxRotation( BoxLengthMC, BoxLengthInverseMC, BoxVolumeMC )
              END IF
            END IF
          ! Retrieve old properties of the system configuration
          ELSE
            BoxVolumeMC        = OldBoxVolume        ! Retrieve old box volume
            BoxLengthMC        = OldBoxLength        ! Retrieve old box length
            BoxLengthInverseMC = OldBoxLengthInverse ! Retrieve old box length (inverse)
            pPositionMC        = pPositionSaveMC     ! Retrieve old position of particles
            imPositionMC       = imPositionSaveMC    ! Retrieve old position of images of particles
            cPositionMC        = cPositionSaveMC     ! Retrieve old position of cylinders
            imcPositionMC      = imcPositionSaveMC   ! Retrieve old position of cylindrical images
          END IF

        ! Retrieve old properties of the simulation box
        ELSE

          BoxVolumeMC        = OldBoxVolume        ! Retrieve old box volume
          BoxLengthMC        = OldBoxLength        ! Retrieve old box length
          BoxLengthInverseMC = OldBoxLengthInverse ! Retrieve old box length (inverse)

        END IF ! Enthalpy criterion

      END IF ! Box distortion criterion

    END IF

    ! Adjustment of maximum displacements
    IF( iCycle <= nEquilibrationCycles ) THEN ! During equilibration only

      ! Acceptance ratio (translation)
      IF( MOD( iCycle, nAdjustmentMovementFrequency ) == 0 .AND. nMovementTranslationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          MaxTranslationalDisplacement  = 0.95D0 * MaxTranslationalDisplacement
        ELSE
          MaxTranslationalDisplacement  = 1.05D0 * MaxTranslationalDisplacement
        END IF
        ! Ratio data (volume scaling)
        WRITE( 30, "(7G0)" ) iCycle, ",", Ratio, ",", MaxTranslationalDisplacement, ",", AcceptanceRatioTranslation
        FLUSH( 30 )
        ! Reset counter
        nAcceptanceTranslation      = 0
        nMovementTranslationCounter = 0
        ! Avoid large translations
        BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( BoxLengthMC(1:3), BoxLengthMC(1:3) ) )
        BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( BoxLengthMC(4:6), BoxLengthMC(4:6) ) )
        BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( BoxLengthMC(7:9), BoxLengthMC(7:9) ) )
        IF( MaxTranslationalDisplacement > 2.D0 * MAXVAL( BoxEdgeLength ) ) THEN
          MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxEdgeLength )
        END IF
      END IF

      ! Adjustment of maximum displacement (rotation)
      IF( MOD( iCycle, nAdjustmentMovementFrequency ) == 0 .AND. nMovementRotationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotation adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
        ELSE
          MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
        END IF
        ! Ratio data
        WRITE( 40, "(7G0)" ) iCycle, ",", Ratio, ",", MaxAngularDisplacement, ",", AcceptanceRatioRotation
        FLUSH( 40 )
        ! Reset counter
        nAcceptanceRotation      = 0
        nMovementRotationCounter = 0
        ! Avoid multiple turns
        IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
          MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
        END IF
      END IF

      ! Adjustment of maximum displacement (isotropic volume scaling)
      IF( MOD( iCycle, nAdjustmentVolumeFrequency ) == 0 .AND. nMovementIsoVolumeChangeCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceIsotropicVolumeChange ) / DBLE( nMovementIsoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioIsoVolumeChange ) THEN
          MaxIsoVolumetricDisplacement = 0.95D0 * MaxIsoVolumetricDisplacement
        ELSE
          MaxIsoVolumetricDisplacement = 1.05D0 * MaxIsoVolumetricDisplacement
        END IF
        ! Ratio data
        WRITE( 50, "(9G0)" ) iCycle, ",", Ratio, ",", MaxIsoVolumetricDisplacement, ",", AcceptanceRatioIsoVolumeChange, &
        &                    ",", VolumeType
        FLUSH( 50 )
        ! Reset counter
        nAcceptanceIsotropicVolumeChange = 0
        nMovementIsoVolumeChangeCounter  = 0
      END IF

      ! Adjustment of maximum displacement (anisotropic volume scaling)
      IF( MOD( iCycle, nAdjustmentVolumeFrequency ) == 0 .AND. nMovementAnisoVolumeChangeCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceAnisotropicVolumeChange ) / DBLE( nMovementAnisoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioAnisoVolumeChange ) THEN
          MaxAnisoVolumetricDisplacement = 0.95D0 * MaxAnisoVolumetricDisplacement
        ELSE
          MaxAnisoVolumetricDisplacement = 1.05D0 * MaxAnisoVolumetricDisplacement
        END IF
        ! Ratio data
        WRITE( 50, "(9G0)" ) iCycle, ",", Ratio, ",", MaxAnisoVolumetricDisplacement, ",", AcceptanceRatioAnisoVolumeChange, &
        &                    ",", VolumeType
        FLUSH( 50 )
        ! Reset counter
        nAcceptanceAnisotropicVolumeChange = 0
        nMovementAnisoVolumeChangeCounter  = 0
      END IF

    END IF

    ! Trajectory data
    IF( TrajectoryLogical ) THEN
      IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
        WRITE( 20, "(G0)" ) nParticles * 4 + nParticles * nImages * 4
        ! Descriptor string
        DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
        WRITE( 20, DescriptorString ) 'Lattice="', (1.D0 + 2.D0 * DBLE( nLayers) ) * BoxLengthMC(1:9), '" Origin="', -0.5D0 * &
        &                             (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLengthMC(1) + BoxLengthMC(4) + BoxLengthMC(7) ), &
        &                             -0.5D0 * (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLengthMC(2) + BoxLengthMC(5) + &
        &                             BoxLengthMC(8) ), -0.5D0 * (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLengthMC(3) + &
        &                             BoxLengthMC(6) + BoxLengthMC(9) ), '" ', &
        &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
        DO pParticle = 1, nParticles
          ! Position of cylinders
          DO cCylinder = 1, 4
            TrajectoryPosition(1) = cPositionMC(1,cCylinder,pParticle)
            TrajectoryPosition(2) = cPositionMC(2,cCylinder,pParticle)
            TrajectoryPosition(3) = cPositionMC(3,cCylinder,pParticle)
            WRITE( 20, "(11(G0,1X))" ) Atom(1), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), &
            &                          pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), &
            &                          pQuaternionMC(0,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
            DO pImage = 1, nImages
              TrajectoryPosition(1) = imcPositionMC(1,cCylinder,pImage,pParticle)
              TrajectoryPosition(2) = imcPositionMC(2,cCylinder,pImage,pParticle)
              TrajectoryPosition(3) = imcPositionMC(3,cCylinder,pImage,pParticle)
              WRITE( 20, "(11(G0,1X))" ) Atom(2), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), &
              &                          imQuaternionMC(1,pImage,pParticle), imQuaternionMC(2,pImage,pParticle), &
              &                          imQuaternionMC(3,pImage,pParticle), imQuaternionMC(0,pImage,pParticle), &
              &                          0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
            END DO
          END DO
        END DO
        FLUSH( 20 )
      END IF
    END IF

    ! Results data (packing fraction, number density, box volume, and pressure)
    IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
      WRITE( 60, "(9G0)" ) iCycle, ",", PackingFraction, ",", TotalNumberDensity, ",", BoxVolumeMC, ",", ReducedPressure(rPressure)
      FLUSH( 60 )
    END IF

    ! Order parameter data
    IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
      ! Nematic order parameter (Q-tensor method)
      CALL UniaxialNematicOrderParameter( OrderParameterMC, pOrientationMC )
      WRITE( 80, "(3G0)" ) iCycle, ",", OrderParameterMC
      FLUSH( 80 )
    END IF

    ! Floppy-box data
    IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
      REWIND( 70 )
      WRITE( 70, "(2G0)" ) "Arrangement_type: ", Arrangement
      WRITE( 70, "(2G0)" ) "Cylindrical_diameter: ", cDiameter
      WRITE( 70, "(2G0)" ) "Cylindrical_length: ", cLength
      WRITE( 70, "(2G0)" ) "Cylindrical_aspect_ratio: ", cAspectRatio
      WRITE( 70, "(18G0)" ) "Box_length: ", BoxLengthMC(1), ", ", BoxLengthMC(2), ", ", BoxLengthMC(3), ", ", BoxLengthMC(4), &
      &                     ", ", BoxLengthMC(5), ", ", BoxLengthMC(6), ", ", BoxLengthMC(7), ", ", BoxLengthMC(8), ", ", &
      &                     BoxLengthMC(9)
      WRITE( 70, "(2G0)" ) "Box_volume: ", BoxVolumeMC
      WRITE( 70, "(2G0)" ) "Packing_fraction: ", PackingFraction
      WRITE( 70, "(2G0)" ) "Number_density: ", TotalNumberDensity
      WRITE( 70, "(2G0)" ) "Number_of_particles: ", nParticles
      WRITE( 70, "(2G0)" ) "Nematic_order_parameter: ", OrderParameterMC
      WRITE( 70, "(G0)" ) " "
      WRITE( 70, "(G0)" ) "Positions_and_orientations: "
      DO pParticle = 1, nParticles
          TrajectoryPosition(1) = pPositionMC(1,pParticle)
          TrajectoryPosition(2) = pPositionMC(2,pParticle)
          TrajectoryPosition(3) = pPositionMC(3,pParticle)
          WRITE( 70, "(21G0)" ) "R", ", ", TrajectoryPosition(1), ", ", TrajectoryPosition(2), ", ", TrajectoryPosition(3), ", ", &
          &                     pQuaternionMC(0,pParticle), ", ", pQuaternionMC(1,pParticle), ", ", pQuaternionMC(2,pParticle), &
          &                     ", ", pQuaternionMC(3,pParticle)
      END DO
      FLUSH( 70 )
    END IF

  END DO

  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "

END DO

! Final calculations
fBoxLength       = BoxLengthMC
fBoxVolume       = BoxVolumeMC
fPackingFraction = PackingFraction
fNumberDensity   = TotalNumberDensity

! Nematic order parameter (Q-tensor method)
CALL UniaxialNematicOrderParameter( OrderParameterMC, pOrientationMC )

! Last box parameters
REWIND( 70 )
WRITE( 70, "(2G0)" ) "Arrangement_type: ", Arrangement
WRITE( 70, "(2G0)" ) "Cylindrical_diameter: ", cDiameter
WRITE( 70, "(2G0)" ) "Cylindrical_length: ", cLength
WRITE( 70, "(2G0)" ) "Cylindrical_aspect_ratio: ", cAspectRatio
WRITE( 70, "(18G0)" ) "Box_length: ", fBoxLength(1), ", ", fBoxLength(2), ", ", fBoxLength(3), ", ", fBoxLength(4), &
&                     ", ", fBoxLength(5), ", ", fBoxLength(6), ", ", fBoxLength(7), ", ", fBoxLength(8), ", ", &
&                     fBoxLength(9)
WRITE( 70, "(2G0)" ) "Box_volume: ", fBoxVolume
WRITE( 70, "(2G0)" ) "Packing_fraction: ", fPackingFraction
WRITE( 70, "(2G0)" ) "Number_density: ", fNumberDensity
WRITE( 70, "(2G0)" ) "Number_of_particles: ", nParticles
WRITE( 70, "(2G0)" ) "Nematic_order_parameter: ", OrderParameterMC
WRITE( 70, "(G0)" ) " "
WRITE( 70, "(G0)" ) "Positions_and_orientations: "
DO pParticle = 1, nParticles
    TrajectoryPosition(1) = pPositionMC(1,pParticle)
    TrajectoryPosition(2) = pPositionMC(2,pParticle)
    TrajectoryPosition(3) = pPositionMC(3,pParticle)
    WRITE( 70, "(21G0)" ) "R", ", ", TrajectoryPosition(1), ", ", TrajectoryPosition(2), ", ", TrajectoryPosition(3), ", ", &
    &                     pQuaternionMC(0,pParticle), ", ", pQuaternionMC(1,pParticle), ", ", pQuaternionMC(2,pParticle), &
    &                     ", ", pQuaternionMC(3,pParticle)
END DO
FLUSH( 70 )

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"FINAL PROPERTIES"//REPEAT( " ", 20 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,9(G0.5,1X),G0)" ) "Final box length: [ ", fBoxLength, "]"
WRITE( *, "(G0,G0.5)" ) "Final box volume: ", fBoxVolume
WRITE( *, "(G0,G0.5)" ) "Final packing fraction: ", fPackingFraction
WRITE( *, "(G0,G0.5)" ) "Final number density: ", fNumberDensity
WRITE( *, "(G0,G0.5)" ) "Final nematic order parameter: ", OrderParameterMC
WRITE( *, "(G0)" ) " "

! Output units
IF( TrajectoryLogical ) THEN
  CLOSE( 20 )
END IF
CLOSE( 30 )
CLOSE( 40 )
CLOSE( 50 )
CLOSE( 60 )
CLOSE( 70 )

! Finish simulation timer
CALL CPU_Time( StopTimer )

! Deallocation
DEALLOCATE( pQuaternionMC, imQuaternionMC, pPositionMC, imPositionMC, pPositionSaveMC, imPositionSaveMC, pOrientationMC )
DEALLOCATE( imOrientationMC, cPositionMC, imcPositionMC, cPositionSaveMC, imcPositionSaveMC )
DEALLOCATE( pQuaternion, imQuaternion, pPosition, imPosition, pOrientation, imOrientation, cPosition, imcPosition )

! Execution time
WRITE( *, "(G0,G0.5,G0)" ) "Execution time: ", StopTimer - StartTimer, "s"
WRITE( *, "(G0)" ) " "

END PROGRAM FloppyBoxMonteCarlo