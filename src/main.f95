! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!     The Perram-Wertheim Method (J. W. Perram & M. S. Wertheim, 1985) is used to search for      !
!             molecular overlaps between ellipsoids of revolution after a trial move.             !
!              The Vega-Lago Method (C. Vega & S. Lago, 1993) is used to search for               !
!                 molecular overlaps between spherocylinders after a trial move.                  !
!  The algorithm of Lopes et al. (J. Lopes, F. Romano, E. Grelet, L. Franco, A. Giacometti, 2021) !
!         is used to search for molecular overlaps between cylinders after a trial move.          !
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
! Main References:                 J. W. Perram, M. S. Wertheim                                   !
!                               J. Comput. Phys 15, 409-416 (1985)                                !
!                                DOI: 10.1016/0021-9991(85)90171-8                                !
!                             --------------------------------------                              !
!                                   J. W. Perram, J. Rasmussen                                    !
!                                  Phys. Rev. E 54, 6565 (1996)                                   !
!                                 DOI: 10.1103/PhysRevE.54.6565                                   !
!                             --------------------------------------                              !
!                                        C. Vega, S. Lago                                         !
!                                 Computers Chem. 18, 55-59 (1993)                                !
!                                DOI: 10.1016/0097-8485(94)80023-5                                !
!                             --------------------------------------                              !
!                    J. Lopes, F. Romano, E. Grelet, L. Franco, A. Giacometti                     !
!                                 Chem. Phys. 154, 104902 (2021)                                  !
!                                     DOI: 10.1063/5.0040942                                      !
!                             --------------------------------------                              !
!                                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                      N. Metropolis et al.                                       !
!                                 J. Chem. Phys. 21, 1087 (1953)                                  !
!                                     DOI: 10.1063/1.1699114                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
PROGRAM Main

! Used seven modules: global variables, variable initialization, initial configuration, overlap check, potential energy, linked lists, and directory creator
USE GlobalVar
USE InitializeVariables
USE InitialSystemConfiguration, ONLY: GeometrySelection, InitialConfigurationSelection, InitialConfigurationStructure
USE OverlapCheck
USE ComputePotential
USE LinkedLists, ONLY: FinalizeList, FinalizeListPotential, MakeList, MakeListPotential, BoxCheckNPT, ParticleTranslationNVT
USE Folders, ONLY: InquireFolders, DateFolders, RangeFolders, BackupFolder

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER                 :: SeedSize ! Seed array size
INTEGER                 :: iSeed    ! Seed array component
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: RandomSeed ! Random seed array

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: iParticle, pParticle               ! Counter (particles)
INTEGER( Kind= Int64 )                 :: bEdge                              ! Counter (box edges)
INTEGER( Kind= Int64 )                 :: pCycle                             ! Counter of cycles
INTEGER( Kind= Int64 )                 :: iComponent, cComponent             ! Counter (component)
INTEGER( Kind= Int64 )                 :: rRange                             ! Counter (potential range)
INTEGER( Kind= Int64 )                 :: FirstCycle                         ! Initial cycle
INTEGER( Kind= Int64 )                 :: BoxMatrixComponent                 ! Box matrix component
INTEGER( Kind= Int64 )                 :: iCycle                             ! Counter (simulation cycles)
INTEGER( Kind= Int64 )                 :: nAcceptanceTranslation             ! Move acceptance counter: Translation
INTEGER( Kind= Int64 )                 :: nAcceptanceRotation                ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 )                 :: nAcceptanceIsotropicVolumeChange   ! Move acceptance counter: Volume scaling
INTEGER( Kind= Int64 )                 :: nAcceptanceAnisotropicVolumeChange ! Move acceptance counter: Volume scaling
INTEGER( Kind= Int64 )                 :: nMovementTranslationCounter        ! Move counter (Translation)
INTEGER( Kind= Int64 )                 :: nMovementRotationCounter           ! Move counter (Rotation)
INTEGER( Kind= Int64 )                 :: nMovementIsoVolumeChangeCounter    ! Move counter (Volume scaling)
INTEGER( Kind= Int64 )                 :: nMovementAnisoVolumeChangeCounter  ! Move counter (Volume scaling)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: LastLine                           ! Last line of a file

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: StartTimer                               ! Start timer of Monte Carlo simulation
REAL( Kind= Real64 )                   :: StopTimer                                ! Stop timer of Monte Carlo simulation
REAL( Kind= Real64 )                   :: ContactDistance                          ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: Ratio                                    ! Acceptance ratio (simulation)
REAL( Kind= Real64 )                   :: OrderParameterMC                         ! Nematic order parameter
REAL( Kind= Real64 )                   :: OldBoxVolume, NewBoxVolume               ! Volume of simulation box (before/after a trial move)
REAL( Kind= Real64 )                   :: EnthalpyChange                           ! Enthalpy criterion (reduced)
REAL( Kind= Real64 )                   :: VolumeScalingFactor                      ! Scaling factor
REAL( Kind= Real64 )                   :: BoxDistortionMC                          ! Box distortion
REAL( Kind= Real64 )                   :: ExecutionTime                            ! Execution time
REAL( Kind= Real64 )                   :: BoxVolumeMC                              ! Reduced volume of the simulation box
REAL( Kind= Real64 )                   :: MaxTranslationalDisplacement             ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                   :: MaxAngularDisplacement                   ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                   :: MaxIsoVolumetricDisplacement             ! Maximum displacement [+/-] (Isotropic volume scaling)
REAL( Kind= Real64 )                   :: MaxAnisoVolumetricDisplacement           ! Maximum displacement [+/-] (Anisotropic volume scaling)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxCutoff                                ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox                   ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOldOrientation, iNewOrientation         ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOldPosition, iNewPosition               ! Position (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: OldBoxLength, NewBoxLength               ! Length of simulation box (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: OldBoxLengthInverse, NewBoxLengthInverse ! Length of simulation box (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthMC                              ! Length (x,y,z) of the simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthInverseMC                       ! Length (x,y,z) of simulation box (inverse)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxVectorAngle                           ! Cossine of angle between box vectors
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxEdgeLength                            ! Length of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxEdgeRatio                             ! Length ratio of box edges
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iOldQuaternion, iNewQuaternion           ! Quaternion (before/after a trial move)

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: TotalPotentialEnergyMC                   ! Potential energy array
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: iOldPotentialEnergy, iNewPotentialEnergy ! Perturbed potential energy (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: iPotentialEnergyDifference               ! Perturbed potential energy difference between two microstates m and n
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: CoefficientTPT1, CoefficientTPT2         ! First- and second-order TPT coefficients
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: rFreeEnergy                              ! Full Helmholtz free energy of the perturbed system
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: CoefficientTPT1Deviation                 ! First-order TPT coefficient (standard deviation)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: CoefficientTPT2Deviation                 ! Second-order TPT coefficient (standard deviation)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: rFreeEnergyDeviation                     ! Full Helmholtz free energy of the perturbed system (standard deviation)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: PositionSaveMC                           ! Old position of particles

! *********************************************************************************************** !
! CHARACTER STRINGS (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
CHARACTER( LEN= 02 ) :: VolumeType             ! Expansion/Compression type
CHARACTER( LEN= 32 ) :: FormatHour             ! String format for output file
CHARACTER( LEN= 32 ) :: FormatDate             ! String format for output folder
CHARACTER( LEN= 14 ) :: DescriptorBackupFile   ! Descriptor for output file
CHARACTER( LEN= 20 ) :: DescriptorBackupString ! Descriptor for backup variable

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap                          ! Detects overlap between two particles (PW) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: MovementRotationLogical          ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical       ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementIsoVolumeChangeLogical   ! Volumetric movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementAnisoVolumeChangeLogical ! Volumetric movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: CheckBoxDistortion               ! Detects if a box deformation is valid or not : TRUE = ignore box deformation; FALSE = consider box deformation
LOGICAL :: LatticeReductionLogical          ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: FlagLogical                      ! Generic TRUE/FALSE flag
LOGICAL :: FileExist                        ! Checks whether a file exists or not
LOGICAL :: FolderExist                      ! Checks whether a folder exists or not

! Title
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 6 )//"NVT/NPT-MONTE CARLO SIMULATION OF MIXTURES"//REPEAT( " ", 7 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR

! Restoration from backup file
OPEN( Unit= 10, File= "ini_control.ini", Action= "READ" )
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, RestoreBackupFileInquiry
CALL ToUpper( RestoreBackupFileInquiry, LEN_TRIM( RestoreBackupFileInquiry ), RestoreBackupFileInquiry )
! Transforms characters into logical variables
IF( RestoreBackupFileInquiry == "Y" ) THEN
  RestoreBackupFileLogical = .TRUE.
  WRITE( *, "(G0)" ) "User has chosen to restore a backup file from a previous simulation. Type the 14-digit backup descriptor: "
  READ( *, * ) DescriptorBackupFile
  WRITE( *, "(G0)" ) " "
  INQUIRE( File= "Backup/", Exist= FolderExist )
  IF( .NOT. FolderExist ) THEN
    WRITE( *, "(G0)" ) "Backup folder not found! Exiting..."
    CALL Sleep( 1 )
    CALL Exit(  )
  END IF
  INQUIRE( File= "Backup/"//TRIM( DescriptorBackupFile )//"_simulation.backup", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    INQUIRE( File= "Backup/"//TRIM( DescriptorBackupFile )//"_variables.backup", Exist= FileExist )
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0)" ) "Variable and simulation backup files not found! Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    ELSE
      WRITE( *, "(G0)" ) "Simulation backup file not found! Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
  ELSE
    INQUIRE( File= "Backup/"//TRIM( DescriptorBackupFile )//"_variables.backup", Exist= FileExist )
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0)" ) "Variable backup file not found! Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    ELSE
      WRITE( *, "(G0)" ) "Backup files found! Restoring previous simulation..."
      CALL Sleep( 1 )
      WRITE( *, "(G0)" ) " "
    END IF
  END IF
ELSE
  RestoreBackupFileLogical = .FALSE.
END IF
CLOSE( 10 )

! Initialize variables if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  ! Molecular geometry selection (see 'InitialConfiguration' module)
  CALL GeometrySelection(  )
  ! Molecular configuration selection (see 'InitialConfiguration' module)
  CALL InitialConfigurationSelection(  )
  ! Initialization of Monte Carlo parameters (see 'InitializeVariables' module)
  CALL MonteCarloVariables(  )
  ! Initialization of Inquiry/Control variables (see 'InitializeVariables' module)
  CALL VariablesInquest(  )
  ! Initialization of common variables (see 'InitializeVariables' module)
  CALL GeneralVariables(  )
  ! Initialization of potential variables (see 'InitializeVariables' module)
  CALL ForceFieldVariables(  )
  ! Pseudorandom number generator seed
  IF( FixedSeedLogical ) THEN
    SeedValue   = 123456789
    InitialSeed = SeedValue
  ELSE IF( .NOT. FixedSeedLogical ) THEN
    CALL RANDOM_SEED( Size= SeedSize )
    ALLOCATE( RandomSeed(SeedSize) )
    CALL RANDOM_SEED( Get= RandomSeed )
    CALL RANDOM_NUMBER( RandomNumber )
    iSeed       = INT( RandomNumber * DBLE( SeedSize ) ) + 1
    SeedValue   = ABS( RandomSeed(iSeed) )
    InitialSeed = SeedValue
  END IF
! Initialize variables from backup file
ELSE IF( RestoreBackupFileLogical ) THEN
  ! Summary
  WRITE( *, "(G0)", Advance= "NO" ) "Restoring common variables from a previous run..."
  ! Restore backup variables
  CALL RestoreBackupVariables( DescriptorBackupFile, DateTime )
  CALL Sleep( 1 )
  ! Allocation
  IF( PotentialTypeLogical(2) ) THEN
    ALLOCATE( CoefficientTPT1(nRange), CoefficientTPT2(nRange), rFreeEnergy(nRange) )
    ALLOCATE( CoefficientTPT1Deviation(nRange), CoefficientTPT2Deviation(nRange), rFreeEnergyDeviation(nRange) )
    ALLOCATE( TotalPotentialEnergy(nRange), TotalPotentialEnergyMC(nRange) )
    ALLOCATE( iNewPotentialEnergy(nRange), iOldPotentialEnergy(nRange), iPotentialEnergyDifference(nRange) )
    ALLOCATE( cPotentialRange(nComponents,nRange) )
  END IF
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Allocation
IF( .NOT. RestoreBackupFileLogical ) THEN
  ALLOCATE( pQuaternion(0:3,nParticles), pQuaternionMC(0:3,nParticles) )
  ALLOCATE( pPosition(3,nParticles), pPositionMC(3,nParticles) )
  ALLOCATE( pOrientation(3,nParticles), pOrientationMC(3,nParticles) )
  IF( PotentialTypeLogical(2) ) THEN
    ALLOCATE( CoefficientTPT1(nRange), CoefficientTPT2(nRange), rFreeEnergy(nRange) )
    ALLOCATE( CoefficientTPT1Deviation(nRange), CoefficientTPT2Deviation(nRange), rFreeEnergyDeviation(nRange) )
    ALLOCATE( TotalPotentialEnergy(nRange), TotalPotentialEnergyMC(nRange) )
    ALLOCATE( iNewPotentialEnergy(nRange), iOldPotentialEnergy(nRange), iPotentialEnergyDifference(nRange) )
    ALLOCATE( cPotentialRange(nComponents,nRange) )
  END IF
END IF

! Diameter of circumscribing sphere
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      IF( cAspectRatio(cComponent) > 0.D0 .AND. cAspectRatio(cComponent) < 1.D0 ) THEN
        cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
      ELSE IF( cAspectRatio(cComponent) > 1.D0 ) THEN
        cCircumscribingSphereDiameter(cComponent) = cLength(cComponent)
      END IF
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent) + cLength(cComponent)
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent) + cLength(cComponent)
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
END IF

! Hard-core volumetric relation (EOR/SPC/HC and SPHERES)
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ) ! Ellipsoids of revolution
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( 1.D0 + ( 1.5D0 * cAspectRatio(cComponent) ) ) ** &
      &                                       ( 1.D0 / 3.D0 ) ) ! Spherocylinders
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ) ! Cylinders
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
END IF

! Effective range of attraction
IF( PotentialTypeLogical(2) ) THEN
  DO cComponent = 1, nComponents
    cPotentialRange(cComponent,:) = PotentialRange(:) * cDiameterEquivalentSphere(cComponent)
  END DO
END IF

! Diameter of the largest circumscribing sphere
cLargestSphereDiameter = 0.D0
IF( CellListLogical ) THEN
  ! Geometrical reference
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        IF( cAspectRatio(cComponent) > 0.D0 .AND. cAspectRatio(cComponent) < 1.D0 ) THEN
          IF( cDiameter(cComponent) >= cLargestSphereDiameter ) cLargestSphereDiameter = cDiameter(cComponent)
        ELSE IF( cAspectRatio(cComponent) > 1.D0 ) THEN
          IF( cLength(cComponent) >= cLargestSphereDiameter ) cLargestSphereDiameter = cLength(cComponent)
        END IF
      ELSE
        IF( cDiameter(cComponent) >= cLargestSphereDiameter ) cLargestSphereDiameter = cDiameter(cComponent)
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        IF( ( cDiameter(cComponent) + cLength(cComponent) ) >= cLargestSphereDiameter ) cLargestSphereDiameter = &
        &     cDiameter(cComponent) + cLength(cComponent)
      ELSE
        IF( cDiameter(cComponent) >= cLargestSphereDiameter ) cLargestSphereDiameter = cDiameter(cComponent)
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        IF( ( cDiameter(cComponent) + cLength(cComponent) ) >= cLargestSphereDiameter ) cLargestSphereDiameter = &
        &     cDiameter(cComponent) + cLength(cComponent)
      ELSE
        IF( cDiameter(cComponent) >= cLargestSphereDiameter ) cLargestSphereDiameter = cDiameter(cComponent)
      END IF
    END DO
  END IF
END IF

! Diameter of the largest spherical well
cLargestSphericalWell = 0.D0
IF( CellListLogical ) THEN
  ! Potential reference
  IF( PotentialTypeLogical(2) ) THEN
    DO cComponent = 1, nComponents
      DO rRange = 1, nRange
        IF( cPotentialRange(cComponent,rRange) >= cLargestSphericalWell ) cLargestSphericalWell = &
        &   cPotentialRange(cComponent,rRange)
      END DO
    END DO
  END IF
END IF

! Atom ID (Required in some visualization and analysis software)
ALLOCATE( cIndex(nComponents) )
DO cComponent = 1, nComponents
  cIndex(cComponent) = cComponent
END DO

! Initialize descriptors if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  ! CPU Clock
  CALL DATE_AND_TIME( VALUES= DateTime )
  ! Date format (YYYY/MM/DD)
  FormatDate = "(I4,2I2.2)"
  WRITE( DescriptorDate, FormatDate ) DateTime(1), DateTime(2), DateTime(3)
  ! Time format (HH:MM:SS)
  FormatHour = "(3I2.2)"
  ! Hour descriptor
  WRITE( DescriptorHour, FormatHour ) DateTime(5), DateTime(6), DateTime(7)
  ! Output file descriptors (p. fraction/pressure [1], # of components [2], and geometry [3])
  FormatFileThermoVariable = "(F0.5)"
  IF( EnsembleMC == "NVT" ) THEN
    WRITE( DescriptorFileThermoVariable, FormatFileThermoVariable ) PackingFraction
  ELSE IF( EnsembleMC == "NPT" ) THEN
    WRITE( DescriptorFileThermoVariable, FormatFileThermoVariable ) ReducedPressure
  END IF
  FormatFileComponents = "(I0.3)"
  WRITE( DescriptorFileComponents, FormatFileComponents ) nComponents
  FormatFileGeometry = "(A3)"
  WRITE( DescriptorFileGeometry, FormatFileGeometry ) GeometryAcronym
END IF

! Initialize initial structure if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  ! Initial configuration (see 'InitialConfiguration' module)
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
  WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"INITIAL CONFIGURATION"//REPEAT( " ", 17 )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
  WRITE( *, "(G0)" ) "Setting up initial configuration folder..."
  CALL Sleep( 1 )
  WRITE( *, "(G0)" ) " "
  ! Initial configuration structure
  CALL InitialConfigurationStructure(  )
END IF

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"FOLDER ORGANIZER"//REPEAT( " ", 20 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Setting up folders..."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Create output directories (see 'Folders' module)
CALL InquireFolders(  )

! Create date subfolders (see 'Folders' module)
CALL DateFolders(  )

! Create backup folder (see 'Folders' module)
IF( BackupFileLogical ) THEN
  CALL BackupFolder(  )
  ! Date format (YYYY/MM/DD)
  FormatDate = "(I4,2I2.2)"
  ! Time format (HH:MM:SS)
  FormatHour = "(3I2.2)"
  ! Output file descriptors (p. fraction/pressure [1], # of components [2], and geometry [3])
  FormatFileThermoVariable = "(F0.5)"
  FormatFileComponents     = "(I0.3)"
  FormatFileGeometry       = "(A3)"
END IF

! Initialization of the attractive range subfolders (see 'Folders' module)
IF( PotentialTypeLogical(2) ) THEN
  CALL RangeFolders(  )
END IF

! Re-initialize cell arrays
IF( CellListLogical .AND. .NOT. RestoreBackupFileLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLength(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLength(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLength(9)
  CALL MakeList( BoxCutoff, pPosition, BoxLengthInverse )
  IF( PotentialTypeLogical(2) ) THEN
    BoxCutoff(1) = cLargestSphericalWell / BoxLength(1)
    BoxCutoff(2) = cLargestSphericalWell / BoxLength(5)
    BoxCutoff(3) = cLargestSphericalWell / BoxLength(9)
    CALL MakeListPotential( BoxCutoff, pPosition, BoxLengthInverse )
  END IF
END IF

! Search for overlapping configurations if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  ! Active transformation (orientation of particles)
  DO pParticle = 1, nParticles
    CALL ActiveTransformation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
  END DO
  ! Status
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
  WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 6 )//"OVERLAP CHECK FOR THE INITIAL CONFIGURATION"//REPEAT( " ", 6 )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
  WRITE( *, "(G0)" ) "Verifying initial configuration..."
  CALL Sleep( 1 )
  WRITE( *, "(G0)" ) " "
  ! Overlap check (initial configuration)
  IF( .NOT. CellListControl ) THEN
    ! Whole system
    CALL OverlapCheckInitialConfiguration(  )
  ELSE
    ! Linked lists
    CALL FullListOverlapCheckInitialConfiguration( ContactDistance, BoxLength, BoxLengthInverse, HalfNeighboursControl )
  END IF
  ! Status
  WRITE( *, "(G0)" ) "No overlaps found in the initial configuration. Resuming..."
  CALL Sleep( 1 )
  WRITE( *, "(G0)" ) " "
END IF

! Start simulation timer
CALL CPU_Time( StartTimer )

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"FILE ORGANIZER"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)", Advance= "NO" ) "Creating files..."

! Output file units
CALL FileHandler( LastLine )

! Status
CALL Sleep( 1 )
WRITE( *, "(G0)", Advance= "YES" ) " Done!"
WRITE( *, "(G0)" ) " "

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 16 )//"MONTE CARLO SIMULATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( EnsembleMC == "NVT" ) THEN
  WRITE( *, "(G0)" ) "Starting up the NVT-Monte Carlo simulation..."
ELSE IF( EnsembleMC == "NPT" ) THEN
  WRITE( *, "(G0)" ) "Starting up the NPT-Monte Carlo simulation..."
END IF
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Computation of total potential energy (initial configuration) if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  IF( PotentialTypeLogical(2) ) THEN
    WRITE( *, "(G0)", Advance= "NO" ) "Computing total potential energy of the initial configuration..."
    IF( .NOT. CellListControl ) THEN
      ! Whole system
      CALL ComputeTotalPotentialEnergy(  )
    ELSE
      ! Linked lists
      CALL ListComputeTotalPotentialEnergyInitialConfiguration( BoxLength, BoxLengthInverse, HalfNeighboursControl )
    END IF
    ! Status
    CALL Sleep( 1 )
    WRITE( *, "(G0)", Advance= "YES" ) " Done!"
    WRITE( *, "(G0)" ) " "
  END IF
END IF

! Pseudorandom number generator seed (if restoration is not selected)
IF( .NOT. RestoreBackupFileLogical ) THEN
  IF( FixedSeedLogical ) THEN
    SeedValue   = 123456789
    InitialSeed = SeedValue
  ELSE IF( .NOT. FixedSeedLogical ) THEN
    CALL RANDOM_NUMBER( RandomNumber )
    iSeed       = INT( RandomNumber * DBLE( SeedSize ) ) + 1
    SeedValue   = ABS( RandomSeed(iSeed) )
    InitialSeed = SeedValue
  END IF
END IF

! Initialize box distortion parameter if restoration is not selected
IF( .NOT. RestoreBackupFileLogical ) THEN
  CALL LatticeReduction( BoxLength, BoxDistortionMC, LatticeReductionLogical )
END IF

! Backup file
IF( BackupFileLogical ) THEN
  CALL WriteBackupVariables( DateTime )
END IF

! *********************************************************************************************** !
! Monte Carlo Parameters                                                                          !
! *********************************************************************************************** !
IF( .NOT. RestoreBackupFileLogical ) THEN
  FirstCycle                         = 0                                  ! First cycle of the simulation                    (initial value)
  MovementTranslationLogical         = .FALSE.                            ! Translational move selector                      (initial value)
  MovementRotationLogical            = .FALSE.                            ! Rotational move selector                         (initial value)
  MovementIsoVolumeChangeLogical     = .FALSE.                            ! Volume move selector (Isotropic)                 (initial value)
  MovementAnisoVolumeChangeLogical   = .FALSE.                            ! Volume move selector (Anisotropic)               (initial value)
  nAcceptanceTranslation             = 0                                  ! Translational move acceptance counter            (initial value)
  nAcceptanceRotation                = 0                                  ! Rotational move acceptance counter               (initial value)
  nAcceptanceIsotropicVolumeChange   = 0                                  ! Volumetric move acceptance counter (Isotropic)   (initial value)
  nAcceptanceAnisotropicVolumeChange = 0                                  ! Volumetric move acceptance counter (Anisotropic) (initial value)
  nMovementTranslationCounter        = 0                                  ! Translational move counter                       (initial value)
  nMovementRotationCounter           = 0                                  ! Rotational move counter                          (initial value)
  nMovementIsoVolumeChangeCounter    = 0                                  ! Volume scaling counter (Isotropic)               (initial value)
  nMovementAnisoVolumeChangeCounter  = 0                                  ! Volume scaling counter (Anisotropic)             (initial value)
  pQuaternionMC                      = pQuaternion                        ! Rotation quaternions                             (initial value)
  pPositionMC                        = pPosition                          ! Position of particles                            (initial value)
  pOrientationMC                     = pOrientation                       ! Orientation of particles                         (initial value)
  IF( PotentialTypeLogical(2) ) THEN
    TotalPotentialEnergyMC           = TotalPotentialEnergy               ! Total potential energy                           (initial value)
  END IF
  BoxLengthMC                        = BoxLength                          ! Box length                                       (initial value)
  BoxLengthInverseMC                 = BoxLengthInverse                   ! Box length (inverse)                             (initial value)
  BoxVolumeMC                        = BoxVolume                          ! Box volume                                       (initial value)
  VolumeScalingFactor                = 1.D0                               ! Scaling factor                                   (initial value)
  LatticeReductionLogical            = .FALSE.                            ! Lattice reduction                                (initial value)
  MaxTranslationalDisplacement       = UserMaxTranslationalDisplacement   ! Maximum translational displacement               (initial value)
  MaxAngularDisplacement             = UserMaxRotationalDisplacement      ! Maximum rotational displacement                  (initial value)
  MaxIsoVolumetricDisplacement       = UserMaxIsoVolumetricDisplacement   ! Maximum isovolumetric displacement               (initial value)
  MaxAnisoVolumetricDisplacement     = UserMaxAnisoVolumetricDisplacement ! Maximum anisovolumetric displacement             (initial value)
ELSE
  WRITE( *, "(G0)", Advance= "NO" ) "Restoring simulation variables from a previous run..."
  OPEN( Unit= 105, File= "Backup/"//TRIM( DescriptorBackupFile )//"_simulation.backup", Action= "READ" )
  READ( 105, * ) Dummy, FirstCycle
  IF( CellListLogical ) THEN
    READ( 105, * ) Dummy, CellListControl
    IF( PotentialTypeLogical(2) ) READ( 105, * ) Dummy, CellListControlPotential
    READ( 105, * ) Dummy, pCells
    IF( PotentialTypeLogical(2) ) READ( 105, * ) Dummy, pCellsPotential
    ! Allocation
    ALLOCATE( pCellHead(0:(pCells(1)-1),0:(pCells(2)-1),0:(pCells(3)-1)) )
    IF( PotentialTypeLogical(2) ) THEN
      ALLOCATE( pCellHeadPotential(0:(pCellsPotential(1)-1),0:(pCellsPotential(2)-1),0:(pCellsPotential(3)-1)) )
    END IF
    READ( 105, * ) Dummy, pCellHead
    IF( PotentialTypeLogical(2) ) READ( 105, * ) Dummy, pCellHeadPotential
    READ( 105, * ) Dummy, pCellList
    READ( 105, * ) Dummy, pCellIndex
    IF( PotentialTypeLogical(2) ) THEN
      READ( 105, * ) Dummy, pCellListPotential
      READ( 105, * ) Dummy, pCellIndexPotential
    END IF
  END IF
  READ( 105, * ) Dummy, pPositionMC
  READ( 105, * ) Dummy, pOrientationMC
  READ( 105, * ) Dummy, pQuaternionMC
  IF( PotentialTypeLogical(2) ) THEN
    READ( 105, * ) Dummy, TotalPotentialEnergyMC
  END IF
  IF( CellListLogical ) THEN
    READ( 105, * ) Dummy, OldBoxLength
  END IF
  READ( 105, * ) Dummy, BoxLengthMC
  READ( 105, * ) Dummy, BoxLengthInverseMC
  READ( 105, * ) Dummy, BoxVolumeMC
  READ( 105, * ) Dummy, BoxDistortionMC
  READ( 105, * ) Dummy, nAcceptanceTranslation
  READ( 105, * ) Dummy, nAcceptanceRotation
  READ( 105, * ) Dummy, nAcceptanceIsotropicVolumeChange
  READ( 105, * ) Dummy, nAcceptanceAnisotropicVolumeChange
  READ( 105, * ) Dummy, nMovementTranslationCounter
  READ( 105, * ) Dummy, nMovementRotationCounter
  READ( 105, * ) Dummy, nMovementIsoVolumeChangeCounter
  READ( 105, * ) Dummy, nMovementAnisoVolumeChangeCounter
  READ( 105, * ) Dummy, MaxTranslationalDisplacement
  READ( 105, * ) Dummy, MaxAngularDisplacement
  READ( 105, * ) Dummy, MaxIsoVolumetricDisplacement
  READ( 105, * ) Dummy, MaxAnisoVolumetricDisplacement
  READ( 105, * ) Dummy, SeedValue
  CLOSE( 105 )
  CALL Sleep( 1 )
  WRITE( *, "(G0)", Advance= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Deallocation
DEALLOCATE( pQuaternion, pPosition, pOrientation, SphericalComponentInquiry )
IF( PotentialTypeLogical(2) ) THEN
  DEALLOCATE( TotalPotentialEnergy )
END IF

! Allocation
ALLOCATE( PositionSaveMC(3,nParticles) )

! Metropolis Algorithm - Importance Sampling
WRITE( *, "(G0)" ) "Running Metropolis algorithm..."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Initialize cell list
IF( CellListLogical .AND. .NOT. RestoreBackupFileLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLengthMC(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLengthMC(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLengthMC(9)
  CALL MakeList( BoxCutoff, pPositionMC, BoxLengthInverseMC )
  IF( PotentialTypeLogical(2) ) THEN
    BoxCutoff(1) = cLargestSphericalWell / BoxLengthMC(1)
    BoxCutoff(2) = cLargestSphericalWell / BoxLengthMC(5)
    BoxCutoff(3) = cLargestSphericalWell / BoxLengthMC(9)
    CALL MakeListPotential( BoxCutoff, pPositionMC, BoxLengthInverseMC )
  END IF
  OldBoxLength = BoxLengthMC
END IF

! *********************************************************************************************** !
! Simulation cycles                                                                               !
! *********************************************************************************************** !
!  A 'cycle' is characterized by N trial moves (rotation or translation) of a random particle or  !
!  by a single change of the simulation box volume.                                               !
! *********************************************************************************************** !
DO iCycle = FirstCycle + 1, MaxSimulationCycles

  ! Simulation progress
  CALL ProgressBarMC( iCycle, MaxSimulationCycles, EnsembleMC )

  ! Prepare simulation box for single-particle moves
  IF( CellListLogical .AND. (MovementAnisoVolumeChangeLogical .OR. MovementIsoVolumeChangeLogical) ) THEN ! Check cells only after a volume change
    CALL BoxCheckNPT( pPositionMC, OldBoxLength, BoxLengthMC, BoxLengthInverseMC )
  END IF

  ! Generates a random number for the NPT-simulation
  IF( EnsembleMC == "NPT" ) THEN
    CALL RandomNumberGenLCG(  )
  END IF

  ! Movement (Translation or Rotation)
  IF( RandomNumber <= MovementProbability .OR. EnsembleMC == "NVT" ) THEN

    ! Disable volume scaling
    MovementIsoVolumeChangeLogical   = .FALSE.
    MovementAnisoVolumeChangeLogical = .FALSE.

    ! Particle loop
    DO pCycle = 1, nParticles

      ! Component index
      IF( nComponents > 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        CALL RandomNumberGenLCG(  )
        ! Get component index
        iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
      ELSE
        ! Get component index
        iComponent = 1
      END IF

      ! Avoid components with molar fraction of 0
      DO WHILE( cParticles(iComponent) < 1 )
        ! Component index
        IF( nComponents > 1 ) THEN
          ! Pseudorandom number generator (uniform distribution)
          CALL RandomNumberGenLCG(  )
          ! Get component index
          iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
        ELSE
          ! Get component index
          iComponent = 1
        END IF
      END DO

      ! Forbid rotation if component is spherical
      IF( SphericalComponentLogical(iComponent) ) THEN
        MovementTranslationLogical  = .TRUE.  ! Enable translation
        MovementRotationLogical     = .FALSE. ! Disable rotation
        nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
      ! Allow rotation if component is nonspherical
      ELSE
        ! Pseudorandom number generator (uniform distribution)
        CALL RandomNumberGenLCG(  )
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
      END IF

      ! Pseudorandom number generator (uniform distribution)
      CALL RandomNumberGenLCG(  )
      ! Random selection of particles of component cComponent
      iParticle = SUM( cParticles(0:(iComponent-1)) ) + INT( RandomNumber * DBLE( cParticles(iComponent) ) ) + 1

      ! Assignment of previous configuration (Microstate m)
      iOldPosition(:)    = pPositionMC(:,iParticle)    ! Position
      iOldQuaternion(:)  = pQuaternionMC(:,iParticle)  ! Quaternion
      iOldOrientation(:) = pOrientationMC(:,iParticle) ! Orientation

      ! Translational movement
      IF( MovementTranslationLogical ) THEN
        ! Random translation along x-axis
        CALL RandomNumberGenLCG(  )
        iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along y-axis
        CALL RandomNumberGenLCG(  )
        iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along z-axis
        CALL RandomNumberGenLCG(  )
        iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
      ! No translation
      ELSE IF( .NOT. MovementTranslationLogical ) THEN
        iNewPosition = iOldPosition
      END IF

      ! Rotational movement
      IF( MovementRotationLogical ) THEN
        ! Random quaternion
        CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
        ! Active transformation
        CALL ActiveTransformation( zAxis, iNewQuaternion, iNewOrientation )
      ! No rotation
      ELSE IF( .NOT. MovementRotationLogical ) THEN
        iNewQuaternion  = iOldQuaternion
        iNewOrientation = iOldOrientation
      END IF

      ! Overlap check after displacement of a particle
      IF( .NOT. CellListControl ) THEN
        ! Whole system
        CALL ParticleOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
        &                          BoxLengthMC, BoxLengthInverseMC, Overlap )
      ELSE
        ! Linked lists
        CALL ListOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
        &                      BoxLengthMC, BoxLengthInverseMC, Overlap, .FALSE. )
      END IF

      ! Acceptance criterion
      IF( .NOT. Overlap ) THEN
        ! System configuration update
        pPositionMC(:,iParticle)    = iNewPosition(:)    ! Update position
        pQuaternionMC(:,iParticle)  = iNewQuaternion(:)  ! Update quaternion
        pOrientationMC(:,iParticle) = iNewOrientation(:) ! Update orientation
        ! Update total potential energy
        IF( PotentialTypeLogical(2) ) THEN
          IF( .NOT. CellListControl ) THEN ! Whole system
            ! Computation of potential energy of particle i (microstate m)
            CALL ComputeParticlePotentialEnergy( iComponent, iParticle, iOldPosition, iOldPotentialEnergy, BoxLengthMC, &
            &                                    BoxLengthInverseMC )
            ! Computation of potential energy of particle i (microstate n)
            CALL ComputeParticlePotentialEnergy( iComponent, iParticle, iNewPosition, iNewPotentialEnergy, BoxLengthMC, &
            &                                    BoxLengthInverseMC )
          ELSE ! Linked lists
            ! Computation of potential energy of particle i (microstate m)
            CALL ListComputeParticlePotentialEnergy( iComponent, iParticle, iOldPosition, iOldPotentialEnergy, BoxLengthMC, &
            &                                        BoxLengthInverseMC, .FALSE. )
            ! Computation of potential energy of particle i (microstate n)
            CALL ListComputeParticlePotentialEnergy( iComponent, iParticle, iNewPosition, iNewPotentialEnergy, BoxLengthMC, &
            &                                        BoxLengthInverseMC, .FALSE. )
          END IF
          ! Computation of energy difference of microstates n and m
          iPotentialEnergyDifference = iNewPotentialEnergy - iOldPotentialEnergy
          ! System energy update
          TotalPotentialEnergyMC = TotalPotentialEnergyMC + iPotentialEnergyDifference
        END IF
        ! Displacement counter update
        IF( MovementTranslationLogical ) THEN
          IF( CellListControl ) CALL ParticleTranslationNVT( iParticle, ScalingDistanceUnitBox ) ! Update cell
          nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
        ELSE IF ( MovementRotationLogical ) THEN
          nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
        END IF
      ELSE
        ! Retrieve old configuration
        pPositionMC(:,iParticle)    = iOldPosition(:)    ! Retrieve old position
        pQuaternionMC(:,iParticle)  = iOldQuaternion(:)  ! Retrieve old quaternion
        pOrientationMC(:,iParticle) = iOldOrientation(:) ! Retrieve old orientation
      END IF

    END DO

  ! Volume scaling
  ELSE IF( RandomNumber > MovementProbability .AND. EnsembleMC == "NPT" ) THEN

    ! Disable translation and rotation
    MovementTranslationLogical = .FALSE.
    MovementRotationLogical    = .FALSE.

    ! Assignment of previous configuration (box)
    OldBoxLength        = BoxLengthMC        ! Box length
    OldBoxLengthInverse = BoxLengthInverseMC ! Box length (inverse)
    OldBoxVolume        = BoxVolumeMC        ! Box volume

    ! Expansion/compression type
    CALL RandomNumberGenLCG(  )

    ! Isotropic volume scaling
    IF( RandomNumber <= IsoVolumetricProbability ) THEN
      ! Random walk in the logarithm of the volume
      CALL RandomNumberGenLCG(  )
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
      MovementIsoVolumeChangeLogical   = .TRUE.  ! Enable isotropic volume scaling
      MovementAnisoVolumeChangeLogical = .FALSE. ! Disable anisotropic volume scaling
    ! Anisotropic volume scaling
    ELSE IF( RandomNumber > IsoVolumetricProbability ) THEN
      ! Random box component
      CALL RandomNumberGenLCG(  )
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
      ! Random deformation of the box
      CALL RandomNumberGenLCG(  )
      NewBoxLength(BoxMatrixComponent) = OldBoxLength(BoxMatrixComponent) + MaxAnisoVolumetricDisplacement * (RandomNumber - 0.5D0)
      VolumeType = "AN"
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
      ! Movement counter
      nMovementAnisoVolumeChangeCounter = nMovementAnisoVolumeChangeCounter + 1
      ! Movement type
      MovementAnisoVolumeChangeLogical = .TRUE.  ! Enable anisotropic volume scaling
      MovementIsoVolumeChangeLogical   = .FALSE. ! Disable isotropic volume scaling
    END IF

    ! Reset condition of anisotropic volume scaling (we must be careful not to induce a bias in the system)
    CheckBoxDistortion = .FALSE.

    ! Condition of anisotropic volume scaling (box distortion)
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
      IF( MovementIsoVolumeChangeLogical )   EnthalpyChange = ( ReducedPressure * ( NewBoxVolume - OldBoxVolume ) ) - &
      &                                                       ( DBLE( nParticles + 1 ) * DLOG( NewBoxVolume / OldBoxVolume ) )
      IF( MovementAnisoVolumeChangeLogical ) EnthalpyChange = ( ReducedPressure * ( NewBoxVolume - OldBoxVolume ) ) - &
      &                                                       ( DBLE( nParticles ) * DLOG( NewBoxVolume / OldBoxVolume ) )

      ! Random number
      CALL RandomNumberGenLCG(  )

      ! Enthalpy change criterion
      IF( DEXP( - EnthalpyChange ) >= RandomNumber ) THEN

        ! System configuration update
        PositionSaveMC = pPositionMC ! Old configuration

        ! Isotropic volume scaling
        IF( MovementIsoVolumeChangeLogical ) THEN
          ! Rescale positions of all particles accordingly
          DO pParticle = 1, nParticles
            pPositionMC(:,pParticle) = pPositionMC(:,pParticle) * VolumeScalingFactor
          END DO
        ! Anisotropic volume scaling
        ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
          ! Rescale positions of all particles accordingly
          DO pParticle = 1, nParticles
            ! Scaling coordinates using the old box length
            CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
            ! New real coordinates using the new box length
            CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceUnitBox, pPositionMC(:,pParticle) )
          END DO
        END IF

        ! Overlap check after expansion/compression of the simulation box
        IF( .NOT. CellListControl ) THEN
          ! Whole system
          CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        ELSE
          ! Linked lists
          CALL FullListOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap, HalfNeighboursControl )
          ! In case the number of cells in one direction (x, y, or z) becomes less than 3
          IF( .NOT. CellListControl ) CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        END IF

        ! Acceptance criterion
        IF( .NOT. Overlap ) THEN
          ! System configuration update
          BoxVolumeMC        = NewBoxVolume        ! Update volume
          BoxLengthMC        = NewBoxLength        ! Update box length
          BoxLengthInverseMC = NewBoxLengthInverse ! Update box length (inverse)
          ! Displacement counter update
          IF( MovementIsoVolumeChangeLogical ) THEN
            nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Isotropic move counter
          ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
            nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Anisotropic move counter
          END IF
          ! Update packing fraction and reduced number density
          PackingFraction    = ( TotalParticleVolume / NewBoxVolume )
          TotalNumberDensity = ( DBLE( nParticles ) / NewBoxVolume )
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
          pPositionMC        = PositionSaveMC      ! Retrieve old position of particles
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

    ! Adjustment of maximum displacement (translation)
    IF( DBLE( nMovementTranslationCounter ) >= DBLE( nAdjustmentFrequency * nParticles ) * AcceptanceRatioTranslation ) THEN

      ! Acceptance ratio (translation)
      IF( nMovementTranslationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
        ELSE
          MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
        END IF
        ! Ratio data (translation)
        IF( iCycle > LastLine(1) ) THEN
          WRITE( 30, "(7G0)" ) iCycle, ",", Ratio, ",", MaxTranslationalDisplacement, ",", AcceptanceRatioTranslation
          FLUSH( 30 )
        END IF
        ! Reset counter
        nAcceptanceTranslation      = 0
        nMovementTranslationCounter = 0
      END IF

      ! Avoid large translations
      BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( BoxLengthMC(1:3), BoxLengthMC(1:3) ) )
      BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( BoxLengthMC(4:6), BoxLengthMC(4:6) ) )
      BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( BoxLengthMC(7:9), BoxLengthMC(7:9) ) )
      IF( MaxTranslationalDisplacement > 2.D0 * MAXVAL( BoxEdgeLength ) ) THEN
        MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxEdgeLength )
      END IF

    END IF

    ! Adjustment of maximum displacement (rotation)
    IF( DBLE( nMovementRotationCounter ) >= DBLE( nAdjustmentFrequency * nParticles ) * AcceptanceRatioRotation ) THEN

      ! Acceptance ratio (rotation)
      IF( nMovementRotationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotational adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
        ELSE
          MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
        END IF
        ! Ratio data (rotation)
        IF( iCycle > LastLine(2) ) THEN
          WRITE( 40, "(7G0)" ) iCycle, ",", Ratio, ",", MaxAngularDisplacement, ",", AcceptanceRatioRotation
          FLUSH( 40 )
        END IF
        ! Reset counter
        nAcceptanceRotation      = 0
        nMovementRotationCounter = 0
      END IF

      ! Avoid 4π-rotations
      IF( MaxAngularDisplacement > 4.D0 * cPi ) THEN
        MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
      END IF

    END IF

    ! Adjustment of maximum displacement (isotropic volume scaling)
    IF( DBLE( nMovementIsoVolumeChangeCounter ) >= DBLE( nAdjustmentFrequency ) * AcceptanceRatioIsoVolumeChange ) THEN

      IF( nMovementIsoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (isotropic volume scaling)
        Ratio = DBLE( nAcceptanceIsotropicVolumeChange ) / DBLE( nMovementIsoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioIsoVolumeChange ) THEN
          MaxIsoVolumetricDisplacement = 0.95D0 * MaxIsoVolumetricDisplacement
        ELSE
          MaxIsoVolumetricDisplacement = 1.05D0 * MaxIsoVolumetricDisplacement
        END IF
        ! Ratio data (volume scaling)
        IF( EnsembleMC == "NPT" ) THEN
          IF( iCycle > LastLine(3) ) THEN
            WRITE( 50, "(9G0)" ) iCycle, ",", Ratio, ",", MaxIsoVolumetricDisplacement, ",", AcceptanceRatioIsoVolumeChange, &
            &                    ",", VolumeType
            FLUSH( 50 )
          END IF
        END IF
        ! Reset counter
        nAcceptanceIsotropicVolumeChange = 0
        nMovementIsoVolumeChangeCounter  = 0
      END IF

    END IF

    ! Adjustment of maximum displacement (anisotropic volume scaling)
    IF( DBLE( nMovementAnisoVolumeChangeCounter ) >= DBLE( nAdjustmentFrequency ) * AcceptanceRatioAnisoVolumeChange ) THEN

      IF( nMovementAnisoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (isotropic volume scaling)
        Ratio = DBLE( nAcceptanceAnisotropicVolumeChange ) / DBLE( nMovementAnisoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioAnisoVolumeChange ) THEN
          MaxAnisoVolumetricDisplacement = 0.95D0 * MaxAnisoVolumetricDisplacement
        ELSE
          MaxAnisoVolumetricDisplacement = 1.05D0 * MaxAnisoVolumetricDisplacement
        END IF
        ! Ratio data (volume scaling)
        IF( EnsembleMC == "NPT" ) THEN
          IF( iCycle > LastLine(3) ) THEN
            WRITE( 50, "(9G0)" ) iCycle, ",", Ratio, ",", MaxAnisoVolumetricDisplacement, ",", AcceptanceRatioAnisoVolumeChange, &
            &                    ",", VolumeType
            FLUSH( 50 )
          END IF
        END IF
        ! Reset counter
        nAcceptanceAnisotropicVolumeChange = 0
        nMovementAnisoVolumeChangeCounter  = 0
      END IF

    END IF

  END IF

  ! Order parameter data
  IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
    ! Nematic order parameter (Q-tensor method)
    CALL UniaxialNematicOrderParameter( OrderParameterMC, pOrientationMC )
    WRITE( 60, "(3G0)" ) iCycle, ",", OrderParameterMC
    FLUSH( 60 )
  END IF

  ! Trajectory data
  IF( TrajectoryLogical ) THEN
    IF( MOD ( iCycle, nSavingFrequency ) == 0 ) THEN
      WRITE( 20, "(G0)" ) nParticles
      WRITE( 20, * ) " "
      IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
        DO cComponent = 1, nComponents
          IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cLength(cComponent)
            END DO
          ELSE
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
            END DO
          END IF
        END DO
      ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
        DO cComponent = 1, nComponents
          IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), cLength(cComponent)
            END DO
          ELSE
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
            END DO
          END IF
        END DO
      ELSE IF( GeometryType(3) ) THEN ! Cylinders
        DO cComponent = 1, nComponents
          IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), cLength(cComponent)
            END DO
          ELSE
            DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
              WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPositionMC(1,pParticle), pPositionMC(2,pParticle), &
              &                          pPositionMC(3,pParticle), pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), &
              &                          pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), 0.5D0 * cDiameter(cComponent), &
              &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
            END DO
          END IF
        END DO
      END IF
      FLUSH( 20 )
    END IF
  END IF

  ! Results data (packing fraction, number density, box volume, and pressure)
  IF( MOD( iCycle, nSavingFrequency ) == 0 .AND. EnsembleMC == "NPT" ) THEN
    WRITE( 70, "(9G0)" ) iCycle, ",", PackingFraction, ",", TotalNumberDensity, ",", BoxVolumeMC, ",", ReducedPressure
    FLUSH( 70 )
  END IF

  ! Box properties data
  IF( MOD( iCycle, nSavingFrequency ) == 0 .AND. EnsembleMC == "NPT" ) THEN
    WRITE( 55, "(23G0)" ) iCycle, ",", BoxDistortionMC, ",", BoxVolumeMC, ",", BoxLengthMC(1), ",", BoxLengthMC(2), ",", &
    &                     BoxLengthMC(3), ",", BoxLengthMC(4), ",", BoxLengthMC(5), ",", BoxLengthMC(6), ",", &
    &                     BoxLengthMC(7), ",", BoxLengthMC(8), ",", BoxLengthMC(9)
    FLUSH( 55 )
  END IF

  ! Potential data (equilibration and production)
  IF( PotentialTypeLogical(2) ) THEN
    IF( .NOT. PotentialEnergyLogical ) THEN
      IF( MOD( iCycle, nSavingFrequency ) == 0 ) THEN
        DO rRange = 1, nRange
          WRITE( (80 + rRange), "(3G0)" ) iCycle, ",", TotalPotentialEnergyMC(rRange)
          FLUSH( (80 + rRange) )
        END DO
      END IF
    ! Potential data (production-related only)
    ELSE IF( PotentialEnergyLogical ) THEN
      IF( iCycle > nEquilibrationCycles .AND. MOD( iCycle, nSavingFrequency ) == 0 ) THEN
        DO rRange = 1, nRange
          WRITE( (80 + rRange), "(3G0)" ) iCycle, ",", TotalPotentialEnergyMC(rRange)
          FLUSH( (80 + rRange) )
        END DO
      END IF
    END IF
  END IF

  ! Backup file
  IF( BackupFileLogical .AND. (MOD( iCycle, nSavingFrequency ) == 0 .OR. iCycle == FirstCycle + 1) ) THEN
    OPEN( Unit= 105, File= "Backup/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_simulation.backup" )
    WRITE( 105, "(2G0)" ) "Current_Cycle: ", iCycle
    IF( CellListLogical ) THEN
      WRITE( 105, "(2G0)" ) "Cell_List_Control: ", CellListControl
      IF( PotentialTypeLogical(2) ) WRITE( 105, "(2G0)" ) "Cell_List_Control_Potential: ", CellListControlPotential
      WRITE( 105, "(G0,3(G0,1X))" ) "Current_Number_of_Cells: ", pCells
      IF( PotentialTypeLogical(2) ) WRITE( 105, "(G0,3(G0,1X))" ) "Current_Number_of_Cells_Potential: ", pCellsPotential
      WRITE( DescriptorBackupString, "(G0)" ) pCells(1) * pCells(2) * pCells(3)
      DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
      WRITE( 105, DescriptorBackupString ) "Cell_Head: ", pCellHead
      IF( PotentialTypeLogical(2) ) WRITE( 105, DescriptorBackupString ) "Cell_Head_Potential: ", pCellHeadPotential
      WRITE( DescriptorBackupString, "(G0)" ) nParticles
      DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
      WRITE( 105, DescriptorBackupString ) "Cell_List: ", pCellList
      WRITE( DescriptorBackupString, "(G0)" ) nParticles * 3
      DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
      WRITE( 105, DescriptorBackupString ) "Cell_Index: ", pCellIndex
      IF( PotentialTypeLogical(2) ) THEN
        WRITE( DescriptorBackupString, "(G0)" ) nParticles
        DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
        WRITE( 105, DescriptorBackupString ) "Cell_List_Potential: ", pCellListPotential
        WRITE( DescriptorBackupString, "(G0)" ) nParticles * 3
        DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
        WRITE( 105, DescriptorBackupString ) "Cell_Index_Potential: ", pCellIndexPotential
      END IF
    END IF
    WRITE( DescriptorBackupString, "(G0)" ) nParticles * 3
    DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
    WRITE( 105, DescriptorBackupString ) "Particle_Position: ", pPositionMC
    WRITE( 105, DescriptorBackupString ) "Particle_Orientation: ", pOrientationMC
    WRITE( DescriptorBackupString, "(G0)" ) nParticles * 4
    DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
    WRITE( 105, DescriptorBackupString ) "Particle_Quaternion: ", pQuaternionMC
    IF( PotentialTypeLogical(2) ) THEN
      WRITE( DescriptorBackupString, "(G0)" ) nRange
      DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
      WRITE( 105, DescriptorBackupString ) "Total_Potential_MC: ", TotalPotentialEnergyMC
    END IF
    IF( CellListLogical ) THEN
      WRITE( 105, "(G0,9(G0.11,1X))" ) "Old_Box_Length: ", OldBoxLength
    END IF
    WRITE( 105, "(G0,9(G0.11,1X))" ) "Box_Length: ", BoxLengthMC
    WRITE( 105, "(G0,9(G0.11,1X))" ) "Box_Length_Inverse: ", BoxLengthInverseMC
    WRITE( 105, "(G0,G0.11)" ) "Box_Volume: ", BoxVolumeMC
    WRITE( 105, "(G0,G0.11)" ) "Box_Distortion: ", BoxDistortionMC
    WRITE( 105, "(2G0)" ) "Acceptance_Translation: ", nAcceptanceTranslation
    WRITE( 105, "(2G0)" ) "Acceptance_Rotation: ", nAcceptanceRotation
    WRITE( 105, "(2G0)" ) "Acceptance_Isotropic_Volume_Change: ", nAcceptanceIsotropicVolumeChange
    WRITE( 105, "(2G0)" ) "Acceptance_Anisotropic_Volume_Change: ", nAcceptanceAnisotropicVolumeChange
    WRITE( 105, "(2G0)" ) "Movement_Translation_Counter: ", nMovementTranslationCounter
    WRITE( 105, "(2G0)" ) "Movement_Rotation_Counter: ", nMovementRotationCounter
    WRITE( 105, "(2G0)" ) "Movement_Iso_Volume_Change_Counter: ", nMovementIsoVolumeChangeCounter
    WRITE( 105, "(2G0)" ) "Movement_Aniso_Volume_Change_Counter: ", nMovementAnisoVolumeChangeCounter
    WRITE( 105, "(2G0)" ) "Max_Translational_Displacement: ", MaxTranslationalDisplacement
    WRITE( 105, "(2G0)" ) "Max_Angular_Displacement: ", MaxAngularDisplacement
    WRITE( 105, "(2G0)" ) "Max_Iso_Volumetric_Displacement: ", MaxIsoVolumetricDisplacement
    WRITE( 105, "(2G0)" ) "Max_Aniso_Volumetric_Displacement: ", MaxAnisoVolumetricDisplacement
    WRITE( 105, "(2G0)" ) "Seed_Value: ", SeedValue
    CLOSE( 105 )
  END IF

END DO

! End of Metropolis algorithm
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Monte Carlo simulation finished successfully! See directories for results."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Output unit
IF( TrajectoryLogical ) THEN
  CLOSE( 20 )
END IF
CLOSE( 30 )
CLOSE( 40 )
CLOSE( 60 )
IF( EnsembleMC == "NPT" ) THEN
  CLOSE( 50 )
  CLOSE( 55 )
  CLOSE( 70 )
END IF
IF( PotentialTypeLogical(2) ) THEN
  DO rRange = 1, nRange
    CLOSE ( 80 + rRange )
  END DO
END IF

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"SIMULATION LENGTH"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR

! End simulation timer
CALL CPU_Time( StopTimer )
WRITE( *, "(G0,G0.5,G0)" ) "Elapsed Time: ", (StopTimer - StartTimer), "s."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Deallocation of arrays
DEALLOCATE( pQuaternionMC, pPositionMC, pOrientationMC, PositionSaveMC )
IF( PotentialTypeLogical(2) ) THEN
  DEALLOCATE( TotalPotentialEnergyMC, iNewPotentialEnergy, iOldPotentialEnergy, iPotentialEnergyDifference )
  DEALLOCATE( cPotentialRange )
END IF
IF( CellListLogical ) CALL FinalizeList(  )
IF( CellListLogical .AND. PotentialTypeLogical(2) ) CALL FinalizeListPotential(  )

! Calculation of the first- and second-order TPT coefficients
IF( PerturbationCoefficientLogical .AND. EnsembleMC == "NVT" .AND. PotentialTypeLogical(2) ) THEN
  ! Status
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
  WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"TPT PARAMETERS"//REPEAT( " ", 21 )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
  WRITE( *, "(G0)" ) "Initializing calculation of the TPT parameters..."
  WRITE( *, "(G0)" ) " "
  ! Calculation of the TPT coefficients for every potential range
  BlockAveragePotentialRange: DO rRange = 1, nRange
    ! Condition
    IF( ( MinBlocks > MaxBlocks ) ) THEN
      WRITE( *, "(5G0)" ) "Minimum number of blocks (", MinBlocks, ") greater than the maximum number of blocks (", MaxBlocks, &
      &                   "). Skipping..."
      WRITE( *, "(G0)" ) " "
      ! Terminate calculation if min. block greater than max. blocks
      EXIT BlockAveragePotentialRange
    END IF
    ! Status
    WRITE( *, "(G0,G0.5,G0)", Advance= "NO" ) "Current Potential Range: ", PotentialRange(rRange), ". Calculating..."
    WRITE( DescriptorRange, FormatRange ) PotentialRange(rRange)
    ! Initialization
    FlagLogical = .FALSE.
    ! Block-averaging method
    CALL BlockAverage( FlagLogical, CoefficientTPT1(rRange), CoefficientTPT2(rRange), rFreeEnergy(rRange), &
    &                  CoefficientTPT1Deviation(rRange), CoefficientTPT2Deviation(rRange), rFreeEnergyDeviation(rRange), &
    &                  ExecutionTime )
    ! Stop criterion
    IF( FlagLogical ) THEN
      WRITE( *, "(G0)", Advance= "YES" ) " "
      WRITE( *, "(5G0)" ) "Number of blocks (", MaxBlocks - MinBlocks, ") exceed the maximum number of blocks (", 10000, &
      &                   "). Skipping..."
      WRITE( *, "(G0)" ) " "
      ! Terminate calculation if min. block greater than max. blocks
      EXIT BlockAveragePotentialRange
    END IF
    ! Status
    WRITE( *, "(G0,G0.5,G0)", Advance= "YES" ) " Done in ", ExecutionTime, "s."
  END DO BlockAveragePotentialRange
  ! Subroutine returns no errors
  IF( .NOT. FlagLogical .AND. ( MinBlocks <= MaxBlocks ) ) THEN
    DO rRange = 1, nRange
      WRITE ( DescriptorRange, FormatRange ) PotentialRange(rRange)
      OPEN( Unit= 150, File= "Perturbed_Coefficient/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/" &
      &                      //TRIM( DescriptorHour )//"_TPT_coefficients_η"//TRIM( DescriptorFileThermoVariable )//"_C" &
      &                      //TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".csv")
      WRITE( 150, "(G0,G0)" ) "# of Particles: ", nParticles
      WRITE( 150, "(G0,G0.5)" ) "Reduced Temperature: ", ReducedTemperature
      WRITE( 150, "(G0,G0.5)" ) "Reduced Number Density: ", TotalNumberDensity
      WRITE( 150, "(G0,G0.5)" ) "Aspect Ratio (L/D): ", cAspectRatio
      WRITE( 150, "(G0,G0.5)" ) "Attractive Range: ", PotentialRange(rRange)
      WRITE( 150, "(G0)" ) " "
      WRITE( 150, "(11G0)" ) "'1st_Order_TPTCoefficient'", ",", "'1st_Order_TPTCoefficient_STD'", ",", &
      &                      "'2nd_Order_TPTCoefficient'", ",", "'2nd_Order_TPTCoefficient_STD'", ",", &
      &                      "'Perturbed_Helmholtz_FEnergy'", ",", "'Perturbed_Helmholtz_FEnergy_STD'"
      WRITE( 150, "(G0.7,5(G0,G0.7))" ) CoefficientTPT1(rRange), ",", CoefficientTPT1Deviation(rRange), ",", &
      &                                 CoefficientTPT2(rRange), ",", CoefficientTPT2Deviation(rRange), ",", &
      &                                 rFreeEnergy(rRange), ",", rFreeEnergyDeviation(rRange)
      CLOSE( 150 )
    END DO
    WRITE( *, "(G0)", Advance= "YES" ) " "
  END IF
END IF

! Write down results
WRITE( *, "(G0)", Advance= "NO" ) "Writing simulation log..."
CALL Sleep( 1 )
CALL SimulationLog( BoxVolumeMC, BoxLengthMC, StopTimer - StartTimer, DateTime )

! Status
WRITE( *, "(G0)", Advance= "YES" ) " Done!"

END PROGRAM Main
