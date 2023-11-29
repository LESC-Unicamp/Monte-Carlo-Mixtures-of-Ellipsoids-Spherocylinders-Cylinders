! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!         to generate a log file containg all pertinent information about the simulation.         !
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
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                       This subroutine generates a log for the simulation                        !
! *********************************************************************************************** !
SUBROUTINE SimulationLog( BoxVolumeMC, BoxLengthMC, ExecutionTime, DateTime )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iString    ! Counter (string)
INTEGER( Kind= Int64 ) :: cComponent ! Counter (component)
INTEGER( Kind= Int64 ) :: rRange     ! Counter (potential range)
INTEGER( Kind= Int64 ) :: Remainder  ! Remainder of division

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: ExecutionTime ! Execution time
REAL( Kind= Real64 )                 :: BoxVolumeMC   ! Reduced volume of the simulation box
REAL( Kind= Real64 ), DIMENSION( 9 ) :: BoxLengthMC   ! Length (x,y,z) of the simulation box

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: StringsSizeHeader   ! String size (header)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: StringsSizeSubtitle ! String size (subtitle)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: StringsSizeText     ! String size (text)

! *********************************************************************************************** !
! CHARACTER STRINGS (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
CHARACTER( LEN= 64 ) :: FormatClock ! String format (clock)

! *********************************************************************************************** !
! CHARACTER STRINGS (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
CHARACTER( LEN= 68 ), DIMENSION( : ), ALLOCATABLE     :: StringsHeader   ! Strings used in the simulation log file (header)
CHARACTER( LEN= 68 ), DIMENSION( : ), ALLOCATABLE     :: StringsSubtitle ! Strings used in the simulation log file (subtitle)
CHARACTER( LEN= 68 ), DIMENSION( :, : ), ALLOCATABLE  :: StringsText     ! Strings used in the simulation log file (text)
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: StringsPRange   ! Strings used in the simulation log file (range)
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: CharLabel       ! Simulation log
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: CharLabelPRange ! Simulation log

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FileExist ! Checks whether a file exists or not

! Allocation
ALLOCATE( CharLabel(73,nComponents) )
ALLOCATE( CharLabelPRange(1,nRange ) )

! Simulation log descriptors
WRITE( CharLabel(1,1), "(G0)"    ) nComponents
WRITE( CharLabel(2,1), "(G0.5)"  ) PackingFraction
DO cComponent = 1, nComponents
  WRITE( CharLabel(3,cComponent), "(G0.5)"  ) cDiameter(cComponent)
  WRITE( CharLabel(4,cComponent), "(G0.5)"  ) cLength(cComponent)
  WRITE( CharLabel(5,cComponent), "(G0.5)"  ) cAspectRatio(cComponent)
  WRITE( CharLabel(6,cComponent), "(G0.5)"  ) cMolarFraction(cComponent)
  WRITE( CharLabel(7,cComponent), "(G0)"    ) cParticles(cComponent)
  WRITE( CharLabel(8,cComponent), "(G0.5)"  ) cMolecularVolume(cComponent)
  WRITE( CharLabel(9,cComponent), "(G0.5)"  ) cNumberDensity(cComponent)
END DO
WRITE( CharLabel(10,1), "(G0)"   ) nParticles
WRITE( CharLabel(11,1), "(G0.5)" ) TotalParticleVolume
WRITE( CharLabel(12,1), "(G0.5)" ) TotalNumberDensity
WRITE( CharLabel(13,1), "(G0.5)" ) AbsoluteTemperature
WRITE( CharLabel(14,1), "(G0.5)" ) ReducedPressure
WRITE( CharLabel(15,1), "(G0.5)" ) BoxVolumeMC
WRITE( CharLabel(16,1), "(2(E11.4,G0),E11.4)" ) BoxLengthMC(1), ", ", BoxLengthMC(2), ", ", BoxLengthMC(3)
WRITE( CharLabel(17,1), "(2(E11.4,G0),E11.4)" ) BoxLengthMC(4), ", ", BoxLengthMC(5), ", ", BoxLengthMC(6)
WRITE( CharLabel(18,1), "(2(E11.4,G0),E11.4)" ) BoxLengthMC(7), ", ", BoxLengthMC(8), ", ", BoxLengthMC(9)
WRITE( CharLabel(19,1), "(G0)"   ) MaxSimulationCycles
WRITE( CharLabel(20,1), "(G0)"   ) nEquilibrationCycles
WRITE( CharLabel(21,1), "(G0)"   ) MaxSimulationCycles - nEquilibrationCycles
WRITE( CharLabel(22,1), "(G0)"   ) nSavingFrequency
WRITE( CharLabel(23,1), "(G0)"   ) nAdjustmentFrequency
WRITE( CharLabel(24,1), "(G0.5)" ) AcceptanceRatioTranslation
WRITE( CharLabel(25,1), "(G0.5)" ) AcceptanceRatioRotation
WRITE( CharLabel(26,1), "(G0.5)" ) AcceptanceRatioIsoVolumeChange
WRITE( CharLabel(54,1), "(G0.5)" ) AcceptanceRatioAnisoVolumeChange
WRITE( CharLabel(60,1), "(G0.5)" ) MovementProbability
WRITE( CharLabel(27,1), "(G0.5)" ) TranslationalProbability
WRITE( CharLabel(28,1), "(G0.5)" ) RotationalProbability
WRITE( CharLabel(29,1), "(G0.5)" ) VolumeChangeProbability
WRITE( CharLabel(30,1), "(G0.5)" ) IsoVolumetricProbability
WRITE( CharLabel(55,1), "(G0.5)" ) AnisoVolumetricProbability
WRITE( CharLabel(31,1), "(G0.5)" ) UserMaxTranslationalDisplacement
WRITE( CharLabel(32,1), "(G0.5)" ) UserMaxRotationalDisplacement
WRITE( CharLabel(33,1), "(G0.5)" ) UserMaxIsoVolumetricDisplacement
WRITE( CharLabel(56,1), "(G0.5)" ) UserMaxAnisoVolumetricDisplacement
WRITE( CharLabel(34,1), "(G0)"   ) BodyFixedAxisInquiry
WRITE( CharLabel(35,1), "(G0.5)" ) QuaternionAngle
WRITE( CharLabel(60,1), "(G0.5)" ) MaxBoxDistortion
WRITE( CharLabel(61,1), "(G0.5)" ) BoxEdgeMaxRatio
WRITE( CharLabel(62,1), "(G0.5)" ) BoxVectorMaxAngle * 180.D0 / cPi
IF( LatticeReductionTypeLogical(1) ) THEN
  WRITE( CharLabel(63,1), "(G0)" ) "Gottwald"
ELSE IF( LatticeReductionTypeLogical(2) ) THEN
  WRITE( CharLabel(63,1), "(G0)" ) "Lenstra-Lenstra-Lovász"
END IF
WRITE( CharLabel(36,1), "(G0)"   ) InitialConfiguration
IF( ( ExecutionTime ) < 60.D0 ) THEN
  WRITE( CharLabel(37,1), "(2G0)" ) FLOOR( ExecutionTime ), "s"
ELSE IF( ( ExecutionTime ) < 3600.D0 ) THEN
  WRITE( CharLabel(37,1), "(4G0)" ) FLOOR( ( ExecutionTime ) / 60.D0 ), "m:", &
  &                                 FLOOR( ( ExecutionTime ) - 60.D0 * DBLE( FLOOR( ( ExecutionTime ) / &
  &                                 60.D0 ) ) ), "s"
ELSE IF( ( ExecutionTime ) < 86400.D0 ) THEN
  WRITE( CharLabel(37,1), "(6G0)" ) FLOOR( ( ExecutionTime ) / 3.6d3 ), "h:", &
  &                                 FLOOR( DBLE( FLOOR( ( ExecutionTime ) / 60.D0 ) ) - 60.D0 * &
  &                                 DBLE( FLOOR( ( ExecutionTime ) / 3.6d3 ) ) ), "m:", &
  &                                 FLOOR( ( ExecutionTime ) - 3600.D0 * DBLE( FLOOR( ( ExecutionTime ) / &
  &                                 3.6d3 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( ExecutionTime ) / 60.D0 ) ) - &
  &                                 60.D0 * DBLE( FLOOR( ( ExecutionTime ) / 3.6d3 ) ) ) ) ), "s"
ELSE IF( ( ExecutionTime ) >= 86400.D0 ) THEN
  WRITE( CharLabel(37,1), "(8G0)" ) FLOOR( ( ExecutionTime ) / 864.D2 ), "d:", &
  &                                 FLOOR( DBLE( FLOOR( ( ExecutionTime ) / 3600.D0 ) ) - 24.D0 * &
  &                                 DBLE( FLOOR( ( ExecutionTime ) / 864.D2 ) ) ), "h:", &
  &                                 FLOOR( DBLE( FLOOR( ( ExecutionTime ) / 60.D0 ) ) - 1440.D0 * &
  &                                 DBLE( FLOOR( ( ExecutionTime ) / 864.D2 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( &
  &                                 ExecutionTime ) / 3600.D0 ) ) - 24.D0 * DBLE( FLOOR( ( ExecutionTime ) / 864.D2 ) ) ) ) ), &
  &                                 "m:", FLOOR( ( ExecutionTime ) - 864.D2 * DBLE( FLOOR( ( ExecutionTime ) / &
  &                                 864.D2 ) ) - 3600.D0 * DBLE( FLOOR( DBLE( FLOOR( ( ExecutionTime ) / 3600.D0 ) ) - &
  &                                 24.D0 * DBLE( FLOOR( ( ExecutionTime ) / 864.D2 ) ) ) ) - 60.D0 * &
  &                                 DBLE( FLOOR( ( ( ExecutionTime ) / 60.D0 ) - 1440.D0 * DBLE( FLOOR( ( ExecutionTime ) / &
  &                                 864.D2 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( ExecutionTime ) / &
  &                                 3600.D0 ) ) - 24.D0 * DBLE( FLOOR( ( ExecutionTime ) / 864.D2 ) ) ) ) ) ) ), "s"
END IF
FormatClock = "(I4,G0,I2.2,G0,I2.2,G0,I2.2,G0,I2.2,G0,I2.2)"
WRITE( CharLabel(38,1), FormatClock ) DateTime(1), "/", DateTime(2), "/", DateTime(3), " ", &
&                                     DateTime(5), ":", DateTime(6), ":", DateTime(7)
DO cComponent = 1, nComponents
  WRITE( CharLabel(39,cComponent), "(G0)" ) cComponent
END DO
IF( ConfigurationSelection(4) ) THEN
  WRITE( CharLabel(40,1), "(G0.5)" ) PackingFractionInitialConfiguration
  WRITE( CharLabel(42,1), "(G0.5)" ) PressureRandomConfig
  WRITE( CharLabel(43,1), "(G0)"   ) nAdjustmentRandomConfig
  WRITE( CharLabel(44,1), "(G0.5)" ) MaxTranslationalDisplacementRandomConfig
  WRITE( CharLabel(45,1), "(G0.5)" ) MaxAngularDisplacementRandomConfig
  WRITE( CharLabel(46,1), "(G0.5)" ) MaxIsoVolumetricDisplacementRandomConfig
  WRITE( CharLabel(57,1), "(G0.5)" ) MaxAnisoVolumetricDisplacementRandomConfig
  WRITE( CharLabel(47,1), "(G0.5)" ) MinVolumetricDisplacementRandomConfig
  WRITE( CharLabel(48,1), "(G0.5)" ) MovementProbabilityRandomConfig
  WRITE( CharLabel(49,1), "(G0.5)" ) VolumeChangeProbabilityRandomConfig
  WRITE( CharLabel(50,1), "(G0.5)" ) TranslationalProbabilityRandomConfig
  WRITE( CharLabel(51,1), "(G0.5)" ) RotationalProbabilityRandomConfig
  WRITE( CharLabel(58,1), "(G0.5)" ) IsoVolumetricProbabilityRandomConfig
  WRITE( CharLabel(59,1), "(G0.5)" ) AnisoVolumetricProbabilityRandomConfig
END IF
IF( PresetInitialConfiguration ) THEN
  WRITE( CharLabel(52,1), "(G0)" ) "[YES]"
ELSE IF( .NOT. PresetInitialConfiguration ) THEN
  WRITE( CharLabel(52,1), "(G0)" ) "[NO]"
END IF
IF( FixedSeedLogical ) THEN
  WRITE( CharLabel(53,1), "(G0)" ) "[YES]"
ELSE IF( .NOT. FixedSeedLogical ) THEN
  WRITE( CharLabel(53,1), "(G0)" ) "[NO]"
END IF
IF( PotentialTypeLogical(1) ) THEN
  WRITE( CharLabel(64,1), "(G0)" ) "Purely Repulsive Hard-Core"
ELSE IF( PotentialTypeLogical(2) ) THEN
  WRITE( CharLabel(64,1), "(G0)" ) "Spherical Square-Well"
END IF
IF( PotentialTypeLogical(2) ) THEN
  WRITE( CharLabel(65,1), "(G0)"   ) nRange
  DO rRange = 1, nRange
    WRITE( CharLabelPRange(1,rRange), "(G0.5)" ) PotentialRange(rRange)
  END DO
  WRITE( CharLabel(67,1), "(G0.5)" ) ReducedTemperature
  IF( PotentialEnergyLogical ) THEN
    WRITE( CharLabel(68,1), "(G0)" ) "Only Production Cycles"
  ELSE IF( .NOT. PotentialEnergyLogical ) THEN
    WRITE( CharLabel(68,1), "(G0)" ) "Equilibration and Production Cycles"
  END IF
  WRITE( CharLabel(69,1), "(G0)"   ) MinBlocks
  WRITE( CharLabel(70,1), "(G0)"   ) MaxBlocks
END IF
DO cComponent = 1, nComponents
  WRITE( CharLabel(71,cComponent), "(G0)"   ) SphericalComponentLogical(cComponent)
END DO
WRITE( CharLabel(72,1), "(G0)"   ) InitialSeed
IF( CellListLogical ) THEN
  WRITE( CharLabel(73,1), "(G0)" ) "[ENABLED]"
ELSE
  WRITE( CharLabel(73,1), "(G0)" ) "[DISABLED]"
END IF

! Log strings
ALLOCATE( StringsHeader(6), StringsText(79,nComponents), StringsSubtitle(6) )
ALLOCATE( StringsPRange( 1,( INT( DBLE( nRange ) / 4.D0 ) + 1 ) ) )
ALLOCATE( StringsSizeHeader(6), StringsSizeText(79,nComponents), StringsSizeSubtitle(6) )

! String name (header)
StringsHeader(1) = "MONTE CARLO SIMULATION LOG"
StringsHeader(2) = "NVT/NPT-Monte Carlo algorithm for cylindrycally-symmetric molecules"
StringsHeader(3) = "Nathan Barros de Souza"
StringsHeader(4) = "University of Campinas"
StringsHeader(5) = "School of Chemical Engineering"
StringsHeader(6) = "Supervisor: Luis Fernando Mercier Franco"

! String name (text)
StringsText(1,1) = "Execution Date: "//TRIM( CharLabel(38,1) )
StringsText(2,1) = "Ensemble: "//TRIM( EnsembleMC )
StringsText(3,1) = "Molecular Shape: "//TRIM( MolecularGeometry )
StringsText(4,1) = "Number of Components: "//TRIM( CharLabel(1,1) )
StringsText(5,1) = "Packing Fraction: "//TRIM( CharLabel(2,1) )
DO cComponent = 1, nComponents
  StringsText(6,cComponent)  = "● COMPONENT #"//TRIM( CharLabel(39,cComponent) )
  IF( SphericalComponentLogical(cComponent) ) THEN
    StringsText(77,cComponent) = "Spherical Component"
  ELSE
    StringsText(77,cComponent) = "Non-Spherical Component"
  END IF
  StringsText(7,cComponent)  = "Diameter: "//TRIM( CharLabel(3,cComponent) )//"Å"
  StringsText(8,cComponent)  = "Length: "//TRIM( CharLabel(4,cComponent) )//"Å"
  StringsText(9,cComponent)  = "Elongation: "//TRIM( CharLabel(5,cComponent) )
  StringsText(10,cComponent) = "Molar Fraction: "//TRIM( CharLabel(6,cComponent) )
  StringsText(11,cComponent) = "Number of Particles: "//TRIM( CharLabel(7,cComponent) )
  StringsText(12,cComponent) = "Particle Volume: "//TRIM( CharLabel(8,cComponent) )//"Å³"
  StringsText(13,cComponent) = "Number Density: "//TRIM( CharLabel(9,cComponent) )//"Å⁻³"
END DO
StringsText(14,1) = "Total Number of Particles: "//TRIM( CharLabel(10,1) )
StringsText(15,1) = "Total Particle Volume: "//TRIM( CharLabel(11,1) )//"Å³"
StringsText(16,1) = "Total Number Density: "//TRIM( CharLabel(12,1) )//"Å⁻³"
StringsText(17,1) = "Absolute Temperature: "//TRIM( CharLabel(13,1) )//"K"
StringsText(18,1) = "Reduced Pressure: "//TRIM( CharLabel(14,1) )
StringsText(19,1) = "Box Volume: "//TRIM( CharLabel(15,1) )//"Å³"
StringsText(20,1) = "Box Length (Å):"
StringsText(21,1) = TRIM( CharLabel(16,1) )
StringsText(22,1) = TRIM( CharLabel(17,1) )
StringsText(23,1) = TRIM( CharLabel(18,1) )
StringsText(24,1) = "Total Number of Cycles: "//TRIM( CharLabel(19,1) )
StringsText(25,1) = "Equilibration Cycles: "//TRIM( CharLabel(20,1) )
StringsText(26,1) = "Production Cycles: "//TRIM( CharLabel(21,1) )
StringsText(27,1) = "Saving Frequency (Cycles): "//TRIM( CharLabel(22,1) )
StringsText(28,1) = "Adjustment Frequency (Cycles): "//TRIM( CharLabel(23,1) )
StringsText(29,1) = "Ratio Threshold (Translation): "//TRIM( CharLabel(24,1) )
StringsText(30,1) = "Ratio Threshold (Rotation): "//TRIM( CharLabel(25,1) )
StringsText(31,1) = "Ratio Threshold (Isotropic): "//TRIM( CharLabel(26,1) )
StringsText(58,1) = "Ratio Threshold (Anisotropic): "//TRIM( CharLabel(54,1) )
StringsText(59,1) = "Probability (Movement): "//TRIM( CharLabel(60,1) )
StringsText(32,1) = "Probability (Translation): "//TRIM( CharLabel(27,1) )
StringsText(33,1) = "Probability (Rotation): "//TRIM( CharLabel(28,1) )
StringsText(34,1) = "Probability (Volume): "//TRIM( CharLabel(29,1) )
StringsText(35,1) = "Probability (Isotropic): "//TRIM( CharLabel(30,1) )
StringsText(60,1) = "Probability (Anisotropic): "//TRIM( CharLabel(55,1) )
StringsText(36,1) = "Initial Maximum Displacement (Translation): "//TRIM( CharLabel(31,1) )
StringsText(37,1) = "Initial Maximum Displacement (Rotation): "//TRIM( CharLabel(32,1) )
StringsText(38,1) = "Initial Maximum Displacement (Isotropic): "//TRIM( CharLabel(33,1) )
StringsText(61,1) = "Initial Maximum Displacement (Anisotropic): "//TRIM( CharLabel(56,1) )
StringsText(65,1) = "Maximum Box Distortion: "//TRIM( CharLabel(60,1) )
StringsText(66,1) = "Maximum Box Length Distortion: "//TRIM( CharLabel(61,1) )
StringsText(67,1) = "Maximum Box Angle Distortion: "//TRIM( CharLabel(62,1) )//"°"
StringsText(68,1) = "Lattice Reduction Algorithm: "//TRIM( CharLabel(63,1) )
StringsText(39,1) = "Unrotated Axis (Initial Configuration): "//TRIM( CharLabel(34,1) )
StringsText(40,1) = "Quaternion Angle (Initial Configuration): "//TRIM( CharLabel(35,1) )
StringsText(41,1) = "● Configuration: "//TRIM( CharLabel(36,1) )//" Structure"
StringsText(42,1) = "Initial Packing Fraction (Random Configuration): "//TRIM( CharLabel(40,1) )
StringsText(44,1) = "Target Reduced Pressure (Random Configuration): "//TRIM( CharLabel(42,1) )
StringsText(45,1) = "Adjustment Frequency (Random Configuration): Every "//TRIM( CharLabel(43,1) )//" Cycle(s)"
StringsText(46,1) = "Maximum Rotation (Random Configuration): "//TRIM( CharLabel(44,1) )
StringsText(47,1) = "Maximum Translation (Random Configuration): "//TRIM( CharLabel(45,1) )
StringsText(48,1) = "Maximum Isotropic Volume Scaling (Random Configuration): "//TRIM( CharLabel(46,1) )
StringsText(62,1) = "Maximum Anisotropic Volume Scaling (Random Configuration): "//TRIM( CharLabel(57,1) )
StringsText(49,1) = "Minimum Volume Scaling (Random Configuration): "//TRIM( CharLabel(47,1) )
StringsText(50,1) = "Movement Probability (Random Configuration): "//TRIM( CharLabel(48,1) )
StringsText(51,1) = "Volume Scaling Probability (Random Configuration): "//TRIM( CharLabel(49,1) )
StringsText(52,1) = "Translation Probability (Random Configuration): "//TRIM( CharLabel(50,1) )
StringsText(53,1) = "Rotation Probability (Random Configuration): "//TRIM( CharLabel(51,1) )
StringsText(63,1) = "Isotropic Probability (Random Configuration): "//TRIM( CharLabel(58,1) )
StringsText(64,1) = "Anisotropic Probability (Random Configuration): "//TRIM( CharLabel(59,1) )
StringsText(54,1) = "Preset Initial Configuration: "//TRIM( CharLabel(52,1) )
IF( TrajectoryLogical ) THEN
  StringsText(55,1) = "Trajectory of particles computed."
ELSE IF( .NOT. TrajectoryLogical ) THEN
  StringsText(55,1) = "Trajectory of particles not computed."
END IF
StringsText(56,1) = "Fixed seed: "//TRIM( CharLabel(53,1) )
StringsText(57,1) = "Simulation length: "//TRIM( CharLabel(37,1) )
StringsText(69,1) = "● Potential Type: "//TRIM( CharLabel(64,1) )//" Potential"
StringsText(70,1) = "Number of Attractive Range Points: "//TRIM( CharLabel(65,1) )
IF( nRange >= 4 ) THEN
  Remainder = MOD( INT( nRange ), 4 )
ELSE
  Remainder = nRange
END IF
IF( nRange >= 4 .AND. Remainder == 0 ) THEN
  DO rRange = 1, INT( DBLE( nRange ) / 4.D0 ) - 1
    StringsPRange(1,rRange) = TRIM( CharLabelPRange(1,4*rRange-3) )//", "//TRIM( CharLabelPRange(1,4*rRange-2) )//", "// &
    &                         TRIM( CharLabelPRange(1,4*rRange-1) )//", "//TRIM( CharLabelPRange(1,4*rRange-0) )//", "
  END DO
  rRange = INT( DBLE( nRange ) / 4.D0 )
  StringsPRange(1,rRange) = TRIM( CharLabelPRange(1,4*rRange-3) )//", "//TRIM( CharLabelPRange(1,4*rRange-2) )//", "// &
  &                         TRIM( CharLabelPRange(1,4*rRange-1) )//", "//TRIM( CharLabelPRange(1,4*rRange-0) )
ELSE IF( nRange >= 4 .AND. Remainder /= 0 ) THEN
  DO rRange = 1, INT( DBLE( nRange ) / 4.D0 )
    StringsPRange(1,rRange) = TRIM( CharLabelPRange(1,4*rRange-3) )//", "//TRIM( CharLabelPRange(1,4*rRange-2) )//", "// &
    &                         TRIM( CharLabelPRange(1,4*rRange-1) )//", "//TRIM( CharLabelPRange(1,4*rRange-0) )//", "
  END DO
END IF
IF( nRange < 4 .OR. (nRange >= 4 .AND. Remainder /= 0) ) THEN
  IF( Remainder == 1 ) THEN
    StringsPRange(1,INT( DBLE( nRange ) / 4.D0 ) + 1) = TRIM( CharLabelPRange(1,nRange-Remainder+1) )
  ELSE IF( Remainder == 2 ) THEN
    StringsPRange(1,INT( DBLE( nRange ) / 4.D0 ) + 1) = TRIM( CharLabelPRange(1,nRange-Remainder+1) )//", "// &
    &                                                   TRIM( CharLabelPRange(1,nRange-Remainder+2) )
  ELSE IF( Remainder == 3 ) THEN
    StringsPRange(1,INT( DBLE( nRange ) / 4.D0 ) + 1) = TRIM( CharLabelPRange(1,nRange-Remainder+1) )//", "// &
    &                                                   TRIM( CharLabelPRange(1,nRange-Remainder+2) )//", "// &
    &                                                   TRIM( CharLabelPRange(1,nRange-Remainder+3) )
  END IF
END IF
StringsText(71,1) = "Attractive Range Points:"
StringsText(72,1) = "Reduced Temperature: "//TRIM( CharLabel(67,1) )
StringsText(73,1) = "Potential Files: "//TRIM( CharLabel(68,1) )
IF( PerturbationCoefficientLogical ) THEN
  StringsText(74,1) = "Perturbation coefficients computed."
ELSE IF( .NOT. PerturbationCoefficientLogical ) THEN
  StringsText(74,1) = "Perturbation coefficients not computed."
END IF
StringsText(75,1) = "Minimum Number of Blocks (Block Averaging): "//TRIM( CharLabel(69,1) )
StringsText(76,1) = "Maximum Number of Blocks (Block Averaging): "//TRIM( CharLabel(70,1) )
StringsText(78,1) = "Seed Value: "//TRIM( CharLabel(72,1) )
StringsText(79,1) = "Linked Lists: "//TRIM( CharLabel(73,1) )

! String name (subtitle)
IF( nComponents > 1 ) THEN
  StringsSubtitle(1) = "PARAMETERS OF THE COMPONENTS"
ELSE IF( nComponents == 1 ) THEN
  StringsSubtitle(1) = "PARAMETERS OF THE COMPONENT"
END IF
StringsSubtitle(2) = "PARAMETERS OF THE SYSTEM"
StringsSubtitle(3) = "PARAMETERS OF THE SIMULATION"
StringsSubtitle(4) = "PARAMETERS OF THE INITIAL CONFIGURATION"
StringsSubtitle(5) = "OTHER PARAMETERS"
StringsSubtitle(6) = "PARAMETERS OF THE POTENTIAL"

! String size
DO iString = 1, 6
  StringsSizeHeader(iString) = ( 70.D0 - DBLE( LEN( TRIM( StringsHeader(iString) ) ) ) ) * 0.5D0
END DO
DO iString = 1, 79
  DO cComponent = 1, nComponents
    StringsSizeText(iString,cComponent) = ( 69 - LEN( TRIM( StringsText(iString,cComponent) ) ) )
  END DO
END DO
DO iString = 1, 6
  StringsSizeSubtitle(iString) = ( 70.D0 - DBLE( LEN( TRIM( StringsSubtitle(iString) ) ) ) ) * 0.5D0
END DO

! File inquest
INQUIRE( File= "Simulation_Log.txt", EXIST= FileExist )

! Simulation log
IF( .NOT. FileExist ) THEN
  OPEN( Unit= 95, File= "Simulation_Log.txt" )
  WRITE( 95, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(1) ) )//TRIM( StringsHeader(1) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(2) ) )//TRIM( StringsHeader(2) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(2) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(3) ) )//TRIM( StringsHeader(3) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(3) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(4) ) )//TRIM( StringsHeader(4) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(4) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(5) ) )//TRIM( StringsHeader(5) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(5) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeHeader(6) ) )//TRIM( StringsHeader(6) )//REPEAT( " ", &
  &                                       CEILING( StringsSizeHeader(6) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VL//REPEAT( CH_HS, 70 )//CH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(1,1) )//REPEAT( " ", NINT( StringsSizeText(1,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(2,1) )//REPEAT( " ", NINT( StringsSizeText(2,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(3,1) )//REPEAT( " ", NINT( StringsSizeText(3,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(4,1) )//REPEAT( " ", NINT( StringsSizeText(4,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(5,1) )//REPEAT( " ", NINT( StringsSizeText(5,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                                       LEN( TRIM( StringsSubtitle(1) ) ) )//SS_UR//REPEAT( " ", &
  &                                       CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
  IF( nComponents > 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_VL//TRIM( StringsSubtitle(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_BL//REPEAT( SS_HS, &
    &                   LEN( TRIM( StringsSubtitle(1) ) ) )// SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    DO cComponent = 1, nComponents
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(6,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(6,cComponent) ) + 2 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(77,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(77,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(7,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(7,cComponent) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(8,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(8,cComponent) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(9,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(9,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(10,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(10,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(11,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(11,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(12,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(12,cComponent) ) + 2 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(13,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(13,cComponent) ) + 4 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    END DO
  ELSE IF( nComponents == 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_VL//TRIM( StringsSubtitle(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_BL//REPEAT( SS_HS, &
    &                   LEN( TRIM( StringsSubtitle(1) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(77,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(77,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(7,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(7,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(8,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(8,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(9,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(9,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(10,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(10,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(11,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(11,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(12,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(12,1) ) + 2 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(13,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(13,1) ) + 4 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(2) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_VL//TRIM( StringsSubtitle(2) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(2) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(2) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(14,1) )//REPEAT( " ", NINT( StringsSizeText(14,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(15,1) )//REPEAT( " ", NINT( StringsSizeText(15,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(16,1) )//REPEAT( " ", NINT( StringsSizeText(16,1) ) + 4 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(17,1) )//REPEAT( " ", NINT( StringsSizeText(17,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(18,1) )//REPEAT( " ", NINT( StringsSizeText(18,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(19,1) )//REPEAT( " ", NINT( StringsSizeText(19,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )// &
  &                   SS_UL//REPEAT( " ", LEN( TRIM( StringsText(21,1) ) ) )// &
  &                   SS_UR//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )//SS_VS//TRIM( StringsText(21,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(20,1) )//REPEAT( " ", 1 )//SS_VS//TRIM( StringsText(22,1) )// &
  &                   SS_VS//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )//SS_VS//TRIM( StringsText(23,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )// &
  &                   SS_BL//REPEAT( " ", LEN( TRIM( StringsText(21,1) ) ) )// &
  &                   SS_BR//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(3) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_VL//TRIM( StringsSubtitle(3) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(3) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(3) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(24,1) )//REPEAT( " ", NINT( StringsSizeText(24,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(25,1) )//REPEAT( " ", NINT( StringsSizeText(25,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(26,1) )//REPEAT( " ", NINT( StringsSizeText(26,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(27,1) )//REPEAT( " ", NINT( StringsSizeText(27,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(28,1) )//REPEAT( " ", NINT( StringsSizeText(28,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(29,1) )//REPEAT( " ", NINT( StringsSizeText(29,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(30,1) )//REPEAT( " ", NINT( StringsSizeText(30,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(31,1) )//REPEAT( " ", NINT( StringsSizeText(31,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(58,1) )//REPEAT( " ", NINT( StringsSizeText(58,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(59,1) )//REPEAT( " ", NINT( StringsSizeText(59,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(32,1) )//REPEAT( " ", NINT( StringsSizeText(32,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(33,1) )//REPEAT( " ", NINT( StringsSizeText(33,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(34,1) )//REPEAT( " ", NINT( StringsSizeText(34,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(35,1) )//REPEAT( " ", NINT( StringsSizeText(35,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(60,1) )//REPEAT( " ", NINT( StringsSizeText(60,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(36,1) )//REPEAT( " ", NINT( StringsSizeText(36,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(37,1) )//REPEAT( " ", NINT( StringsSizeText(37,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(38,1) )//REPEAT( " ", NINT( StringsSizeText(38,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(61,1) )//REPEAT( " ", NINT( StringsSizeText(61,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(65,1) )//REPEAT( " ", NINT( StringsSizeText(65,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(66,1) )//REPEAT( " ", NINT( StringsSizeText(66,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(67,1) )//REPEAT( " ", NINT( StringsSizeText(67,1) ) + 1 )//CH_VS
  IF( LatticeReductionTypeLogical(1) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(68,1) )//REPEAT( " ", NINT( StringsSizeText(68,1) ) )//CH_VS
  ELSE IF( LatticeReductionTypeLogical(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(68,1) )//REPEAT( " ", NINT( StringsSizeText(68,1) ) + 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(4) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_VL//TRIM( StringsSubtitle(4) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(4) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(4) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(39,1) )//REPEAT( " ", NINT( StringsSizeText(39,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(40,1) )//REPEAT( " ", NINT( StringsSizeText(40,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(41,1) )//REPEAT( " ", NINT( StringsSizeText(41,1) ) + 2 )//CH_VS
  IF( ConfigurationSelection(4) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(42,1) )//REPEAT( " ", NINT( StringsSizeText(42,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(44,1) )//REPEAT( " ", NINT( StringsSizeText(44,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(45,1) )//REPEAT( " ", NINT( StringsSizeText(45,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(46,1) )//REPEAT( " ", NINT( StringsSizeText(46,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(47,1) )//REPEAT( " ", NINT( StringsSizeText(47,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(48,1) )//REPEAT( " ", NINT( StringsSizeText(48,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(62,1) )//REPEAT( " ", NINT( StringsSizeText(62,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(49,1) )//REPEAT( " ", NINT( StringsSizeText(49,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(50,1) )//REPEAT( " ", NINT( StringsSizeText(50,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(51,1) )//REPEAT( " ", NINT( StringsSizeText(51,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(52,1) )//REPEAT( " ", NINT( StringsSizeText(52,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(53,1) )//REPEAT( " ", NINT( StringsSizeText(53,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(63,1) )//REPEAT( " ", NINT( StringsSizeText(63,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(64,1) )//REPEAT( " ", NINT( StringsSizeText(64,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(54,1) )//REPEAT( " ", NINT( StringsSizeText(54,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(6) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_VL//TRIM( StringsSubtitle(6) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(6) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(6) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(69,1) )//REPEAT( " ", NINT( StringsSizeText(69,1) ) + 2 )//CH_VS
  IF( PotentialTypeLogical(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(70,1) )//REPEAT( " ", NINT( StringsSizeText(70,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(71,1) )//REPEAT( " ", NINT( StringsSizeText(71,1) ) - 1 )//CH_VS
    IF( nRange == 4 ) THEN
      DO rRange = 1, INT( DBLE( nRange ) / 4.D0 )
        WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 )//CH_VS
      END DO
    ELSE IF( nRange /= 4 ) THEN
      DO rRange = 1, INT( DBLE( nRange ) / 4.D0 )
        IF( rRange == INT( DBLE( nRange ) / 4.D0 ) .AND. Remainder == 0 ) THEN
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 )//CH_VS
        ELSE
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 35 )//CH_VS
        END IF
      END DO
    END IF
    IF( Remainder /= 0 ) THEN
      rRange = INT( DBLE( nRange ) / 4.D0 ) + 1
      Remainder = 4 - Remainder
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 + (8 * Remainder) )//CH_VS
    END IF
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(72,1) )//REPEAT( " ", NINT( StringsSizeText(72,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(73,1) )//REPEAT( " ", NINT( StringsSizeText(73,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(74,1) )//REPEAT( " ", NINT( StringsSizeText(74,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(75,1) )//REPEAT( " ", NINT( StringsSizeText(75,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(76,1) )//REPEAT( " ", NINT( StringsSizeText(76,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(5) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_VL//TRIM( StringsSubtitle(5) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(5) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(5) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(55,1) )//REPEAT( " ", NINT( StringsSizeText(55,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(56,1) )//REPEAT( " ", NINT( StringsSizeText(56,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(78,1) )//REPEAT( " ", NINT( StringsSizeText(78,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(79,1) )//REPEAT( " ", NINT( StringsSizeText(79,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(57,1) )//REPEAT( " ", NINT( StringsSizeText(57,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  CLOSE( 95 )
! Simulation log (appending)
ELSE IF( FileExist ) THEN
  OPEN( Unit= 95, File= "Simulation_Log.txt", Position= "APPEND" )
  WRITE( 95, "(G0)" ) REPEAT( " ", 1 )//REPEAT( C_FUL, 70 )//REPEAT( " ", 1 )
  WRITE( 95, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(1,1) )//REPEAT( " ", NINT( StringsSizeText(1,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(2,1) )//REPEAT( " ", NINT( StringsSizeText(2,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(3,1) )//REPEAT( " ", NINT( StringsSizeText(3,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(4,1) )//REPEAT( " ", NINT( StringsSizeText(4,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(5,1) )//REPEAT( " ", NINT( StringsSizeText(5,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(1) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
  IF( nComponents > 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_VL//TRIM( StringsSubtitle(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_BL//REPEAT( SS_HS, &
    &                   LEN( TRIM( StringsSubtitle(1) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    DO cComponent = 1, nComponents
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(6,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(6,cComponent) ) + 2 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(77,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(77,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(7,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(7,cComponent) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(8,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(8,cComponent) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(9,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(9,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(10,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(10,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(11,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(11,cComponent) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(12,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(12,cComponent) ) + 2 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(13,cComponent) )//REPEAT( " ", &
      &                   NINT( StringsSizeText(13,cComponent) ) + 4 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    END DO
  ELSE IF( nComponents == 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_VL//TRIM( StringsSubtitle(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(1) - 1 ) )//SS_BL//REPEAT( SS_HS, &
    &                   LEN( TRIM( StringsSubtitle(1) ) ) )// SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(77,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(77,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(7,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(7,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(8,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(8,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(9,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(9,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(10,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(10,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(11,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(11,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(12,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(12,1) ) + 2 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(13,1) )//REPEAT( " ", &
    &                   NINT( StringsSizeText(13,1) ) + 4 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(2) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_VL//TRIM( StringsSubtitle(2) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(2) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(2) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(2) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(14,1) )//REPEAT( " ", NINT( StringsSizeText(14,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(15,1) )//REPEAT( " ", NINT( StringsSizeText(15,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(16,1) )//REPEAT( " ", NINT( StringsSizeText(16,1) ) + 4 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(17,1) )//REPEAT( " ", NINT( StringsSizeText(17,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(18,1) )//REPEAT( " ", NINT( StringsSizeText(18,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(19,1) )//REPEAT( " ", NINT( StringsSizeText(19,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )// &
  &                   SS_UL//REPEAT( " ", LEN( TRIM( StringsText(21,1) ) ) )// &
  &                   SS_UR//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )//SS_VS//TRIM( StringsText(21,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(20,1) )//REPEAT( " ", 1 )//SS_VS//TRIM( StringsText(22,1) )// &
  &                   SS_VS//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )//SS_VS//TRIM( StringsText(23,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( StringsText(20,1) ) ) + 1 )// &
  &                   SS_BL//REPEAT( " ", LEN( TRIM( StringsText(21,1) ) ) )// &
  &                   SS_BR//REPEAT( " ", 67 - ( LEN( TRIM( StringsText(20,1) ) ) + LEN( TRIM( StringsText(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(3) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_VL//TRIM( StringsSubtitle(3) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(3) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(3) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(3) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(24,1) )//REPEAT( " ", NINT( StringsSizeText(24,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(25,1) )//REPEAT( " ", NINT( StringsSizeText(25,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(26,1) )//REPEAT( " ", NINT( StringsSizeText(26,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(27,1) )//REPEAT( " ", NINT( StringsSizeText(27,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(28,1) )//REPEAT( " ", NINT( StringsSizeText(28,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(29,1) )//REPEAT( " ", NINT( StringsSizeText(29,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(30,1) )//REPEAT( " ", NINT( StringsSizeText(30,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(31,1) )//REPEAT( " ", NINT( StringsSizeText(31,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(58,1) )//REPEAT( " ", NINT( StringsSizeText(58,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(59,1) )//REPEAT( " ", NINT( StringsSizeText(59,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(32,1) )//REPEAT( " ", NINT( StringsSizeText(32,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(33,1) )//REPEAT( " ", NINT( StringsSizeText(33,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(34,1) )//REPEAT( " ", NINT( StringsSizeText(34,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(35,1) )//REPEAT( " ", NINT( StringsSizeText(35,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(60,1) )//REPEAT( " ", NINT( StringsSizeText(60,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(36,1) )//REPEAT( " ", NINT( StringsSizeText(36,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(37,1) )//REPEAT( " ", NINT( StringsSizeText(37,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(38,1) )//REPEAT( " ", NINT( StringsSizeText(38,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(61,1) )//REPEAT( " ", NINT( StringsSizeText(61,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(65,1) )//REPEAT( " ", NINT( StringsSizeText(65,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(66,1) )//REPEAT( " ", NINT( StringsSizeText(66,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(67,1) )//REPEAT( " ", NINT( StringsSizeText(67,1) ) + 1 )//CH_VS
  IF( LatticeReductionTypeLogical(1) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(68,1) )//REPEAT( " ", NINT( StringsSizeText(68,1) ) )//CH_VS
  ELSE IF( LatticeReductionTypeLogical(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(68,1) )//REPEAT( " ", NINT( StringsSizeText(68,1) ) + 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(4) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_VL//TRIM( StringsSubtitle(4) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(4) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(4) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(4) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(39,1) )//REPEAT( " ", NINT( StringsSizeText(39,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(40,1) )//REPEAT( " ", NINT( StringsSizeText(40,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(41,1) )//REPEAT( " ", NINT( StringsSizeText(41,1) ) + 2 )//CH_VS
  IF( ConfigurationSelection(4) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(42,1) )//REPEAT( " ", NINT( StringsSizeText(42,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(44,1) )//REPEAT( " ", NINT( StringsSizeText(44,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(45,1) )//REPEAT( " ", NINT( StringsSizeText(45,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(46,1) )//REPEAT( " ", NINT( StringsSizeText(46,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(47,1) )//REPEAT( " ", NINT( StringsSizeText(47,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(48,1) )//REPEAT( " ", NINT( StringsSizeText(48,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(62,1) )//REPEAT( " ", NINT( StringsSizeText(62,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(49,1) )//REPEAT( " ", NINT( StringsSizeText(49,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(50,1) )//REPEAT( " ", NINT( StringsSizeText(50,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(51,1) )//REPEAT( " ", NINT( StringsSizeText(51,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(52,1) )//REPEAT( " ", NINT( StringsSizeText(52,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(53,1) )//REPEAT( " ", NINT( StringsSizeText(53,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(63,1) )//REPEAT( " ", NINT( StringsSizeText(63,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(64,1) )//REPEAT( " ", NINT( StringsSizeText(64,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(54,1) )//REPEAT( " ", NINT( StringsSizeText(54,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(6) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_VL//TRIM( StringsSubtitle(6) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(6) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(6) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(6) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(69,1) )//REPEAT( " ", NINT( StringsSizeText(69,1) ) + 2 )//CH_VS
  IF( PotentialTypeLogical(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(70,1) )//REPEAT( " ", NINT( StringsSizeText(70,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(71,1) )//REPEAT( " ", NINT( StringsSizeText(71,1) ) - 1 )//CH_VS
    IF( nRange == 4 ) THEN
      DO rRange = 1, INT( DBLE( nRange ) / 4.D0 )
        WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 )//CH_VS
      END DO
    ELSE IF( nRange /= 4 ) THEN
      DO rRange = 1, INT( DBLE( nRange ) / 4.D0 )
        IF( rRange == INT( DBLE( nRange ) / 4.D0 ) .AND. Remainder == 0 ) THEN
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 )//CH_VS
        ELSE
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 35 )//CH_VS
        END IF
      END DO
    END IF
    IF( Remainder /= 0 ) THEN
      rRange = INT( DBLE( nRange ) / 4.D0 ) + 1
      Remainder = 4 - Remainder
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( StringsPRange(1,rRange) )//REPEAT( " ", 36 + (8 * Remainder) )//CH_VS
    END IF
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(72,1) )//REPEAT( " ", NINT( StringsSizeText(72,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(73,1) )//REPEAT( " ", NINT( StringsSizeText(73,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(74,1) )//REPEAT( " ", NINT( StringsSizeText(74,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(75,1) )//REPEAT( " ", NINT( StringsSizeText(75,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( StringsText(76,1) )//REPEAT( " ", NINT( StringsSizeText(76,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_UL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(5) ) ) )//SS_UR//REPEAT( " ", CEILING( StringsSizeSubtitle(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_VL//TRIM( StringsSubtitle(5) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( StringsSizeSubtitle(5) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( StringsSizeSubtitle(5) - 1 ) )//SS_BL//REPEAT( SS_HS, &
  &                   LEN( TRIM( StringsSubtitle(5) ) ) )//SS_BR//REPEAT( " ", CEILING( StringsSizeSubtitle(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(55,1) )//REPEAT( " ", NINT( StringsSizeText(55,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(56,1) )//REPEAT( " ", NINT( StringsSizeText(56,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(78,1) )//REPEAT( " ", NINT( StringsSizeText(78,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(79,1) )//REPEAT( " ", NINT( StringsSizeText(79,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( StringsText(57,1) )//REPEAT( " ", NINT( StringsSizeText(57,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  CLOSE( 95 )
END IF

RETURN

END SUBROUTINE SimulationLog
