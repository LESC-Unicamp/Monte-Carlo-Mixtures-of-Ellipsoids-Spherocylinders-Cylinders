! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!  This module initialize common system variables (number of particles, reduced number density,   !
!     reduced temperature etc.), molecular properties (geometry and dimensions), Monte Carlo      !
!  parameters (ensemble type, total number of cycles, number of equilibration cycles etc.), and   !
!                                 potential parameters (if any).                                  !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
! which results will be written out in external files and to enable post-processing subroutines.  !
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
MODULE InitializeVariables

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
!                                     VARIABLE INITIALIZATION                                     !
!              All variables should be first specified in an input file (*.ini file)              !
!                       We provided an example .ini file to guide the user.                       !
!                     Please also check our README file for more information.                     !
! *********************************************************************************************** !

CONTAINS

! *********************************************************************************************** !
!                            Initialization of Monte Carlo parameters                             !
! *********************************************************************************************** !
SUBROUTINE MonteCarloVariables(  )

IMPLICIT NONE

! Simulation parameters
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )

! Total number of cycles
READ( 10, * ) Dummy, MaxSimulationCycles
! Condition
IF( MaxSimulationCycles < 1 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of cycles [", MaxSimulationCycles, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Number of equilibration cycles
READ( 10, * ) Dummy, nEquilibrationCycles
! Condition 1
IF( nEquilibrationCycles >= MaxSimulationCycles ) THEN
  WRITE( *, "(6G0)" ) "The number of equilibration cycles [", nEquilibrationCycles, "] cannot be greater than or equal to the ", &
  &                   "maximum number of cycles [", MaxSimulationCycles, "]. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( nEquilibrationCycles < 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of equilibration cycles [", nEquilibrationCycles, "] cannot be a negative integer. Exiting... "
  CALL Exit(  )
END IF

! Saving frequency
READ( 10, * ) Dummy, nSavingFrequency
! Condition
IF( nSavingFrequency < 1 ) THEN
  WRITE( *, "(3G0)" ) "The saving frequency [", nSavingFrequency, "] cannot be a negative integer nor zero. Exiting... "
  CALL Exit(  )
END IF

! Adjustment frequency
READ( 10, * ) Dummy, nAdjustmentFrequency
! Condition
IF( nAdjustmentFrequency < 1 ) THEN
  WRITE( *, "(4G0)" ) "The adjustment frequency of the simulation [", nAdjustmentFrequency, "] cannot be a negative integer ", &
  &                   "nor zero. Exiting... "
  CALL Exit(  )
END IF

! Adjustment frequency (random configuration)
READ( 10, * ) Dummy, nAdjustmentRandomConfig
! Condition
IF( nAdjustmentRandomConfig < 1 ) THEN
  WRITE( *, "(4G0)" ) "The adjustment frequency of the random configuration [", nAdjustmentRandomConfig, "] cannot be a ", &
  &                   "negative integer nor zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum translational displacement
READ( 10, * ) Dummy, UserMaxTranslationalDisplacement
! Condition
IF( DABS( UserMaxTranslationalDisplacement - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum translational displacement of the simulation [", UserMaxTranslationalDisplacement, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum translational displacement (random configuration)
READ( 10, * ) Dummy, MaxTranslationalDisplacementRandomConfig
! Condition
IF( DABS( MaxTranslationalDisplacementRandomConfig - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum translational displacement of the random configuration [", &
  &                          MaxTranslationalDisplacementRandomConfig, "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum rotational displacement
READ( 10, * ) Dummy, UserMaxRotationalDisplacement
! Condition
IF( DABS( UserMaxRotationalDisplacement - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum rotational displacement of the simulation [", UserMaxRotationalDisplacement, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum rotational displacement (random configuration)
READ( 10, * ) Dummy, MaxAngularDisplacementRandomConfig
! Condition
IF( DABS( MaxAngularDisplacementRandomConfig - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum rotational displacement of the random configuration [", &
  &                          MaxAngularDisplacementRandomConfig, "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum isotropic volume scaling
READ( 10, * ) Dummy, UserMaxIsoVolumetricDisplacement
! Condition
IF( DABS( UserMaxIsoVolumetricDisplacement - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum isotropic volume scaling of the simulation [", UserMaxIsoVolumetricDisplacement, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum anisotropic volume scaling
READ( 10, * ) Dummy, UserMaxAnisoVolumetricDisplacement
! Condition
IF( DABS( UserMaxAnisoVolumetricDisplacement - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum anisotropic volume scaling of the simulation [", UserMaxAnisoVolumetricDisplacement, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum isotropic volume scaling (random configuration)
READ( 10, * ) Dummy, MaxIsoVolumetricDisplacementRandomConfig
! Condition
IF( DABS( MaxIsoVolumetricDisplacementRandomConfig - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum isotropic volume scaling of the random configuration [", &
  &                          MaxIsoVolumetricDisplacementRandomConfig, "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum anisotropic volume scaling (random configuration)
READ( 10, * ) Dummy, MaxAnisoVolumetricDisplacementRandomConfig
! Condition
IF( DABS( MaxAnisoVolumetricDisplacementRandomConfig - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum anisotropic volume scaling of the random configuration [", &
  &                          MaxAnisoVolumetricDisplacementRandomConfig, "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Minimum volumetric displacement (random configuration)
READ( 10, * ) Dummy, MinVolumetricDisplacementRandomConfig
! Condition
IF( DABS( MinVolumetricDisplacementRandomConfig ) >= DABS( MaxIsoVolumetricDisplacementRandomConfig ) .OR. &
&   DABS( MinVolumetricDisplacementRandomConfig ) >= DABS( MaxAnisoVolumetricDisplacementRandomConfig ) ) THEN
  WRITE( *, "(2(G0,G0.5,G0),G0.5,G0)" ) "The absolute minimum volume scaling of the random configuration [", &
  &                                     MinVolumetricDisplacementRandomConfig, "] cannot be greater than or equal to the ", &
  &                                     "absolute maximum isotropic [", MaxIsoVolumetricDisplacementRandomConfig, &
  &                                     "] or anisotropic [", MaxAnisoVolumetricDisplacementRandomConfig, &
  &                                     "] volume scaling. Exiting... "
  CALL Exit(  )
END IF

! Maximum box distortion
READ( 10 , * ) Dummy, MaxBoxDistortion
! Condition
IF( MaxBoxDistortion < 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum box distortion [", MaxBoxDistortion, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Maximum length ratio (box distortion)
READ( 10, * ) Dummy, BoxEdgeMaxRatio
! Condition
IF( BoxEdgeMaxRatio <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum linear distortion of the box [", BoxEdgeMaxRatio, "] cannot be less than ", &
  &                           "or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Maximum angle (box distortion)
READ( 10, * ) Dummy, BoxVectorMaxAngle
! Condition
IF( BoxVectorMaxAngle <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum angular distortion of the box [", BoxVectorMaxAngle, "] cannot be less than or ", &
  &                           "equal to 0°. Exiting... "
  CALL Exit(  )
END IF
! Convert to radians
BoxVectorMaxAngle = BoxVectorMaxAngle * cPi / 180.D0

! Lattice reduction method
READ( 10 , * ) Dummy, LatticeReductionType
CALL ToUpper( LatticeReductionType, LEN_TRIM( LatticeReductionType ), LatticeReductionType )
LatticeReductionTypeLogical = .FALSE.
! Lattice reduction: Gottwald method
IF( LatticeReductionType == "FBM" ) THEN
  LatticeReductionTypeLogical(1) = .TRUE.
! Lattice reduction: Lenstra-Lenstra-Lovász method
ELSE IF( LatticeReductionType == "LLL" ) THEN
  LatticeReductionTypeLogical(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(4G0)" ) "The user-defined [", TRIM( LatticeReductionType ), "] is not an available lattice reduction method. ", &
  &                   "Exiting... "
  CALL Exit(  )
END IF

! Simulation ensemble
READ( 10, * ) Dummy, EnsembleMC
CALL ToUpper( EnsembleMC, LEN_TRIM( EnsembleMC ), EnsembleMC )
! Condition
IF( EnsembleMC /= "NVT" .AND. EnsembleMC /= "NPT" ) THEN
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( EnsembleMC ), "] is not an available simulation ensemble. Exiting... "
  CALL Exit(  )
END IF
WRITE( *, "(3G0)" ) "Monte Carlo ensemble: [", EnsembleMC, "]"
WRITE( *, "(G0)" ) " "

CLOSE( 10 )

RETURN

END SUBROUTINE MonteCarloVariables

! *********************************************************************************************** !
!                           Initialization of Inquiry/Control variables                           !
! *********************************************************************************************** !
SUBROUTINE VariablesInquest(  )

IMPLICIT NONE

! Inquiry variables
OPEN( Unit= 10, File= "ini_control.ini", Action= "READ" )

! Trajectory inquiry
READ( 10, * ) Dummy, TrajectoryInquiry
CALL ToUpper( TrajectoryInquiry, LEN_TRIM( TrajectoryInquiry ), TrajectoryInquiry )
! Transforms characters into logical variables
IF( TrajectoryInquiry == "Y" ) THEN
  TrajectoryLogical = .TRUE.
ELSE
  TrajectoryLogical = .FALSE.
END IF

! Fixed/Random seed inquiry
READ( 10, * ) Dummy, SeedTypeInquiry
CALL ToUpper( SeedTypeInquiry, LEN_TRIM( SeedTypeInquiry ), SeedTypeInquiry )
! Transforms characters into logical variables
IF( SeedTypeInquiry == "Y" ) THEN
  FixedSeedLogical = .TRUE.
ELSE
  FixedSeedLogical = .FALSE.
END IF

! Backup file inquiry
READ( 10, * ) Dummy, BackupFileInquiry
CALL ToUpper( BackupFileInquiry, LEN_TRIM( BackupFileInquiry ), BackupFileInquiry )
! Transforms characters into logical variables
IF( BackupFileInquiry == "Y" ) THEN
  BackupFileLogical = .TRUE.
ELSE
  BackupFileLogical = .FALSE.
END IF

! Skip
READ( 10, * ) Dummy, Dummy

! Potential type
READ( 10, * ) Dummy, PotentialType
CALL ToUpper( PotentialType, LEN_TRIM( PotentialType ), PotentialType )
! Transforms characters into logical variables
PotentialTypeLogical = .FALSE.
IF( PotentialType == "HARDCORE" ) THEN
  PotentialTypeLogical(1) = .TRUE.
ELSE IF( PotentialType == "SQUAREWELL" ) THEN
  PotentialTypeLogical(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( PotentialType ), "] is not an available force field. Exiting..."
  CALL Exit(  )
END IF

! Cell list inquiry
READ( 10, * ) Dummy, CellListLogical
IF( .NOT. CellListLogical ) CellListControl = .FALSE.
IF( .NOT. CellListLogical ) CellListControlPotential = .FALSE.

! Potential inquiry
IF( .NOT. PotentialTypeLogical(1) ) THEN
  READ( 10, * ) Dummy, PotentialEnergyInquiry
  CALL ToUpper( PotentialEnergyInquiry, LEN_TRIM( PotentialEnergyInquiry ), PotentialEnergyInquiry )
  ! Transforms characters into logical variables
  IF( PotentialEnergyInquiry == "Y" ) THEN
    PotentialEnergyLogical = .TRUE.
  ELSE
    PotentialEnergyLogical = .FALSE.
  END IF
END IF

! TPT coefficients inquiry
IF( .NOT. PotentialTypeLogical(1) ) THEN
  READ( 10, * ) Dummy, PerturbationCoefficientInquiry
  CALL ToUpper( PerturbationCoefficientInquiry, LEN_TRIM( PerturbationCoefficientInquiry ), PerturbationCoefficientInquiry )
  ! Transforms characters into logical variables
  IF( PerturbationCoefficientInquiry == "Y" ) THEN
    PerturbationCoefficientLogical = .TRUE.
  ELSE
    PerturbationCoefficientLogical = .FALSE.
  END IF
END IF

CLOSE( 10 )

RETURN

END SUBROUTINE VariablesInquest

! *********************************************************************************************** !
!                               Initialization of common variables                                !
! *********************************************************************************************** !
SUBROUTINE GeneralVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: cComponent, pParticle ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: TotalMolarFraction ! Total molar fraction
REAL( Kind= Real64 ) :: CubicRootFloat     ! Checks the cube root for the number of particles, N, in each crystalline structure

! System variables
OPEN( Unit= 100, File= "ini_system.ini", Action= "READ" )

! Packing fraction
READ( 100, * ) Dummy, PackingFraction
! Condition 1
IF( PackingFraction <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction [", PackingFraction, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( PackingFraction > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction [", PackingFraction, "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF

! Number of components
READ( 100, * ) Dummy, nComponents
! Condition
IF( nComponents < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of components [", nComponents, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Allocation
ALLOCATE( cDiameter(nComponents), cLength(nComponents), cAspectRatio(nComponents), cMolecularVolume(nComponents) )
ALLOCATE( cParticles(0:nComponents), cMolarFraction(nComponents), cNumberDensity(nComponents) )
ALLOCATE( cCircumscribingSphereDiameter(nComponents) )
ALLOCATE( cDiameterEquivalentSphere(nComponents) )
ALLOCATE( SphericalComponentInquiry(nComponents) )
ALLOCATE( SphericalComponentLogical(nComponents) )

! Component sphericity
READ( 100, * ) Dummy, SphericalComponentInquiry
DO cComponent = 1, nComponents
  CALL ToUpper( SphericalComponentInquiry(cComponent), LEN_TRIM( SphericalComponentInquiry(cComponent) ), &
  &             SphericalComponentInquiry(cComponent) )
  IF( SphericalComponentInquiry(cComponent) == "T" ) THEN
    SphericalComponentLogical(cComponent) = .TRUE.
  ELSE
    SphericalComponentLogical(cComponent) = .FALSE.
  END IF
END DO

! Diameter of component i
READ( 100, * ) Dummy, cDiameter
! Condition
DO cComponent = 1, nComponents
  IF( cDiameter(cComponent) <= 0.D0 ) THEN
    WRITE( *, "(3G0,G0.5,2G0)" ) "The diameter of component #", cComponent, " [", cDiameter(cComponent), "] cannot be less ", &
    &                            "than or equal to 0. Exiting... "
    CALL Exit(  )
  END IF
END DO

! Length of component i
READ( 100, * ) Dummy, cLength
! Condition
DO cComponent = 1, nComponents
  IF( GeometryType(1) .OR. GeometryType(3) ) THEN
    IF( cLength(cComponent) <= 0.D0 ) THEN
      WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", cComponent, " [", cLength(cComponent), "] cannot be less than ", &
      &                            "or equal to 0. Exiting... "
      CALL Exit(  )
    END IF
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    IF( cLength(cComponent) < 0.D0 ) THEN
      WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", cComponent, " [", cLength(cComponent), "] cannot be less than ", &
      &                            "0. Exiting... "
      CALL Exit(  )
    END IF
  END IF
END DO

! Check if component is spherical
DO cComponent = 1, nComponents
  IF( GeometryType(1) .AND. DABS( cLength(cComponent) - cDiameter(cComponent) ) >= EPSILON( 1.D0 ) .AND. &
  &   SphericalComponentLogical(cComponent) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", cComponent, " [", cLength(cComponent), "] cannot be ", &
    &                               "different from its diameter [", cDiameter(cComponent), "] when the component is ", &
    &                               "considered spherical under the ellipsoidal framework. Exiting..."
    CALL Exit(  )
  ELSE IF( GeometryType(2) .AND. DABS( cLength(cComponent) - 0.D0 ) >= EPSILON( 1.D0 ) .AND. &
    &      SphericalComponentLogical(cComponent) ) THEN
    WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", cComponent, " [", cLength(cComponent), "] must be 0 when the ", &
    &                            "component is considered spherical under the spherocylindrical framework. Exiting..."
    CALL Exit(  )
  ELSE IF( GeometryType(3) .AND. DABS( cLength(cComponent) - cDiameter(cComponent) ) >= EPSILON( 1.D0 ) .AND. &
    &      SphericalComponentLogical(cComponent) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", cComponent, " [", cLength(cComponent), "] cannot be ", &
    &                               "different from its diameter [", cDiameter(cComponent), "] when the component is ", &
    &                               "considered spherical under the cylindrical framework. Exiting..."
    CALL Exit(  )
  END IF
  IF( GeometryType(1) .AND. DABS( cLength(cComponent) - cDiameter(cComponent) ) < EPSILON( 1.D0 ) .AND. &
  &   .NOT. SphericalComponentLogical(cComponent) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", cComponent, " [", cLength(cComponent), "] is equal to its ", &
    &                               "diameter [", cDiameter(cComponent), "] even though the component is not considered ", &
    &                               "spherical under the ellipsoidal framework. Exiting..."
    CALL Exit(  )
  ELSE IF( GeometryType(2) .AND. DABS( cLength(cComponent) - 0.D0 ) < EPSILON( 1.D0 ) .AND. &
    &      .NOT. SphericalComponentLogical(cComponent) ) THEN
    WRITE( *, "(3G0,G0.5,3G0)" ) "The length of component #", cComponent, " [", cLength(cComponent), "] is equal to 0 even ", &
    &                            "though the component is not considered spherical under the spherocylindrical framework. ", &
    &                            "Exiting..."
    CALL Exit(  )
  END IF
END DO

! Molar fraction of component i
READ( 100, * ) Dummy, cMolarFraction
! Condition
DO cComponent = 1, nComponents
  IF( cMolarFraction(cComponent) < 0.D0 ) THEN
    WRITE( *, "(3G0,G0.5,G0)" ) "The molar fraction of component #", cComponent, " [", cMolarFraction(cComponent), &
    &                           "] cannot be less than 0. Exiting... "
    CALL Exit(  )
  END IF
END DO

! Aspect ratio
DO cComponent = 1, nComponents
  cAspectRatio(cComponent) = cLength(cComponent) / cDiameter(cComponent)
END DO

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"COMPONENT DETAILS"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Number of Components: ", nComponents
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  IF( SphericalComponentLogical( cComponent ) ) THEN
    WRITE( *, "(G0,G0,G0)" ) "Is Component #", cComponent, " Spherical? [YES]"
  ELSE
    WRITE( *, "(G0,G0,G0)" ) "Is Component #", cComponent, " Spherical? [NO]"
  END IF
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Diameter of Component #", cComponent, ": ", cDiameter(cComponent), "Å"
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Length of Component #", cComponent, ": ", cLength(cComponent), "Å"
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Aspect Ratio of Component #", cComponent, ": ", cAspectRatio(cComponent)
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Molar Fraction of Component #", cComponent, ": ", cMolarFraction(cComponent)
END DO
TotalMolarFraction = 0.D0
DO cComponent = 1, nComponents
  TotalMolarFraction = TotalMolarFraction + cMolarFraction(cComponent)
END DO
WRITE( *, "(G0)" ) " "
IF( DABS( TotalMolarFraction - 1.D0 ) >= EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0)" ) "Total molar fraction not equal to unit! Recalculating... "
  WRITE( *, "(G0)" ) " "
  DO cComponent = 1, nComponents
    cMolarFraction(cComponent) = cMolarFraction(cComponent) / TotalMolarFraction
    WRITE( *, "(G0,G0,G0,G0.5)" ) "Normalized Molar Fraction of Component #", cComponent, ": ", cMolarFraction(cComponent)
  END DO
  WRITE( *, "(G0)" ) " "
  TotalMolarFraction = 0.D0
  DO cComponent = 1, nComponents
    TotalMolarFraction = TotalMolarFraction + cMolarFraction(cComponent)
  END DO
END IF

! Number of particles
READ( 100, * ) Dummy, nParticles
! Condition
IF( nParticles < 2 ) THEN
  WRITE( *, "(3G0)" ) "The number of particles [", nParticles, "] cannot be less than 2. Exiting... "
  CALL Exit(  )
END IF
! Allocation
ALLOCATE( pComponents(nParticles) )

! Number of particles of component i
cParticles = 0
DO cComponent = 1, nComponents
  cParticles(cComponent) = NINT( cMolarFraction(cComponent) * DBLE(nParticles) )
END DO

! Component index of a particle
DO cComponent = 1, nComponents
  DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
    pComponents(pParticle) = cComponent
  END DO
END DO

! Particle volume
DO cComponent = 1, nComponents
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cMolecularVolume(cComponent) = (cPi / 6.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cLength(cComponent)
    ELSE
      cMolecularVolume(cComponent) = (cPi / 6.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cDiameter(cComponent)
    END IF
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cMolecularVolume(cComponent) = (cPi / 6.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cDiameter(cComponent) + &
      &                              (cPi / 4.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cLength(cComponent)
    ELSE
      cMolecularVolume(cComponent) = (cPi / 6.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cDiameter(cComponent)
    END IF
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cMolecularVolume(cComponent) = (cPi / 4.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cLength(cComponent)
    ELSE
      cMolecularVolume(cComponent) = (cPi / 6.D0) * cDiameter(cComponent) * cDiameter(cComponent) * cDiameter(cComponent)
    END IF
  END IF
END DO

! Total particle volume
TotalParticleVolume = 0.D0
DO cComponent = 1, nComponents
  TotalParticleVolume = TotalParticleVolume + DBLE( cParticles(cComponent) ) * cMolecularVolume(cComponent)
END DO

! Box volume
BoxVolume = TotalParticleVolume / PackingFraction

! Number density
TotalNumberDensity = 0.D0
DO cComponent = 1, nComponents
  cNumberDensity(cComponent) = DBLE( cParticles(cComponent) ) / BoxVolume
  TotalNumberDensity = TotalNumberDensity + cNumberDensity(cComponent)
END DO

! Summary
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Number of Particles of Component #", cComponent, ": ", cParticles(cComponent)
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Molecular Volume of Component #", cComponent, ": ", cMolecularVolume(cComponent), "Å³"
END DO
WRITE( *, "(G0)" ) " "
DO cComponent = 1, nComponents
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Number Density of Component #", cComponent, ": ", cNumberDensity(cComponent), "Å⁻³"
END DO
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"GLOBAL DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( SUM( cParticles ) == nParticles ) THEN
  WRITE( *, "(G0,G0)" ) "Number of Particles: ", nParticles
ELSE IF( SUM( cParticles ) /= nParticles ) THEN
  WRITE( *, "(G0,G0,G0,G0,G0,G0,G0)" ) "Molar-based number of particles (", SUM( cParticles ), ") not equal to the ", &
  &                                    "user-defined number of particles (", nParticles, "). Overwriting... "
  WRITE( *, "(G0)" ) " "
  nParticles = SUM( cParticles )
  WRITE( *, "(G0,G0)" ) "Overwritten Number of Particles: ", nParticles
END IF
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Total Number Density: ", TotalNumberDensity, "Å⁻³"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Total Molecular Volume: ", TotalParticleVolume, "Å³"

! Cube root check
pLoop: DO
  IF ( ConfigurationSelection(1) ) THEN
    CubicRootFloat = DBLE( nParticles ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CubicRootFloat - DBLE( DNINT( CubicRootFloat ) ) ) <= EPSILON( 1.D0 ) ) THEN
      EXIT pLoop
    ELSE
      WRITE( *, "(G0)" ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( InitialConfiguration )//" configuration."
      WRITE( *, "(G0)") "The total number of particles must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL Exit(  )
    END IF
  ELSE IF ( ConfigurationSelection(2) ) THEN
    CubicRootFloat = ( 0.5D0 * DBLE( nParticles ) ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CubicRootFloat - DBLE( DNINT( CubicRootFloat ) ) ) <= EPSILON( 1.D0 ) ) THEN
      EXIT pLoop
    ELSE
      WRITE( *, "(G0)" ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( InitialConfiguration )//" configuration."
      WRITE( *, "(G0)") "The total number of particles divided by 2 must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL Exit(  )
    END IF
  ELSE IF ( ConfigurationSelection(3) ) THEN
    CubicRootFloat = ( 0.25D0 * DBLE( nParticles ) ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CubicRootFloat - DBLE( DNINT( CubicRootFloat ) ) ) <= EPSILON( 1.D0 ) ) THEN
      EXIT pLoop
    ELSE
      WRITE( *, "(G0)" ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( InitialConfiguration )//" configuration."
      WRITE( *, "(G0)") "The total number of particles divided by 4 must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL Exit(  )
    END IF
  ELSE IF ( ConfigurationSelection(4) .OR. ConfigurationSelection(5) ) THEN
    EXIT pLoop
  END IF
END DO pLoop
WRITE( *, "(G0)" ) " "

! Absolute temperature
READ( 100, * ) Dummy, AbsoluteTemperature
! Condition
IF( AbsoluteTemperature <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The absolute temperature [", AbsoluteTemperature, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Reduced pressure
READ( 100, * ) Dummy, ReducedPressure
! Condition
IF( ReducedPressure <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The reduced pressure [", ReducedPressure, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"SYSTEM DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( EnsembleMC == "NVT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Packing Fraction: ", PackingFraction
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Volume of the Simulation Box: ", BoxVolume, "Å³"
  WRITE( *, "(G0)" ) " "
ELSE IF( EnsembleMC == "NPT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Initial Packing Fraction: ", PackingFraction
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Initial Volume of the Simulation Box: ", BoxVolume, "Å³"
  WRITE( *, "(G0)" ) " "
END IF
WRITE( *, "(G0,G0.5,G0)" ) "Absolute Temperature: ", AbsoluteTemperature, "K"
WRITE( *, "(G0)" ) " "
IF( EnsembleMC == "NPT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Reduced Pressure: ", ReducedPressure
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Real Pressure: ", (ReducedPressure * cBoltzmann * AbsoluteTemperature) / 1.D-30 / 1.D6, "MPa"
  WRITE( *, "(G0)" ) " "
END IF

! Initial configuration variables
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )

! Skip
READ( 100, * ) Dummy, Dummy
READ( 100, * ) Dummy, Dummy

! Unrotated reference axis (initial configuration)
READ( 100, * ) Dummy, BodyFixedAxisInquiry
CALL ToUpper( BodyFixedAxisInquiry, LEN_TRIM( BodyFixedAxisInquiry ), BodyFixedAxisInquiry )
! Transforms characters into logical variables
AxisSelection = .FALSE.
IF( BodyFixedAxisInquiry == "X" ) THEN
  AxisSelection(1) = .TRUE.
ELSE IF( BodyFixedAxisInquiry == "Y" ) THEN
  AxisSelection(2) = .TRUE.
ELSE IF( BodyFixedAxisInquiry == "Z" ) THEN
  AxisSelection(3) = .TRUE.
END IF
! Condition
IF( BodyFixedAxisInquiry /= "X" .AND. BodyFixedAxisInquiry /= "Y" .AND. BodyFixedAxisInquiry /= "Z" ) THEN
  WRITE( *, "(G0)" ) "The unrotated reference axis can only be X, Y, or Z. Exiting... "
  CALL Exit(  )
END IF

! Quaternion angle (initial configuration)
READ( 100, * ) Dummy, QuaternionAngle

! Random configuration
IF( ConfigurationSelection(4) ) THEN
  ! Initial packing fraction for the random configuration
  READ( 100, * ) Dummy, PackingFractionInitialConfiguration
  ! Condition 1
  IF( PackingFractionInitialConfiguration <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction of the random configuration [", PackingFractionInitialConfiguration, &
    &                          "] cannot be less than or equal to 0. Exiting... "
    CALL Exit(  )
  END IF
  ! Condition 2
  IF( PackingFractionInitialConfiguration > 1.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction of the random configuration [", PackingFractionInitialConfiguration, &
    &                          "] cannot be greater than 1. Exiting... "
    CALL Exit(  )
  END IF
  ! Initial pressure for the random configuration
  READ( 100, * ) Dummy, PressureRandomConfig
  ! Condition
  IF( PressureRandomConfig <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The reduced pressure of the random configuration [", PressureRandomConfig, &
    &                          "] cannot be less than or equal to 0. Exiting... "
    CALL Exit(  )
  END IF
ELSE IF( .NOT. ConfigurationSelection(4) ) THEN
  READ( 100, * ) Dummy, Dummy
  READ( 100, * ) Dummy, Dummy
END IF

! Preset initial configuration
READ( 100, * ) Dummy, PresetInitialConfiguration

CLOSE( 100 )

! Simulation variables
OPEN( Unit= 100, File= "ini_ratios.ini", Action= "READ" )

! Acceptance ratio (translation)
READ( 100, * ) Dummy, AcceptanceRatioTranslation
! Condition 1
IF( AcceptanceRatioTranslation <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The translational acceptance ratio [", AcceptanceRatioTranslation, &
  &                          "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( AcceptanceRatioTranslation > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The translational acceptance ratio [", AcceptanceRatioTranslation, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF

! Acceptance ratio (rotation)
READ( 100, * ) Dummy, AcceptanceRatioRotation
! Condition 1
IF( AcceptanceRatioRotation <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The rotational acceptance ratio [", AcceptanceRatioRotation, &
  &                          "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( AcceptanceRatioRotation > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The rotational acceptance ratio [", AcceptanceRatioRotation, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF

! Acceptance ratio (isotropic volume scaling)
READ( 100, * ) Dummy, AcceptanceRatioIsoVolumeChange
! Condition 1
IF( AcceptanceRatioIsoVolumeChange <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (isotropic) [", AcceptanceRatioIsoVolumeChange, &
  &                          "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( AcceptanceRatioIsoVolumeChange > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (isotropic) [", AcceptanceRatioIsoVolumeChange, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF

! Acceptance ratio (anisotropic volume scaling)
READ( 100, * ) Dummy, AcceptanceRatioAnisoVolumeChange
! Condition 1
IF( AcceptanceRatioAnisoVolumeChange <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (anisotropic) [", AcceptanceRatioAnisoVolumeChange, &
  &                          "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( AcceptanceRatioAnisoVolumeChange > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (anisotropic) [", AcceptanceRatioAnisoVolumeChange, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF

CLOSE( 100 )

! Simulation variables
OPEN( Unit= 100, File= "ini_probabilities.ini", Action= "READ" )

! Movement/Volume scaling probability
READ( 100, * ) Dummy, MovementProbability
! Condition 1
IF( MovementProbability < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement [", MovementProbability, "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( MovementProbability > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement [", MovementProbability, "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
VolumeChangeProbability = 1.D0 - MovementProbability
IF( DABS( VolumeChangeProbability ) - 0.D0 <= EPSILON( 1.D0 ) .AND. EnsembleMC == "NPT" ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The probability of volume scaling [", VolumeChangeProbability, "] cannot be 0 when ", &
  &                           "the NPT ensemble is selected. Exiting... "
  CALL Exit(  )
END IF

! Movement/Volume scaling probability (initial configuration)
READ( 100, * ) Dummy, MovementProbabilityRandomConfig
! Condition 1
IF( MovementProbabilityRandomConfig < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement of the random configuration [", MovementProbabilityRandomConfig, &
  &                          "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( MovementProbabilityRandomConfig > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement of the random configuration [", MovementProbabilityRandomConfig, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
VolumeChangeProbabilityRandomConfig = 1.D0 - MovementProbabilityRandomConfig

! Translational/Rotational movement probability
READ( 100, * ) Dummy, TranslationalProbability
! Condition 1
IF( TranslationalProbability < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements [", TranslationalProbability, &
  &                          "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( TranslationalProbability > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements [", TranslationalProbability, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
RotationalProbability = 1.D0 - TranslationalProbability

! Translational/Rotational movement probability (initial configuration)
READ( 100, * ) Dummy, TranslationalProbabilityRandomConfig
! Condition 1
IF( TranslationalProbabilityRandomConfig < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements of the random configuration [", &
  &                          TranslationalProbabilityRandomConfig, "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( TranslationalProbabilityRandomConfig > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements of the random configuration [", &
  &                          TranslationalProbabilityRandomConfig, "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
RotationalProbabilityRandomConfig = 1.D0 - TranslationalProbabilityRandomConfig

! Isotropic/Anisotropic volume scaling
READ( 100, * ) Dummy, IsoVolumetricProbability
! Condition 1
IF( IsoVolumetricProbability < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes [", IsoVolumetricProbability, &
  &                          "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( IsoVolumetricProbability > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes [", IsoVolumetricProbability, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
AnisoVolumetricProbability = 1.D0 - IsoVolumetricProbability

! Isotropic/Anisotropic volume scaling (initial configuration)
READ( 100, * ) Dummy, IsoVolumetricProbabilityRandomConfig
! Condition 1
IF( IsoVolumetricProbabilityRandomConfig < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes of the random configuration [", &
  &                          IsoVolumetricProbabilityRandomConfig, "] cannot be less than 0. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( IsoVolumetricProbabilityRandomConfig > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes of the random configuration [", &
  &                          IsoVolumetricProbabilityRandomConfig, "] cannot be greater than 1. Exiting... "
  CALL Exit(  )
END IF
AnisoVolumetricProbabilityRandomConfig = 1.D0 - IsoVolumetricProbabilityRandomConfig

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"INITIAL CONFIGURATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Fixed-Body Axis: ", BodyFixedAxisInquiry
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Quaternion Angle: ", QuaternionAngle, "°"
WRITE( *, "(G0)" ) " "
IF( ConfigurationSelection(4) ) THEN
  WRITE( *, "(G0,G0.5)" ) "Initial Packing Fraction (Random Configuration): ", PackingFractionInitialConfiguration
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5)" ) "Target Reduced Pressure (Random Configuration): ", PressureRandomConfig
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0,G0)" ) "Adjustment Frequency (Random Configuration): Every ", nAdjustmentRandomConfig, " Cycle(s)"
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5)" ) "Maximum Rotation (Random Configuration): ", MaxTranslationalDisplacementRandomConfig
  WRITE( *, "(G0,G0.5)" ) "Maximum Translation (Random Configuration): ", MaxAngularDisplacementRandomConfig
  WRITE( *, "(G0,G0.5)" ) "Maximum Isotropic Volume Scaling (Random Configuration): ", MaxIsoVolumetricDisplacementRandomConfig
  WRITE( *, "(G0,G0.5)" ) "Maximum Anisotropic Volume Scaling (Random Configuration): ", MaxAnisoVolumetricDisplacementRandomConfig
  WRITE( *, "(G0,G0.5)" ) "Minimum Volume Scaling (Random Configuration): ", MinVolumetricDisplacementRandomConfig
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Random Configuration): ", MovementProbabilityRandomConfig * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Volume Scaling Probability (Random Configuration): ", VolumeChangeProbabilityRandomConfig * 100.D0, &
  &                          "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Translation Probability (Random Configuration): ", TranslationalProbabilityRandomConfig * 100.D0, &
  &                          "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Rotation Probability (Random Configuration): ", RotationalProbabilityRandomConfig * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Isotropic Volume Scaling Probability (Random Configuration): ", &
  &                          IsoVolumetricProbabilityRandomConfig * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Anisotropic Volume Scaling Probability (Random Configuration): ", &
  &                          AnisoVolumetricProbabilityRandomConfig * 100.D0, "%"
  WRITE( *, "(G0)" ) " "
END IF
IF( PresetInitialConfiguration ) THEN
  WRITE( *, "(G0)" ) "Preset Initial Configuration: [YES]"
ELSE IF( .NOT. PresetInitialConfiguration ) THEN
  WRITE( *, "(G0)" ) "Preset Initial Configuration: [NO]"
END IF
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 18 )//"MONTE CARLO DETAILS"//REPEAT( " ", 18 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Number of Cycles: ", MaxSimulationCycles
WRITE( *, "(G0,G0)" ) "Number of Equilibration Cycles: ", nEquilibrationCycles
WRITE( *, "(G0,G0)" ) "Number of Production Cycles: ", MaxSimulationCycles - nEquilibrationCycles
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0,G0)" ) "Saving Frequency: Every ", nSavingFrequency, " Cycle(s)"
WRITE( *, "(G0,G0,G0)" ) "Adjustment Frequency: Every ", nAdjustmentFrequency, " Equilibration Cycle(s)"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Translation): ", AcceptanceRatioTranslation
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Rotation): ", AcceptanceRatioRotation
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Isotropic Volume Scaling): ", AcceptanceRatioIsoVolumeChange
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Anisotropic Volume Scaling): ", AcceptanceRatioAnisoVolumeChange
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Rotation): ", UserMaxTranslationalDisplacement
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Translation): ", UserMaxRotationalDisplacement
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Isotropic Volume Scaling): ", UserMaxIsoVolumetricDisplacement
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Anisotropic Volume Scaling): ", UserMaxAnisoVolumetricDisplacement
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Movement): ", MovementProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Volume Scaling): ", VolumeChangeProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Translation): ", TranslationalProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Rotation): ", RotationalProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Scaling Probability (Isotropic): ", IsoVolumetricProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Scaling Probability (Anisotropic): ", AnisoVolumetricProbability * 100.D0, "%"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Box Distortion: ", MaxBoxDistortion
WRITE( *, "(G0,G0.5)" ) "Maximum Box Length Distortion: ", BoxEdgeMaxRatio
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Box Angle Distortion: ", BoxVectorMaxAngle * 180.D0 / cPi, "°"
IF( LatticeReductionTypeLogical(1) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Algorithm: Gottwald"
ELSE IF( LatticeReductionTypeLogical(2) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Algorithm: Lenstra-Lenstra-Lovász"
END IF
WRITE( *, "(G0)" ) " "
IF( TrajectoryLogical ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [YES]"
ELSE IF( .NOT. TrajectoryLogical ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [NO]"
END IF
WRITE( *, "(G0)" ) " "
IF( BackupFileLogical ) THEN
  WRITE( *, "(G0)" ) "Backup File: [YES]"
ELSE IF( .NOT. BackupFileLogical ) THEN
  WRITE( *, "(G0)" ) "Backup File: [NO]"
END IF
IF( RestoreBackupFileLogical ) THEN
  WRITE( *, "(G0)" ) "Restore Backup File? [YES]"
ELSE IF( .NOT. RestoreBackupFileLogical ) THEN
  WRITE( *, "(G0)" ) "Restore Backup File? [NO]"
END IF
WRITE( *, "(G0)" ) " "
IF( FixedSeedLogical ) THEN
  WRITE( *, "(G0)" ) "Seed Type: [FIXED]"
ELSE IF( .NOT. FixedSeedLogical ) THEN
  WRITE( *, "(G0)" ) "Seed Type: [RANDOM]"
END IF
WRITE( *, "(G0)" ) " "
IF( CellListLogical ) THEN
  WRITE( *, "(G0)" ) "Linked Lists: [ENABLED]"
ELSE IF( .NOT. CellListLogical ) THEN
  WRITE( *, "(G0)" ) "Linked Lists: [DISABLED]"
END IF
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE GeneralVariables

! *********************************************************************************************** !
!                             Initialization of force field variables                             !
! *********************************************************************************************** !
SUBROUTINE ForceFieldVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: rRange ! Counter

! Potential variables
OPEN( Unit= 10, File= "ini_potential.ini", Action= "READ" )

! Square-well potential
IF( PotentialTypeLogical(2) ) THEN
  ! Number of attractive range points
  READ( 10, * ) Dummy, nRange
  ! Condition
  IF( nRange < 1 ) THEN
    WRITE( *, "(3G0)" ) "The number of potential ranges [", nRange, "] cannot be less than 1. Exiting... "
    CALL Exit(  )
  END IF
  ! Allocation
  ALLOCATE( PotentialRange(nRange) )
  ! Attractive range (λ)
  READ( 10, * ) Dummy, PotentialRange
  ! Condition
  DO rRange = 1, nRange
    IF( PotentialRange(rRange) <= 1.D0 ) THEN
      WRITE( *, "(3G0,G0.5,G0)" ) "The attractive range #", rRange, " [", PotentialRange(rRange), &
      &                           "] cannot be less than or equal to 1. Exiting... "
      CALL Exit(  )
    END IF
  END DO
  ! Reduced Temperature
  READ( 10, * ) Dummy, ReducedTemperature
  ! Condition
  IF( ReducedTemperature <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The reduced temperature [", ReducedTemperature, "] cannot be less than or equal to 0. Exiting... "
    CALL Exit(  )
  END IF
  ! Minimum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
  READ( 10, * ) Dummy, MinBlocks
  ! Maximum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
  READ( 10, * ) Dummy, MaxBlocks
  ! Condition 1
  IF( MinBlocks < 1 ) THEN
    WRITE( *, "(3G0)" ) "The minimum number of blocks [", MinBlocks, "] cannot be less than 1. Exiting... "
    CALL Exit(  )
  END IF
  ! Condition 2
  IF( MaxBlocks < 1 ) THEN
    WRITE( *, "(3G0)" ) "The maximum number of blocks [", MaxBlocks, "] cannot be less than 1. Exiting... "
    CALL Exit(  )
  END IF
  ! Condition 3
  IF( MinBlocks >= MaxBlocks ) THEN
    WRITE( *, "(6G0)" ) "The minimum number of blocks [", MinBlocks, "] cannot be greater than or equal to the maximum ", &
    &                   "number of blocks [", MaxBlocks, "]. Exiting... "
    CALL Exit(  )
  END IF
END IF

CLOSE( 10 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"POTENTIAL DETAILS"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( PotentialTypeLogical(1) ) THEN
  WRITE( *, "(G0)" ) "Potential Type: [HARDCORE]"
ELSE IF( PotentialTypeLogical(2) ) THEN
  WRITE( *, "(G0)" ) "Potential Type: [SQUAREWELL]"
END IF
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(2) ) THEN
  WRITE( *, "(G0,G0)" ) "Number of Attractive Range Points: ", nRange
  WRITE( *, "(G0)", Advance= "NO" ) "Attractive Range Values: ["
  DO rRange = 1, nRange
    WRITE( *, "(G0.5)", Advance= "NO" ) PotentialRange(rRange)
    IF( rRange /= nRange ) THEN
      WRITE( *, "(G0)", Advance= "NO" ) ", "
    ELSE
      WRITE( *, "(G0)", Advance= "YES" ) "]"
    END IF
  END DO
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5)" ) "Reduced Temperature: ", ReducedTemperature
  WRITE( *, "(G0)" ) " "
  IF( PotentialEnergyLogical ) THEN
    WRITE( *, "(G0)" ) "Potential Files: Production-Only"
  ELSE IF( .NOT. PotentialEnergyLogical ) THEN
    WRITE( *, "(G0)" ) "Potential Files: Equilibration and Production"
  END IF
  WRITE( *, "(G0)" ) " "
  IF( PerturbationCoefficientLogical .AND. EnsembleMC == "NVT" ) THEN
    WRITE( *, "(G0)" ) "Perturbation Coefficients: [YES]"
    WRITE( *, "(2G0)" ) "Minimum Number of Blocks (Block Averaging): ", MinBlocks
    WRITE( *, "(2G0)" ) "Maximum Number of Blocks (Block Averaging): ", MaxBlocks
  ELSE IF( .NOT. PerturbationCoefficientLogical .AND. EnsembleMC == "NVT" ) THEN
    WRITE( *, "(G0)" ) "Perturbation Coefficients: [NO]"
  END IF
  WRITE( *, "(G0)" ) " "
END IF

! Summary
WRITE( *, "(G0)" ) "RESUME? [Y/N]"
READ( *, * ) Dummy
CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
IF( Dummy /= "Y" ) THEN
  CALL Exit(  )
END IF
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE ForceFieldVariables

! *********************************************************************************************** !
!                           Writes backup variables in an external file                           !
! *********************************************************************************************** !
SUBROUTINE WriteBackupVariables( DateTime )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 20 ) :: DescriptorBackupString ! Descriptor for backup variable

! Backup file of common variables
OPEN( Unit= 105, File= "Backup/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_variables.backup" )
WRITE( 105, "(2G0)" ) "Number_of_Particles: ", nParticles
WRITE( 105, "(2G0)" ) "Number_of_Components: ", nComponents
WRITE( 105, "(2G0)" ) "Cell_List_Inquiry: ", CellListLogical
WRITE( 105, "(G0,5(G0,1X))" ) "Potential_Type_Logical: ", PotentialTypeLogical
IF( PotentialTypeLogical(2) ) THEN
  WRITE( 105, "(2G0)" ) "Number_of_Potential_Ranges: ", nRange
  WRITE( DescriptorBackupString, "(G0)" ) nRange
  DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
  WRITE( 105, DescriptorBackupString ) "Potential_Ranges: ", PotentialRange
  WRITE( 105, "(2G0)" ) "Reduced_Temperature: ", ReducedTemperature
  WRITE( 105, "(2G0)" ) "Min_Blocks: ", MinBlocks
  WRITE( 105, "(2G0)" ) "Max_Blocks: ", MaxBlocks
END IF
WRITE( 105, "(2G0)" ) "Geometry_Inquiry: ", "'"//TRIM( GeometryInquiry )//"'"
WRITE( 105, "(2G0)" ) "Molecular_Geometry: ", "'"//TRIM( MolecularGeometry )//"'"
WRITE( 105, "(2G0)" ) "Geometry_Acronym: ", "'"//TRIM( GeometryAcronym )//"'"
WRITE( 105, "(G0,3(G0,1X))" ) "Geometry_Type: ", GeometryType
WRITE( 105, "(2G0)" ) "Configuration_Inquiry: ", "'"//TRIM( ConfigurationInquiry )//"'"
WRITE( 105, "(2G0)" ) "Initial_Configuration: ", "'"//TRIM( InitialConfiguration )//"'"
WRITE( 105, "(G0,5(G0,1X))" ) "Configuration_Selection: ", ConfigurationSelection
WRITE( 105, "(2G0)" ) "Max_Simulation_Cycles: ", MaxSimulationCycles
WRITE( 105, "(2G0)" ) "Number_of_Equilibration_Cycles: ", nEquilibrationCycles
WRITE( 105, "(2G0)" ) "Saving_Frequency: ", nSavingFrequency
WRITE( 105, "(2G0)" ) "Adjustment_Frequency: ", nAdjustmentFrequency
WRITE( 105, "(2G0)" ) "Adjustment_Frequency_Random_Config: ", nAdjustmentRandomConfig
WRITE( 105, "(2G0)" ) "User_Max_Translational_Displacement: ", UserMaxTranslationalDisplacement
WRITE( 105, "(2G0)" ) "Max_Translational_Displacement_Random_Config: ", MaxTranslationalDisplacementRandomConfig
WRITE( 105, "(2G0)" ) "User_Max_Rotational_Displacement: ", UserMaxRotationalDisplacement
WRITE( 105, "(2G0)" ) "Max_Angular_Displacement_Random_Config: ", MaxAngularDisplacementRandomConfig
WRITE( 105, "(2G0)" ) "User_Max_Iso_Volumetric_Displacement: ", UserMaxIsoVolumetricDisplacement
WRITE( 105, "(2G0)" ) "User_Max_Aniso_Volumetric_Displacement: ", UserMaxAnisoVolumetricDisplacement
WRITE( 105, "(2G0)" ) "Max_Iso_Volumetric_Displacement_Random_Config: ", MaxIsoVolumetricDisplacementRandomConfig
WRITE( 105, "(2G0)" ) "Max_Aniso_Volumetric_Displacement_Random_Config: ", MaxAnisoVolumetricDisplacementRandomConfig
WRITE( 105, "(2G0)" ) "Min_Volumetric_Displacement_Random_Config: ", MinVolumetricDisplacementRandomConfig
WRITE( 105, "(2G0)" ) "Max_Box_Distortion: ", MaxBoxDistortion
WRITE( 105, "(2G0)" ) "Box_Edge_Max_Ratio: ", BoxEdgeMaxRatio
WRITE( 105, "(2G0)" ) "Box_Vector_Max_Angle: ", BoxVectorMaxAngle
WRITE( 105, "(2G0)" ) "Lattice_Reduction_Type: ", "'"//TRIM( LatticeReductionType )//"'"
WRITE( 105, "(G0,2(G0,1X))" ) "Lattice_Reduction_Type_Logical: ", LatticeReductionTypeLogical
WRITE( 105, "(2G0)" ) "Ensemble_MC: ", "'"//TRIM( EnsembleMC )//"'"
WRITE( 105, "(2G0)" ) "Trajectory_Inquiry: ", "'"//TRIM( TrajectoryInquiry )//"'"
WRITE( 105, "(2G0)" ) "Trajectory_Logical: ", TrajectoryLogical
WRITE( 105, "(2G0)" ) "Seed_Type_Inquiry: ", "'"//TRIM( SeedTypeInquiry )//"'"
WRITE( 105, "(2G0)" ) "Fixed_Seed_Logical: ", FixedSeedLogical
WRITE( 105, "(2G0)" ) "Backup_File_Inquiry: ", "'"//TRIM( BackupFileInquiry )//"'"
WRITE( 105, "(2G0)" ) "Backup_File_Logical: ", BackupFileLogical
WRITE( 105, "(2G0)" ) "Potential_Type: ", "'"//TRIM( PotentialType )//"'"
IF( PotentialTypeLogical(2) ) THEN
  WRITE( 105, "(2G0)" ) "Potential_Energy_Inquiry: ", "'"//TRIM( PotentialEnergyInquiry )//"'"
  WRITE( 105, "(2G0)" ) "Potential_Energy_Logical: ", PotentialEnergyLogical
  WRITE( 105, "(2G0)" ) "Perturbation_Coefficient_Inquiry: ", "'"//TRIM( PerturbationCoefficientInquiry )//"'"
  WRITE( 105, "(2G0)" ) "Perturbation_Coefficient_Logical: ", PerturbationCoefficientLogical
END IF
WRITE( DescriptorBackupString, "(G0)" ) nComponents
DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
WRITE( 105, DescriptorBackupString ) "Spherical_Component_Logical: ", SphericalComponentLogical
WRITE( 105, DescriptorBackupString ) "Diameter_of_Component: ", cDiameter
WRITE( 105, DescriptorBackupString ) "Length_of_Component: ", cLength
WRITE( 105, DescriptorBackupString ) "Molar_Fraction_of_Component: ", cMolarFraction
WRITE( 105, DescriptorBackupString ) "Aspect_Ratio_of_Component: ", cAspectRatio
WRITE( DescriptorBackupString, "(G0)" ) nComponents + 1
DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
WRITE( 105, DescriptorBackupString ) "Number_of_Particles_of_Component: ", cParticles
WRITE( DescriptorBackupString, "(G0)" ) nComponents
DescriptorBackupString = "(G0,"//TRIM( DescriptorBackupString )//"(G0,1X)"//")"
WRITE( 105, DescriptorBackupString ) "Molecular_Volume_of_Component: ", cMolecularVolume
WRITE( 105, "(2G0)" ) "Total_Particle_Volume: ", TotalParticleVolume
WRITE( 105, DescriptorBackupString ) "Number_Density_of_Component: ", cNumberDensity
WRITE( 105, "(2G0)" ) "Absolute_Temperature: ", AbsoluteTemperature
WRITE( 105, "(2G0)" ) "Reduced_Pressure: ", ReducedPressure
WRITE( 105, "(2G0)" ) "Body_Fixed_Axis_Inquiry: ", "'"//TRIM( BodyFixedAxisInquiry )//"'"
WRITE( 105, "(G0,3(G0,1X))" ) "Axis_Selection: ", AxisSelection
WRITE( 105, "(2G0)" ) "Quaternion_Angle: ", QuaternionAngle
WRITE( 105, "(2G0)" ) "Packing_Fraction_Initial_Configuration: ", PackingFractionInitialConfiguration
WRITE( 105, "(2G0)" ) "Pressure_Random_Config: ", PressureRandomConfig
WRITE( 105, "(2G0)" ) "Preset_Initial_Configuration: ", PresetInitialConfiguration
WRITE( 105, "(2G0)" ) "Acceptance_Ratio_Translation: ", AcceptanceRatioTranslation
WRITE( 105, "(2G0)" ) "Acceptance_Ratio_Rotation: ", AcceptanceRatioRotation
WRITE( 105, "(2G0)" ) "Acceptance_Ratio_Iso_Volume_Change: ", AcceptanceRatioIsoVolumeChange
WRITE( 105, "(2G0)" ) "Acceptance_Ratio_Aniso_Volume_Change: ", AcceptanceRatioAnisoVolumeChange
WRITE( 105, "(2G0)" ) "Movement_Probability: ", MovementProbability
WRITE( 105, "(2G0)" ) "Volume_Change_Probability: ", VolumeChangeProbability
WRITE( 105, "(2G0)" ) "Movement_Probability_Random_Config: ", MovementProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "Volume_Change_Probability_Random_Config: ", VolumeChangeProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "Translational_Probability: ", TranslationalProbability
WRITE( 105, "(2G0)" ) "Rotational_Probability: ", RotationalProbability
WRITE( 105, "(2G0)" ) "Translational_Probability_Random_Config: ", TranslationalProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "Rotational_Probability_Random_Config: ", RotationalProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "IsoVolumetric_Probability: ", IsoVolumetricProbability
WRITE( 105, "(2G0)" ) "Aniso_Volumetric_Probability: ", AnisoVolumetricProbability
WRITE( 105, "(2G0)" ) "Iso_Volumetric_Probability_Random_Config: ", IsoVolumetricProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "Aniso_Volumetric_Probability_Random_Config: ", AnisoVolumetricProbabilityRandomConfig
WRITE( 105, "(2G0)" ) "Initial_Seed: ", InitialSeed
WRITE( 105, "(G0,8(G0,1X))" ) "Date_Time: ", DateTime
WRITE( 105, "(2G0)" ) "Descriptor_Date: ", "'"//TRIM( DescriptorDate )//"'"
WRITE( 105, "(2G0)" ) "Descriptor_Hour: ", "'"//TRIM( DescriptorHour )//"'"
WRITE( 105, "(2G0)" ) "Descriptor_File_Thermo_Variable: ", "'"//TRIM( DescriptorFileThermoVariable )//"'"
WRITE( 105, "(2G0)" ) "Descriptor_File_Components: ", "'"//TRIM( DescriptorFileComponents )//"'"
WRITE( 105, "(2G0)" ) "Descriptor_File_Geometry: ", "'"//TRIM( DescriptorFileGeometry )//"'"
CLOSE( 105 )

RETURN

END SUBROUTINE WriteBackupVariables

! *********************************************************************************************** !
!                          Restore backup variables from a previous run                           !
! *********************************************************************************************** !
SUBROUTINE RestoreBackupVariables( DescriptorBackupFile, DateTime )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: cComponent, pParticle ! Counter

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 14 ) :: DescriptorBackupFile ! Descriptor for output file

! Backup file of common variables from previous simulations
OPEN( Unit= 105, File= "Backup/"//TRIM( DescriptorBackupFile )//"_variables.backup", Action= "READ" )
READ( 105, * ) Dummy, nParticles
READ( 105, * ) Dummy, nComponents
READ( 105, * ) Dummy, CellListLogical
READ( 105, * ) Dummy, PotentialTypeLogical
! Allocation
ALLOCATE( pQuaternion(0:3,nParticles), pQuaternionMC(0:3,nParticles) )
ALLOCATE( pPosition(3,nParticles), pPositionMC(3,nParticles) )
ALLOCATE( pOrientation(3,nParticles), pOrientationMC(3,nParticles) )
ALLOCATE( cDiameter(nComponents), cLength(nComponents), cAspectRatio(nComponents), cMolecularVolume(nComponents) )
ALLOCATE( pComponents(nParticles) )
ALLOCATE( cParticles(0:nComponents), cMolarFraction(nComponents), cNumberDensity(nComponents) )
ALLOCATE( cCircumscribingSphereDiameter(nComponents) )
ALLOCATE( cDiameterEquivalentSphere(nComponents) )
ALLOCATE( SphericalComponentInquiry(nComponents) )
ALLOCATE( SphericalComponentLogical(nComponents) )
IF( CellListLogical ) ALLOCATE( pCellList(nParticles), pCellIndex(3,nParticles) )
IF( CellListLogical .AND. PotentialTypeLogical(2) ) ALLOCATE( pCellListPotential(nParticles), pCellIndexPotential(3,nParticles) )
IF( PotentialTypeLogical(2) ) THEN
  READ( 105, * ) Dummy, nRange
  ! Allocation
  ALLOCATE( PotentialRange(nRange) )
  READ( 105, * ) Dummy, PotentialRange
  READ( 105, * ) Dummy, ReducedTemperature
  READ( 105, * ) Dummy, MinBlocks
  READ( 105, * ) Dummy, MaxBlocks
END IF
READ( 105, * ) Dummy, GeometryInquiry
READ( 105, * ) Dummy, MolecularGeometry
READ( 105, * ) Dummy, GeometryAcronym
READ( 105, * ) Dummy, GeometryType
READ( 105, * ) Dummy, ConfigurationInquiry
READ( 105, * ) Dummy, InitialConfiguration
READ( 105, * ) Dummy, ConfigurationSelection
READ( 105, * ) Dummy, MaxSimulationCycles
READ( 105, * ) Dummy, nEquilibrationCycles
READ( 105, * ) Dummy, nSavingFrequency
READ( 105, * ) Dummy, nAdjustmentFrequency
READ( 105, * ) Dummy, nAdjustmentRandomConfig
READ( 105, * ) Dummy, UserMaxTranslationalDisplacement
READ( 105, * ) Dummy, MaxTranslationalDisplacementRandomConfig
READ( 105, * ) Dummy, UserMaxRotationalDisplacement
READ( 105, * ) Dummy, MaxAngularDisplacementRandomConfig
READ( 105, * ) Dummy, UserMaxIsoVolumetricDisplacement
READ( 105, * ) Dummy, UserMaxAnisoVolumetricDisplacement
READ( 105, * ) Dummy, MaxIsoVolumetricDisplacementRandomConfig
READ( 105, * ) Dummy, MaxAnisoVolumetricDisplacementRandomConfig
READ( 105, * ) Dummy, MinVolumetricDisplacementRandomConfig
READ( 105, * ) Dummy, MaxBoxDistortion
READ( 105, * ) Dummy, BoxEdgeMaxRatio
READ( 105, * ) Dummy, BoxVectorMaxAngle
READ( 105, * ) Dummy, LatticeReductionType
READ( 105, * ) Dummy, LatticeReductionTypeLogical
READ( 105, * ) Dummy, EnsembleMC
READ( 105, * ) Dummy, TrajectoryInquiry
READ( 105, * ) Dummy, TrajectoryLogical
READ( 105, * ) Dummy, SeedTypeInquiry
READ( 105, * ) Dummy, FixedSeedLogical
READ( 105, * ) Dummy, BackupFileInquiry
READ( 105, * ) Dummy, BackupFileLogical
READ( 105, * ) Dummy, PotentialType
IF( PotentialTypeLogical(2) ) THEN
  READ( 105, * ) Dummy, PotentialEnergyInquiry
  READ( 105, * ) Dummy, PotentialEnergyLogical
  READ( 105, * ) Dummy, PerturbationCoefficientInquiry
  READ( 105, * ) Dummy, PerturbationCoefficientLogical
END IF
READ( 105, * ) Dummy, SphericalComponentLogical
READ( 105, * ) Dummy, cDiameter
READ( 105, * ) Dummy, cLength
READ( 105, * ) Dummy, cMolarFraction
READ( 105, * ) Dummy, cAspectRatio
READ( 105, * ) Dummy, cParticles
DO cComponent = 1, nComponents
  DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
    pComponents(pParticle) = cComponent
  END DO
END DO
READ( 105, * ) Dummy, cMolecularVolume
READ( 105, * ) Dummy, TotalParticleVolume
READ( 105, * ) Dummy, cNumberDensity
READ( 105, * ) Dummy, AbsoluteTemperature
READ( 105, * ) Dummy, ReducedPressure
READ( 105, * ) Dummy, BodyFixedAxisInquiry
READ( 105, * ) Dummy, AxisSelection
READ( 105, * ) Dummy, QuaternionAngle
READ( 105, * ) Dummy, PackingFractionInitialConfiguration
READ( 105, * ) Dummy, PressureRandomConfig
READ( 105, * ) Dummy, PresetInitialConfiguration
READ( 105, * ) Dummy, AcceptanceRatioTranslation
READ( 105, * ) Dummy, AcceptanceRatioRotation
READ( 105, * ) Dummy, AcceptanceRatioIsoVolumeChange
READ( 105, * ) Dummy, AcceptanceRatioAnisoVolumeChange
READ( 105, * ) Dummy, MovementProbability
READ( 105, * ) Dummy, VolumeChangeProbability
READ( 105, * ) Dummy, MovementProbabilityRandomConfig
READ( 105, * ) Dummy, VolumeChangeProbabilityRandomConfig
READ( 105, * ) Dummy, TranslationalProbability
READ( 105, * ) Dummy, RotationalProbability
READ( 105, * ) Dummy, TranslationalProbabilityRandomConfig
READ( 105, * ) Dummy, RotationalProbabilityRandomConfig
READ( 105, * ) Dummy, IsoVolumetricProbability
READ( 105, * ) Dummy, AnisoVolumetricProbability
READ( 105, * ) Dummy, IsoVolumetricProbabilityRandomConfig
READ( 105, * ) Dummy, AnisoVolumetricProbabilityRandomConfig
READ( 105, * ) Dummy, InitialSeed
READ( 105, * ) Dummy, DateTime
READ( 105, * ) Dummy, DescriptorDate
READ( 105, * ) Dummy, DescriptorHour
READ( 105, * ) Dummy, DescriptorFileThermoVariable
READ( 105, * ) Dummy, DescriptorFileComponents
READ( 105, * ) Dummy, DescriptorFileGeometry
CLOSE( 105 )

RETURN

END SUBROUTINE RestoreBackupVariables

END MODULE InitializeVariables
