! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!      This module initialize common variables (number of particles, reduced number density,      !
!       reduced temperature etc.), Monte Carlo parameters (total number of cycles, number of      !
! equilibration cycles etc.), and Block Averaging parameters (maximum/minimum number of blocks).  !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
!     which results are written out in external files and enable post-processing subroutines.     !
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
! Main Reference:                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the USE of this code.        !
! ############################################################################################### !
MODULE VariableInitialization

! Uses one modules: global variables
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
!                               Initialization of common variables                                !
! *********************************************************************************************** !
SUBROUTINE CommonVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: kLayer    ! Counter (layers)
INTEGER( Kind= Int64 ) :: rPressure ! Counter (pressure)
INTEGER( Kind= Int64 ) :: iLine     ! Counter (lines)
INTEGER( Kind= Int64 ) :: nLines    ! Number of lines (format)
INTEGER( Kind= Int64 ) :: Remainder ! Remainder of division

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: SurfaceArea ! Surface area of the considered molecular geometry

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( Len= 30 ) :: DescriptorFormat ! Descriptor of formats

! *********************************************************************************************** !
! System Properties                                                                               !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_system.ini", Action= "READ" )

! Number of particles
READ( 100, * ) Dummy, nParticles
! Condition 1
IF( nParticles < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of particles [", nParticles, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF
! Condition 2
IF( nParticles > 12 ) THEN
  WRITE( *, "(3G0)" ) "The number of particles [", nParticles, "] cannot be greater than 12. Exiting... "
  CALL Exit(  )
END IF

! Diameter (cylinders)
READ( 100, * ) Dummy, cDiameter ! Å
! Condition
IF( cDiameter <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The diameter of the cylinders [", cDiameter, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Length (cylinders)
READ( 100, * ) Dummy, cLength ! Å
! Condition
IF( cLength <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The length of the cylinders [", cLength, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Aspect ratio (cylinders)
cAspectRatio = cLength / cDiameter

! Number of layers
READ( 100, * ) Dummy, nLayers
! Condition
IF( nLayers < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of layers [", nLayers, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Allocation
ALLOCATE( nImagesLayer(nLayers) )

! Number of images per layer
DO kLayer = 1, nLayers
  nImagesLayer(kLayer) = (24 * kLayer * kLayer) + 2
END DO

! Total number of images
nImages = SUM( nImagesLayer )

! Absolute temperature
READ( 100, * ) Dummy, AbsoluteTemperature ! K
! Condition
IF( AbsoluteTemperature <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The absolute temperature [", AbsoluteTemperature, "] cannot be less than or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Cross-sectional area
CALL CrossSectionArea( SurfaceArea )

! Molecular volume
pMolecularVolume = SurfaceArea * cLength ! Å³

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"PARTICLE DETAILS"//REPEAT( " ", 20 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(3G0)" ) "Geometry Type: [", Arrangement, "]"
WRITE( *, "(G0)" ) " "
WRITE( *, "(2G0)" ) "Number of Cylinders: ", 4 * nParticles
WRITE( *, "(G0,G0.5,G0)" ) "Diameter of Cylinders: ", cDiameter, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Length of Cylinders: ", cLength, "Å"
WRITE( *, "(G0,G0.5)" ) "Aspect Ratio of Cylinders: ", cAspectRatio
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Surface Area: ", SurfaceArea, "Å²"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Particle Volume: ", pMolecularVolume, "Å³"
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"GLOBAL DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(2G0)" ) "Number of Particles: ", nParticles
WRITE( *, "(G0)" ) " "
WRITE( *, "(2G0)" ) "Layers of Periodic Images: ", nLayers
WRITE( *, "(2G0)" ) "Total Number of Periodic Images / Particle: ", nImages
WRITE( *, "(2G0)" ) "Total Number of Periodic Images of Cylinders / Particle: ", nImages * 4
WRITE( *, "(G0)" ) " "

! *********************************************************************************************** !
! System Properties                                                                               !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_pramp.ini", Action= "READ" )

! Pressure points
READ( 100, * ) Dummy, nPressure
! Condition
IF( nPressure < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of pressure points [", nPressure, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Allocation
ALLOCATE( ReducedPressure(nPressure) )

! Setpoint pressure (ramp)
READ( 100, * ) Dummy, ReducedPressure ! P* = Pσ₀³/(kT), σ₀ = 1Å
! Condition
DO rPressure = 1, nPressure
  IF( ReducedPressure(rPressure) <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The reduced pressure point [", ReducedPressure(rPressure), &
    &                          "] cannot be less than or equal to 0. Exiting... "
    CALL Exit(  )
  END IF
END DO
! Sort array in ascending order
CALL RealArrayAscendingOrder( ReducedPressure, nPressure )

! Pressure index after which anisotropic volume changes are allowed
READ( 100, * ) Dummy, AnisotropicVolumeChangeIndex
! Condition
IF( AnisotropicVolumeChangeIndex < 1 .OR. AnisotropicVolumeChangeIndex > nPressure ) THEN
  WRITE( *, "(6G0)" ) "The index of the pressure point [", AnisotropicVolumeChangeIndex, "] cannot be less than 1 or greater ", &
  &                   "than the total number of pressure points [", nPressure, "]. Exiting... "
  CALL Exit(  )
END IF

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"SYSTEM DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(2G0)" ) "Number of Reduced Pressure Points: ", nPressure
Remainder = MOD( nPressure, 5_Int64 )
IF( Remainder == 0 ) THEN
  nLines = nPressure / 5
ELSE
  nLines = nPressure / 5 + 1
END IF
IF( Remainder == 0 .AND. nLines == 1 ) THEN
  WRITE( DescriptorFormat, "(G0)" ) nPressure
  WRITE( *, "(G0,"//TRIM( DescriptorFormat )//"(G0.5,1X),G0)" ) "Reduced Pressure Points: [ ", ReducedPressure, "]"
ELSE IF( Remainder == 0 .AND. nLines > 1 ) THEN
  WRITE( *, "(G0,5(G0.5,1X),G0)" ) "Reduced Pressure Points: [ ", ReducedPressure(1:5)
  DO iLine = 1, nLines - 1
    IF( iLine == nLines - 1 ) THEN
      WRITE( *, "(G0,5(G0.5,1X),G0)" ) REPEAT( " ", 27 ), ReducedPressure(5*iLine+1:5*(iLine+1)), "]"
    ELSE
      WRITE( *, "(G0,5(G0.5,1X))" ) REPEAT( " ", 27 ), ReducedPressure(5*iLine+1:5*(iLine+1))
    END IF
  END DO
ELSE IF( Remainder /= 0 .AND. nLines == 1 ) THEN
  WRITE( DescriptorFormat, "(G0)" ) Remainder
  WRITE( *, "(G0,"//TRIM( DescriptorFormat )//"(G0.5,1X),G0)" ) "Reduced Pressure Points: [ ", ReducedPressure, "]"
ELSE IF( Remainder /= 0 .AND. nLines > 1 ) THEN
  WRITE( *, "(G0,5(G0.5,1X))" ) "Reduced Pressure Points: [ ", ReducedPressure(1:5)
  DO iLine = 1, nLines - 2
    WRITE( *, "(G0,5(G0.5,1X))" ) REPEAT( " ", 27 ), ReducedPressure(5*iLine+1:5*(iLine+1))
  END DO
  WRITE( DescriptorFormat, "(G0)" ) Remainder
  WRITE( *, "(G0,"//TRIM( DescriptorFormat )//"(G0.5,1X),G0)" ) REPEAT( " ", 27 ), &
  &         ReducedPressure(5*(nLines-1)+1:5*(nLines-1)+Remainder), "]"
END IF
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Setpoint pressure to allow anisotropic volume changes: ", ReducedPressure(AnisotropicVolumeChangeIndex)
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Absolute Temperature: ", AbsoluteTemperature, "K"
WRITE( *, "(G0)" ) " "

! *********************************************************************************************** !
! Initial Configuration Properties                                                                !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )

! Skip
READ( 100 , * ) Dummy, Dummy
READ( 100 , * ) Dummy, Dummy

! Orientational angle (initial configuration only)
READ( 100 , * ) Dummy, QuaternionAngle ! Degrees

! Unrotated reference axis (initial configuration only)
READ( 100 , * ) Dummy, BodyFixedAxisInquiry
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

! Total number of cycles (NVT-simulation)
READ( 100 , * ) Dummy, InitialConfigurationMaxCycles
! Condition
IF( InitialConfigurationMaxCycles < 1 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of cycles for the initial configuration [", InitialConfigurationMaxCycles, &
  &                   "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"INITIAL CONFIGURATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(3G0)" ) "Initial Configuration: [", TRIM( ConfigurationInquiry ), "]"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0)" ) "Fixed-Body Axis: ", BodyFixedAxisInquiry
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Quaternion Angle: ", QuaternionAngle, "°"
WRITE( *, "(G0)" ) " "
WRITE( *, "(2G0)" ) "Total number of cycles (NVT-simulation): ", InitialConfigurationMaxCycles
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE CommonVariables

! *********************************************************************************************** !
!                             Initialization of Monte Carlo variables                             !
! *********************************************************************************************** !
SUBROUTINE MonteCarloVariables(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Simulation Properties                                                                           !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_montecarlo.ini", Action= "READ" )

! Total number of cycles
READ( 100, * ) Dummy, MaxSimulationCycles
! Condition
IF( MaxSimulationCycles < 1 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of cycles [", MaxSimulationCycles, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Number of equilibration cycles
READ( 100, * ) Dummy, nEquilibrationCycles
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
READ( 100 , * ) Dummy, nSavingFrequency
! Condition
IF( nSavingFrequency < 1 ) THEN
  WRITE( *, "(3G0)" ) "The saving frequency [", nSavingFrequency, "] cannot be a negative integer nor zero. Exiting... "
  CALL Exit(  )
END IF

! Adjustment frequency (movement)
READ( 100 , * ) Dummy, nAdjustmentMovementFrequency
! Condition
IF( nAdjustmentMovementFrequency < 1 ) THEN
  WRITE( *, "(4G0)" ) "The movement adjustment frequency of the simulation [", nAdjustmentMovementFrequency, "] cannot be a ", &
  &                   "negative integer nor zero. Exiting... "
  CALL Exit(  )
END IF

! Adjustment frequency (volume change)
READ( 100 , * ) Dummy, nAdjustmentVolumeFrequency
! Condition
IF( nAdjustmentVolumeFrequency < 1 ) THEN
  WRITE( *, "(4G0)" ) "The volumetric adjustment frequency of the simulation [", nAdjustmentVolumeFrequency, "] cannot be a ", &
  &                   "negative integer nor zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum translational displacement
READ( 100 , * ) Dummy, UserMaxTranslation
! Condition
IF( DABS( UserMaxTranslation - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum translational displacement of the simulation [", UserMaxTranslation, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum rotational displacement
READ( 100 , * ) Dummy, UserMaxRotation
! Condition
IF( DABS( UserMaxRotation - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum rotational displacement of the simulation [", UserMaxRotation, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum volumetric change
READ( 100 , * ) Dummy, UserMaxIsotropicVol
! Condition
IF( DABS( UserMaxIsotropicVol - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum isotropic volume scaling of the simulation [", UserMaxIsotropicVol, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum volumetric change
READ( 100 , * ) Dummy, UserMaxAnisotropicVol
! Condition
IF( DABS( UserMaxAnisotropicVol - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum anisotropic volume scaling of the simulation [", UserMaxAnisotropicVol, &
  &                          "] cannot be zero. Exiting... "
  CALL Exit(  )
END IF

! Maximum box distortion
READ( 100 , * ) Dummy, MaxBoxDistortion
! Condition
IF( MaxBoxDistortion < 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum box distortion [", MaxBoxDistortion, "] cannot be less than 1. Exiting... "
  CALL Exit(  )
END IF

! Lattice reduction method
READ( 100 , * ) Dummy, LatticeReductionType
CALL ToUpper( LatticeReductionType, LEN_TRIM( LatticeReductionType ), LatticeReductionType )
LatticeReductionTypeLogical = .FALSE.
! Lattice reduction: Gottwald method
IF( LatticeReductionType == "FBM" ) THEN
  LatticeReductionTypeLogical(1) = .TRUE.
! lattice reduction: Lenstra-Lenstra-Lovász method
ELSE IF( LatticeReductionType == "LLL" ) THEN
  LatticeReductionTypeLogical(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(4G0)" ) "The user-defined [", TRIM( LatticeReductionType ), "] is not an available lattice reduction method. ", &
  &                   "Exiting... "
  CALL Exit(  )
END IF

CLOSE( 100 )

! *********************************************************************************************** !
! Simulation Properties                                                                           !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_probabilities.ini", Action= "READ" )

! Displacement probability
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
! Volume change probability
VolumeChangeProbability = 1.D0 - MovementProbability

! Translational probability
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
! Rotational probability
RotationalProbability = 1.D0 - TranslationalProbability

! Isotropic probability
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
! Anisotropic probability
AnisoVolumetricProbability = 1.D0 - IsoVolumetricProbability

CLOSE( 100 )

! *********************************************************************************************** !
! Simulation Properties                                                                           !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_ratios.ini", Action= "READ" )

! Translational acceptance ratio threshold
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

! Rotational acceptance ratio threshold
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

! Volumetric acceptance ratio threshold (isotropic)
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

! Volumetric acceptance ratio threshold (anisotropic)
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

RETURN

END SUBROUTINE MonteCarloVariables

! *********************************************************************************************** !
!                               Initialization of Control Variables                               !
! *********************************************************************************************** !
SUBROUTINE VariablesInquest(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Control Variables                                                                               !
! *********************************************************************************************** !
OPEN( Unit= 100, File= "ini_control.ini", Action= "READ" )

! Trajectory inquiry
READ( 100 , * ) Dummy, TrajectoryInquiry
CALL ToUpper( TrajectoryInquiry, LEN_TRIM( TrajectoryInquiry ), TrajectoryInquiry )
! Transforms characters into logical variables
IF( TrajectoryInquiry == "Y" ) THEN
  TrajectoryLogical = .TRUE.
ELSE
  TrajectoryLogical = .FALSE.
END IF

! Maximum length ratio during anisotropic volume changes
READ( 100, * ) Dummy, BoxEdgeMaxRatio
! Condition
IF( BoxEdgeMaxRatio <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum linear distortion of the box [", BoxEdgeMaxRatio, "] cannot be less than ", &
  &                           "or equal to 0. Exiting... "
  CALL Exit(  )
END IF

! Maximum angle between box vectors during anisotropic volume changes
READ( 100, * ) Dummy, BoxVectorMaxAngle
! Condition
IF( BoxVectorMaxAngle <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum angular distortion of the box [", BoxVectorMaxAngle, "] cannot be less than or ", &
  &                           "equal to 0°. Exiting... "
  CALL Exit(  )
END IF
! Convert to radians
BoxVectorMaxAngle = BoxVectorMaxAngle * cPi / 180.D0

! Random number generator inquiry
RNGeneratorLogical = .FALSE.
READ( 100, * ) Dummy, RNGeneratorInquiry
CALL ToUpper( RNGeneratorInquiry, LEN_TRIM( RNGeneratorInquiry ), RNGeneratorInquiry )
! Transforms characters into logical variables
IF( RNGeneratorInquiry == "FORTRAN" ) THEN
  RNGeneratorLogical(1) = .TRUE.
ELSE IF( RNGeneratorInquiry == "BITWISE" ) THEN
  RNGeneratorLogical(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( RNGeneratorInquiry ), "] is not an available random number generator. Exiting..."
  CALL Exit(  )
END IF

CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 18 )//"MONTE CARLO DETAILS"//REPEAT( " ", 18 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Number of Cycles: ", MaxSimulationCycles
WRITE( *, "(G0,G0)" ) "Number of Equilibration Cycles: ", nEquilibrationCycles
WRITE( *, "(G0,G0)" ) "Number of Production Cycles: ", MaxSimulationCycles - nEquilibrationCycles
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0,G0)" ) "Saving Frequency: Every ", nSavingFrequency, " Cycle(s)"
WRITE( *, "(G0,G0,G0)" ) "Movement Adjustment Frequency: Every ", nAdjustmentMovementFrequency, " Equilibration Cycle(s)"
WRITE( *, "(G0,G0,G0)" ) "Volumetric Adjustment Frequency: Every ", nAdjustmentVolumeFrequency, " Equilibration Cycle(s)"
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Translation): ", AcceptanceRatioTranslation
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Rotation): ", AcceptanceRatioRotation
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Isotropic Volume Scaling): ", AcceptanceRatioIsoVolumeChange
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Anisotropic Volume Scaling): ", AcceptanceRatioAnisoVolumeChange
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Rotation): ", UserMaxTranslation
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Translation): ", UserMaxRotation
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Isotropic Volume Scaling): ", UserMaxIsotropicVol
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Anisotropic Volume Scaling): ", UserMaxAnisotropicVol
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
IF( RNGeneratorLogical(1) ) THEN
  WRITE( *, "(G0)" ) "Random Number Generator: [FORTRAN]"
ELSE IF( RNGeneratorLogical(2) ) THEN
  WRITE( *, "(G0)" ) "Random Number Generator: [BITWISE]"
END IF
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Resume? [Y/N]"
READ( *, * ) Dummy
CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
IF( Dummy /= "Y" ) CALL Exit(  )

RETURN

END SUBROUTINE VariablesInquest

END MODULE VariableInitialization
