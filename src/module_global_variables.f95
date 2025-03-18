! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!           This module defines the variables used by the main program and most of the            !
!         subroutines and functions. A brief description is presented for each variable.          !
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
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE GlobalVar

! Use kind Real64
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nParticles                   ! Number of particles
INTEGER( Kind= Int64 ) :: Arrangement                  ! Cross-sectional geometry
INTEGER( Kind= Int64 ) :: nLayers                      ! Number of layers of periodic images
INTEGER( Kind= Int64 ) :: nImages                      ! Total number of periodic images
INTEGER( Kind= Int64 ) :: nPressure                    ! Total number of pressure points
INTEGER( Kind= Int64 ) :: AnisotropicVolumeChangeIndex ! Index in the setpoint pressure array after which anisotropic volume changes are allowed

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL, ALLOCATABLE)                                                        !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE :: SeedValue    ! Random number generator seed
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE :: nImagesLayer ! Number of images per layer

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: MaxSimulationCycles           ! Total number of cycles
INTEGER( Kind= Int64 ) :: InitialConfigurationMaxCycles ! Total number of cycles (initial configuration only)
INTEGER( Kind= Int64 ) :: nEquilibrationCycles          ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nSavingFrequency              ! Saving frequency
INTEGER( Kind= Int64 ) :: nAdjustmentMovementFrequency  ! Adjustment frequency
INTEGER( Kind= Int64 ) :: nAdjustmentVolumeFrequency    ! Adjustment frequency

! *********************************************************************************************** !
! REAL VARIABLES (GENERAL)                                                                        !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: QuaternionAngle                    ! Quaternion angle [real part, W] (initial configuration only)
REAL( Kind= Real64 )                    :: BoxVolume                          ! Volume of simulation box
REAL( Kind= Real64 )                    :: iBoxVolume, fBoxVolume             ! Volume of simulation box (initial and final)
REAL( Kind= Real64 )                    :: TotalNumberDensity                 ! Number density
REAL( Kind= Real64 )                    :: iNumberDensity, fNumberDensity     ! Number density (initial and final)
REAL( Kind= Real64 )                    :: PackingFraction                    ! Packing fraction
REAL( Kind= Real64 )                    :: iPackingFraction, fPackingFraction ! Packing fraction (initial and final)
REAL( Kind= Real64 )                    :: AbsoluteTemperature                ! Temperature
REAL( Kind= Real64 )                    :: pMolecularVolume                   ! Molecular volume
REAL( Kind= Real64 )                    :: cAspectRatio                       ! Length-to-diameter (L/D) aspect ratio
REAL( Kind= Real64 )                    :: cDiameter                          ! Diameter
REAL( Kind= Real64 )                    :: cLength                            ! Length
REAL( Kind= Real64 )                    :: RandomNumber                       ! Pseudorandom number
REAL( Kind= Real64 )                    :: StartTimer                         ! Start timestamp
REAL( Kind= Real64 )                    :: StopTimer                          ! End timestamp
REAL( Kind= Real64 )                    :: MovementProbability                ! Probability of molecular displacement
REAL( Kind= Real64 )                    :: TranslationalProbability           ! Probability of movement (translation)
REAL( Kind= Real64 )                    :: RotationalProbability              ! Probability of movement (rotation)
REAL( Kind= Real64 )                    :: VolumeChangeProbability            ! Probability of volumetric change
REAL( Kind= Real64 )                    :: IsoVolumetricProbability           ! Probability of volume change (isotropic)
REAL( Kind= Real64 )                    :: AnisoVolumetricProbability         ! Probability of volume change (anisotropic)
REAL( Kind= Real64 )                    :: BoxEdgeMaxRatio                    ! Maximum length ratio of simulation box during anisotropic volume changes
REAL( Kind= Real64 )                    :: BoxVectorMaxAngle                  ! Maximum angle between box vectors during anisotropic volume changes
REAL( Kind= Real64 )                    :: pSquaredSphereDistance             ! Squared cutoff diameter (larger sphere)
REAL( Kind= Real64 )                    :: pSquaredSpheroCylinderDistance     ! Squared diameter (spherocylinder)
REAL( Kind= Real64 )                    :: cSquaredSphereDistance             ! Squared cutoff diameter (smaller sphere)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLength                          ! Length (x,y,z) of simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverse                   ! Length (x,y,z) of simulation box (inverse)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: iBoxLength, fBoxLength             ! Length (x,y,z) of simulation box (initial and final)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: xAxis                              ! Body-fixed axis of rotation (x) (see Allen and Tildesley, page 111)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: yAxis                              ! Body-fixed axis of rotation (y) (see Allen and Tildesley, page 111)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: zAxis                              ! Body-fixed axis of rotation (z) (see Allen and Tildesley, page 111)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cReferencePosition                 ! Body-fixed position

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS)                                                         !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: UserMaxTranslation               ! User maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: UserMaxRotation                  ! User maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: UserMaxIsotropicVol              ! User maximum displacement [+/-] (Isotropic)
REAL( Kind= Real64 ) :: UserMaxAnisotropicVol            ! User maximum displacement [+/-] (Anistropic)
REAL( Kind= Real64 ) :: AcceptanceRatioTranslation       ! Acceptance ratio threshold (Translation)
REAL( Kind= Real64 ) :: AcceptanceRatioRotation          ! Acceptance ratio threshold (Rotation)
REAL( Kind= Real64 ) :: AcceptanceRatioIsoVolumeChange   ! Acceptance ratio threshold (Isotropic volume change)
REAL( Kind= Real64 ) :: AcceptanceRatioAnisoVolumeChange ! Acceptance ratio threshold (Anisotropic volume change)
REAL( Kind= Real64 ) :: MaxBoxDistortion                 ! Maximum box distortion before lattice reduction

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS, ALLOCATABLE)                                            !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE          :: ReducedPressure                ! Reduced pressure (setpoint)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: pQuaternion, pQuaternionMC     ! Quaternion array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: pPosition, pPositionMC         ! Position array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: pOrientation, pOrientationMC   ! Orientation array
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: cPosition, cPositionMC         ! Position array (cylinders)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: imQuaternion, imQuaternionMC   ! Quaternion array (images)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: imPosition, imPositionMC       ! Position array (images)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: imOrientation, imOrientationMC ! Orientation array (images)
REAL( Kind= Real64 ), DIMENSION( :, :, :, : ), ALLOCATABLE :: imcPosition, imcPositionMC     ! Position array (cylinder images)

! *********************************************************************************************** !
! REAL PARAMETERS (CONSTANTS)                                                                     !
! *********************************************************************************************** !
REAL( Kind= Real64 ), PARAMETER :: cPi = 4.D0 * DATAN( 1.D0 ) ! π

! *********************************************************************************************** !
! CHARACTER STRINGS (GENERAL)                                                                     !
! *********************************************************************************************** !
CHARACTER( LEN= 01 )                 :: TrajectoryInquiry    ! Trajectory output inquiry
CHARACTER( LEN= 03 )                 :: ConfigurationInquiry ! Molecular configuration inquiry
CHARACTER( LEN= 01 ), DIMENSION( 2 ) :: Atom                 ! Atom ID

! *********************************************************************************************** !
! CHARACTER PARAMETERS                                                                            !
! *********************************************************************************************** !
CHARACTER( LEN= 3 ), PARAMETER :: CH_HS = "═" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VS = "║" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_UL = "╔" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_BL = "╚" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_UR = "╗" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_BR = "╝" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VL = "╠" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: CH_VR = "╣" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SH_VL = "╟" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SH_VR = "╢" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_HS = "─" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VS = "│" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VR = "├" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_VL = "┤" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_UL = "┌" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_BL = "└" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_UR = "┐" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: SS_BR = "┘" ! Box drawing symbol
CHARACTER( LEN= 3 ), PARAMETER :: C_FUL = "█" ! Box drawing symbol

! *********************************************************************************************** !
! CHARACTER STRINGS (FILE/FOLDER ORGANIZER)                                                       !
! *********************************************************************************************** !
CHARACTER( LEN= 10 ) :: DescriptorFileParticles ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorFileGeometry  ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorHour          ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorDate          ! Descriptor for output folder
CHARACTER( LEN= 32 ) :: FormatFileParticles     ! String format for output file
CHARACTER( LEN= 32 ) :: FormatFileGeometry      ! String format for output file
CHARACTER( LEN= 32 ) :: FormatHour              ! String format for output file
CHARACTER( LEN= 32 ) :: FormatDate              ! String format for output folder
CHARACTER( LEN= 01 ) :: BodyFixedAxisInquiry    ! Unrotated axis inquiry
CHARACTER( LEN= 03 ) :: LatticeReductionType    ! Lattice reduction type
CHARACTER( LEN= 32 ) :: InitialConfiguration    ! Molecular structure
CHARACTER( LEN= 07 ) :: RNGeneratorInquiry      ! Random number generator inquiry
CHARACTER( LEN= 01 ) :: Dummy                   ! Dummy (Character)

! *********************************************************************************************** !
! LOGICAL VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
LOGICAL                 :: TrajectoryLogical           ! Trajectory output selection
LOGICAL, DIMENSION( 1 ) :: ConfigurationSelection      ! Molecular configuration selection
LOGICAL, DIMENSION( 2 ) :: PetalArrangementLogical     ! Cross-sectional geometry selection
LOGICAL, DIMENSION( 2 ) :: RNGeneratorLogical          ! Checks which random number generator will be used
LOGICAL, DIMENSION( 2 ) :: LatticeReductionTypeLogical ! Lattice reduction selection
LOGICAL, DIMENSION( 3 ) :: AxisSelection               ! Unrotated reference axis selection

END MODULE GlobalVar