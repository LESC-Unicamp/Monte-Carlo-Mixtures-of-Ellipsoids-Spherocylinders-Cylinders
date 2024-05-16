! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!           This module defines the variables used by the main program and most of the            !
!         subroutines and functions. A brief description is presented for each variable.          !
!                                                                                                 !
! Version number: 1.4.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                         May 15th, 2024                                          !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE GlobalVar

! Use kind Real64 and Int64
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64

! OpenMP API
#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER PARAMETER                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER :: HalfNeighboursControl = 1 ! Checks whether a cell and its 26 surrounding cells are searched or just its 13 neighbour cells

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: nParticles      ! Number of particles
INTEGER( Kind= Int64 )                 :: nComponents     ! Number of components
INTEGER( Kind= Int64 )                 :: nRange          ! Number of attractive range parameters
INTEGER( Kind= Int64 )                 :: nThreads        ! Number of threads (OpenMP)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: pCells          ! Number of cells
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: pCellsPotential ! Number of cells (for the perturbed potential)

! *********************************************************************************************** !
! INTEGER VARIABLES (BLOCK-AVERAGING PARAMETERS)                                                  !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: MinBlocks ! Minimum number of blocks
INTEGER( Kind= Int64 ) :: MaxBlocks ! Maximum number of blocks

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: MaxSimulationCycles             ! Total number of simulation cycles
INTEGER( Kind= Int64 ) :: nEquilibrationCycles            ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nSavingFrequency                ! Saving frequency
INTEGER( Kind= Int64 ) :: nAdjustmentMovementFrequency    ! Movement adjustment frequency
INTEGER( Kind= Int64 ) :: nAdjustmentVolumeFrequency      ! Volumetric adjustment frequency
INTEGER( Kind= Int64 ) :: nAdjustmentMovementRandomConfig ! Movement adjustment frequency (initial configuration)
INTEGER( Kind= Int64 ) :: nAdjustmentVolumeRandomConfig   ! Volumetric adjustment frequency (initial configuration)

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: SeedValue           ! Random number generator seed
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: InitialSeed         ! Initial seed
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: cParticles          ! Number of particles of a component
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: pComponents         ! Component index of a particle
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: cIndex              ! Component index
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: pCellList           ! Cell list
INTEGER( Kind= Int64 ), DIMENSION( : ),       ALLOCATABLE :: pCellListPotential  ! Cell list (for the perturbed potential)
INTEGER( Kind= Int64 ), DIMENSION( :, : ),    ALLOCATABLE :: pCellIndex          ! 3D-cell index of each particle
INTEGER( Kind= Int64 ), DIMENSION( :, : ),    ALLOCATABLE :: pCellIndexPotential ! 3D-cell index of each particle (for the perturbed potential)
INTEGER( Kind= Int64 ), DIMENSION( :, :, : ), ALLOCATABLE :: pCellHead           ! Cell head
INTEGER( Kind= Int64 ), DIMENSION( :, :, : ), ALLOCATABLE :: pCellHeadPotential  ! Cell head (for the perturbed potential)

! *********************************************************************************************** !
! REAL VARIABLES (GENERAL)                                                                        !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: QuaternionAngle                     ! Quaternion angle [real part, W] (for the initial configuration only)
REAL( Kind= Real64 )                 :: PressureRandomConfig                ! Reduced pressure (for the initial configuration only)
REAL( Kind= Real64 )                 :: TotalNumberDensity                  ! Total number density
REAL( Kind= Real64 )                 :: TotalParticleVolume                 ! Total particle volume
REAL( Kind= Real64 )                 :: AbsoluteTemperature                 ! Absolute temperature
REAL( Kind= Real64 )                 :: ReducedTemperature                  ! Reduced temperature
REAL( Kind= Real64 )                 :: ReducedPressure                     ! Reduced pressure
REAL( Kind= Real64 )                 :: cLargestSphereDiameter              ! Diameter of the largest circumscribing sphere
REAL( Kind= Real64 )                 :: cLargestSphericalWell               ! Diameter of the largest spherical, attractive well
REAL( Kind= Real64 )                 :: BoxVolume                           ! Volume of the simulation box
REAL( Kind= Real64 )                 :: RandomNumber                        ! Random number from a pseudorandom number generator subroutine
REAL( Kind= Real64 )                 :: PackingFraction                     ! Packing fraction
REAL( Kind= Real64 )                 :: BoxEdgeMaxRatio                     ! Maximum length ratio of simulation box during anisotropic volume changes
REAL( Kind= Real64 )                 :: BoxVectorMaxAngle                   ! Maximum angle between box vectors during anisotropic volume changes
REAL( Kind= Real64 )                 :: PackingFractionInitialConfiguration ! Packing fraction (for the initial configuration only)
REAL( Kind= Real64 ), DIMENSION( 9 ) :: BoxLength                           ! Length (x,y,z) of the simulation box
REAL( Kind= Real64 ), DIMENSION( 9 ) :: BoxLengthInverse                    ! Length (x,y,z) of simulation box (inverse)

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS)                                                         !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: UserMaxTranslationalDisplacement           ! User maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: UserMaxRotationalDisplacement              ! User maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: UserMaxIsoVolumetricDisplacement           ! User maximum displacement [+/-] (Isotropic Volume Scaling)
REAL( Kind= Real64 ) :: UserMaxAnisoVolumetricDisplacement         ! User maximum displacement [+/-] (Anisotropic Volume Scaling)
REAL( Kind= Real64 ) :: MovementProbability                        ! Movement probability
REAL( Kind= Real64 ) :: TranslationalProbability                   ! Movement probability (Translation)
REAL( Kind= Real64 ) :: RotationalProbability                      ! Movement probability (Rotation)
REAL( Kind= Real64 ) :: VolumeChangeProbability                    ! Volume scaling probability
REAL( Kind= Real64 ) :: IsoVolumetricProbability                   ! Volume scaling probability (Isotropic)
REAL( Kind= Real64 ) :: AnisoVolumetricProbability                 ! Volume scaling probability (Anisotropic)
REAL( Kind= Real64 ) :: MovementProbabilityRandomConfig            ! Movement probability for the initial configuration
REAL( Kind= Real64 ) :: TranslationalProbabilityRandomConfig       ! Movement probability for the initial configuration (Translation)
REAL( Kind= Real64 ) :: RotationalProbabilityRandomConfig          ! Movement probability for the initial configuration (Rotation)
REAL( Kind= Real64 ) :: VolumeChangeProbabilityRandomConfig        ! Volume scaling probability for the initial configuration
REAL( Kind= Real64 ) :: IsoVolumetricProbabilityRandomConfig       ! Volume scaling probability for the initial configuration (Isotropic)
REAL( Kind= Real64 ) :: AnisoVolumetricProbabilityRandomConfig     ! Volume scaling probability for the initial configuration (Anisotropic)
REAL( Kind= Real64 ) :: MaxTranslationalDisplacementRandomConfig   ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: MaxAngularDisplacementRandomConfig         ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: MaxIsoVolumetricDisplacementRandomConfig   ! Maximum displacement [+/-] (Isotropic volume scaling)
REAL( Kind= Real64 ) :: MaxAnisoVolumetricDisplacementRandomConfig ! Maximum displacement [+/-] (Anisotropic volume scaling)
REAL( Kind= Real64 ) :: MinVolumetricDisplacementRandomConfig      ! Minimum displacement [+]   (Volume)
REAL( Kind= Real64 ) :: AcceptanceRatioTranslation                 ! Acceptance ratio threshold (Translation)
REAL( Kind= Real64 ) :: AcceptanceRatioRotation                    ! Acceptance ratio threshold (Rotation)
REAL( Kind= Real64 ) :: AcceptanceRatioIsoVolumeChange             ! Acceptance ratio threshold (Isotropic volume scaling)
REAL( Kind= Real64 ) :: AcceptanceRatioAnisoVolumeChange           ! Acceptance ratio threshold (Anisotropic volume scaling)
REAL( Kind= Real64 ) :: MaxBoxDistortion                           ! Maximum box distortion before performing a lattice reduction

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cDiameter                     ! Diameter of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cLength                       ! Length of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cMolarFraction                ! Molar fraction of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cMolecularVolume              ! Particle volume of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cNumberDensity                ! Number density of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cAspectRatio                  ! Aspect ratio of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cCircumscribingSphereDiameter ! Cutoff diameter
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: cDiameterEquivalentSphere     ! Diameter of a sphere with same volume of the nonspherical geometry of component c
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: TotalPotentialEnergy          ! Total potential energy
REAL( Kind= Real64 ), DIMENSION( : ),    ALLOCATABLE :: PotentialRange                ! Attractive range parameter
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pQuaternion, pQuaternionMC    ! Quaternion array of particles
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pPosition, pPositionMC        ! Position array of particles
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pOrientation, pOrientationMC  ! Orientation array of particles
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: cPotentialRange               ! Effective range of attraction of component c

! *********************************************************************************************** !
! REAL PARAMETERS                                                                                 !
! *********************************************************************************************** !
REAL( Kind= Real64 ), PARAMETER                 :: cPi = 4.D0 * DATAN( 1.D0 ) ! Pi number
REAL( Kind= Real64 ), PARAMETER                 :: cBoltzmann = 1.380649D-23  ! Boltzmann constant
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: xAxis = [1.D0,0.D0,0.D0]   ! Body-fixed axis of rotation along x-direction
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: yAxis = [0.D0,1.D0,0.D0]   ! Body-fixed axis of rotation along y-direction
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER :: zAxis = [0.D0,0.D0,1.D0]   ! Body-fixed axis of rotation along z-direction

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 10 ) :: DescriptorFileThermoVariable   ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorFileComponents       ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorFileGeometry         ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorHour                 ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DescriptorDate                 ! Descriptor for output folder
CHARACTER( LEN= 10 ) :: DescriptorRange                ! Descriptor for output folder
CHARACTER( LEN= 32 ) :: FormatFileThermoVariable       ! String format for output file
CHARACTER( LEN= 32 ) :: FormatFileComponents           ! String format for output file
CHARACTER( LEN= 32 ) :: FormatFileGeometry             ! String format for output file
CHARACTER( LEN= 32 ) :: FormatRange                    ! String format for output folder
CHARACTER( LEN= 01 ) :: TrajectoryInquiry              ! Trajectory output inquiry
CHARACTER( LEN= 03 ) :: ConfigurationInquiry           ! Molecular structure inquiry
CHARACTER( LEN= 03 ) :: GeometryInquiry                ! Molecular geometry inquiry
CHARACTER( LEN= 01 ) :: PotentialEnergyInquiry         ! Potential output inquiry
CHARACTER( LEN= 01 ) :: PerturbationCoefficientInquiry ! Coefficient output inquiry
CHARACTER( LEN= 01 ) :: BackupFileInquiry              ! Backup file output inquiry
CHARACTER( LEN= 01 ) :: RestoreBackupFileInquiry       ! Backup file restoration inquiry
CHARACTER( LEN= 32 ) :: InitialConfiguration           ! Molecular structure
CHARACTER( LEN= 03 ) :: LatticeReductionType           ! Lattice reduction type
CHARACTER( LEN= 32 ) :: PerturbedPotentialType         ! Potential type (perturbed system)
CHARACTER( LEN= 32 ) :: FullPotentialType              ! Potential type (reference system)
CHARACTER( LEN= 32 ) :: MolecularGeometry              ! Molecular geometry (extended)
CHARACTER( LEN= 03 ) :: GeometryAcronym                ! Molecular geometry (acronym)
CHARACTER( LEN= 01 ) :: BodyFixedAxisInquiry           ! Body-fixed axis inquiry
CHARACTER( LEN= 03 ) :: EnsembleMC                     ! Monte Carlo ensemble
CHARACTER( LEN= 01 ) :: SeedTypeInquiry                ! Fixed/random seed inquiry
CHARACTER( LEN= 07 ) :: RNGeneratorInquiry             ! Random number generator inquiry
CHARACTER( LEN= 01 ) :: Dummy                          ! Dummy (character)

! *********************************************************************************************** !
! CHARACTER STRINGS (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
CHARACTER( LEN= 01 ), DIMENSION( : ), ALLOCATABLE :: SphericalComponentInquiry ! Spherical component inquiry

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
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL                 :: TrajectoryLogical              ! Checks whether the trajectory file will be written out
LOGICAL                 :: PotentialEnergyLogical         ! Checks whether the potential file will be written out for production-related cycles or for all cycles
LOGICAL                 :: PerturbationCoefficientLogical ! Checks whether the perturbation coefficients will be calculated
LOGICAL                 :: BackupFileLogical              ! Checks whether a backup file will be generated
LOGICAL                 :: RestoreBackupFileLogical       ! Checks whether a backup file will be used to restore previous simulations
LOGICAL                 :: PresetInitialConfiguration     ! Checks whether the initial configuration will be read from an external file
LOGICAL                 :: FixedSeedLogical               ! Checks whether the seed for the random number generator will be fixed or random
LOGICAL                 :: CellListLogical                ! Checks whether a cell list will be used to compute the potential (NVT only) or search for overlaps
LOGICAL                 :: CellListControl                ! Toggles control of cell list
LOGICAL                 :: CellListControlPotential       ! Toggles control of cell list (for the perturbed potential)
LOGICAL, DIMENSION( 5 ) :: ConfigurationSelection         ! Checks the selected molecular configuration
LOGICAL, DIMENSION( 3 ) :: GeometryType                   ! Checks the selected molecular geometry
LOGICAL, DIMENSION( 2 ) :: LatticeReductionTypeLogical    ! Lattice reduction selection
LOGICAL, DIMENSION( 3 ) :: AxisSelection                  ! Checks the selected unrotated reference axis
LOGICAL, DIMENSION( 2 ) :: RNGeneratorLogical             ! Checks which random number generator will be used
LOGICAL, DIMENSION( 5 ) :: PerturbedPotentialTypeLogical  ! Checks the selected potential type (perturbed system)
LOGICAL, DIMENSION( 5 ) :: FullPotentialTypeLogical       ! Checks the selected potential type (reference system)

! *********************************************************************************************** !
! LOGICAL VARIABLES (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
LOGICAL, DIMENSION( : ), ALLOCATABLE :: SphericalComponentLogical ! Checks whether the component is spherical or not

END MODULE GlobalVar
