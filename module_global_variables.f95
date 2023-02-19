! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
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
!                                       February 15th, 2023                                       !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                 C. R. A. Abreu, F. A. Escobedo                                  !
!                                J. Chem. Phys. 124, 054116 (2006)                                !
!                                     DOI: 10.1063/1.2165188                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

MODULE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER, DIMENSION( 8 ) :: DATE_TIME ! Computer clock (date and time)

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER*8 :: SEED        ! Random number generator seed
INTEGER*8 :: N_PARTICLES ! Number of particles
INTEGER*8 :: COMPONENTS  ! Number of components

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER*8 :: MAX_CYCLES    ! Total number of cycles
INTEGER*8 :: N_EQUIL       ! Number of equilibration cycles
INTEGER*8 :: N_SAVE        ! Saving frequency
INTEGER*8 :: N_ADJUST      ! Adjustment frequency
INTEGER*8 :: N_ADJUST_INIT ! Adjustment frequency (initial configuration)
INTEGER*8 :: MAX_ATTEMPTS  ! Maximum number of attempts for the hit-and-miss algorithm

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
INTEGER*8, DIMENSION( : ), ALLOCATABLE :: N_COMPONENT ! Number of particles of a component
INTEGER*8, DIMENSION( : ), ALLOCATABLE :: INDEX_P     ! Molecule index

! *********************************************************************************************** !
! REAL VARIABLES (GENERAL)                                                                        !
! *********************************************************************************************** !
REAL*8                 :: QUATERNION_ANGLE   ! Quaternion angle [real part, W] (for the initial configuration only)
REAL*8                 :: ETA_INI            ! Packing fraction (for the initial configuration only)
REAL*8                 :: PRESS_RND          ! Reduced pressure (for the initial configuration only)
REAL*8                 :: TOTAL_RHO          ! Total number density (reduced)
REAL*8                 :: TOTAL_VP           ! Total particle volume (reduced)
REAL*8                 :: TEMP               ! Absolute temperature
REAL*8                 :: PRESS              ! Reduced pressure
REAL*8                 :: BOX_VOLUME, BOXVMC ! Reduced volume of the simulation box
REAL*8                 :: CHECK_ROOT         ! Checks the cube root for the number of particles, N, in each crystalline structure
REAL*8                 :: RANDOM_N           ! Random number from a pseudorandom number generator subroutine
REAL*8                 :: RANDOM_ANGLE       ! Random angle
REAL*8                 :: TOTAL_MOLAR_F      ! Total molar fraction
REAL*8                 :: START_TIMER        ! Start timer of Monte Carlo simulation
REAL*8                 :: STOP_TIMER         ! Stop timer of Monte Carlo simulation
REAL*8                 :: PACKING_F          ! Packing fraction
REAL*8                 :: ANYNUMBER          ! Dummy (number)
REAL*8, DIMENSION( 9 ) :: BOX_LENGTH         ! Length (x,y,z) of the simulation box
REAL*8, DIMENSION( 9 ) :: BOX_LENGTH_I       ! Length (x,y,z) of simulation box (inverse)
REAL*8, DIMENSION( 3 ) :: AXISX              ! Body-fixed axis of rotation along x-direction
REAL*8, DIMENSION( 3 ) :: AXISY              ! Body-fixed axis of rotation along y-direction
REAL*8, DIMENSION( 3 ) :: AXISZ              ! Body-fixed axis of rotation along z-direction
REAL*8, DIMENSION( 3 ) :: AXISN              ! Body-fixed axis of rotation (for the initial configuration only)

! *********************************************************************************************** !
! REAL VARIABLES (MONTE CARLO PARAMETERS)                                                         !
! *********************************************************************************************** !
REAL*8 :: MAX_TRANS       ! User maximum displacement [+/-] (Translation)
REAL*8 :: MAX_ROT         ! User maximum displacement [+/-] (Rotation)
REAL*8 :: MAX_VOL         ! User maximum displacement [+/-] (Volume)
REAL*8 :: DRMAX           ! Maximum displacement [+/-] (Translation)
REAL*8 :: ANGMAX          ! Maximum displacement [+/-] (Rotation)
REAL*8 :: DVMAX           ! Maximum displacement [+/-] (Volume)
REAL*8 :: PROB_MOV        ! Movement probability
REAL*8 :: PROB_TRANS      ! Movement probability (Translation)
REAL*8 :: PROB_ROT        ! Movement probability (Rotation)
REAL*8 :: PROB_VOL        ! Volume change probability
REAL*8 :: PROB_VOL_ISO    ! Volume change probability (Isotropic)
REAL*8 :: PROB_VOL_ANISO  ! Volume change probability (Anisotropic)
REAL*8 :: PROB_MOV_INIT   ! Movement probability for the initial configuration
REAL*8 :: PROB_TRANS_INIT ! Movement probability for the initial configuration (Translation)
REAL*8 :: PROB_ROT_INIT   ! Movement probability for the initial configuration (Rotation)
REAL*8 :: PROB_VOL_INIT   ! Volume change probability for the initial configuration (Isotropic)
REAL*8 :: DRMAX_INIT      ! Maximum displacement [+/-] (Translation)
REAL*8 :: ANGMAX_INIT     ! Maximum displacement [+/-] (Rotation)
REAL*8 :: DVMAX_INIT      ! Maximum displacement [+/-] (Volume)
REAL*8 :: DVMIN_INIT      ! Minimum displacement [+]   (Volume)
REAL*8 :: R_ACC_T         ! Acceptance ratio threshold (Translation)
REAL*8 :: R_ACC_R         ! Acceptance ratio threshold (Rotation)
REAL*8 :: R_ACC_V         ! Acceptance ratio threshold (Volume)

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL*8, DIMENSION( : ),    ALLOCATABLE :: DIAMETER      ! Diameter
REAL*8, DIMENSION( : ),    ALLOCATABLE :: LENGTH        ! Length
REAL*8, DIMENSION( : ),    ALLOCATABLE :: MOLAR_F       ! Molar fraction
REAL*8, DIMENSION( : ),    ALLOCATABLE :: PARTICLE_VOL  ! Particle volume of component c
REAL*8, DIMENSION( : ),    ALLOCATABLE :: RHO_PARTICLE  ! Number density of component c
REAL*8, DIMENSION( : ),    ALLOCATABLE :: ASPECT_RATIO  ! Aspect ratio of component c
REAL*8, DIMENSION( :, : ), ALLOCATABLE :: Q, QMC        ! Rotation quaternion array
REAL*8, DIMENSION( :, : ), ALLOCATABLE :: R, RMC        ! Position array
REAL*8, DIMENSION( :, : ), ALLOCATABLE :: E, EMC        ! Orientation array

! *********************************************************************************************** !
! REAL PARAMETERS                                                                                 !
! *********************************************************************************************** !
REAL*8, PARAMETER :: PI = 4.D0 * DATAN ( 1.D0 ) ! Pi number
REAL*8, PARAMETER :: C_BOLTZMANN = 1.380649D-23 ! Boltzmann constant

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 14 ) :: CODE_DESCRIPTOR  ! Descriptor for output file (initial configuration)
CHARACTER( LEN= 02 ) :: VOL_TYPE         ! Expansion/Compression type
CHARACTER( LEN= 10 ) :: DESCRIPTOR_FILE1 ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DESCRIPTOR_FILE2 ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DESCRIPTOR_FILE3 ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DESCRIPTOR_HOUR  ! Descriptor for output file
CHARACTER( LEN= 10 ) :: DESCRIPTOR_DATE  ! Descriptor for output folder
CHARACTER( LEN= 32 ) :: FORMAT_FILE1     ! String format for output file
CHARACTER( LEN= 32 ) :: FORMAT_FILE2     ! String format for output file
CHARACTER( LEN= 32 ) :: FORMAT_FILE3     ! String format for output file
CHARACTER( LEN= 32 ) :: FORMAT_HOUR      ! String format for output file
CHARACTER( LEN= 32 ) :: FORMAT_DATE      ! String format for output folder
CHARACTER( LEN= 64 ) :: FORMAT_SELF      ! String format (general)
CHARACTER( LEN= 01 ) :: TRAJ_INQ         ! Trajectory output inquiry
CHARACTER( LEN= 03 ) :: CONFIG_INQ       ! Molecular structure inquiry
CHARACTER( LEN= 03 ) :: GEOM_INQ         ! Molecular geometry inquiry
CHARACTER( LEN= 32 ) :: GET              ! Variable names (.ini input files)
CHARACTER( LEN= 32 ) :: CONFIGURATION    ! Molecular structure
CHARACTER( LEN= 32 ) :: GEOMETRY         ! Molecular geometry (extended)
CHARACTER( LEN= 03 ) :: GEO_ACRONYM      ! Molecular geometry (acronym)
CHARACTER( LEN= 01 ) :: UNROT_AXIS       ! Unrotated axis inquiry
CHARACTER( LEN= 03 ) :: MC_ENSEMBLE      ! Monte Carlo Ensemble
CHARACTER( LEN= 01 ) :: SEED_INQ         ! Fixed/random seed inquiry
CHARACTER( LEN= 01 ) :: DUMMY            ! Dummy (character)

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
LOGICAL                 :: FILE_EXIST   ! Checks whether a file exists or not
LOGICAL                 :: TRAJ_CHECK   ! Checks whether the trajectory file will be written out
LOGICAL                 :: INIT_CONF    ! Checks whether the initial configuration will be read from an external file
LOGICAL                 :: FSEED        ! Checks whether the seed for the random number generator will be fixed or random
LOGICAL, DIMENSION( 5 ) :: CONFIG_SELEC ! Checks the selected molecular configuration
LOGICAL, DIMENSION( 3 ) :: GEOM_SELEC   ! Checks the selected molecular geometry
LOGICAL, DIMENSION( 3 ) :: AXIS_SELEC   ! Checks the selected unrotated reference axis

END MODULE GLOBALVAR