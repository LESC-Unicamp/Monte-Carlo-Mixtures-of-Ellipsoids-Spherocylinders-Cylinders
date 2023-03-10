! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!      This module initialize common variables (number of particles, reduced number density,      !
!     reduced temperature etc.) and Monte Carlo parameters (total number of cycles, number of     !
!                                   equilibration cycles etc.).                                   !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
! which results will be written out in external files and to enable post-processing subroutines.  !
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
! Main Reference:                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

MODULE INITVAR

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
!                                     VARIABLE INITIALIZATION                                     !
!       Most variables should be first specified in an input file ('mchgo_input.ini' file)        !
!              Some variables are on-the-fly (OTF) and must be specified by the user              !
!                               at the beginning of the simulation.                               !
!                       We provided an example .ini file to guide the user.                       !
!                     Please also check our README file for more information.                     !
! *********************************************************************************************** !

CONTAINS

! *********************************************************************************************** !
!                            Initialization of Monte Carlo parameters                             !
! *********************************************************************************************** !
SUBROUTINE MONTECARLO_VAR(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Simulation parameters                                                                           !
! *********************************************************************************************** !
OPEN( UNIT= 10, FILE= "ini_montecarlo.ini" )

! *********************************************************************************************** !
! Total number of cycles                                                                          !
! *********************************************************************************************** !
READ( 10, * ) GET, MAX_CYCLES

! *********************************************************************************************** !
! Number of equilibration cycles                                                                  !
! *********************************************************************************************** !
READ( 10, * ) GET, N_EQUIL

! Condition
IF( N_EQUIL >= MAX_CYCLES ) THEN
  WRITE( *, "(G0)" ) "Number of equilibration cycles greater than or equal to the maximum number of cycles. Exiting... "
  CALL EXIT(  )
END IF

! *********************************************************************************************** !
! Saving frequency                                                                                !
! *********************************************************************************************** !
!  [1    = highest frequency, results will be written out for each cycle ]                        !
!  [1000 = lower frequency, results will be written out every 1000 cycles]                        !
! *********************************************************************************************** !
READ( 10, * ) GET, N_SAVE

! *********************************************************************************************** !
! Maximum displacement adjustment frequency                                                       !
! *********************************************************************************************** !
!  [1    = highest frequency, maximum displacement will be adjusted for each cycle   ]            !
!  [200  = moderate frequency, maximum displacement will be adjusted every 200 cycles]            !
! *********************************************************************************************** !
READ( 10, * ) GET, N_ADJUST
READ( 10, * ) GET, N_ADJUST_INIT

! *********************************************************************************************** !
! Maximum translational displacement                                                              !
! *********************************************************************************************** !
READ( 10, * ) GET, MAX_TRANS
READ( 10, * ) GET, DRMAX_INIT

! *********************************************************************************************** !
! Maximum rotational displacement                                                                 !
! *********************************************************************************************** !
READ( 10, * ) GET, MAX_ROT
READ( 10, * ) GET, ANGMAX_INIT

! *********************************************************************************************** !
! Maximum volumetric displacement                                                                 !
! *********************************************************************************************** !
READ( 10, * ) GET, MAX_VOL
READ( 10, * ) GET, DVMAX_INIT

! *********************************************************************************************** !
! Minimum volumetric displacement (random configuration only)                                     !
! *********************************************************************************************** !
READ( 10, * ) GET, DVMIN_INIT

! *********************************************************************************************** !
! Simulation ensemble                                                                             !
! *********************************************************************************************** !
READ( 10, * ) GET, MC_ENSEMBLE
CALL TO_UPPER( MC_ENSEMBLE, LEN_TRIM( MC_ENSEMBLE ), MC_ENSEMBLE )
WRITE( *, "(G0,G0)" ) "Monte Carlo ensemble: ", MC_ENSEMBLE
WRITE( *, "(G0)" ) " "

CLOSE( 10 )

RETURN

END SUBROUTINE MONTECARLO_VAR

! *********************************************************************************************** !
!                           Initialization of Inquiry/Control variables                           !
! *********************************************************************************************** !
SUBROUTINE INQUERY_VAR(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Inquiry variables                                                                               !
! *********************************************************************************************** !
OPEN( UNIT= 10, FILE= "ini_control.ini" )

! *********************************************************************************************** !
! Trajectory inquiry                                                                              !
! *********************************************************************************************** !
!  Inquires whether a trajectory of particles shall be written out.                               !
!  Answer Y (Yes) to write out trajectories or                                                    !
!         N (No) to ignore it                                                                     !
! *********************************************************************************************** !
READ( 10, * ) GET, TRAJ_INQ
CALL TO_UPPER( TRAJ_INQ, LEN_TRIM( TRAJ_INQ ), TRAJ_INQ )

! Transforms characters into logical variables
IF ( TRAJ_INQ == "Y" ) THEN
  TRAJ_CHECK = .TRUE.
ELSE IF ( TRAJ_INQ == "N" ) THEN
  TRAJ_CHECK = .FALSE.
END IF

! *********************************************************************************************** !
! Random/fixed seed inquiry                                                                       !
! *********************************************************************************************** !
!  Inquires whether the seed used in the pseudorandom number generator subroutine will be         !
!  randomly defined (required for repeatability assays) or fixed (no repeatability).              !
!  Answer Y (Yes) to fix the seed of the pseudorandom number generator subroutine or              !
!         N (No) to randomly generate a seed every time the program is run                        !
! *********************************************************************************************** !
READ( 10, * ) GET, SEED_INQ
CALL TO_UPPER( SEED_INQ, LEN_TRIM( SEED_INQ ), SEED_INQ )

! Transforms characters into logical variables
IF ( SEED_INQ == "Y" ) THEN
  FSEED = .TRUE.
ELSE IF ( SEED_INQ == "N" ) THEN
  FSEED = .FALSE.
END IF

CLOSE( 10 )

RETURN

END SUBROUTINE INQUERY_VAR

! *********************************************************************************************** !
!                               Initialization of common variables                                !
! *********************************************************************************************** !
SUBROUTINE COMMON_VAR(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: C ! Counter

! *********************************************************************************************** !
! System variables                                                                                !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_system.ini" )
! Packing fraction
READ( 100, * ) GET, PACKING_F
! Number of components
READ( 100, * ) GET, COMPONENTS
ALLOCATE( DIAMETER(COMPONENTS), LENGTH(COMPONENTS), MOLAR_F(COMPONENTS) )
ALLOCATE( N_COMPONENT(0:COMPONENTS), PARTICLE_VOL(COMPONENTS), RHO_PARTICLE(COMPONENTS) )
ALLOCATE( ASPECT_RATIO(COMPONENTS) )
! Diameter of component i
READ( 100, * ) GET, DIAMETER
! Length of component i
READ( 100, * ) GET, LENGTH
! Molar fraction of component i
READ( 100, * ) GET, MOLAR_F

! Aspect ratio
DO C = 1, COMPONENTS
  ASPECT_RATIO(C) = LENGTH(C) / DIAMETER(C)
END DO

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"COMPONENT DETAILS"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Number of Components: ", COMPONENTS
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Diameter of Component #", C, ": ", DIAMETER(C), "Å"
END DO
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Length of Component #", C, ": ", LENGTH(C), "Å"
END DO
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Aspect Ratio of Component #", C, ": ", ASPECT_RATIO(C)
END DO
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Molar Fraction of Component #", C, ": ", MOLAR_F(C)
END DO
TOTAL_MOLAR_F = 0.D0
DO C = 1, COMPONENTS
  TOTAL_MOLAR_F = TOTAL_MOLAR_F + MOLAR_F(C)
END DO
WRITE( *, * ) " "
IF( DABS( TOTAL_MOLAR_F - 1.D0 ) >= EPSILON(1.D0) ) THEN
  WRITE( *, "(G0)" ) "Total molar fraction not equal to unit! Recalculating... "
  WRITE( *, * ) " "
  DO C = 1, COMPONENTS
    MOLAR_F(C) = MOLAR_F(C) / TOTAL_MOLAR_F
    WRITE( *, "(G0,G0,G0,G0.5)" ) "Normalized Molar Fraction of Component #", C, ": ", MOLAR_F(C)
  END DO
  WRITE( *, * ) " "
  TOTAL_MOLAR_F = 0.D0
  DO C = 1, COMPONENTS
    TOTAL_MOLAR_F = TOTAL_MOLAR_F + MOLAR_F(C)
  END DO
END IF

! *********************************************************************************************** !
! System variables                                                                                !
! *********************************************************************************************** !
! Number of particles
READ( 100, * ) GET, N_PARTICLES
! Number of particles of component i
N_COMPONENT(:) = 0
DO C = 1, COMPONENTS
  N_COMPONENT(C) = NINT( MOLAR_F(C) * DBLE(N_PARTICLES) )
END DO

! Particle volume
DO C = 1, COMPONENTS
  IF( GEOM_SELEC(1) ) THEN
    PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
  ELSE IF( GEOM_SELEC(2) ) THEN
    PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * DIAMETER(C) + (PI / 4.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
  ELSE IF( GEOM_SELEC(3) ) THEN
    PARTICLE_VOL(C) = (PI / 4.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
  END IF
END DO

! Total particle volume
TOTAL_VP = 0.D0
DO C = 1, COMPONENTS
  TOTAL_VP = TOTAL_VP + DBLE( N_COMPONENT(C) ) * PARTICLE_VOL(C)
END DO

! Box volume
BOX_VOLUME = TOTAL_VP / PACKING_F

! Number density
TOTAL_RHO = 0.D0
DO C = 1, COMPONENTS
  RHO_PARTICLE(C) = DBLE( N_COMPONENT(C) ) / BOX_VOLUME
  TOTAL_RHO = TOTAL_RHO + RHO_PARTICLE(C)
END DO

! Summary
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5)" ) "Number of Particles of Component #", C, ": ", N_COMPONENT(C)
END DO
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Molecular Volume of Component #", C, ": ", PARTICLE_VOL(C), "Å³"
END DO
WRITE( *, * ) " "
DO C = 1, COMPONENTS
  WRITE( *, "(G0,G0,G0,G0.5,G0)" ) "Number Density of Component #", C, ": ", RHO_PARTICLE(C), "Å⁻³"
END DO
WRITE( *, * ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"GLOBAL DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( SUM( N_COMPONENT ) == N_PARTICLES ) THEN
  WRITE( *, "(G0,G0)" ) "Number of Particles: ", N_PARTICLES
ELSE IF( SUM( N_COMPONENT ) /= N_PARTICLES ) THEN
  WRITE( *, "(G0,G0,G0,G0,G0,G0,G0)" ) "Molar-based number of particles (", SUM( N_COMPONENT ), ") not equal the ", &
  &                                    "user-defined number of particles (", N_PARTICLES, "). Overwriting... "
  WRITE( *, * ) " "
  N_PARTICLES = SUM( N_COMPONENT )
  WRITE( *, "(G0,G0)" ) "Overwritten Number of Particles: ", N_PARTICLES
END IF
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Total Number Density: ", TOTAL_RHO, "Å⁻³"
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Total Molecular Volume: ", TOTAL_VP, "Å³"

! *********************************************************************************************** !
! Cube root check                                                                                 !
! *********************************************************************************************** !
!  Checks whether the number of particles entered by the user is valid based on each molecular    !
!  cubic structure.                                                                               !
! *********************************************************************************************** !
PARTICLE_LOOP: DO
  IF ( CONFIG_SELEC(1) ) THEN
    CHECK_ROOT = DBLE( N_PARTICLES ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CHECK_ROOT - DNINT( CHECK_ROOT ) ) <= 1.D-10 ) THEN
      EXIT PARTICLE_LOOP
    ELSE
      WRITE( *, * ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( CONFIGURATION )//" configuration."
      WRITE( *, "(G0)") "The total number of particles must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL EXIT(  )
    END IF
  ELSE IF ( CONFIG_SELEC(2) ) THEN
    CHECK_ROOT = ( 0.5D0 * DBLE( N_PARTICLES ) ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CHECK_ROOT - DNINT( CHECK_ROOT ) ) <= 1.D-10 ) THEN
      EXIT PARTICLE_LOOP
    ELSE
      WRITE( *, * ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( CONFIGURATION )//" configuration."
      WRITE( *, "(G0)") "The total number of particles divided by 2 must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL EXIT(  )
    END IF
  ELSE IF ( CONFIG_SELEC(3) ) THEN
    CHECK_ROOT = ( 0.25D0 * DBLE( N_PARTICLES ) ) ** ( 1.D0 / 3.D0 )
    IF ( DABS( CHECK_ROOT - DNINT( CHECK_ROOT ) ) <= 1.D-10 ) THEN
      EXIT PARTICLE_LOOP
    ELSE
      WRITE( *, * ) " "
      WRITE( *, "(G0)") "Invalid number of particles for the "//TRIM( CONFIGURATION )//" configuration."
      WRITE( *, "(G0)") "The total number of particles divided by 4 must be a perfect cube root."
      WRITE( *, "(G0)") "Exiting..."
      CALL EXIT(  )
    END IF
  ELSE IF ( CONFIG_SELEC(4) .OR. CONFIG_SELEC(5) ) THEN
    EXIT PARTICLE_LOOP
  END IF
END DO PARTICLE_LOOP

WRITE( *, * ) " "

! *********************************************************************************************** !
! System variables                                                                                !
! *********************************************************************************************** !
! Absolute Temperature
READ( 100, * ) GET, TEMP
! Reduced Pressure
READ( 100, * ) GET, PRESS
CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"SYSTEM DETAILS"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( MC_ENSEMBLE == "NVT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Packing Fraction: ", PACKING_F
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Volume of the Simulation Box: ", BOX_VOLUME, "Å³"
  WRITE( *, * ) " "
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Initial Packing Fraction: ", PACKING_F
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Initial Volume of the Simulation Box: ", BOX_VOLUME, "Å³"
  WRITE( *, * ) " "
END IF
WRITE( *, "(G0,G0.5,G0)" ) "Absolute Temperature: ", TEMP, "K"
WRITE( *, * ) " "
IF( MC_ENSEMBLE == "NPT" ) THEN
  WRITE( *, "(G0,G0.5)" ) "Reduced Pressure: ", PRESS
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Real Pressure: ", (PRESS * C_BOLTZMANN * TEMP) / 1.D-30 / 1.D6, "MPa"
END IF
WRITE( *, * ) " "

! *********************************************************************************************** !
! Initial configuration variables                                                                 !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_config.ini" )
! Dummy
READ( 100, * ) GET, DUMMY
READ( 100, * ) GET, DUMMY
! Unrotated reference axis (initial configuration)
READ( 100, * ) GET, UNROT_AXIS
CALL TO_UPPER( UNROT_AXIS, LEN_TRIM( UNROT_AXIS ), UNROT_AXIS )
! Quaternion angle (initial configuration)
READ( 100, * ) GET, QUATERNION_ANGLE
IF( CONFIG_SELEC(4) ) THEN
  ! Maximum number of attempts to distribute particles in a random configuration
  READ( 100, * ) GET, MAX_ATTEMPTS
  ! Initial packing fraction for the random configuration
  READ( 100, * ) GET, ETA_INI
  ! Initial pressure for the random configuration
  READ( 100, * ) GET, PRESS_RND
ELSE IF( .NOT. CONFIG_SELEC(4) ) THEN
  READ( 100, * ) GET, DUMMY
  READ( 100, * ) GET, DUMMY
  READ( 100, * ) GET, DUMMY
END IF
! Last Frame
READ( 100, * ) GET, INIT_CONF
CLOSE( 100 )

! Transforms characters into logical variables
AXIS_SELEC(:) = .FALSE.
IF( UNROT_AXIS == "X" ) THEN
  AXIS_SELEC(1) = .TRUE.
ELSE IF( UNROT_AXIS == "Y" ) THEN
  AXIS_SELEC(2) = .TRUE.
ELSE IF( UNROT_AXIS == "Z" ) THEN
  AXIS_SELEC(3) = .TRUE.
END IF

! *********************************************************************************************** !
! Initial configuration variables (NVT/NPT-Monte Carlo)                                           !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_ratios.ini" )
READ( 100, * ) GET, R_ACC_T
READ( 100, * ) GET, R_ACC_R
READ( 100, * ) GET, R_ACC_V
CLOSE( 100 )
OPEN( UNIT= 100, FILE= "ini_probabilities.ini" )
READ( 100, * ) GET, PROB_MOV
PROB_VOL = 1.D0 - PROB_MOV
READ( 100, * ) GET, PROB_MOV_INIT
PROB_VOL_INIT = 1.D0 - PROB_MOV_INIT
READ( 100, * ) GET, PROB_TRANS
PROB_ROT = 1.D0 - PROB_TRANS
READ( 100, * ) GET, PROB_TRANS_INIT
PROB_ROT_INIT = 1.D0 - PROB_TRANS_INIT
READ( 100, * ) GET, PROB_VOL_ISO
PROB_VOL_ANISO = 1.D0 - PROB_VOL_ISO
CLOSE( 100 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"INITIAL CONFIGURATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Fixed-Body Axis: ", UNROT_AXIS
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Quaternion Angle: ", QUATERNION_ANGLE, "°"
WRITE( *, * ) " "
IF( CONFIG_SELEC(4) ) THEN
  WRITE( *, "(G0,G0.5)" ) "Initial Packing Fraction (Random Configuration): ", ETA_INI
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5)" ) "Maximum Attempts (Random Configuration): ", MAX_ATTEMPTS
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5)" ) "Target Reduced Pressure (Random Configuration): ", PRESS_RND
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0,G0)" ) "Adjustment Frequency (Random Configuration): Every ", N_ADJUST_INIT, " Cycle(s)"
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5)" ) "Maximum Rotation (Random Configuration): ", DRMAX_INIT
  WRITE( *, "(G0,G0.5)" ) "Maximum Translation (Random Configuration): ", ANGMAX_INIT
  WRITE( *, "(G0,G0.5)" ) "Maximum Volume Change (Random Configuration): ", DVMAX_INIT
  WRITE( *, "(G0,G0.5)" ) "Minimum Volume Change (Random Configuration): ", DVMIN_INIT
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Random Configuration): ", PROB_MOV_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Random Configuration): ", PROB_VOL_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Translation Probability (Random Configuration): ", PROB_TRANS_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Rotation Probability (Random Configuration): ", PROB_ROT_INIT * 100.D0, "%"
  WRITE( *, * ) " "
END IF
IF( INIT_CONF ) THEN
  WRITE( *, "(G0)" ) "Preset Initial Configuration: [YES]"
ELSE IF( .NOT. INIT_CONF ) THEN
  WRITE( *, "(G0)" ) "Preset Initial Configuration: [NO]"
END IF
WRITE( *, * ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 18 )//"MONTE CARLO DETAILS"//REPEAT( " ", 18 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0,G0)" ) "Number of Cycles: ", MAX_CYCLES
WRITE( *, "(G0,G0)" ) "Number of Equilibration Cycles: ", N_EQUIL
WRITE( *, "(G0,G0)" ) "Number of Production Cycles: ", MAX_CYCLES - N_EQUIL
WRITE( *, * ) " "
WRITE( *, "(G0,G0,G0)" ) "Saving Frequency: Every ", N_SAVE, " Cycle(s)"
WRITE( *, "(G0,G0,G0)" ) "Adjustment Frequency: Every ", N_ADJUST, " Equilibration Cycle(s)"
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Translation): ", R_ACC_T
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Rotation): ", R_ACC_R
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Volume Change): ", R_ACC_V
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Rotation): ", MAX_TRANS
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Translation): ", MAX_ROT
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Volume Change): ", MAX_VOL
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Movement): ", PROB_MOV * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Volume Change): ", PROB_VOL * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Translation): ", PROB_TRANS * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Rotation): ", PROB_ROT * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Isotropic): ", PROB_VOL_ISO * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Anisotropic): ", PROB_VOL_ANISO * 100.D0, "%"
WRITE( *, * ) " "
IF( TRAJ_CHECK ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [YES]"
ELSE IF( .NOT. TRAJ_CHECK ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [NO]"
END IF
WRITE( *, * ) " "

! Summary
WRITE( *, "(G0)" ) "RESUME? [Y/N]"
READ( *, * ) DUMMY
CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
IF( DUMMY /= "Y" ) THEN
  CALL EXIT(  )
END IF
WRITE( *, * ) " "

RETURN

END SUBROUTINE COMMON_VAR

END MODULE INITVAR