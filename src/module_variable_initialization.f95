! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!  This module initialize common system variables (number of particles, reduced number density,   !
!     reduced temperature etc.), molecular properties (geometry and dimensions), Monte Carlo      !
!  parameters (ensemble type, total number of cycles, number of equilibration cycles etc.), and   !
!                                 potential parameters (if any).                                  !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
! which results will be written out in external files and to enable post-processing subroutines.  !
!                                                                                                 !
! Version number: 1.1.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                        August 25th, 2023                                        !
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
!              All variables should be first specified in an input file (*.ini file)              !
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
OPEN( UNIT= 10, FILE= "ini_montecarlo.ini", ACTION= "READ" )

! Total number of cycles
READ( 10, * ) GET, MAX_CYCLES
! Condition
IF( MAX_CYCLES < 1 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of cycles [", MAX_CYCLES, "] cannot be less than 1. Exiting... "
  CALL EXIT(  )
END IF

! Number of equilibration cycles
READ( 10, * ) GET, N_EQUIL
! Condition 1
IF( N_EQUIL >= MAX_CYCLES ) THEN
  WRITE( *, "(6G0)" ) "The number of equilibration cycles [", N_EQUIL, "] cannot be greater than or equal to the maximum number ", &
  &                   "of cycles [", MAX_CYCLES, "]. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( N_EQUIL < 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of equilibration cycles [", N_EQUIL, "] cannot be a negative integer. Exiting... "
  CALL EXIT(  )
END IF

! Saving frequency
READ( 10, * ) GET, N_SAVE
! Condition
IF( N_SAVE < 1 ) THEN
  WRITE( *, "(3G0)" ) "The saving frequency [", N_SAVE, "] cannot be a negative integer nor zero. Exiting... "
  CALL EXIT(  )
END IF

! Adjustment frequency
READ( 10, * ) GET, N_ADJUST
! Condition
IF( N_ADJUST < 1 ) THEN
  WRITE( *, "(4G0)" ) "The adjustment frequency of the simulation [", N_ADJUST, "] cannot be a negative integer nor zero. ", &
  &                   "Exiting... "
  CALL EXIT(  )
END IF

! Adjustment frequency (random configuration)
READ( 10, * ) GET, N_ADJUST_INIT
! Condition
IF( N_ADJUST_INIT < 1 ) THEN
  WRITE( *, "(4G0)" ) "The adjustment frequency of the random configuration [", N_ADJUST_INIT, "] cannot be a negative integer ", &
  &                   "nor zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum translational displacement
READ( 10, * ) GET, MAX_TRANS
! Condition
IF( DABS( MAX_TRANS - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum translational displacement of the simulation [", MAX_TRANS, "] cannot be zero. ", &
  &                           "Exiting... "
  CALL EXIT(  )
END IF

! Maximum translational displacement (random configuration)
READ( 10, * ) GET, DRMAX_INIT
! Condition
IF( DABS( DRMAX_INIT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum translational displacement of the random configuration [", DRMAX_INIT, "] cannot be ", &
  &                           "zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum rotational displacement
READ( 10, * ) GET, MAX_ROT
! Condition
IF( DABS( MAX_ROT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum rotational displacement of the simulation [", MAX_ROT, "] cannot be zero. Exiting... "
  CALL EXIT(  )
END IF


! Maximum rotational displacement (random configuration)
READ( 10, * ) GET, ANGMAX_INIT
! Condition
IF( DABS( ANGMAX_INIT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum rotational displacement of the random configuration [", ANGMAX_INIT, "] cannot be ", &
  &                           "zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum isotropic volume change
READ( 10, * ) GET, MAX_VOLI
! Condition
IF( DABS( MAX_VOLI - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum isotropic volume change of the simulation [", MAX_VOLI, "] cannot be zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum anisotropic volume change
READ( 10, * ) GET, MAX_VOLA
! Condition
IF( DABS( MAX_VOLA - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum anisotropic volume change of the simulation [", MAX_VOLA, "] cannot be zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum isotropic volume change (random configuration)
READ( 10, * ) GET, DVMAXISO_INIT
! Condition
IF( DABS( DVMAXISO_INIT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum isotropic volume change of the random configuration [", DVMAXISO_INIT, &
  &                          "] cannot be zero. Exiting... "
  CALL EXIT(  )
END IF

! Maximum anisotropic volume change (random configuration)
READ( 10, * ) GET, DVMAXANISO_INIT
! Condition
IF( DABS( DVMAXANISO_INIT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum anisotropic volume change of the random configuration [", DVMAXANISO_INIT, &
  &                          "] cannot be zero. Exiting... "
  CALL EXIT(  )
END IF

! Minimum volumetric displacement (random configuration)
READ( 10, * ) GET, DVMIN_INIT
! Condition
IF( DVMIN_INIT >= DVMAXISO_INIT .OR. DVMIN_INIT >= DVMAXANISO_INIT ) THEN
  WRITE( *, "(3(G0,G0.5),G0)" ) "The minimum volume change of the random configuration [", DVMIN_INIT, &
  &                             "] cannot be greater than or equal to the maximum isotropic [", DVMAXISO_INIT, &
  &                             "] or anisotropic [", DVMAXANISO_INIT, "] volume change. Exiting... "
  CALL EXIT(  )
END IF

! Maximum box distortion
READ( 10 , * ) GET, BOX_DIST
! Condition
IF( BOX_DIST < 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The maximum box distortion [", BOX_DIST, "] cannot be less than 1. Exiting... "
  CALL EXIT(  )
END IF

! Lattice reduction method
READ( 10 , * ) GET, LRTYPE
CALL TO_UPPER( LRTYPE, LEN_TRIM( LRTYPE ), LRTYPE )
LRED_SELEC(:) = .FALSE.
! Lattice reduction: Gottwald method
IF( LRTYPE == "FBM" ) THEN
  LRED_SELEC(1) = .TRUE.
! Lattice reduction: Lenstra-Lenstra-Lovász method
ELSE IF( LRTYPE == "LLL" ) THEN
  LRED_SELEC(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( LRTYPE ), "] is not an available lattice reduction method. Exiting... "
  CALL EXIT(  )
END IF

! Simulation ensemble
READ( 10, * ) GET, MC_ENSEMBLE
CALL TO_UPPER( MC_ENSEMBLE, LEN_TRIM( MC_ENSEMBLE ), MC_ENSEMBLE )
! Condition
IF( MC_ENSEMBLE /= "NVT" .AND. MC_ENSEMBLE /= "NPT" ) THEN
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( MC_ENSEMBLE ), "] is not an available simulation ensemble. Exiting... "
  CALL EXIT(  )
END IF
WRITE( *, "(3G0)" ) "Monte Carlo ensemble: [", MC_ENSEMBLE, "]"
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
OPEN( UNIT= 10, FILE= "ini_control.ini", ACTION= "READ" )

! Trajectory inquiry
READ( 10, * ) GET, TRAJ_INQ
CALL TO_UPPER( TRAJ_INQ, LEN_TRIM( TRAJ_INQ ), TRAJ_INQ )
! Transforms characters into logical variables
IF( TRAJ_INQ == "Y" ) THEN
  TRAJ_CHECK = .TRUE.
ELSE
  TRAJ_CHECK = .FALSE.
END IF

! Fixed/Random seed inquiry
READ( 10, * ) GET, SEED_INQ
CALL TO_UPPER( SEED_INQ, LEN_TRIM( SEED_INQ ), SEED_INQ )
! Transforms characters into logical variables
IF( SEED_INQ == "Y" ) THEN
  FSEED = .TRUE.
ELSE
  FSEED = .FALSE.
END IF

! Potential type
READ( 10, * ) GET, POTENTIAL_TYPE
CALL TO_UPPER( POTENTIAL_TYPE, LEN_TRIM( POTENTIAL_TYPE ), POTENTIAL_TYPE )
! Transforms characters into logical variables
POTENTIAL_SELEC(:) = .FALSE.
IF( POTENTIAL_TYPE == "HARDCORE" ) THEN
  POTENTIAL_SELEC(1) = .TRUE.
ELSE IF( POTENTIAL_TYPE == "SQUAREWELL" ) THEN
  POTENTIAL_SELEC(2) = .TRUE.
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined [", TRIM( POTENTIAL_TYPE ), "] is not an available force field. Exiting..."
  CALL EXIT(  )
END IF

! Potential inquiry
IF( .NOT. POTENTIAL_SELEC(1) ) THEN
  READ( 10, * ) GET, POT_INQ
  CALL TO_UPPER( POT_INQ, LEN_TRIM( POT_INQ ), POT_INQ )
  ! Transforms characters into logical variables
  IF( POT_INQ == "Y" ) THEN
    POTENTIAL_CHECK = .TRUE.
  ELSE
    POTENTIAL_CHECK = .FALSE.
  END IF
END IF

! TPT coefficients inquiry
IF( .NOT. POTENTIAL_SELEC(1) ) THEN
  READ( 10, * ) GET, COEF_INQ
  CALL TO_UPPER( COEF_INQ, LEN_TRIM( COEF_INQ ), COEF_INQ )
  ! Transforms characters into logical variables
  IF( COEF_INQ == "Y" ) THEN
    COEF_CHECK = .TRUE.
  ELSE
    COEF_CHECK = .FALSE.
  END IF
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
INTEGER( KIND= INT64 ) :: C ! Counter

! *********************************************************************************************** !
! System variables                                                                                !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_system.ini", ACTION= "READ" )

! Packing fraction
READ( 100, * ) GET, PACKING_F
! Condition 1
IF( PACKING_F <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction [", PACKING_F, "] cannot be less than or equal to 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PACKING_F > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The packing fraction [", PACKING_F, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF

! Number of components
READ( 100, * ) GET, COMPONENTS
! Condition
IF( COMPONENTS < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of components [", COMPONENTS, "] cannot be less than 1. Exiting... "
  CALL EXIT(  )
END IF

! Allocation
ALLOCATE( DIAMETER(COMPONENTS), LENGTH(COMPONENTS), MOLAR_F(COMPONENTS) )
ALLOCATE( N_COMPONENT(0:COMPONENTS), PARTICLE_VOL(COMPONENTS), RHO_PARTICLE(COMPONENTS) )
ALLOCATE( ASPECT_RATIO(COMPONENTS) )
ALLOCATE( SIGSPHERE(COMPONENTS) )
ALLOCATE( SPHERCOMP(COMPONENTS) )
ALLOCATE( SPHCOMP_INQ(COMPONENTS) )

! Component sphericity
READ( 100, * ) GET, SPHCOMP_INQ
DO C = 1, COMPONENTS
  CALL TO_UPPER( SPHCOMP_INQ(C), LEN_TRIM( SPHCOMP_INQ(C) ), SPHCOMP_INQ(C) )
  IF( SPHCOMP_INQ(C) == "T" ) THEN
    SPHERCOMP(C) = .TRUE.
  ELSE
    SPHERCOMP(C) = .FALSE.
  END IF
END DO

! Diameter of component i
READ( 100, * ) GET, DIAMETER
! Condition
DO C = 1, COMPONENTS
  IF( DIAMETER(C) <= 0.D0 ) THEN
    WRITE( *, "(3G0,G0.5,2G0)" ) "The diameter of component #", C, " [", DIAMETER(C), "] cannot be less than or equal to 0. ", &
    &                            "Exiting... "
    CALL EXIT(  )
  END IF
END DO

! Length of component i
READ( 100, * ) GET, LENGTH
! Condition
DO C = 1, COMPONENTS
  IF( GEOM_SELEC(1) .OR. GEOM_SELEC(3) ) THEN
    IF( LENGTH(C) <= 0.D0 ) THEN
      WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", C, " [", LENGTH(C), "] cannot be less than or equal to 0. ", &
      &                            "Exiting... "
      CALL EXIT(  )
    END IF
  ELSE IF( GEOM_SELEC(2) ) THEN
    IF( LENGTH(C) < 0.D0 ) THEN
      WRITE( *, "(3G0,G0.5,G0)" ) "The length of component #", C, " [", LENGTH(C), "] cannot be less than 0. Exiting... "
      CALL EXIT(  )
    END IF
  END IF
END DO

! Check is component is spherical
DO C = 1, COMPONENTS
  IF( GEOM_SELEC(1) .AND. DABS( LENGTH(C) - DIAMETER(C) ) >= EPSILON( 1.D0 ) .AND. SPHERCOMP(C) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", C, " [", LENGTH(C), "] cannot be different from its diameter ", &
    &                               "[", DIAMETER(C), "] when the component is considered spherical under the ellipsoidal ", &
    &                               "framework. Exiting..."
    CALL EXIT(  )
  ELSE IF( GEOM_SELEC(2) .AND. DABS( LENGTH(C) - 0.D0 ) >= EPSILON( 1.D0 ) .AND. SPHERCOMP(C) ) THEN
    WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", C, " [", LENGTH(C), "] must be 0 when the component is ", &
    &                            "considered spherical under the spherocylindrical framework. Exiting..."
    CALL EXIT(  )
  ELSE IF( GEOM_SELEC(3) .AND. DABS( LENGTH(C) - DIAMETER(C) ) >= EPSILON( 1.D0 ) .AND. SPHERCOMP(C) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", C, " [", LENGTH(C), "] cannot be different from its diameter ", &
    &                               "[", DIAMETER(C), "] when the component is considered spherical under the cylindrical ", &
    &                               "framework. Exiting..."
    CALL EXIT(  )
  END IF
  IF( GEOM_SELEC(1) .AND. DABS( LENGTH(C) - DIAMETER(C) ) < EPSILON( 1.D0 ) .AND. .NOT. SPHERCOMP(C) ) THEN
    WRITE( *, "(3G0,2(G0.5,2G0))" ) "The length of component #", C, " [", LENGTH(C), "] is equal to its diameter ", &
    &                               "[", DIAMETER(C), "] even though the component is not considered spherical under the ", &
    &                               "ellipsoidal framework. Exiting..."
    CALL EXIT(  )
  ELSE IF( GEOM_SELEC(2) .AND. DABS( LENGTH(C) - 0.D0 ) < EPSILON( 1.D0 ) .AND. .NOT. SPHERCOMP(C) ) THEN
    WRITE( *, "(3G0,G0.5,2G0)" ) "The length of component #", C, " [", LENGTH(C), "] is equal to 0 even though the ", &
    &                            "component is not considered spherical under the spherocylindrical framework. Exiting..."
    CALL EXIT(  )
  END IF
END DO

! Molar fraction of component i
READ( 100, * ) GET, MOLAR_F
! Condition
DO C = 1, COMPONENTS
  IF( MOLAR_F(C) < 0.D0 ) THEN
    WRITE( *, "(3G0,G0.5,G0)" ) "The molar fraction of component #", C, " [", MOLAR_F(C), "] cannot be less than 0. Exiting... "
    CALL EXIT(  )
  END IF
END DO

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

! Number of particles
READ( 100, * ) GET, N_PARTICLES
! Condition
IF( N_PARTICLES < 2 ) THEN
  WRITE( *, "(3G0)" ) "The number of particles [", N_PARTICLES, "] cannot be less than 2. Exiting... "
  CALL EXIT(  )
END IF

! Number of particles of component i
N_COMPONENT(:) = 0
DO C = 1, COMPONENTS
  N_COMPONENT(C) = NINT( MOLAR_F(C) * DBLE(N_PARTICLES) )
END DO

! Particle volume
DO C = 1, COMPONENTS
  IF( GEOM_SELEC(1) ) THEN
    IF( .NOT. SPHERCOMP(C) ) THEN
      PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
    ELSE
      PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * DIAMETER(C)
    END IF
  ELSE IF( GEOM_SELEC(2) ) THEN
    IF( .NOT. SPHERCOMP(C) ) THEN
      PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * DIAMETER(C) + (PI / 4.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
    ELSE
      PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * DIAMETER(C)
    END IF
  ELSE IF( GEOM_SELEC(3) ) THEN
    IF( .NOT. SPHERCOMP(C) ) THEN
      PARTICLE_VOL(C) = (PI / 4.D0) * DIAMETER(C) * DIAMETER(C) * LENGTH(C)
    ELSE
      PARTICLE_VOL(C) = (PI / 6.D0) * DIAMETER(C) * DIAMETER(C) * DIAMETER(C)
    END IF
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

! Cube root check
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

! Absolute Temperature
READ( 100, * ) GET, TEMP
! Condition
IF( TEMP <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The absolute temperature [", TEMP, "] cannot be less than or equal to 0. Exiting... "
  CALL EXIT(  )
END IF

! Reduced Pressure
READ( 100, * ) GET, PRESS
! Condition
IF( PRESS <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The reduced pressure [", PRESS, "] cannot be less than or equal to 0. Exiting... "
  CALL EXIT(  )
END IF

! Maximum length ratio (box distortion)
READ( 100, * ) GET, MAX_LENGTH_RATIO
! Condition
IF( MAX_LENGTH_RATIO <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum linear distortion of the box [", MAX_LENGTH_RATIO, "] cannot be less than ", &
  &                           "or equal to 0. Exiting... "
  CALL EXIT(  )
END IF

! Maximum angle (box distortion)
READ( 100, * ) GET, MAX_ANGLE
! Condition
IF( MAX_ANGLE <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The maximum angular distortion of the box [", MAX_ANGLE, "] cannot be less than or equal to 0°. ", &
  &                           "Exiting... "
  CALL EXIT(  )
END IF
MAX_ANGLE = MAX_ANGLE * PI / 180.D0

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
OPEN( UNIT= 100, FILE= "ini_config.ini", ACTION= "READ" )

! Skip
READ( 100, * ) GET, DUMMY
READ( 100, * ) GET, DUMMY

! Unrotated reference axis (initial configuration)
READ( 100, * ) GET, UNROT_AXIS
CALL TO_UPPER( UNROT_AXIS, LEN_TRIM( UNROT_AXIS ), UNROT_AXIS )
! Transforms characters into logical variables
AXIS_SELEC(:) = .FALSE.
IF( UNROT_AXIS == "X" ) THEN
  AXIS_SELEC(1) = .TRUE.
ELSE IF( UNROT_AXIS == "Y" ) THEN
  AXIS_SELEC(2) = .TRUE.
ELSE IF( UNROT_AXIS == "Z" ) THEN
  AXIS_SELEC(3) = .TRUE.
END IF
! Condition
IF( UNROT_AXIS /= "X" .AND. UNROT_AXIS /= "Y" .AND. UNROT_AXIS /= "Z" ) THEN
  WRITE( *, "(G0)" ) "The unrotated reference axis can only be X, Y, or Z. Exiting... "
  CALL EXIT(  )
END IF

! Quaternion angle (initial configuration)
READ( 100, * ) GET, QUATERNION_ANGLE

! Random configuration
IF( CONFIG_SELEC(4) ) THEN
  ! Initial packing fraction for the random configuration
  READ( 100, * ) GET, ETA_INI
  ! Condition 1
  IF( ETA_INI <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,2G0)" ) "The packing fraction of the random configuration [", ETA_INI, "] cannot be less than or equal ", &
    &                           "to 0. Exiting... "
    CALL EXIT(  )
  END IF
  ! Condition 2
  IF( ETA_INI > 1.D0 ) THEN
    WRITE( *, "(G0,G0.5,2G0)" ) "The packing fraction of the random configuration [", ETA_INI, "] cannot be greater than 1. ", &
    &                           "Exiting... "
    CALL EXIT(  )
  END IF
  ! Initial pressure for the random configuration
  READ( 100, * ) GET, PRESS_RND
  ! Condition
  IF( PRESS_RND <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,2G0)" ) "The reduced pressure of the random configuration [", PRESS_RND, "] cannot be less than or ", &
    &                           "equal to 0. Exiting... "
    CALL EXIT(  )
  END IF
ELSE IF( .NOT. CONFIG_SELEC(4) ) THEN
  READ( 100, * ) GET, DUMMY
  READ( 100, * ) GET, DUMMY
END IF

! Last Frame
READ( 100, * ) GET, INIT_CONF

CLOSE( 100 )

! *********************************************************************************************** !
! Simulation variables                                                                            !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_ratios.ini", ACTION= "READ" )

! Acceptance ratio (translation)
READ( 100, * ) GET, R_ACC_T
! Condition 1
IF( R_ACC_T <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The translational acceptance ratio [", R_ACC_T, "] cannot be less than or equal to 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( R_ACC_T > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The translational acceptance ratio [", R_ACC_T, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF

! Acceptance ratio (rotation)
READ( 100, * ) GET, R_ACC_R
! Condition 1
IF( R_ACC_R <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The rotational acceptance ratio [", R_ACC_R, "] cannot be less than or equal to 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( R_ACC_R > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The rotational acceptance ratio [", R_ACC_R, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF

! Acceptance ratio (isotropic volume change)
READ( 100, * ) GET, R_ACC_VI
! Condition 1
IF( R_ACC_VI <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The volumetric acceptance ratio (isotropic) [", R_ACC_VI, "] cannot be less than or equal to 0. ", &
  &                           "Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( R_ACC_VI > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (isotropic) [", R_ACC_VI, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF

! Acceptance ratio (anisotropic volume change)
READ( 100, * ) GET, R_ACC_VA
! Condition 1
IF( R_ACC_VA <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "The volumetric acceptance ratio (anisotropic) [", R_ACC_VA, "] cannot be less than or equal to ", &
  &                           "0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( R_ACC_VA > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The volumetric acceptance ratio (anisotropic) [", R_ACC_VA, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF

CLOSE( 100 )

! *********************************************************************************************** !
! Simulation variables                                                                            !
! *********************************************************************************************** !
OPEN( UNIT= 100, FILE= "ini_probabilities.ini", ACTION= "READ" )

! Movement/Volume change probability
READ( 100, * ) GET, PROB_MOV
! Condition 1
IF( PROB_MOV < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement [", PROB_MOV, "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_MOV > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement [", PROB_MOV, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_VOL = 1.D0 - PROB_MOV

! Movement/Volume change probability (initial configuration)
READ( 100, * ) GET, PROB_MOV_INIT
! Condition 1
IF( PROB_MOV_INIT < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement of the random configuration [", PROB_MOV_INIT, &
  &                          "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_MOV_INIT > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of movement of the random configuration [", PROB_MOV_INIT, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_VOL_INIT = 1.D0 - PROB_MOV_INIT

! Translational/Rotational movement probability
READ( 100, * ) GET, PROB_TRANS
! Condition 1
IF( PROB_TRANS < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements [", PROB_TRANS, "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_TRANS > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements [", PROB_TRANS, "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_ROT = 1.D0 - PROB_TRANS

! Translational/Rotational movement probability (initial configuration)
READ( 100, * ) GET, PROB_TRANS_INIT
! Condition 1
IF( PROB_TRANS_INIT < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements of the random configuration [", PROB_TRANS_INIT, &
  &                          "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_TRANS_INIT > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of translational movements of the random configuration [", PROB_TRANS_INIT, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_ROT_INIT = 1.D0 - PROB_TRANS_INIT

! Isotropic/Anisotropic volume change
READ( 100, * ) GET, PROB_VOL_ISO
! Condition 1
IF( PROB_VOL_ISO < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes [", PROB_VOL_ISO, "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_VOL_ISO > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes [", PROB_VOL_ISO, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_VOL_ANISO = 1.D0 - PROB_VOL_ISO

! Isotropic/Anisotropic volume change (initial configuration)
READ( 100, * ) GET, PROB_ISO_INIT
! Condition 1
IF( PROB_ISO_INIT < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes of the random configuration [", PROB_ISO_INIT, &
  &                          "] cannot be less than 0. Exiting... "
  CALL EXIT(  )
END IF
! Condition 2
IF( PROB_ISO_INIT > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "The probability of isotropic volume changes of the random configuration [", PROB_ISO_INIT, &
  &                          "] cannot be greater than 1. Exiting... "
  CALL EXIT(  )
END IF
PROB_ANISO_INIT = 1.D0 - PROB_ISO_INIT

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
  WRITE( *, "(G0,G0.5)" ) "Target Reduced Pressure (Random Configuration): ", PRESS_RND
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0,G0)" ) "Adjustment Frequency (Random Configuration): Every ", N_ADJUST_INIT, " Cycle(s)"
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5)" ) "Maximum Rotation (Random Configuration): ", DRMAX_INIT
  WRITE( *, "(G0,G0.5)" ) "Maximum Translation (Random Configuration): ", ANGMAX_INIT
  WRITE( *, "(G0,G0.5)" ) "Maximum Isotropic Volume Change (Random Configuration): ", DVMAXISO_INIT
  WRITE( *, "(G0,G0.5)" ) "Maximum Anisotropic Volume Change (Random Configuration): ", DVMAXANISO_INIT
  WRITE( *, "(G0,G0.5)" ) "Minimum Volume Change (Random Configuration): ", DVMIN_INIT
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Random Configuration): ", PROB_MOV_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Random Configuration): ", PROB_VOL_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Translation Probability (Random Configuration): ", PROB_TRANS_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Rotation Probability (Random Configuration): ", PROB_ROT_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Isotropic Volume Change Probability (Random Configuration): ", PROB_ISO_INIT * 100.D0, "%"
  WRITE( *, "(G0,G0.5,G0)" ) "Anisotropic Volume Change Probability (Random Configuration): ", PROB_ANISO_INIT * 100.D0, "%"
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
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Isotropic Volume Change): ", R_ACC_VI
WRITE( *, "(G0,G0.5)" ) "Acceptance Ratio (Anisotropic Volume Change): ", R_ACC_VA
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Rotation): ", MAX_TRANS
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Translation): ", MAX_ROT
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Isotropic Volume Change): ", MAX_VOLI
WRITE( *, "(G0,G0.5)" ) "Maximum Displacement (Anisotropic Volume Change): ", MAX_VOLA
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Movement): ", PROB_MOV * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Global Probability (Volume Change): ", PROB_VOL * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Translation): ", PROB_TRANS * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Movement Probability (Rotation): ", PROB_ROT * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Isotropic): ", PROB_VOL_ISO * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Volume Change Probability (Anisotropic): ", PROB_VOL_ANISO * 100.D0, "%"
WRITE( *, * ) " "
WRITE( *, "(G0,G0.5)" ) "Maximum Box Distortion: ", BOX_DIST
WRITE( *, "(G0,G0.5)" ) "Maximum Box Length Distortion: ", MAX_LENGTH_RATIO
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Box Angle Distortion: ", MAX_ANGLE * 180.D0 / PI, "°"
IF( LRED_SELEC(1) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Algorithm: Gottwald"
ELSE IF( LRED_SELEC(2) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Algorithm: Lenstra-Lenstra-Lovász"
END IF
WRITE( *, * ) " "
IF( TRAJ_CHECK ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [YES]"
ELSE IF( .NOT. TRAJ_CHECK ) THEN
  WRITE( *, "(G0)" ) "Trajectory Files: [NO]"
END IF
WRITE( *, * ) " "
IF( FSEED ) THEN
  WRITE( *, "(G0)" ) "Seed Type: [FIXED]"
ELSE IF( .NOT. FSEED ) THEN
  WRITE( *, "(G0)" ) "Seed Type: [RANDOM]"
END IF

WRITE( *, * ) " "

RETURN

END SUBROUTINE COMMON_VAR

SUBROUTINE POTENTIAL_VAR(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: C_LAMB ! Counter

! *********************************************************************************************** !
! Potential variables                                                                             !
! *********************************************************************************************** !
OPEN( UNIT= 10, FILE= "ini_potential.ini", ACTION= "READ" )

! Square-well potential
IF( POTENTIAL_SELEC(2) ) THEN
  ! Number of attractive range points
  READ( 10, * ) GET, N_LAMBDA
  ! Condition
  IF( N_LAMBDA < 1 ) THEN
    WRITE( *, "(3G0)" ) "The number of potential ranges [", N_LAMBDA, "] cannot be less than 1. Exiting... "
    CALL EXIT(  )
  END IF
  ! Allocation
  ALLOCATE( L_RANGE(N_LAMBDA), LFEXIST(N_LAMBDA) )
  ! Attractive range (λ)
  READ( 10, * ) GET, L_RANGE
  ! Condition
  DO C_LAMB = 1, N_LAMBDA
    IF( L_RANGE(C_LAMB) <= 1.D0 ) THEN
      WRITE( *, "(3G0,G0.5,G0)" ) "The attractive range #", C_LAMB, " [", L_RANGE(C_LAMB), &
      &                           "] cannot be less than or equal to 1. Exiting... "
      CALL EXIT(  )
    END IF
  END DO
  ! Reduced Temperature
  READ( 10, * ) GET, RED_TEMP
  ! Condition
  IF( RED_TEMP <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "The reduced temperature [", RED_TEMP, "] cannot be less than or equal to 0. Exiting... "
    CALL EXIT(  )
  END IF
  ! Minimum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
  READ( 10, * ) GET, MIN_BLOCKS
  ! Maximum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
  READ( 10, * ) GET, MAX_BLOCKS
  ! Condition 1
  IF( MIN_BLOCKS < 1 ) THEN
    WRITE( *, "(3G0)" ) "The minimum number of blocks [", MIN_BLOCKS, "] cannot be less than 1. Exiting... "
    CALL EXIT(  )
  END IF
  ! Condition 2
  IF( MAX_BLOCKS < 1 ) THEN
    WRITE( *, "(3G0)" ) "The maximum number of blocks [", MAX_BLOCKS, "] cannot be less than 1. Exiting... "
    CALL EXIT(  )
  END IF
  ! Condition 3
  IF( MIN_BLOCKS >= MAX_BLOCKS ) THEN
    WRITE( *, "(6G0)" ) "The minimum number of blocks [", MIN_BLOCKS, "] cannot be greater than or equal to the maximum ", &
    &                   "number of blocks [", MAX_BLOCKS, "]. Exiting... "
    CALL EXIT(  )
  END IF
END IF

CLOSE( 10 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"POTENTIAL DETAILS"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( POTENTIAL_SELEC(1) ) THEN
  WRITE( *, "(G0)" ) "Potential Type: [HARDCORE]"
ELSE IF( POTENTIAL_SELEC(2) ) THEN
  WRITE( *, "(G0)" ) "Potential Type: [SQUAREWELL]"
END IF
WRITE( *, * ) " "
IF( .NOT. POTENTIAL_SELEC(1) ) THEN
  WRITE( *, "(G0,G0)" ) "Number of Attractive Range Points: ", N_LAMBDA
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Attractive Range Values: ["
  DO C_LAMB = 1, N_LAMBDA
    WRITE( *, "(G0.5)", ADVANCE= "NO" ) L_RANGE(C_LAMB)
    IF( C_LAMB /= N_LAMBDA ) THEN
      WRITE( *, "(G0)", ADVANCE= "NO" ) ", "
    ELSE
      WRITE( *, "(G0)", ADVANCE= "YES" ) "]"
    END IF
  END DO
  WRITE( *, * ) " "
  WRITE( *, "(G0,G0.5)" ) "Reduced Temperature: ", RED_TEMP
  WRITE( *, * ) " "
  IF( POTENTIAL_CHECK ) THEN
    WRITE( *, "(G0)" ) "Potential Files: Production-Only"
  ELSE IF( .NOT. POTENTIAL_CHECK ) THEN
    WRITE( *, "(G0)" ) "Potential Files: Equilibration and Production"
  END IF
  WRITE( *, * ) " "
  IF( COEF_CHECK ) THEN
    WRITE( *, "(G0)" ) "Perturbation Coefficients: [YES]"
    WRITE( *, "(2G0)" ) "Minimum Number of Blocks (Block Averaging): ", MIN_BLOCKS
    WRITE( *, "(2G0)" ) "Maximum Number of Blocks (Block Averaging): ", MAX_BLOCKS
  ELSE IF( .NOT. COEF_CHECK ) THEN
    WRITE( *, "(G0)" ) "Perturbation Coefficients: [NO]"
  END IF
  WRITE( *, * ) " "
END IF

! Summary
WRITE( *, "(G0)" ) "RESUME? [Y/N]"
READ( *, * ) DUMMY
CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
IF( DUMMY /= "Y" ) THEN
  CALL EXIT(  )
END IF
WRITE( *, * ) " "

RETURN

END SUBROUTINE POTENTIAL_VAR

END MODULE INITVAR
