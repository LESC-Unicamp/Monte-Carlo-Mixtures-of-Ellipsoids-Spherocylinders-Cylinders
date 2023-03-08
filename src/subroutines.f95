! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains all subroutines used in the main program.                   !
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
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                               Linear congruential generator (LCG)                               !
!    This subroutine generates a random number from the uniform distribution over the range       !
!   0 ≤ x < 1. It does not take any arguments. The number generator seed has an in/out intent,    !
!             i. e., its value is changed every time the ranf() subroutine is called.             !
!                   This seed is the heart of the pseudorandom number generator                   !
!                        and also responsible for ensuring repeatability.                         !
!         See Allen and Tildesley, 2nd Edition (2017), Appendix E, for more information.          !
! *********************************************************************************************** !
SUBROUTINE RANF(  )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER PARAMETERS                                                                              !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ), PARAMETER :: L = 1029
INTEGER( KIND= REAL64 ), PARAMETER :: C = 221591
INTEGER( KIND= REAL64 ), PARAMETER :: G = 1048576

! *********************************************************************************************** !
! Finite modulus arithmetic                                                                       !
! *********************************************************************************************** !
SEED     = MOD( ( (SEED * L) + C ), G ) ! Number generator seed
RANDOM_N = DBLE( SEED ) / DBLE( G )     ! Pseudorandom number

RETURN

END SUBROUTINE RANF

! *********************************************************************************************** !
!        This subroutine takes the body-fixed orientation and the rotation quaternion and         !
!                 generates a space-fixed orientation using a 3D rotation matrix.                 !
!        See Allen and Tildesley, 2nd Edition (2017), pages 106-111 for more information.         !
! *********************************************************************************************** !
SUBROUTINE ACTIVE_TRANSFORMATION( EFIXED, QP, EP )

! Uses one module: global variables
use globalvar

implicit none

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: I, J ! Counters

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EFIXED ! Body-fixed orientation
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EP     ! Space-fixed orientation
REAL( KIND= REAL64 ), DIMENSION( 0:3 )  :: QP     ! Rotation quaternion
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: A      ! Rotation matrix
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: AT     ! Transpose of rotation matrix

! *********************************************************************************************** !
! Rotation Matrix A - Allen and Tildesley, 2nd edition, page 110                                  !
! *********************************************************************************************** !
! First row
A(1,1) = ( QP(0) * QP(0) ) + ( QP(1) * QP(1) ) - ( QP(2) * QP(2) ) - ( QP(3) * QP(3) )
A(1,2) = 2.D0 * ( ( QP(1) * QP(2) ) + ( QP(0) * QP(3) ) )
A(1,3) = 2.D0 * ( ( QP(1) * QP(3) ) - ( QP(0) * QP(2) ) )
! Second row
A(2,1) = 2.D0 * ( ( QP(1) * QP(2) ) - ( QP(0) * QP(3) ) )
A(2,2) = ( QP(0) * QP(0) ) - ( QP(1) * QP(1) ) + ( QP(2) * QP(2) ) - ( QP(3) * QP(3) )
A(2,3) = 2.D0 * ( ( QP(2) * QP(3) ) + ( QP(0) * QP(1) ) )
! Third row
A(3,1) = 2.D0 * ( ( QP(1) * QP(3) ) + ( QP(0) * QP(2) ) )
A(3,2) = 2.D0 * ( ( QP(2) * QP(3) ) - ( QP(0) * QP(1) ) )
A(3,3) = ( QP(0) * QP(0) ) - ( QP(1) * QP(1) ) - ( QP(2) * QP(2) ) + ( QP(3) * QP(3) )

! *********************************************************************************************** !
! Transpose of rotation matrix                                                                    !
! *********************************************************************************************** !
DO I = 1, 3
  DO J = 1, 3
    AT(I,J) = A(J,I)
  END DO
END DO

! *********************************************************************************************** !
! Active tranformation (dot product of body-fixed orientation and transpose of rotation matrix)   !
! *********************************************************************************************** !
EP(1) = ( EFIXED(1) * AT(1,1) ) + ( EFIXED(2) * AT(1,2) ) + ( EFIXED(3) * AT(1,3) )
EP(2) = ( EFIXED(1) * AT(2,1) ) + ( EFIXED(2) * AT(2,2) ) + ( EFIXED(3) * AT(2,3) )
EP(3) = ( EFIXED(1) * AT(3,1) ) + ( EFIXED(2) * AT(3,2) ) + ( EFIXED(3) * AT(3,3) )

RETURN

END SUBROUTINE ACTIVE_TRANSFORMATION

! *********************************************************************************************** !
!        This subroutine takes a rotation quaternion (qm) and combines it with a randomly         !
!        generated quaternion (qr) through quaternion multiplication, creating a randomly         !
!                                composed rotation quaternion (qn)                                !
! *********************************************************************************************** !
SUBROUTINE COMPOSED_QUATERNION( QM, QN, ANGLE )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                   :: ANGLE ! Maximum angular displacement
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QM    ! Reference rotation quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QR    ! Random quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QN    ! Composed rotation quaternion

! *********************************************************************************************** !
! Random quaternion generator                                                                     !
! *********************************************************************************************** !
CALL RANDOM_QUATERNION( QR, ANGLE )

! *********************************************************************************************** !
! Quaternion multiplication (composed rotation)                                                   !
! *********************************************************************************************** !
CALL MULTIPLY_QUATERNIONS( QR, QM, QN )

RETURN

END SUBROUTINE COMPOSED_QUATERNION

! *********************************************************************************************** !
!        This subroutine generates a random quaternion from a random angle and random axis        !
! *********************************************************************************************** !
SUBROUTINE RANDOM_QUATERNION( QR, ANGLE )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                   :: ANGLE ! Maximum angular displacement
REAL( KIND= REAL64 ), DIMENSION( 3 )   :: SR    ! Random vector
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QR    ! Random quaternion

! *********************************************************************************************** !
! Random rotation angle                                                                           !
! *********************************************************************************************** !
CALL RANF(  )
RANDOM_ANGLE = ( (2.D0 * RANDOM_N) - 1.D0 ) * ANGLE ! Range [-angmax,angmax]

! *********************************************************************************************** !
! Random vector generator                                                                         !
! *********************************************************************************************** !
CALL RANDOM_VECTOR( SR )

! *********************************************************************************************** !
! Quaternion algebra                                                                              !
! *********************************************************************************************** !
QR(0) = DCOS( RANDOM_ANGLE * 0.5D0 )         ! Real part
QR(1) = DSIN( RANDOM_ANGLE * 0.5D0 ) * SR(1) ! Imaginary part (Vector)
QR(2) = DSIN( RANDOM_ANGLE * 0.5D0 ) * SR(2) ! Imaginary part (Vector)
QR(3) = DSIN( RANDOM_ANGLE * 0.5D0 ) * SR(3) ! Imaginary part (Vector)

RETURN

END SUBROUTINE RANDOM_QUATERNION

! *********************************************************************************************** !
!            This subroutine generates a random vector on the surface of a unit sphere            !
!           See Allen and Tildesley, 2nd Edition (2017), page 514 for more information.           !
!        (Routine 'maths_module.f90' of Code A.1. in Marsaglia, Ann. Math. Statist., 1972)        !
! *********************************************************************************************** !   
SUBROUTINE RANDOM_VECTOR( SR )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: CSI_1, CSI_2, ZETA ! Random numbers
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: SR                 ! Random vector

! *********************************************************************************************** !
! Marsaglia's routine                                                                             !
! *********************************************************************************************** !
ZETA_LOOP: DO
  ! Uniform random number, ξ₁
  CALL RANF(  )
  CSI_1 = (2.D0 * RANDOM_N) - 1.D0
  ! Uniform random number, ξ₂
  CALL RANF(  )
  CSI_2 = (2.D0 * RANDOM_N) - 1.D0
  ! Sum of squares
  ZETA  = (CSI_1 * CSI_1) + (CSI_2 * CSI_2)
  ! Marseglia's criterion
  IF ( ZETA < 1.D0 ) THEN
    EXIT ZETA_LOOP
  END IF
END DO ZETA_LOOP

! *********************************************************************************************** !
! Random vector                                                                                   !
! *********************************************************************************************** !
SR(1) = (2.D0 * CSI_1) * DSQRT(1.D0 - ZETA)
SR(2) = (2.D0 * CSI_2) * DSQRT(1.D0 - ZETA)
SR(3) = 1.D0 - (2.D0 * ZETA)

RETURN

END SUBROUTINE RANDOM_VECTOR

! *********************************************************************************************** !
!             This subroutine creates a composed rotation quaternion by multiplying a             !
!                          reference quaternion and a random quaternion                           !
! *********************************************************************************************** !
SUBROUTINE MULTIPLY_QUATERNIONS( QR, QM, QN )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QR ! Random quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QM ! Reference quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 ) :: QN ! Composed quaternion

! *********************************************************************************************** !
! Cross product of quaternions (qr × qm)                                                          !
! *********************************************************************************************** !
QN(0) = ( QR(0) * QM(0) ) - ( QR(1) * QM(1) ) - ( QR(2) * QM(2) ) - ( QR(3) * QM(3) )
QN(1) = ( QR(1) * QM(0) ) + ( QR(0) * QM(1) ) - ( QR(3) * QM(2) ) + ( QR(2) * QM(3) )
QN(2) = ( QR(2) * QM(0) ) + ( QR(3) * QM(1) ) + ( QR(0) * QM(2) ) - ( QR(1) * QM(3) )
QN(3) = ( QR(3) * QM(0) ) - ( QR(2) * QM(1) ) + ( QR(1) * QM(2) ) + ( QR(0) * QM(3) )

RETURN

END SUBROUTINE MULTIPLY_QUATERNIONS

! *********************************************************************************************** !
! This subroutine checks if a random displacement (translation or rotation) of a fixed particle i !
!                           causes any overlaps with other particles j                            !
! *********************************************************************************************** !
SUBROUTINE CHECK_OVERLAP( CI, I, QI, EI, RI, CONTACT_D, BL, BLI, OVERLAP )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: I, J      ! Counters
INTEGER( KIND= REAL64 ) :: C, CI, CJ ! Component indexes

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                          :: RIJSQ     ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                          :: CONTACT_D ! Contact distance (Perram-Wertheim or Vega-Lago Methods)
REAL( KIND= REAL64 )                          :: CUTOFF_D  ! Cutoff distance
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: DLAMBDAEI ! Auxiliar vector (cylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: DMUEJ     ! Auxiliar vector (cylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RIJ       ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RI, RJ    ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: EI, EJ    ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: S12       ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BL        ! Box length
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BLI       ! Box length (inverse)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )        :: QI, QJ    ! Quaternions of particles i and j
REAL( KIND= REAL64 ), DIMENSION( COMPONENTS ) :: CUTOFF    ! Cutoff diameter

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP     ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_HER ! Detects overlap between two ellipsoids of revolution: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_SPC ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_CYL ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PARALLEL    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
OVERLAP     = .FALSE.
OVERLAP_HER = .FALSE.
OVERLAP_CYL = .FALSE.

! Diameter of circumscribing sphere
IF( GEOM_SELEC(1) ) THEN
  DO C = 1, COMPONENTS
    IF( ASPECT_RATIO(C) > 0.D0 .AND. ASPECT_RATIO(C) <= 1.D0 ) THEN
      CUTOFF(C) = DIAMETER(C)
    ELSE IF( ASPECT_RATIO(C) > 1.D0 ) THEN
      CUTOFF(C) = LENGTH(C)
    END IF
  END DO
ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
  DO C = 1, COMPONENTS
    CUTOFF(C) = DIAMETER(C) + LENGTH(C)
  END DO
END IF

! *********************************************************************************************** !
! Component and Particle Loops (Component index less than Ci)                                     !
! *********************************************************************************************** !
DO CJ = 1, CI - 1
  ! Unique loop takes only particles whose component indexes are less than Ci
  DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, SUM( N_COMPONENT(0:CJ) )
    ! Position of particle j
    RJ(1)  = RMC(1,J)
    RJ(2)  = RMC(2,J)
    RJ(3)  = RMC(3,J)
    ! Orientation of particle j
    EJ(1)  = EMC(1,J)
    EJ(2)  = EMC(2,J)
    EJ(3)  = EMC(3,J)
    ! Quaternion of particle j
    QJ(0)  = QMC(0,J)
    QJ(1)  = QMC(1,J)
    QJ(2)  = QMC(2,J)
    QJ(3)  = QMC(3,J)
    ! Vector distance between particles i and j
    RIJ(1) = RJ(1) - RI(1)
    RIJ(2) = RJ(2) - RI(2)
    RIJ(3) = RJ(3) - RI(3)
    ! Minimum Image Convention
    CALL MULTI_MATRIX( BLI, RIJ, S12 )
    S12 = S12 - ANINT(S12)
    CALL MULTI_MATRIX( BL, S12, RIJ )
    ! Magnitude of the vector distance (squared)
    RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
    ! Cutoff distance (squared)
    CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
    CUTOFF_D = CUTOFF_D * CUTOFF_D
    ! Preliminary test (circumscribing spheres)
    IF( RIJSQ <= CUTOFF_D ) THEN
      ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
      IF( GEOM_SELEC(1) ) THEN
        CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CONTACT_D, OVERLAP_HER )
        ! Overlap criterion
        IF( OVERLAP_HER ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      ! Overlap test for spherocylinders (Vega-Lago Method)
      ELSE IF( GEOM_SELEC(2) ) THEN
        CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
        ! Overlap criterion
        IF( OVERLAP_SPC ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      ! Overlap test for cylinders (Lopes et al. Method)
      ELSE IF( GEOM_SELEC(3) ) THEN
        ! Preliminary test (circumscribing spherocylinders)
        OVERLAP_SPC = .FALSE.
        CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
        ! Overlap criterion
        IF( OVERLAP_SPC ) THEN
          RJ(1) = RI(1) + RIJ(1)
          RJ(2) = RI(2) + RIJ(2)
          RJ(3) = RI(3) + RIJ(3)
          CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, DLAMBDAEI, DMUEJ, PARALLEL, OVERLAP_CYL )
          ! Overlap criterion
          IF( OVERLAP_CYL ) THEN
            ! Overlap detected
            OVERLAP = .TRUE.
            RETURN
          END IF
        END IF
      END IF
    END IF
  END DO
END DO

! *********************************************************************************************** !
! Component and Particle Loops (Component index greater than Ci)                                  !
! *********************************************************************************************** !
DO CJ = CI + 1, COMPONENTS
  ! Unique loop takes only particles whose component indexes are greater than Ci
  DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, SUM( N_COMPONENT(0:CJ) )
    ! Position of particle j
    RJ(1)  = RMC(1,J)
    RJ(2)  = RMC(2,J)
    RJ(3)  = RMC(3,J)
    ! Orientation of particle j
    EJ(1)  = EMC(1,J)
    EJ(2)  = EMC(2,J)
    EJ(3)  = EMC(3,J)
    ! Quaternion of particle j
    QJ(0)  = QMC(0,J)
    QJ(1)  = QMC(1,J)
    QJ(2)  = QMC(2,J)
    QJ(3)  = QMC(3,J)
    ! Vector distance between particles i and j
    RIJ(1) = RJ(1) - RI(1)
    RIJ(2) = RJ(2) - RI(2)
    RIJ(3) = RJ(3) - RI(3)
    ! Minimum Image Convention
    CALL MULTI_MATRIX( BLI, RIJ, S12 )
    S12 = S12 - ANINT(S12)
    CALL MULTI_MATRIX( BL, S12, RIJ )
    ! Magnitude of the vector distance (squared)
    RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
    ! Cutoff distance (squared)
    CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
    CUTOFF_D = CUTOFF_D * CUTOFF_D
    ! Preliminary test (circumscribing spheres)
    IF( RIJSQ <= CUTOFF_D ) THEN
      ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
      IF( GEOM_SELEC(1) ) THEN
        CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CONTACT_D, OVERLAP_HER )
        ! Overlap criterion
        IF( OVERLAP_HER ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      ! Overlap test for spherocylinders (Vega-Lago Method)
      ELSE IF( GEOM_SELEC(2) ) THEN
        CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
        ! Overlap criterion
        IF( OVERLAP_SPC ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      ! Overlap test for cylinders (Lopes et al. Method)
      ELSE IF( GEOM_SELEC(3) ) THEN
        ! Preliminary test (circumscribing spherocylinders)
        OVERLAP_SPC = .FALSE.
        CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
        ! Overlap criterion
        IF( OVERLAP_SPC ) THEN
          RJ(1) = RI(1) + RIJ(1)
          RJ(2) = RI(2) + RIJ(2)
          RJ(3) = RI(3) + RIJ(3)
          CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, DLAMBDAEI, DMUEJ, PARALLEL, OVERLAP_CYL )
          ! Overlap criterion
          IF( OVERLAP_CYL ) THEN
            ! Overlap detected
            OVERLAP = .TRUE.
            RETURN
          END IF
        END IF
      END IF
    END IF
  END DO
END DO

! *********************************************************************************************** !
! Component and Particle Loops (Component index equals Ci)                                        !
! *********************************************************************************************** !
CJ = CI
! First loop takes only particles whose j-indexes are below the i-index of the particles of the component Ci
DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, I - 1
  ! Position of particle j
  RJ(1)  = RMC(1,J)
  RJ(2)  = RMC(2,J)
  RJ(3)  = RMC(3,J)
  ! Orientation of particle j
  EJ(1)  = EMC(1,J)
  EJ(2)  = EMC(2,J)
  EJ(3)  = EMC(3,J)
  ! Quaternion of particle j
  QJ(0)  = QMC(0,J)
  QJ(1)  = QMC(1,J)
  QJ(2)  = QMC(2,J)
  QJ(3)  = QMC(3,J)
  ! Vector distance between particles i and j
  RIJ(1) = RJ(1) - RI(1)
  RIJ(2) = RJ(2) - RI(2)
  RIJ(3) = RJ(3) - RI(3)
  ! Minimum Image Convention
  CALL MULTI_MATRIX( BLI, RIJ, S12 )
  S12 = S12 - ANINT(S12)
  CALL MULTI_MATRIX( BL, S12, RIJ )
  ! Magnitude of the vector distance (squared)
  RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
  ! Cutoff distance (squared)
  CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
  CUTOFF_D = CUTOFF_D * CUTOFF_D
  ! Preliminary test (circumscribing spheres)
  IF( RIJSQ <= CUTOFF_D ) THEN
    ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
    IF( GEOM_SELEC(1) ) THEN
      CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CONTACT_D, OVERLAP_HER )
      ! Overlap criterion
      IF( OVERLAP_HER ) THEN
        ! Overlap detected
        OVERLAP = .TRUE.
        RETURN
      END IF
    ! Overlap test for spherocylinders (Vega-Lago Method)
    ELSE IF( GEOM_SELEC(2) ) THEN
      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
      ! Overlap criterion
      IF( OVERLAP_SPC ) THEN
        ! Overlap detected
        OVERLAP = .TRUE.
        RETURN
      END IF
    ! Overlap test for cylinders (Lopes et al. Method)
    ELSE IF( GEOM_SELEC(3) ) THEN
      ! Preliminary test (circumscribing spherocylinders)
      OVERLAP_SPC = .FALSE.
      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
      ! Overlap criterion
      IF( OVERLAP_SPC ) THEN
        RJ(1) = RI(1) + RIJ(1)
        RJ(2) = RI(2) + RIJ(2)
        RJ(3) = RI(3) + RIJ(3)
        CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, DLAMBDAEI, DMUEJ, PARALLEL, OVERLAP_CYL )
        ! Overlap criterion
        IF( OVERLAP_CYL ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      END IF
    END IF
  END IF
END DO
! Second loop takes only particles whose j-indexes are above the i-index of the particles of the component Ci
DO J = I + 1, SUM( N_COMPONENT(0:CJ) )
  ! Position of particle j
  RJ(1)  = RMC(1,J)
  RJ(2)  = RMC(2,J)
  RJ(3)  = RMC(3,J)
  ! Orientation of particle j
  EJ(1)  = EMC(1,J)
  EJ(2)  = EMC(2,J)
  EJ(3)  = EMC(3,J)
  ! Quaternion of particle j
  QJ(0)  = QMC(0,J)
  QJ(1)  = QMC(1,J)
  QJ(2)  = QMC(2,J)
  QJ(3)  = QMC(3,J)
  ! Vector distance between particles i and j
  RIJ(1) = RJ(1) - RI(1)
  RIJ(2) = RJ(2) - RI(2)
  RIJ(3) = RJ(3) - RI(3)
  ! Minimum Image Convention
  CALL MULTI_MATRIX( BLI, RIJ, S12 )
  S12 = S12 - ANINT(S12)
  CALL MULTI_MATRIX( BL, S12, RIJ )
  ! Magnitude of the vector distance (squared)
  RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
  ! Cutoff distance (squared)
  CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
  CUTOFF_D = CUTOFF_D * CUTOFF_D
  ! Preliminary test (circumscribing spheres)
  IF( RIJSQ <= CUTOFF_D ) THEN
    ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
    IF( GEOM_SELEC(1) ) THEN
      CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CONTACT_D, OVERLAP_HER )
      ! Overlap criterion
      IF( OVERLAP_HER ) THEN
        ! Overlap detected
        OVERLAP = .TRUE.
        RETURN
      END IF
    ! Overlap test for spherocylinders (Vega-Lago Method)
    ELSE IF( GEOM_SELEC(2) ) THEN
      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
      ! Overlap criterion
      IF( OVERLAP_SPC ) THEN
        ! Overlap detected
        OVERLAP = .TRUE.
        RETURN
      END IF
    ! Overlap test for cylinders (Lopes et al. Method)
    ELSE IF( GEOM_SELEC(3) ) THEN
      ! Preliminary test (circumscribing spherocylinders)
      OVERLAP_SPC = .FALSE.
      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, CONTACT_D, PARALLEL, OVERLAP_SPC )
      ! Overlap criterion
      IF( OVERLAP_SPC ) THEN
        RJ(1) = RI(1) + RIJ(1)
        RJ(2) = RI(2) + RIJ(2)
        RJ(3) = RI(3) + RIJ(3)
        CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, DLAMBDAEI, DMUEJ, PARALLEL, OVERLAP_CYL )
        ! Overlap criterion
        IF( OVERLAP_CYL ) THEN
          ! Overlap detected
          OVERLAP = .TRUE.
          RETURN
        END IF
      END IF
    END IF
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE CHECK_OVERLAP

! *********************************************************************************************** !
!    This subroutine calculates the order parameter of a nematic phase via the Q-tensor method    !
! *********************************************************************************************** !
SUBROUTINE ORDER_PARAMETER( S, EP )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: I, C        ! Counter
INTEGER( KIND= REAL64 ) :: N_S         ! Non-spherical particles counter
INTEGER( KIND= REAL64 ) :: ALPHA, BETA ! Unit vector specifiers (î, ĵ, k̂)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                              :: S           ! Nematic order parameter
REAL( KIND= REAL64 )                              :: KRONECKER   ! Kronecker delta
REAL( KIND= REAL64 )                              :: M           ! Trace (matrix)
REAL( KIND= REAL64 )                              :: DET         ! Determinant (matrix)
REAL( KIND= REAL64 )                              :: P           ! Sum of squares of elements of a matrix
REAL( KIND= REAL64 )                              :: HQ, PHI, PQ ! Auxiliar variables
REAL( KIND= REAL64 ), DIMENSION( 3 )              :: EIGENVECTOR ! Eigenvector of order tensor Q
REAL( KIND= REAL64 ), DIMENSION( 3, 3 )           :: QABI        ! Order tensor Q of particle i (3 x 3 Matrix)
REAL( KIND= REAL64 ), DIMENSION( 3, 3 )           :: QAB         ! Order tensor Q of particle i (Average)
REAL( KIND= REAL64 ), DIMENSION( 3, N_PARTICLES ) :: EP          ! Order tensor Q of particle i (Average)

! *********************************************************************************************** !
! Initialization of Q matrix                                                                      !
! *********************************************************************************************** !
QABI(:,:) = 0.D0  

! *********************************************************************************************** !
! Matrix construction                                                                             !
! *********************************************************************************************** !
N_S = 0
DO C = 1, COMPONENTS
  ! Skip if component is spherical
  IF( ( GEOM_SELEC(1) .AND. DABS( ASPECT_RATIO(C) - 1.D0 ) < EPSILON( 1.D0 ) ) .OR. &
  &   ( GEOM_SELEC(2) .AND. DABS( ASPECT_RATIO(C) - 0.D0 ) < EPSILON( 1.D0 ) ) ) THEN
    CYCLE
  END IF
  DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
    DO ALPHA = 1, 3
      DO BETA = 1, 3
        ! Dyadic product (Kronecker delta)
        IF ( ALPHA == BETA ) THEN
          KRONECKER = 1.D0
        ELSE
          KRONECKER = 0.D0
        END IF
        ! Second-order Legendre polynomial (P₂)
        QABI(ALPHA,BETA) = QABI(ALPHA,BETA) + ( 1.5D0 * EP(ALPHA,I) * EP(BETA,I) ) - (0.5D0 * KRONECKER)
      END DO
    END DO
    N_S = N_S + 1
  END DO
END DO

! *********************************************************************************************** !
! Arithmetic mean                                                                                 !
! *********************************************************************************************** !
QAB(:,:) = QABI(:,:) / DBLE( N_S )

! All components are spherical
IF( DABS( SUM( QAB ) - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  S = 0.D0
  RETURN
END IF

! *********************************************************************************************** !
! Cardano's trigonometric solution                                                                !
! *********************************************************************************************** !
! Trace of symmetric Q matrix
M   = QAB(1,1) + QAB(2,2) + QAB(3,3)
! One-third of the trace of the symmetric Q matrix
M   = M / 3.D0
! Determinant of (Q - mI), where I is the identity matrix
DET = ( QAB(1,1) - M ) * ( QAB(2,2) - M ) * ( QAB(3,3) - M ) + ( QAB(1,2) * QAB(2,3) * QAB(3,1) ) + &
&     ( QAB(1,3) * QAB(3,2) * QAB(2,1) ) - ( QAB(1,3) * ( QAB(2,2) - M ) * QAB(3,1) ) - &
&     ( QAB(2,3) * QAB(3,2) * ( QAB(1,1) - M ) ) - ( ( QAB(3,3) - M ) * QAB(2,1) * QAB(1,2) )
! Half of the determinant of (Q - mI)
HQ  = 0.5D0 * DET
! Sum of squares of elements of (Q - mI)
P   = ( QAB(1,1) - M ) * ( QAB(1,1) - M ) + ( QAB(1,2) * QAB(1,2) ) + ( QAB(1,3) * QAB(1,3) ) + &
&     ( QAB(2,2) - M ) * ( QAB(2,2) - M ) + ( QAB(2,1) * QAB(2,1) ) + ( QAB(2,3) * QAB(2,3) ) + &
&     ( QAB(3,3) - M ) * ( QAB(3,3) - M ) + ( QAB(3,1) * QAB(3,1) ) + ( QAB(3,2) * QAB(3,2) )
! One-sixth of the sum of squares of elements of (Q - mI)
P   = P / 6.D0
! Test condition (Discriminant)
PQ  = ( P * P * P ) - ( HQ * HQ )
! Real eigenvalues condition (p³ ≥ hq²)                                     
IF ( PQ >= 0.D0 ) THEN
  PHI = DATAN( DSQRT( PQ ) / HQ ) / 3.D0 ! 0 ≤ ϕ ≤ π
ELSE
  PHI = 0.D0
END IF

! *********************************************************************************************** !
! Eigenvector of Q (characteristic roots from Cardano's Formula)                                  !
! *********************************************************************************************** !
EIGENVECTOR(1) = M + ( 2.D0 * DSQRT( P ) * DCOS( PHI ) )
EIGENVECTOR(2) = M - DSQRT( P ) * ( DCOS( PHI ) + ( DSQRT( 3.D0 ) * DSIN( PHI ) ) )
EIGENVECTOR(3) = M - DSQRT( P ) * ( DCOS( PHI ) - ( DSQRT( 3.D0 ) * DSIN( PHI ) ) )

! *********************************************************************************************** !
! Nematic order parameter                                                                         !
! *********************************************************************************************** !
S = MAXVAL( EIGENVECTOR ) ! Largest eigenvalue

RETURN

END SUBROUTINE ORDER_PARAMETER

! *********************************************************************************************** !
!               This subroutine calculates the inverse of a matrix using cofactors                !
! *********************************************************************************************** !
SUBROUTINE INVERSE_COF( A, B, DET )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: DET  ! Determinant (volume)
REAL( KIND= REAL64 )                 :: RDET ! Determinant (inverse)
REAL( KIND= REAL64 ), DIMENSION( 9 ) :: A    ! Box length
REAL( KIND= REAL64 ), DIMENSION( 9 ) :: B    ! Box length (inverse)

! *********************************************************************************************** !
! Tranpose matrix of the matrix of cofactors of matrix A                                          !
! *********************************************************************************************** !
B(1) = A(5) * A(9) - A(6) * A(8)
B(2) = A(3) * A(8) - A(2) * A(9)
B(3) = A(2) * A(6) - A(3) * A(5)
B(4) = A(6) * A(7) - A(4) * A(9)
B(5) = A(1) * A(9) - A(3) * A(7)
B(6) = A(3) * A(4) - A(1) * A(6)
B(7) = A(4) * A(8) - A(5) * A(7)
B(8) = A(2) * A(7) - A(1) * A(8)
B(9) = A(1) * A(5) - A(2) * A(4)

! *********************************************************************************************** !
! Determinant of matrix A                                                                         !
! *********************************************************************************************** !
DET  = A(1) * B(1) + A(4) * B(2) + A(7) * B(3)
RDET = 0.0D0

IF( DABS( DET ) > 0.0D0 ) THEN
  RDET = 1.0D0 / DET
END IF

! *********************************************************************************************** !
! Inverse of matrix A                                                                             !
! *********************************************************************************************** !
B(:) = RDET * B(:)

RETURN

END SUBROUTINE INVERSE_COF

! *********************************************************************************************** !
!                        This subroutine multiplies a matrix and a vector                         !
! *********************************************************************************************** !
SUBROUTINE MULTI_MATRIX( A, B, C )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 9 ) :: A ! Box length
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: B ! Vector distance (original)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: C ! Vector distance (transformation)

! *********************************************************************************************** !
! Multiplication of a matrix and a vector                                                         !
! *********************************************************************************************** !
C(1) = A(1) * B(1) + A(4) * B(2) + A(7) * B(3)
C(2) = A(2) * B(1) + A(5) * B(2) + A(8) * B(3)
C(3) = A(3) * B(1) + A(6) * B(2) + A(9) * B(3)

RETURN

END SUBROUTINE MULTI_MATRIX

! *********************************************************************************************** !
!             This subroutine generates a progress bar for the initial configuration.             !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR( J, TOTAL )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: J     ! Counter
INTEGER( KIND= REAL64 ) :: TOTAL ! Total/Maximum number of cycles

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 16 ) :: BAR ! Progress bar

! *********************************************************************************************** !
! Progress bar (FORMAT)                                                                           !
! *********************************************************************************************** !
BAR = "Progress: ?????%"

! *********************************************************************************************** !
! Progress bar (replace character positions)                                                      !
! *********************************************************************************************** !
WRITE( UNIT= BAR(11:15), FMT= "(F5.1)" ) ( DBLE(J) / DBLE(TOTAL) ) * 100.D0

! *********************************************************************************************** !
! Print progress bar                                                                              !
! *********************************************************************************************** !
WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A16)" , ADVANCE= "NO" ) CHAR(13), BAR

! *********************************************************************************************** !
! Flush standard output unit                                                                      !
! *********************************************************************************************** !
FLUSH( UNIT= OUTPUT_UNIT )

RETURN

END SUBROUTINE PROGRESS_BAR

! *********************************************************************************************** !
!            This subroutine generates a progress bar for the Monte Carlo simulation.             !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR_MC( J, TOTAL, ENS )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= REAL64 ) :: J     ! Counter
INTEGER( KIND= REAL64 ) :: TOTAL ! Total/Maximum number of cycles

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 21 ) :: BAR  ! Progress bar
CHARACTER( LEN= 03 ) :: ENS  ! Ensemble type

! *********************************************************************************************** !
! Progress bar (FORMAT)                                                                           !
! *********************************************************************************************** !
IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 10.D0 ) THEN
  BAR = "Progress("//TRIM( ENS )//"): ???%"
ELSE IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 100.D0 ) THEN
  BAR = "Progress("//TRIM( ENS )//"): ????%"
ELSE
  BAR = "Progress("//TRIM( ENS )//"): ?????%"
END IF

! *********************************************************************************************** !
! Progress bar (replace character positions)                                                      !
! *********************************************************************************************** !
IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 10.D0 ) THEN
  WRITE( UNIT= BAR(16:18), FMT= "(F3.1)" ) ( DBLE(J) / DBLE(TOTAL) ) * 100.D0
ELSE IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 100.D0 ) THEN
  WRITE( UNIT= BAR(16:19), FMT= "(F4.1)" ) ( DBLE(J) / DBLE(TOTAL) ) * 100.D0
ELSE
  WRITE( UNIT= BAR(16:20), FMT= "(F5.1)" ) ( DBLE(J) / DBLE(TOTAL) ) * 100.D0
END IF

! *********************************************************************************************** !
! Print progress bar                                                                              !
! *********************************************************************************************** !
IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 10.D0 ) THEN
  WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A19)" , ADVANCE= "NO" ) CHAR(13), BAR
ELSE IF( ( DBLE(J) / DBLE(TOTAL) ) * 100.D0 < 100.D0 ) THEN
  WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A20)" , ADVANCE= "NO" ) CHAR(13), BAR
ELSE
  WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A21)" , ADVANCE= "NO" ) CHAR(13), BAR
END IF

! *********************************************************************************************** !
! Flush standard output unit                                                                      !
! *********************************************************************************************** !
FLUSH( UNIT= OUTPUT_UNIT )

RETURN

END SUBROUTINE PROGRESS_BAR_MC

! *********************************************************************************************** !
!                 This function converts any string into uppercase (from A to Z)                  !
! *********************************************************************************************** !
SUBROUTINE TO_UPPER( STRIN, LENSTR, STROUT )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER :: LENSTR ! String length
INTEGER :: J      ! ASCII character code
INTEGER :: I      ! Counter

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= LENSTR ) :: STRIN
CHARACTER( LEN= LENSTR ) :: STROUT

! Character positions
DO I = 1, LENSTR
  ! ASCII character code
  J = IACHAR( STRIN(I:I) )
  ! Convert to uppercase (letters only)
  IF( J >= IACHAR( "a" ) .AND. J <= IACHAR( "z" ) ) THEN
    STROUT(I:I) = ACHAR(IACHAR(STRIN(I:I))-32)
  ! Do not convert symbols or numbers (special characters included)
  ELSE
    STROUT(I:I) = STRIN(I:I)
  END IF
END DO

RETURN

END SUBROUTINE TO_UPPER