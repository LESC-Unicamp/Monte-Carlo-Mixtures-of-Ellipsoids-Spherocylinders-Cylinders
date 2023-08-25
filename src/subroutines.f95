! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains all subroutines used in the main program.                   !
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
!            i. e., its value is changed every time the ranf(  ) subroutine is called.            !
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
INTEGER( KIND= INT64 ), PARAMETER :: L = 1029
INTEGER( KIND= INT64 ), PARAMETER :: C = 221591
INTEGER( KIND= INT64 ), PARAMETER :: G = 1048576

! *********************************************************************************************** !
! Finite modulus arithmetic                                                                       !
! *********************************************************************************************** !
SEED     = MOD( ( (SEED * L) + C ), G ) ! Number generator seed
RANDOM_N = DBLE( SEED ) / DBLE( G )     ! Pseudorandom number

RETURN

END SUBROUTINE RANF

! *********************************************************************************************** !
!     This subroutine takes a body-fixed orientation/position and a rotation quaternion and       !
!            generates a space-fixed orientation/position using a 3D-rotation matrix.             !
!        See Allen and Tildesley, 2nd Edition (2017), pages 106-111 for more information.         !
! *********************************************************************************************** !
SUBROUTINE ACTIVE_TRANSFORMATION( EFIXED, QP, EP )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J ! Counters

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EFIXED ! Body-fixed orientation/position
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EP     ! Space-fixed orientation/position
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
! Active tranformation (dot product of body-fixed vector and transpose of rotation matrix)        !
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
!            This subroutine creates a composed rotation quaternion by multiplying the            !
!                         reference quaternion with a random quaternion                           !
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
INTEGER( KIND= INT64 ) :: I, J      ! Counters
INTEGER( KIND= INT64 ) :: C, CI, CJ ! Component indexes

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                          :: RIJSQ     ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                          :: CONTACT_D ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( KIND= REAL64 )                          :: CUTOFF_D  ! Cutoff distance
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: DLAMBDAEI ! Auxiliary vector (cylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: DMUEJ     ! Auxiliary vector (cylinder overlap algorithm)
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
OVERLAP_SPC = .FALSE.
OVERLAP_CYL = .FALSE.
PARALLEL    = .FALSE.

! Diameter of the circumscribing sphere
IF( GEOM_SELEC(1) ) THEN ! Ellipsoids of revolution
  DO C = 1, COMPONENTS
    ! Oblate ellipsoids of revolution or spheres
    IF( ASPECT_RATIO(C) > 0.D0 .AND. ASPECT_RATIO(C) <= 1.D0 ) THEN
      CUTOFF(C) = DIAMETER(C)
    ! Prolate ellipsoids of revolution
    ELSE IF( ASPECT_RATIO(C) > 1.D0 ) THEN
      CUTOFF(C) = LENGTH(C)
    END IF
  END DO
ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN ! Spherocylinders and cylinders
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
    RJ(1) = RMC(1,J)
    RJ(2) = RMC(2,J)
    RJ(3) = RMC(3,J)
    ! Orientation of particle j
    EJ(1) = EMC(1,J)
    EJ(2) = EMC(2,J)
    EJ(3) = EMC(3,J)
    ! Quaternion of particle j
    QJ(0) = QMC(0,J)
    QJ(1) = QMC(1,J)
    QJ(2) = QMC(2,J)
    QJ(3) = QMC(3,J)
    ! Vector distance between particles i and j
    RIJ(1) = RJ(1) - RI(1)
    RIJ(2) = RJ(2) - RI(2)
    RIJ(3) = RJ(3) - RI(3)
    ! Minimum image convention
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
!   This subroutine computes the total potential energy for the initial molecular configuration   !
! *********************************************************************************************** !
SUBROUTINE COMPUTE_TOTAL_ENERGY(  )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J   ! Counters
INTEGER( KIND= INT64 ) :: CI, CJ ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                        :: RIJSQ  ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: RI, RJ ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: RIJ    ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: S12    ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( N_LAMBDA ) :: VIJ    ! Pair potential energy

! Initialization
V(:) = 0.D0

! Anisomorphic molecules (unlike components)
DO CI = 1, COMPONENTS - 1
  DO CJ = CI + 1, COMPONENTS
    ! First loop represents all particles with indexes i of component Ci
    DO I = SUM( N_COMPONENT(0:(CI-1)) ) + 1, SUM( N_COMPONENT(0:CI) )
      ! Second loop represents all particles with indexes j of component Cj
      DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, SUM( N_COMPONENT(0:CJ) )
        ! Position of particle i
        RI(1)  = R(1,I)
        RI(2)  = R(2,I)
        RI(3)  = R(3,I)
        ! Position of particle j
        RJ(1)  = R(1,J)
        RJ(2)  = R(2,J)
        RJ(3)  = R(3,J)
        ! Vector distance between particles i and j
        RIJ(1) = RJ(1) - RI(1)
        RIJ(2) = RJ(2) - RI(2)
        RIJ(3) = RJ(3) - RI(3)
        ! Minimum Image Convention
        CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
        S12 = S12 - ANINT(S12)
        CALL MULTI_MATRIX( BOX_LENGTH, S12, RIJ )
        ! Magnitude of the vector distance (squared)
        RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
        ! Compute pair potential
        IF( POTENTIAL_SELEC(2) ) THEN
          CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
        END IF
        ! Increment total potential energy
        V(:) = V(:) + VIJ(:)
      END DO
    END DO
  END DO
END DO
! Isomorphic molecules (like components)
DO CI = 1, COMPONENTS
  CJ = CI
  ! First loop represents a particle with an index i of component Ci
  DO I = SUM( N_COMPONENT(0:(CI-1)) ) + 1, SUM( N_COMPONENT(0:CI) ) - 1
    ! Second loop represents all other particles with indexes j > i of component Cj = Ci
    DO J = I + 1, SUM( N_COMPONENT(0:CI) )
      ! Position of particle i
      RI(1)  = R(1,I)
      RI(2)  = R(2,I)
      RI(3)  = R(3,I)
      ! Position of particle j
      RJ(1)  = R(1,J)
      RJ(2)  = R(2,J)
      RJ(3)  = R(3,J)
      ! Vector distance between particles i and j
      RIJ(1) = RJ(1) - RI(1)
      RIJ(2) = RJ(2) - RI(2)
      RIJ(3) = RJ(3) - RI(3)
      ! Minimum Image Convention
      CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
      S12 = S12 - ANINT(S12)
      CALL MULTI_MATRIX( BOX_LENGTH, S12, RIJ )
      ! Magnitude of the vector distance (squared)
      RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
      ! Compute pair potential
      IF( POTENTIAL_SELEC(2) ) THEN
        CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
      END IF
      ! Increment total potential energy
      V(:) = V(:) + VIJ(:)
    END DO
  END DO
END DO

RETURN

END SUBROUTINE

! *********************************************************************************************** !
!              This subroutine computes the potential energy of a random particle i               !
! *********************************************************************************************** !
SUBROUTINE COMPUTE_PARTICLE_ENERGY( CI, I, RI, VI, BL, BLI )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J   ! Counters
INTEGER( KIND= INT64 ) :: CI, CJ ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                        :: RIJSQ  ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: RI, RJ ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: RIJ    ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )        :: S12    ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( 9 )        :: BL     ! Box length
REAL( KIND= REAL64 ), DIMENSION( 9 )        :: BLI    ! Box length (inverse)
REAL( KIND= REAL64 ), DIMENSION( N_LAMBDA ) :: VI     ! Potential energy of particle i
REAL( KIND= REAL64 ), DIMENSION( N_LAMBDA ) :: VIJ    ! Pair potential energy

! Initialization
VI(:) = 0.D0

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
    ! Compute pair potential
    IF( POTENTIAL_SELEC(2) ) THEN
      CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
    END IF
    ! Increment total potential energy
    VI(:) = VI(:) + VIJ(:)
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
    ! Compute pair potential
    IF( POTENTIAL_SELEC(2) ) THEN
      CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
    END IF
    ! Increment total potential energy
    VI(:) = VI(:) + VIJ(:)
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
  ! Compute pair potential
  IF( POTENTIAL_SELEC(2) ) THEN
    CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
  END IF
  ! Increment total potential energy
  VI(:) = VI(:) + VIJ(:)
END DO
! Second loop takes only particles whose j-indexes are above the i-index of the particles of the component Ci
DO J = I + 1, SUM( N_COMPONENT(0:CJ) )
  ! Position of particle j
  RJ(1)  = RMC(1,J)
  RJ(2)  = RMC(2,J)
  RJ(3)  = RMC(3,J)
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
  ! Compute pair potential
  IF( POTENTIAL_SELEC(2) ) THEN
    CALL SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )
  END IF
  ! Increment total potential energy
  VI(:) = VI(:) + VIJ(:)
END DO

RETURN

END SUBROUTINE COMPUTE_PARTICLE_ENERGY

! *********************************************************************************************** !
!              This subroutine computes the pair potential between particles i and j              !
!           It applies a discrete square-well potential to compute the pair potential.            !
! *********************************************************************************************** !
SUBROUTINE SW_POTENTIAL( RIJSQ, CI, CJ, VIJ )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: CI, CJ ! Counters (component)
INTEGER( KIND= INT64 ) :: C_LAMB ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                        :: RIJSQ     ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                        :: SWRANGESQ ! Effective range of attraction (squared)
REAL( KIND= REAL64 ), DIMENSION( N_LAMBDA ) :: VIJ       ! Pair potential energy

! Compute pair potential for every attractive range
DO C_LAMB = 1, N_LAMBDA
  ! Effective range of attraction (squared)
  SWRANGESQ = 0.5D0 * ( SWRANGE(CI,C_LAMB) + SWRANGE(CJ,C_LAMB) )
  SWRANGESQ = SWRANGESQ * SWRANGESQ
  ! Overlap in the region of attraction
  IF( RIJSQ <= SWRANGESQ ) THEN
    ! Pair potential (reduced units)
    VIJ(C_LAMB) = -1.D0
  ELSE
    VIJ(C_LAMB) = 0.D0 ! No overlap
  END IF
END DO

RETURN

END SUBROUTINE SW_POTENTIAL

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
INTEGER( KIND= INT64 ) :: I, C        ! Counter
INTEGER( KIND= INT64 ) :: N_S         ! Non-spherical particles counter
INTEGER( KIND= INT64 ) :: ALPHA, BETA ! Unit vector specifiers (î, ĵ, k̂)

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
!     This subroutine uses a block-averaging method to calculate the first- and second-order      !
!         TPT coefficients, the perturbed Helmholtz free energy, and their uncertainties          !
! *********************************************************************************************** !
!                        Original developer: Luis Fernando Mercier Franco                         !
!                     University of Campinas, School of Chemical Engineering                      !
! *********************************************************************************************** !
!         See Allen and Tildesley, 2nd Edition (2017), pages 281-287 for more information         !
! *********************************************************************************************** !
SUBROUTINE BLOCK_AVERAGING( FLAG, A1, A2, APERT, DPA1, DPA2, DPAPERT, EXEC_TIME )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 )            :: EQUIL_SAVE               ! Equilibration data points in the potential file
INTEGER( KIND= INT64 )            :: N_RUN                    ! Production data points in the potential file
INTEGER( KIND= INT64 )            :: N_BLOCK                  ! Counter for the number of blocks
INTEGER( KIND= INT64 )            :: N_DATA                   ! Number of data points in a block
INTEGER( KIND= INT64 )            :: I, J, K, COUNT_AUX, STEP ! Auxiliary counters
INTEGER( KIND= INT64 ), PARAMETER :: MAX_DATA = 1.D4          ! Maximum number of block data for linear regression

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ) :: PERTURBED_MOMENT1              ! First moment (average), <exp(-U*/T*)>
REAL( KIND= REAL64 ) :: PERTURBED_MOMENT2              ! Second moment, <exp(-2U*/T*)>
REAL( KIND= REAL64 ) :: MOMENT1                        ! First moment (average), <U*>
REAL( KIND= REAL64 ) :: MOMENT2                        ! Second moment, <U*²>
REAL( KIND= REAL64 ) :: MOMENT3                        ! Third moment, <U*³>
REAL( KIND= REAL64 ) :: MOMENT4                        ! Fourth moment, <U*4>
REAL( KIND= REAL64 ) :: PERTURBED_VAR1_TOT             ! Variance of the 1st moment, <exp(-2U*/T*)> - <exp(-U*/T*)>²
REAL( KIND= REAL64 ) :: VAR1_TOT                       ! Variance of the 1st moment, <U*²> - <U*>²
REAL( KIND= REAL64 ) :: VAR2_TOT                       ! Variance of the 2nd moment, <U*4> - <U*²>²
REAL( KIND= REAL64 ) :: PERTURBED_AVG_BLK, AVG_BLK     ! Block average
REAL( KIND= REAL64 ) :: PERTURBED_VAR_BLOCK, VAR_BLOCK ! Block variance
REAL( KIND= REAL64 ) :: PERTURBED_VAR_SUM, VAR_SUM     ! Auxiliary summation variable for the block variance
REAL( KIND= REAL64 ) :: A, B                           ! Parameters of linear regression
REAL( KIND= REAL64 ) :: R2, PERTURBED_R2               ! Coefficient of determination
REAL( KIND= REAL64 ) :: PERTURBED_S_RUN, S_RUN         ! Statistical inefficiency
REAL( KIND= REAL64 ) :: PERTURBED_SIGSQ1               ! Variance of <exp(-U*/T*)>
REAL( KIND= REAL64 ) :: SIGSQ1                         ! Variance of <U*>
REAL( KIND= REAL64 ) :: SIGSQ2                         ! Variance of <U*²>
REAL( KIND= REAL64 ) :: COV                            ! Covariance between U* and U*²
REAL( KIND= REAL64 ) :: CORR                           ! Correlation between U* and U*²
REAL( KIND= REAL64 ) :: EXPARGUMENT                    ! Boltzmann factor argument
REAL( KIND= REAL64 ) :: APERT                          ! Perturbed Helmholtz free energy
REAL( KIND= REAL64 ) :: A1, A2                         ! TPT coefficients
REAL( KIND= REAL64 ) :: DPAPERT                        ! Perturbed Helmholtz free energy (standard deviation)
REAL( KIND= REAL64 ) :: DPA1, DPA2                     ! TPT coefficients (standard deviation)
REAL( KIND= REAL64 ) :: INITIAL_TIME                   ! Initial time
REAL( KIND= REAL64 ) :: FINAL_TIME                     ! Final time
REAL( KIND= REAL64 ) :: EXEC_TIME                      ! Execution time

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE :: STEP_BLOCKS      ! Inverse of the number of data points in each block
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE :: PERTURBED_SB, SB ! Statistical inefficiency per number of blocks

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE :: POT ! Potential energy (production only)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FLAG ! Generic true/false flag

! Number of data points
N_RUN      = INT( ( DBLE( MAX_CYCLES ) - DBLE( N_EQUIL ) ) / DBLE( N_SAVE ) ) ! Production points
EQUIL_SAVE = INT( DBLE( N_EQUIL ) / DBLE( N_SAVE ) )                          ! Equilibration points

! Return to main program if number of block data exceeds max_data parameter
IF ( MAX_BLOCKS - MIN_BLOCKS >= MAX_DATA ) THEN
  FLAG = .TRUE.
  RETURN
END IF

! Initialization of variables
PERTURBED_MOMENT1 = 0.D0
PERTURBED_MOMENT2 = 0.D0
MOMENT1 = 0.D0
MOMENT2 = 0.D0
MOMENT3 = 0.D0
MOMENT4 = 0.D0

! Allocation
ALLOCATE( POT(N_RUN) )
ALLOCATE( STEP_BLOCKS(MAX_DATA), PERTURBED_SB(MAX_DATA), SB(MAX_DATA) )

! Start subroutine timer
CALL CPU_TIME( INITIAL_TIME )

! *********************************************************************************************** !
! Initialization of the Largest Exponential Argument                                              !
! *********************************************************************************************** !
!  This corresponds to a method to calculate the natural logarithm of large exponential           !
!  iterations, avoiding mathematical inaccuracies. This is applied to determine the Helmholtz     !
!  free energy of the perturbed system by taking the negative natural logarithm of the canonical  !
!  ensemble average over the configurations visited by the reference fluid.                       !
! *********************************************************************************************** !
!  See Supplementary Material of Abreu and Escobedo (2006) for more information.                  !
! *********************************************************************************************** !
MAXARG = 0.D0

! Largest Exponential Argument and Potential Energy
OPEN( UNIT= 150, FILE= "Potential/"//TRIM( DESCRIPTOR_DATE )//"/Lambda_"//TRIM( DESCRIPTOR_LAMB )//"/"//TRIM( DESCRIPTOR_HOUR )// &
&                      "_thermo_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_"//TRIM( DESCRIPTOR_FILE3 )// &
&                      ".dat" )
! Skip header
READ ( 150, * )
! Skip equilibration data points (if necessary)
IF( .NOT. POTENTIAL_CHECK ) THEN
  DO I = 1, EQUIL_SAVE
    READ ( 150, * )
  END DO
END IF
! Read production data points
DO I = 1, N_RUN
  READ ( 150, * ) STEP, POT(I)
  BETA_ENERGY = -POT(I) / RED_TEMP
  IF( BETA_ENERGY > MAXARG ) THEN
    MAXARG = BETA_ENERGY
  END IF
END DO
CLOSE( 150 )

! Iterative process
DO I = 1, N_RUN
  EXPARGUMENT = ( -POT(I) / RED_TEMP ) - MAXARG
  PERTURBED_MOMENT1 = PERTURBED_MOMENT1 + DEXP( EXPARGUMENT )
  PERTURBED_MOMENT2 = PERTURBED_MOMENT2 + DEXP( 2.D0 * EXPARGUMENT )
  MOMENT1 = MOMENT1 + ( POT(I) )
  MOMENT2 = MOMENT2 + ( POT(I) * POT(I) )
  MOMENT3 = MOMENT3 + ( POT(I) * POT(I) * POT(I) )
  MOMENT4 = MOMENT4 + ( POT(I) * POT(I) * POT(I) * POT(I) )
END DO

! First moments (average)
PERTURBED_MOMENT1 = PERTURBED_MOMENT1 / DBLE( N_RUN )
MOMENT1           = MOMENT1 / DBLE( N_RUN )

! Second moments
PERTURBED_MOMENT2 = PERTURBED_MOMENT2 / DBLE( N_RUN )
MOMENT2           = MOMENT2 / DBLE( N_RUN )

! Third moment
MOMENT3 = MOMENT3 / DBLE( N_RUN )

! Fourth moment
MOMENT4 = MOMENT4 / DBLE( N_RUN )

! Variances of <X>
PERTURBED_VAR1_TOT = PERTURBED_MOMENT2 - ( PERTURBED_MOMENT1 * PERTURBED_MOMENT1 )
VAR1_TOT           = MOMENT2 - ( MOMENT1 * MOMENT1 )

! Variance of <U*²>
VAR2_TOT = MOMENT4 - ( MOMENT2 * MOMENT2 )

! Covariance
COV = MOMENT3 - ( MOMENT1 * MOMENT2 )

! Block Average
DO N_BLOCK = MIN_BLOCKS, MAX_BLOCKS

  ! Steps in each block for a certain number of blocks
  N_DATA = INT( DBLE( N_RUN ) / DBLE( N_BLOCK ) )

  ! Initialization
  PERTURBED_VAR_SUM = 0.D0
  VAR_SUM = 0.D0
  J = 1

  DO K = 1, N_BLOCK
    PERTURBED_AVG_BLK = 0.D0 ! Reset
    AVG_BLK = 0.D0 ! Reset
    DO I = 1, N_DATA
      EXPARGUMENT = ( -POT(J) / RED_TEMP ) - MAXARG
      PERTURBED_AVG_BLK = PERTURBED_AVG_BLK + DEXP( EXPARGUMENT )
      AVG_BLK = AVG_BLK + POT(J)
      ! Increment counter
      J = J + 1
    END DO
    PERTURBED_AVG_BLK = PERTURBED_AVG_BLK / DBLE( N_DATA )
    PERTURBED_VAR_SUM = PERTURBED_VAR_SUM + ( PERTURBED_AVG_BLK - PERTURBED_MOMENT1 ) ** 2.D0
    AVG_BLK = AVG_BLK / DBLE( N_DATA )
    VAR_SUM = VAR_SUM + ( AVG_BLK - MOMENT1 ) ** 2.D0
  END DO

  ! Variance of the block averages
  PERTURBED_VAR_BLOCK = PERTURBED_VAR_SUM / DBLE( N_BLOCK - 1 )
  VAR_BLOCK = VAR_SUM / DBLE( N_BLOCK - 1 )

  ! Auxiliary counter
  COUNT_AUX = N_BLOCK - MIN_BLOCKS + 1

  ! Inverse of block size
  STEP_BLOCKS(COUNT_AUX) = 1.D0 / DBLE( N_DATA )

  ! Statistical inefficiency
  PERTURBED_SB(COUNT_AUX) = DBLE( N_DATA ) * ( PERTURBED_VAR_BLOCK / PERTURBED_VAR1_TOT )
  SB(COUNT_AUX) = DBLE( N_DATA ) * ( VAR_BLOCK / VAR1_TOT )
END DO

! Deallocation
DEALLOCATE( POT )

! Linear regression (perturbed Helmholtz free energy)
CALL LINEAR_FIT( STEP_BLOCKS, PERTURBED_SB, A, B, PERTURBED_R2, COUNT_AUX )

! Statistical inefficiency (perturbed Helmholtz free energy)
PERTURBED_S_RUN = A

! Linear regression (TPT coefficients)
CALL LINEAR_FIT( STEP_BLOCKS, SB, A, B, R2, COUNT_AUX )

! Statistical inefficiency (TPT coefficients)
S_RUN = A

! Deallocation
DEALLOCATE( STEP_BLOCKS, PERTURBED_SB, SB )

! Variances of <X>
PERTURBED_SIGSQ1 = PERTURBED_VAR1_TOT * ( PERTURBED_S_RUN / DBLE( N_RUN ) )
SIGSQ1 = VAR1_TOT * ( S_RUN / DBLE( N_RUN ) )

! Variance of <U*²>
SIGSQ2 = VAR2_TOT * ( S_RUN / DBLE( N_RUN ) )

! Covariance between U* and U*²
COV  = COV * ( S_RUN / DBLE( N_RUN ) )

! Correlation between U* and U*²
CORR = COV / ( DSQRT( SIGSQ1 ) ) / ( DSQRT( SIGSQ2 ) )

! TPT Coefficients
A1 = MOMENT1 / N_PARTICLES
A2 = - 0.5D0 * ( MOMENT2 - ( MOMENT1 * MOMENT1 ) ) / N_PARTICLES

! TPT Coefficients (Uncertainty propagation)
DPA1 = DSQRT( SIGSQ1 ) / N_PARTICLES
DPA2 = 0.5D0 * ( DSQRT( ( 4.D0 * SIGSQ1 * MOMENT1 * MOMENT1 ) + SIGSQ2 - ( 4.D0 * MOMENT1 * COV ) ) ) / N_PARTICLES

! Perturbed Helmholtz free energy
APERT = - ( 1.D0 / N_PARTICLES ) * ( DLOG( PERTURBED_MOMENT1 ) + MAXARG )

! Perturbed Helmholtz free energy (Uncertainty propagation)
DPAPERT = ( 1.D0 / N_PARTICLES ) * ( DSQRT( ( PERTURBED_SIGSQ1 ) / ( PERTURBED_MOMENT1 * PERTURBED_MOMENT1 ) ) )

! Finish subroutine timer
CALL CPU_TIME( FINAL_TIME )

! Execution time
EXEC_TIME = FINAL_TIME - INITIAL_TIME

RETURN

END SUBROUTINE BLOCK_AVERAGING

! *********************************************************************************************** !
!                    This subroutine makes a linear regression of x and y data                    !
! *********************************************************************************************** !
SUBROUTINE LINEAR_FIT( X, Y, A, B, R2, N )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, N

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ) :: SX                ! Sum of X
REAL( KIND= REAL64 ) :: SY                ! Sum of Y
REAL( KIND= REAL64 ) :: SXX               ! Sum of X²
REAL( KIND= REAL64 ) :: SXY               ! Sum of XY
REAL( KIND= REAL64 ) :: SYY               ! Sum of Y²
REAL( KIND= REAL64 ) :: A                 ! Y-intercept
REAL( KIND= REAL64 ) :: B                 ! Slope
REAL( KIND= REAL64 ) :: R2                ! Coefficient of determination
REAL( KIND= REAL64 ), DIMENSION( N ) :: X ! Independent variable
REAL( KIND= REAL64 ), DIMENSION( N ) :: Y ! Dependent variable

! Initialization
SX  = 0.D0
SY  = 0.D0
SXX = 0.D0
SXY = 0.D0
SYY = 0.D0

! Iteration
DO I = 1, N
  SX  = SX + X(I)
  SY  = SY + Y(I)
  SXX = SXX + ( X(I) * X(I) )
  SXY = SXY + ( X(I) * Y(I) )
  SYY = SYY + ( Y(I) * Y(I) )
END DO

! Slope
B  = ( ( DBLE( N ) * SXY ) - ( SX * SY ) ) / ( ( DBLE( N ) * SXX ) - ( SX * SX ) )
! Y-intercept
A  = ( SY - ( B * SX ) ) / DBLE( N )
! Coefficient of correlation
R2 = ( ( DBLE( N ) * SXY ) - ( SX * SY ) ) / DSQRT( ( ( DBLE( N ) * SXX ) - ( SX * SX ) ) * ( ( DBLE( N ) * SYY ) - ( SY * SY ) ) )
! Coefficient of determination
R2 = ( R2 * R2 )

RETURN

END SUBROUTINE LINEAR_FIT

! *********************************************************************************************** !
!           This subroutine calculates the lattice reduction and reshapes the box size            !
! *********************************************************************************************** !
SUBROUTINE LATTICE_REDUCTION( BOXL, DISTORTION, LATTICER )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ), PARAMETER :: DIM = 3 ! Lattice dimension

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                     :: DISTORTION             ! Distortion parameter (surface-to-volume ratio)
REAL( KIND= REAL64 )                     :: MODV1, MODV2, MODV3    ! Magnitude of the box vectors
REAL( KIND= REAL64 )                     :: MODCROSS1              ! Magnitude of the cross product between vectors vx and vy of the simulation box (area of plane xy)
REAL( KIND= REAL64 )                     :: MODCROSS2              ! Magnitude of the cross product between vectors vx and vz of the simulation box (area of plane xz)
REAL( KIND= REAL64 )                     :: MODCROSS3              ! Magnitude of the cross product between vectors vy and vz of the simulation box (area of plane yz)
REAL( KIND= REAL64 )                     :: DOTVC                  ! Dot product of vector vx and the cross product of vectors vy and vz (volume)
REAL( KIND= REAL64 )                     :: SURFAREA               ! Surface area
REAL( KIND= REAL64 ), DIMENSION( 3 )     :: V1, V2, V3             ! Vectors vx, vy, and vz of the simulation box
REAL( KIND= REAL64 ), DIMENSION( 3 )     :: CROSS1, CROSS2, CROSS3 ! Cross product of box vectors
REAL( KIND= REAL64 ), DIMENSION( 9 )     :: BOXL                   ! Box vector matrix

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: LATTICER ! Checks whether a lattice reduction is necessary based on the distortion parameter value

! Box vectors
V1(:) = BOXL(1:3) ! vx
V2(:) = BOXL(4:6) ! vy
V3(:) = BOXL(7:9) ! vz

! Vector moduli
MODV1 = ( V1(1) * V1(1) ) + ( V1(2) * V1(2) ) + ( V1(3) * V1(3) )
MODV1 = DSQRT( MODV1 )
MODV2 = ( V2(1) * V2(1) ) + ( V2(2) * V2(2) ) + ( V2(3) * V2(3) )
MODV2 = DSQRT( MODV2 )
MODV3 = ( V3(1) * V3(1) ) + ( V3(2) * V3(2) ) + ( V3(3) * V3(3) )
MODV3 = DSQRT( MODV3 )

! Cross product of vx and vy
CROSS1(1) = ( V1(2) * V2(3) ) - ( V1(3) * V2(2) )
CROSS1(2) = ( V1(3) * V2(1) ) - ( V1(1) * V2(3) )
CROSS1(3) = ( V1(1) * V2(2) ) - ( V1(2) * V2(1) )

! Cross product of vx and vz
CROSS2(1) = ( V1(2) * V3(3) ) - ( V1(3) * V3(2) )
CROSS2(2) = ( V1(3) * V3(1) ) - ( V1(1) * V3(3) )
CROSS2(3) = ( V1(1) * V3(2) ) - ( V1(2) * V3(1) )

! Cross product of vy and vz
CROSS3(1) = ( V2(2) * V3(3) ) - ( V2(3) * V3(2) )
CROSS3(2) = ( V2(3) * V3(1) ) - ( V2(1) * V3(3) )
CROSS3(3) = ( V2(1) * V3(2) ) - ( V2(2) * V3(1) )

! Dot product of vx and the cross product of vy and vz
DOTVC = ( V1(1) * CROSS3(1) ) + ( V1(2) * CROSS3(2) ) + ( V1(3) * CROSS3(3) )

! Vector moduli
MODCROSS1 = ( CROSS1(1) * CROSS1(1) ) + ( CROSS1(2) * CROSS1(2) ) + ( CROSS1(3) * CROSS1(3) )
MODCROSS1 = DSQRT( MODCROSS1 )
MODCROSS2 = ( CROSS2(1) * CROSS2(1) ) + ( CROSS2(2) * CROSS2(2) ) + ( CROSS2(3) * CROSS2(3) )
MODCROSS2 = DSQRT( MODCROSS2 )
MODCROSS3 = ( CROSS3(1) * CROSS3(1) ) + ( CROSS3(2) * CROSS3(2) ) + ( CROSS3(3) * CROSS3(3) )
MODCROSS3 = DSQRT( MODCROSS3 )

! Lattice distortion
DISTORTION = (MODV1 + MODV2 + MODV3) * (MODCROSS1 + MODCROSS2 + MODCROSS3) / DOTVC
DISTORTION = DISTORTION / 9.D0

! Surface area
SURFAREA = 2.D0 * MODCROSS1 + 2.D0 * MODCROSS2 + 2.D0 * MODCROSS3

! Verification
IF( DISTORTION > BOX_DIST ) THEN
  LATTICER = .TRUE.
  ! Lattice reduction method
  IF( LRED_SELEC(1) ) THEN
    ! Gottwald method
    CALL GOTTWALD( V1, V2, V3, SURFAREA )
    ! Update box vectors
    BOXL(1:3) = V1(:)
    BOXL(4:6) = V2(:)
    BOXL(7:9) = V3(:)
  ELSE IF( LRED_SELEC(2) ) THEN
    ! Lenstra-Lenstra-Lovász method
    CALL LLL( BOXL, DIM )
  END IF
ELSE
  LATTICER = .FALSE.
  RETURN
END IF

RETURN

END SUBROUTINE LATTICE_REDUCTION

! *********************************************************************************************** !
!          This subroutine applies the Gottwald method to orthogonalize the box vectors           !
! *********************************************************************************************** !
SUBROUTINE GOTTWALD( V1, V2, V3, MINSAREA )

USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: COMBINATION ! Lattice combination
INTEGER( KIND= INT64 ) :: MINSAREALOC ! Array location of minimum surface area

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                     :: MODCROSS1              ! Magnitude of the cross product between vectors vx and vy of the simulation box (area of plane xy)
REAL( KIND= REAL64 )                     :: MODCROSS2              ! Magnitude of the cross product between vectors vx and vz of the simulation box (area of plane xz)
REAL( KIND= REAL64 )                     :: MODCROSS3              ! Magnitude of the cross product between vectors vy and vz of the simulation box (area of plane yz)
REAL( KIND= REAL64 )                     :: MINSAREA               ! Minimum surface area
REAL( KIND= REAL64 ), DIMENSION( 3 )     :: V1, V2, V3             ! Vectors vx, vy, and vz of the simulation box
REAL( KIND= REAL64 ), DIMENSION( 3 )     :: NEWV1, NEWV2, NEWV3    ! New vectors vx, vy, and vz of the simulation box
REAL( KIND= REAL64 ), DIMENSION( 3 )     :: CROSS1, CROSS2, CROSS3 ! Cross product of box vectors
REAL( KIND= REAL64 ), DIMENSION( 12 )    :: SURFACEA               ! Surface area of every lattice combination
REAL( KIND= REAL64 ), DIMENSION( 12, 9 ) :: LATTICECOMB            ! Lattice combination

! Gottwald orthogonalization process
LREDUCTION: DO
  ! Lattice combination
  DO COMBINATION = 1, 12
    CALL LATTICE_COMBINATION( COMBINATION, V1, V2, V3, LATTICECOMB(COMBINATION,:) )
    ! New lattice vectors
    NEWV1(:) = LATTICECOMB(COMBINATION,1:3)
    NEWV2(:) = LATTICECOMB(COMBINATION,4:6)
    NEWV3(:) = LATTICECOMB(COMBINATION,7:9)
    ! Cross product of vx and vy
    CROSS1(1) = ( NEWV1(2) * NEWV2(3) ) - ( NEWV1(3) * NEWV2(2) )
    CROSS1(2) = ( NEWV1(3) * NEWV2(1) ) - ( NEWV1(1) * NEWV2(3) )
    CROSS1(3) = ( NEWV1(1) * NEWV2(2) ) - ( NEWV1(2) * NEWV2(1) )
    ! Cross product of vx and vz
    CROSS2(1) = ( NEWV1(2) * NEWV3(3) ) - ( NEWV1(3) * NEWV3(2) )
    CROSS2(2) = ( NEWV1(3) * NEWV3(1) ) - ( NEWV1(1) * NEWV3(3) )
    CROSS2(3) = ( NEWV1(1) * NEWV3(2) ) - ( NEWV1(2) * NEWV3(1) )
    ! Cross product of vy and vz
    CROSS3(1) = ( NEWV2(2) * NEWV3(3) ) - ( NEWV2(3) * NEWV3(2) )
    CROSS3(2) = ( NEWV2(3) * NEWV3(1) ) - ( NEWV2(1) * NEWV3(3) )
    CROSS3(3) = ( NEWV2(1) * NEWV3(2) ) - ( NEWV2(2) * NEWV3(1) )
    ! Vector moduli
    MODCROSS1 = ( CROSS1(1) * CROSS1(1) ) + ( CROSS1(2) * CROSS1(2) ) + ( CROSS1(3) * CROSS1(3) )
    MODCROSS1 = DSQRT( MODCROSS1 )
    MODCROSS2 = ( CROSS2(1) * CROSS2(1) ) + ( CROSS2(2) * CROSS2(2) ) + ( CROSS2(3) * CROSS2(3) )
    MODCROSS2 = DSQRT( MODCROSS2 )
    MODCROSS3 = ( CROSS3(1) * CROSS3(1) ) + ( CROSS3(2) * CROSS3(2) ) + ( CROSS3(3) * CROSS3(3) )
    MODCROSS3 = DSQRT( MODCROSS3 )
    ! Surface area
    SURFACEA(COMBINATION) = 2.D0 * MODCROSS1 + 2.D0 * MODCROSS2 + 2.D0 * MODCROSS3
  END DO
  ! Minimum surface area
  IF( MINVAL( SURFACEA ) <= MINSAREA ) THEN
    MINSAREALOC = MINLOC( SURFACEA, DIM= 1 )
    MINSAREA    = MINVAL( SURFACEA )
    ! Update box vectors
    V1(:) = LATTICECOMB(MINSAREALOC,1:3)
    V2(:) = LATTICECOMB(MINSAREALOC,4:6)
    V3(:) = LATTICECOMB(MINSAREALOC,7:9)
  ELSE
    EXIT LREDUCTION
  END IF
END DO LREDUCTION

RETURN

END SUBROUTINE GOTTWALD

! *********************************************************************************************** !
!                       This subroutine calculates the lattice combinations                       !
! *********************************************************************************************** !
SUBROUTINE LATTICE_COMBINATION( COMBINATION, V1, V2, V3, LATTICECOMB )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: COMBINATION ! Lattice combination

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: V1, V2, V3  ! Vectors vx, vy, and vz of the simulation box
REAL( KIND= REAL64 ), DIMENSION( 9 ) :: LATTICECOMB ! Lattice combination

! Lattice combination (unchanged determinant)
IF( COMBINATION == 1 ) THEN
  LATTICECOMB(1:3) = V1(:) + V2(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 2 ) THEN
  LATTICECOMB(1:3) = V1(:) - V2(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 3 ) THEN
  LATTICECOMB(1:3) = V1(:) + V3(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 4 ) THEN
  LATTICECOMB(1:3) = V1(:) - V3(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 5 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:) + V1(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 6 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:) - V1(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 7 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:) + V3(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 8 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:) - V3(:)
  LATTICECOMB(7:9) = V3(:)
ELSE IF( COMBINATION == 9 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:) + V1(:)
ELSE IF( COMBINATION == 10 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:) - V1(:)
ELSE IF( COMBINATION == 11 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:) + V2(:)
ELSE IF( COMBINATION == 12 ) THEN
  LATTICECOMB(1:3) = V1(:)
  LATTICECOMB(4:6) = V2(:)
  LATTICECOMB(7:9) = V3(:) - V2(:)
END IF

RETURN

END SUBROUTINE LATTICE_COMBINATION

! *********************************************************************************************** !
!  This subroutine applies the Lenstra–Lenstra–Lovász algorithm to orthogonalize the box vectors  !
! *********************************************************************************************** !
SUBROUTINE LLL( BOXV, DIM )

! Uses one module: global variables
USE GLOBALVAR

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 )                 :: I, J, K ! Counters
INTEGER( KIND= INT64 )                 :: DIM     ! Dimension
INTEGER( KIND= INT64 ), DIMENSION( 3 ) :: XYZ     ! XYZ positions
INTEGER( KIND= INT64 ), DIMENSION( 3 ) :: TEMPXYZ ! XYZ positions (temporary)

! *********************************************************************************************** !
! REAL VARIABLES (PARAMETER)                                                                      !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), PARAMETER :: DELTA = 0.75D0 ! Lovász delta

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                         :: LFACTOR       ! Lovász factor
REAL( KIND= REAL64 )                         :: LCONDITION    ! Lovász condition
REAL( KIND= REAL64 ), DIMENSION( DIM * DIM ) :: BOXV          ! Box vectors
REAL( KIND= REAL64 ), DIMENSION( DIM, DIM )  :: LATTICE_BASIS ! Lattice basis
REAL( KIND= REAL64 ), DIMENSION( DIM, DIM )  :: TEMP_BASIS    ! Lattice basis (temporary)
REAL( KIND= REAL64 ), DIMENSION( DIM, DIM )  :: OLD_LATTICE   ! Lattice basis (old)
REAL( KIND= REAL64 ), DIMENSION( DIM, DIM )  :: GS_BASIS      ! Gram-Schmidt orthonormal basis
REAL( KIND= REAL64 ), DIMENSION( DIM, DIM )  :: GS_COEFF      ! Gram-Schmidt coefficients

! Defining the lattice vectors
DO I = 1, DIM
  LATTICE_BASIS(:,I) = BOXV((1+DIM*(I-1)):(DIM*I))
END DO

! X, Y, and Z positions
XYZ(1) = 1
XYZ(2) = 2
XYZ(3) = 3

! Run LLL algorithm until convergence
LOVASZ_LOOP: DO

  ! Store old lattice vectors
  OLD_LATTICE(:,:) = LATTICE_BASIS(:,:)

  ! Initialization (working vector)
  K = 2

  ! Lenstra–Lenstra–Lovász algorithm
  DO WHILE( K <= DIM )
    ! Initialize/reset Gram-Schmidt basis vectors and coefficients
    DO I = 1, K
      ! Initialization
      GS_BASIS(:,I) = LATTICE_BASIS(:,I)
      ! Gram-Schimidt operators
      DO J = 1, I - 1
        ! Calculate the necessary Gram-Schmidt coefficients
        GS_COEFF(I,J) = DOT_PRODUCT( LATTICE_BASIS(:,I), GS_BASIS(:,J) ) / DOT_PRODUCT( GS_BASIS(:,J), GS_BASIS(:,J) )
        ! Calculate the necessary Gram-Schmidt orthonormal basis vectors
        GS_BASIS(:,I) = GS_BASIS(:,I) - GS_COEFF(I,J) * GS_BASIS(:,J)
      END DO
    END DO
    ! Calculate the Gram-Schmidt reduction (K > 1)
    DO J = 1, K - 1
      ! Recalculate the necessary Gram-Schmidt coefficients
      GS_COEFF(K,J) = DOT_PRODUCT( LATTICE_BASIS(:,K), GS_BASIS(:,J) ) / DOT_PRODUCT( GS_BASIS(:,J), GS_BASIS(:,J) )
      ! Compute the reduced lattice vector
      LATTICE_BASIS(:,K) = LATTICE_BASIS(:,K) - NINT( GS_COEFF(K,J) ) * LATTICE_BASIS(:,J)
    END DO
    ! Calculate the Lovász operators
    LFACTOR    = DELTA - ( GS_COEFF(K,K-1) * GS_COEFF(K,K-1) )
    LCONDITION = LFACTOR * DOT_PRODUCT( GS_BASIS(:,K-1), GS_BASIS(:,K-1) )
    ! Check the Lovász condition
    IF( DOT_PRODUCT( GS_BASIS(:,K), GS_BASIS(:,K) ) >= LCONDITION ) THEN
      ! Iteration
      K = K + 1
    ELSE
      ! Swap working vector and preceding lattice vector
      TEMP_BASIS(:,K)      = LATTICE_BASIS(:,K)
      LATTICE_BASIS(:,K)   = LATTICE_BASIS(:,K-1)
      LATTICE_BASIS(:,K-1) = TEMP_BASIS(:,K)
      ! Swap XYZ positions
      TEMPXYZ(K) = XYZ(K)
      XYZ(K)     = XYZ(K-1)
      XYZ(K-1)   = TEMPXYZ(K)
      ! Return one iteration of K until K = 2
      K = MAX( K-1, 2 )
    END IF
  END DO

  ! Compare new lattice vectors with old lattice vectors
  DO I = 1, DIM
    DO J = 1, DIM
      ! Cycle if they differ
      IF( DABS( OLD_LATTICE(J,I) - LATTICE_BASIS(J,I) ) >= EPSILON( 1.D0 ) ) CYCLE LOVASZ_LOOP
    END DO
  END DO

  EXIT LOVASZ_LOOP

END DO LOVASZ_LOOP

! Redefining the box vectors
IF( DIM == 3 ) THEN
  ! Undo LLL swaps (may cause the inversion of the determinant sign)
  BOXV(1:3) = LATTICE_BASIS(:,MINLOC( XYZ, DIM= 1 )) ! X vector
  BOXV(7:9) = LATTICE_BASIS(:,MAXLOC( XYZ, DIM= 1 )) ! Z vector
  DO I = 1, DIM
    IF( I /= MINLOC( XYZ, DIM= 1 ) .AND. I /= MAXLOC( XYZ, DIM= 1 ) ) THEN
      BOXV(4:6) = LATTICE_BASIS(:,I) ! Y vector
    END IF
  END DO
ELSE IF( DIM /= 3 ) THEN
  DO I = 1, DIM
    BOXV((1+DIM*(I-1)):(DIM*I)) = LATTICE_BASIS(:,I)
  END DO
END IF

RETURN

END SUBROUTINE LLL

! *********************************************************************************************** !
!             This subroutine generates a progress bar for the initial configuration.             !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR( J, TOTAL )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: J     ! Counter
INTEGER( KIND= INT64 ) :: TOTAL ! Total/Maximum number of cycles

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
INTEGER( KIND= INT64 ) :: J     ! Counter
INTEGER( KIND= INT64 ) :: TOTAL ! Total/Maximum number of cycles

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
