! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!          This module allows the user to choose one of the molecular configurations for          !
!              isomorphic molecules - Simple Cube (SC), Body-Centered Cube (BCC), or              !
!        Face-Centered Cube (FCC) -, and for anisomorphic molecules - Random Box (RND) or         !
!                                        Packed Box (PB).                                         !
!        Molecules will be then allocated in accordance to the selected crystal structure.        !
!                       See Macpherson et al. (2007) for some information.                        !
!     This module also writes out a file containing all particles' positions and quaternions.     !
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
! Main References:                 M. P. Allen, D. J. Tildesley                                   !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             doi: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                               G. Macpherson, M. K. Borg, J. Reese                               !
!              Molecular Simulation, Taylor & Francis, 2007, 33 (15), pp. 1199-1212.              !
!                                 doi: 10.1080/08927020701730724                                  !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

MODULE INITCONFIG

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
!                                       MOLECULAR GEOMETRY                                        !
!                                 Ellipsoids of revolution = EOR                                  !
!                                      Spherocylinders = SPC                                      !
!                                         Cylinders = CYL                                         !
! *********************************************************************************************** !

! *********************************************************************************************** !
!                                      INITIAL CONFIGURATION                                      !
!                                        Simple cube = SC                                         !
!                                    Face-centered cube = FCC                                     !
!                                    Body-centered cube = BCC                                     !
!                                        Random box = RND                                         !
!                                         Packed box = PB                                         !
! *********************************************************************************************** !

CONTAINS

! *********************************************************************************************** !
!             This subroutine allows the user to choose the geometry of the molecules             !
! *********************************************************************************************** !
SUBROUTINE GEOM_SELECTION(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = EOR                                                                                      !
!  (2) = SPC                                                                                      !
!  (3) = CYL                                                                                      !
! *********************************************************************************************** !
GEOM_SELEC(:) = .FALSE.

OPEN( UNIT= 100, FILE= "ini_config.ini" )
READ( 100, * ) GET, GEOM_INQ
CALL TO_UPPER( GEOM_INQ, LEN_TRIM( GEOM_INQ ), GEOM_INQ )
CLOSE( 100 )

! Extended configuration name
IF( GEOM_INQ == "EOR" ) THEN
  GEOMETRY    = "Ellipsoids of revolution"
  GEO_ACRONYM = "eor"
ELSE IF( GEOM_INQ == "SPC" ) THEN
  GEOMETRY    = "Spherocylinders"
  GEO_ACRONYM = "spc"
ELSE IF( GEOM_INQ == "CYL" ) THEN
  GEOMETRY    = "Cylinders"
  GEO_ACRONYM = "cyl"
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined ", TRIM( GEOM_INQ ), " is not an available molecular geometry. Exiting... "
  CALL EXIT(  )
END IF

! Initial configuration inquiry
WRITE( *, "(G0)") "The molecules are: "//TRIM( GEOMETRY )//". "
IF( GEOM_INQ == "EOR" ) THEN
  GEOM_SELEC(1) = .TRUE.
ELSE IF( GEOM_INQ == "SPC" ) THEN
  GEOM_SELEC(2) = .TRUE.
ELSE IF( GEOM_INQ == "CYL" ) THEN
  GEOM_SELEC(3) = .TRUE.
END IF
WRITE( *, * ) " "

RETURN

END SUBROUTINE GEOM_SELECTION

! *********************************************************************************************** !
!          This subroutine allows the user to choose the initial molecular configuration          !
! *********************************************************************************************** !
SUBROUTINE CONFIG_SELECTION(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = SC                                                                                       !
!  (2) = BCC                                                                                      !
!  (3) = FCC                                                                                      !
!  (4) = RND                                                                                      !
!  (5) = PB                                                                                       !
! *********************************************************************************************** !
CONFIG_SELEC(:) = .FALSE.

OPEN( UNIT= 100, FILE= "ini_config.ini" )
READ( 100, * ) GET, DUMMY
READ( 100, * ) GET, CONFIG_INQ
CALL TO_UPPER( CONFIG_INQ, LEN_TRIM( CONFIG_INQ ), CONFIG_INQ )
CLOSE( 100 )

! Extended configuration name
IF( CONFIG_INQ == "PB" ) THEN
  CONFIGURATION = "Packed Box"
ELSE IF( CONFIG_INQ == "RND" ) THEN
  CONFIGURATION = "Random"
ELSE IF( CONFIG_INQ == "SC" ) THEN
  CONFIGURATION = "Simple Cube"
ELSE IF( CONFIG_INQ == "BCC" ) THEN
  CONFIGURATION = "Body-Centered Cube"
ELSE IF( CONFIG_INQ == "FCC" ) THEN
  CONFIGURATION = "Face-Centered Cube"
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined ", TRIM( CONFIG_INQ ), " is not an available molecular configuration. Exiting... "
  CALL EXIT(  )
END IF

! Initial configuration inquiry
WRITE( *, "(G0)") "Initial configuration is: "//TRIM( CONFIGURATION )//". "
IF( CONFIG_INQ == "SC" ) THEN
  CONFIG_SELEC(1) = .TRUE.
ELSE IF( CONFIG_INQ == "BCC" ) THEN
  CONFIG_SELEC(2) = .TRUE.
ELSE IF( CONFIG_INQ == "FCC" ) THEN
  CONFIG_SELEC(3) = .TRUE.
ELSE IF( CONFIG_INQ == "RND" ) THEN
  CONFIG_SELEC(4) = .TRUE.
ELSE IF( CONFIG_INQ == "PB" ) THEN
  CONFIG_SELEC(5) = .TRUE.
END IF
WRITE( *, * ) " "

RETURN

END SUBROUTINE CONFIG_SELECTION

! *********************************************************************************************** !
!         This subroutine allocates particles according to the SC molecular configuration         !
! *********************************************************************************************** !
SUBROUTINE CONFIG_SC(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K, COUNTER  ! Counters
INTEGER( KIND= INT64 ) :: N_CELLS           ! Number of unit cells

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: CELL_LENGTH ! Length of unit cell (cubic structure)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: S12         ! Scaling factor (unit box)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AXIS_SELEC(1) ) THEN
  AXISN(:) = AXISX(:)
ELSE IF( AXIS_SELEC(2) ) THEN
  AXISN(:) = AXISY(:)
ELSE IF( AXIS_SELEC(3) ) THEN
  AXISN(:) = AXISZ(:)
END IF

! Convert degrees to radians
QUATERNION_ANGLE = QUATERNION_ANGLE * PI / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
Q(0,:) = DCOS( QUATERNION_ANGLE * 0.5D0 )             ! Real part
Q(1,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(1)  ! Imaginary part (Vector)
Q(2,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(2)  ! Imaginary part (Vector)
Q(3,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(3)  ! Imaginary part (Vector)

! Number of unit cells per axis (Simple Cube)
N_CELLS = NINT( DBLE( N_PARTICLES ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (Simple Cube)
CELL_LENGTH = ( 1.D0 / TOTAL_RHO ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BOX_LENGTH(:) = 0.D0
BOX_LENGTH(1) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(5) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(9) = CELL_LENGTH * DBLE( N_CELLS )

! Simulation box length (inverse)
CALL INVERSE_COF( BOX_LENGTH, BOX_LENGTH_I, BOX_VOLUME )

! Position of particles (centers of mass)
COUNTER = 1
DO I = 1, N_CELLS
  DO J = 1, N_CELLS
    DO K = 1, N_CELLS
      ! Particles on the right vertex of unit cell
      R(1,COUNTER) = DBLE( I - 1 ) * CELL_LENGTH
      R(2,COUNTER) = DBLE( J - 1 ) * CELL_LENGTH
      R(3,COUNTER) = DBLE( K - 1 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO I = 1, N_PARTICLES
  CALL MULTI_MATRIX( BOX_LENGTH_I, R(:,I), S12 )
  S12 = S12 - 0.5D0
  CALL MULTI_MATRIX( BOX_LENGTH, S12, R(:,I) )
END DO

RETURN

END SUBROUTINE CONFIG_SC

! *********************************************************************************************** !
!        This subroutine allocates particles according to the BCC molecular configuration         !
! *********************************************************************************************** !
SUBROUTINE CONFIG_BCC(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K, COUNTER  ! Counters
INTEGER( KIND= INT64 ) :: N_CELLS           ! Number of unit cells

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: CELL_LENGTH ! Length of unit cell (cubic structure)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: S12         ! Scaling factor (unit box)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AXIS_SELEC(1) ) THEN
  AXISN(:) = AXISX(:)
ELSE IF( AXIS_SELEC(2) ) THEN
  AXISN(:) = AXISY(:)
ELSE IF( AXIS_SELEC(3) ) THEN
  AXISN(:) = AXISZ(:)
END IF

! Convert degrees to radians
QUATERNION_ANGLE = QUATERNION_ANGLE * PI / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
Q(0,:) = DCOS( QUATERNION_ANGLE * 0.5D0 )             ! Real part
Q(1,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(1)  ! Imaginary part (Vector)
Q(2,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(2)  ! Imaginary part (Vector)
Q(3,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(3)  ! Imaginary part (Vector)

! Number of unit cells per axis (Body-Centered Cube)
N_CELLS = NINT( ( 0.5D0 * DBLE( N_PARTICLES ) ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (Body-Centered Cube)
CELL_LENGTH = ( 2.D0 / TOTAL_RHO ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BOX_LENGTH(:) = 0.D0
BOX_LENGTH(1) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(5) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(9) = CELL_LENGTH * DBLE( N_CELLS )

! Simulation box length (inverse)
CALL INVERSE_COF( BOX_LENGTH, BOX_LENGTH_I, BOX_VOLUME )

! Positioning of particles (centers of mass)
COUNTER = 1
DO I = 1, N_CELLS
  DO J = 1, N_CELLS
    DO K = 1, N_CELLS
      ! Particles on the right vertex of unit cell
      R(1,COUNTER) = DBLE( I - 1 ) * CELL_LENGTH
      R(2,COUNTER) = DBLE( J - 1 ) * CELL_LENGTH
      R(3,COUNTER) = DBLE( K - 1 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
      ! Particles on the center of unit cell
      R(1,COUNTER) = ( DBLE( I ) - 0.5D0 ) * CELL_LENGTH
      R(2,COUNTER) = ( DBLE( J ) - 0.5D0 ) * CELL_LENGTH
      R(3,COUNTER) = ( DBLE( K ) - 0.5D0 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO I = 1, N_PARTICLES
  CALL MULTI_MATRIX( BOX_LENGTH_I, R(:,I), S12 )
  S12 = S12 - 0.5D0
  CALL MULTI_MATRIX( BOX_LENGTH, S12, R(:,I) )
END DO

RETURN

END SUBROUTINE CONFIG_BCC

! *********************************************************************************************** !
!        This subroutine allocates particles according to the FCC molecular configuration         !
! *********************************************************************************************** !
SUBROUTINE CONFIG_FCC(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K, COUNTER  ! Counters
INTEGER( KIND= INT64 ) :: N_CELLS           ! Number of unit cells

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: CELL_LENGTH ! Length of unit cell (cubic structure)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: S12         ! Scaling factor (unit box)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AXIS_SELEC(1) ) THEN
  AXISN(:) = AXISX(:)
ELSE IF( AXIS_SELEC(2) ) THEN
  AXISN(:) = AXISY(:)
ELSE IF( AXIS_SELEC(3) ) THEN
  AXISN(:) = AXISZ(:)
END IF

! Convert degrees to radians
QUATERNION_ANGLE = QUATERNION_ANGLE * PI / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
Q(0,:) = DCOS( QUATERNION_ANGLE * 0.5D0 )             ! Real part
Q(1,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(1)  ! Imaginary part (Vector)
Q(2,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(2)  ! Imaginary part (Vector)
Q(3,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(3)  ! Imaginary part (Vector)

! Number of unit cells per axis (Face-Centered Cube)
N_CELLS = NINT( ( 0.25D0 * DBLE( N_PARTICLES ) ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (Face-Centered Cube)
CELL_LENGTH = ( 4.D0 / TOTAL_RHO ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BOX_LENGTH(:) = 0.D0
BOX_LENGTH(1) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(5) = CELL_LENGTH * DBLE( N_CELLS )
BOX_LENGTH(9) = CELL_LENGTH * DBLE( N_CELLS )

! Simulation box length (inverse)
CALL INVERSE_COF( BOX_LENGTH, BOX_LENGTH_I, BOX_VOLUME )

! Positioning of particles (centers of mass)
COUNTER = 1
DO I = 1, N_CELLS
  DO J = 1, N_CELLS
    DO K = 1, N_CELLS
      ! Particles on the right vertex of unit cell
      R(1,COUNTER) = DBLE( I - 1 ) * CELL_LENGTH
      R(2,COUNTER) = DBLE( J - 1 ) * CELL_LENGTH
      R(3,COUNTER) = DBLE( K - 1 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
      ! Particles on the front face of unit cell
      R(1,COUNTER) = DBLE( I - 1 ) * CELL_LENGTH
      R(2,COUNTER) = ( DBLE( J ) - 0.5D0 ) * CELL_LENGTH
      R(3,COUNTER) = ( DBLE( K ) - 0.5D0 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
      ! Particles on the left face of unit cell
      R(1,COUNTER) = ( DBLE( I ) - 0.5D0 ) * CELL_LENGTH
      R(2,COUNTER) = DBLE( J - 1 ) * CELL_LENGTH
      R(3,COUNTER) = ( DBLE( K ) - 0.5D0 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
      ! Particles on the lower face of unit cell
      R(1,COUNTER) = ( DBLE( I ) - 0.5D0 ) * CELL_LENGTH
      R(2,COUNTER) = ( DBLE( J ) - 0.5D0 ) * CELL_LENGTH
      R(3,COUNTER) = DBLE( K - 1 ) * CELL_LENGTH
      COUNTER = COUNTER + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO I = 1, N_PARTICLES
  CALL MULTI_MATRIX( BOX_LENGTH_I, R(:,I), S12 )
  S12 = S12 - 0.5D0
  CALL MULTI_MATRIX( BOX_LENGTH, S12, R(:,I) )
END DO

RETURN

END SUBROUTINE CONFIG_FCC

! *********************************************************************************************** !
!        This subroutine allocates particles according to a random molecular configuration        !
! *********************************************************************************************** !
SUBROUTINE CONFIG_RND(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K, L ! Counters
INTEGER( KIND= INT64 ) :: C, CI, CJ  ! Component index
INTEGER( KIND= INT64 ) :: ATTEMPTS   ! Counter
INTEGER( KIND= INT64 ) :: CYCLES     ! Counter of cycles
INTEGER( KIND= INT64 ) :: NACCT      ! Move acceptance counter: Translation
INTEGER( KIND= INT64 ) :: NACCR      ! Move acceptance counter: Rotation
INTEGER( KIND= INT64 ) :: NACCVI     ! Move acceptance counter: Isotropic volume change
INTEGER( KIND= INT64 ) :: NACCVA     ! Move acceptance counter: Anistropic volume change
INTEGER( KIND= INT64 ) :: MOVT       ! Move counter (Translation)
INTEGER( KIND= INT64 ) :: MOVR       ! Move counter (Rotation)
INTEGER( KIND= INT64 ) :: MOVVI      ! Move counter (Isotropic volume change)
INTEGER( KIND= INT64 ) :: MOVVA      ! Move counter (Anisotropic volume change)
INTEGER( KIND= INT64 ) :: COMPONENT  ! Box matrix component

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                          :: CD               ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( KIND= REAL64 )                          :: RIJSQ            ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                          :: CUTOFF_D         ! Cutoff distance
REAL( KIND= REAL64 )                          :: BOX_VOLUME_NVT   ! Box volume (NVT Simulation)
REAL( KIND= REAL64 )                          :: P_FRACTION_NVT   ! Packing fraction (NVT Simulation)
REAL( KIND= REAL64 )                          :: BOXVMC_RND       ! Box volume (NPT Simulation)
REAL( KIND= REAL64 )                          :: BOXVM, BOXVN     ! Box volume (before/after a trial move)
REAL( KIND= REAL64 )                          :: SCALE_FACTOR     ! Scale factor of the volume of the simulation box
REAL( KIND= REAL64 )                          :: HNM              ! Enthalpy criterion (reduced)
REAL( KIND= REAL64 )                          :: ETA_NPT          ! Packing fraction (NPT Simulation)
REAL( KIND= REAL64 )                          :: RATIO            ! Acceptance ratio (Simulation)
REAL( KIND= REAL64 )                          :: DISTORTION       ! Box distortion
REAL( KIND= REAL64 )                          :: BOXVROT          ! Volume of simulation box (after undoing box rotation)
REAL( KIND= REAL64 )                          :: THETA            ! Angle between box vector and coordination system
REAL( KIND= REAL64 )                          :: RAXISMAG         ! Magnitude of rotation axis
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: PROJY_XY         ! Projection of the y-vector of the box onto the ZY-plane and the unit vector of the y-axis
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RAXIS            ! Rotation axis
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: AUXV             ! Auxiliary vector
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: V1, V2, V3       ! Box vectors
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOXLROT, BOXIROT ! Length of simulation box (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )        :: QROT             ! Rotation quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 )        :: QAUX             ! Auxiliary quaternion
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOX_LENGTH_NVT   ! Box length (NVT Simulation)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOX_LENGTH_NVT_I ! Inverse of box length (NVT Simulation)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOXLMC           ! Box length (NPT Simulation)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOXLMC_I         ! Inverse of box length (NPT Simulation)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOXLM, BOXLN     ! Box length (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 9 )          :: BOXLM_I, BOXLN_I ! Inverse of box length (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: COSANGLE_VEC     ! Cossine of angle between box vectors
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: LBOX             ! Length of box edges
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: LBOXR            ! Length ratio of box edges
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: S12              ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RM, RN           ! Position (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: EM, EN           ! Orientation (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )        :: QM, QN           ! Quaternion (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RIJ              ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: RI, RJ           ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )          :: EI, EJ           ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 0:3 )        :: QI, QJ           ! Quaternions of particles i and j
REAL( KIND= REAL64 ), DIMENSION( COMPONENTS ) :: CUTOFF           ! Cutoff diameter

! *********************************************************************************************** !
! REAL VARIABLES (Allocatable)                                                                    !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: RMCV  ! Old position of particles
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: RPROT ! Position of the center of mass (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: QPROT ! Quaternion of the center of mass (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: EPROT ! Orientation of the center of mass (after undoing box rotation)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP            ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_PRELIMINAR ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_VALIDATION ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PARALLEL           ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: DISABLE_TRANS_ADJ  ! Disable translational adjustments
LOGICAL :: DISABLE_ROT_ADJ    ! Disable rotational adjustments
LOGICAL :: LATTICER           ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: MOV_ROT            ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_TRANS          ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_VOL_I          ! Isotropic volume change selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_VOL_A          ! Anisotropic volume change selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: IGNORE             ! Detects if a box deformation is valid or not : TRUE = ignore box deformation; FALSE = consider box deformation

! Allocation
ALLOCATE( RMCV(3,N_PARTICLES) )
ALLOCATE( QPROT(0:3,N_PARTICLES) )
ALLOCATE( RPROT(3,N_PARTICLES) )
ALLOCATE( EPROT(3,N_PARTICLES) )

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

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AXIS_SELEC(1) ) THEN
  AXISN(:) = AXISX(:)
ELSE IF( AXIS_SELEC(2) ) THEN
  AXISN(:) = AXISY(:)
ELSE IF( AXIS_SELEC(3) ) THEN
  AXISN(:) = AXISZ(:)
END IF

! *********************************************************************************************** !
! Convert degrees to radians                                                                      !
! *********************************************************************************************** !
QUATERNION_ANGLE = QUATERNION_ANGLE * PI / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
Q(0,:) = DCOS( QUATERNION_ANGLE * 0.5D0 )             ! Real part
Q(1,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(1)  ! Imaginary part (Vector)
Q(2,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(2)  ! Imaginary part (Vector)
Q(3,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(3)  ! Imaginary part (Vector)

! *********************************************************************************************** !
! Box length (cube)                                                                               !
! *********************************************************************************************** !
BOX_LENGTH(:) = 0.D0
BOX_LENGTH(1) = BOX_VOLUME ** (1.D0 / 3.D0)
BOX_LENGTH(5) = BOX_VOLUME ** (1.D0 / 3.D0)
BOX_LENGTH(9) = BOX_VOLUME ** (1.D0 / 3.D0)

! *********************************************************************************************** !
! Simulation box length (inverse)                                                                 !
! *********************************************************************************************** !
CALL INVERSE_COF( BOX_LENGTH, BOX_LENGTH_I, BOX_VOLUME )

! *********************************************************************************************** !
! Active Transformation (Orientation)                                                             !
! *********************************************************************************************** !
DO I = 1, N_PARTICLES
  CALL ACTIVE_TRANSFORMATION( AXISZ, Q(:,I), E(:,I) )
END DO

! *********************************************************************************************** !
! Packing fraction (NVT Simulation)                                                               !
! *********************************************************************************************** !
P_FRACTION_NVT = ETA_INI

! *********************************************************************************************** !
! Box length (NVT Simulation)                                                                     !
! *********************************************************************************************** !
BOX_VOLUME_NVT = 0.D0
DO C = 1, COMPONENTS
  BOX_VOLUME_NVT = BOX_VOLUME_NVT + N_COMPONENT(C) * PARTICLE_VOL(C)
END DO
! Box volume (large cubic box)
BOX_VOLUME_NVT = BOX_VOLUME_NVT / P_FRACTION_NVT
! Box length (large cubic box)
BOX_LENGTH_NVT(:) = 0.D0
BOX_LENGTH_NVT(1) = (BOX_VOLUME_NVT) ** (1.D0 / 3.D0)
BOX_LENGTH_NVT(5) = (BOX_VOLUME_NVT) ** (1.D0 / 3.D0)
BOX_LENGTH_NVT(9) = (BOX_VOLUME_NVT) ** (1.D0 / 3.D0)

! *********************************************************************************************** !
! Inverse of box length (NVT Simulation)                                                          !
! *********************************************************************************************** !
CALL INVERSE_COF( BOX_LENGTH_NVT, BOX_LENGTH_NVT_I, BOX_VOLUME_NVT )

! *********************************************************************************************** !
! Positioning of particles (centers of mass)                                                      !
! *********************************************************************************************** !
R(:,:) = 0.D0

! *********************************************************************************************** !
! Monte Carlo parameters (NVT Simulation)                                                         !
! *********************************************************************************************** !
DISABLE_TRANS_ADJ = .FALSE. ! Translational adjustments             (initial value)
DISABLE_ROT_ADJ   = .FALSE. ! Rotational adjustments                (initial value)
MOV_TRANS         = .FALSE. ! Translational move selector           (initial value)
MOV_ROT           = .FALSE. ! Rotational move selector              (initial value)
IF( GEOM_SELEC(1) ) THEN
  IF( MAXVAL( DIAMETER ) <= MAXVAL( LENGTH ) ) THEN
    DRMAX         = 1.05D0 * MAXVAL( LENGTH )   ! Maximum translational displacement (initial value)
  ELSE
    DRMAX         = 1.05D0 * MAXVAL( DIAMETER ) ! Maximum translational displacement (initial value)
  END IF
ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
  DRMAX           = 1.05D0 * ( MAXVAL( LENGTH ) + MAXVAL( DIAMETER ) ) ! Maximum translational displacement (initial value)
END IF
ANGMAX            = 0.1D0   ! Maximum rotational displacement       (initial value)
NACCT             = 0       ! Translational move acceptance counter (initial value)
NACCR             = 0       ! Rotational move acceptance counter    (initial value)
MOVT              = 0       ! Translational move counter            (initial value)
MOVR              = 0       ! Rotational move counter               (initial value)
QMC(:,:)          = Q(:,:)  ! Quaternion algebra                    (initial value)
RMC(:,:)          = R(:,:)  ! Position of particles                 (initial value)
EMC(:,:)          = E(:,:)  ! Orientation of particles              (initial value)
ATTEMPTS          = 0       ! Number of attempts                    (initial value)

! Summary
WRITE( *, "(G0)" ) "Attempting to randomly distribute particles inside a cubic box. It may take a while..."
WRITE( *, "(G0)" ) " "
CALL SLEEP( 1 )

! *********************************************************************************************** !
! Hit-and-miss (NVT Simulation)                                                                   !
! *********************************************************************************************** !
HIT_AND_MISS_NVT: DO

  DO CYCLES = 1, N_PARTICLES

    ! Component index
    IF( COMPONENTS > 1 ) THEN
      ! Pseudorandom number generator (uniform distribution)
      CALL RANF(  )
      ! Get component index
      CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
    ELSE
      ! Get component index
      CI = 1
    END IF

    ! Avoid components with molar fraction of 0
    DO WHILE( N_COMPONENT(CI) < 1 )
      ! Component index
      IF( COMPONENTS > 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        CALL RANF(  )
        ! Get component index
        CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
      ELSE
        ! Get component index
        CI = 1
      END IF
    END DO

    ! Forbid rotation if component is spherical
    IF( GEOM_SELEC(1) .AND. DABS( ASPECT_RATIO(CI) - 1.D0 ) < EPSILON( 1.D0 ) ) THEN
      MOV_TRANS = .TRUE.   ! Enable translation
      MOV_ROT   = .FALSE.  ! Disable rotation
      MOVT      = MOVT + 1 ! Increment move counter
    ELSE IF( GEOM_SELEC(2) .AND. DABS( ASPECT_RATIO(CI) - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
      MOV_TRANS = .TRUE.   ! Enable translation
      MOV_ROT   = .FALSE.  ! Disable rotation
      MOVT      = MOVT + 1 ! Increment move counter
    ! Allow rotation if component is nonspherical
    ELSE
      ! Pseudorandom number generator (uniform distribution)
      CALL RANF(  )
      ! Translation criterion
      IF( RANDOM_N < PROB_TRANS_INIT ) THEN
        MOV_TRANS = .TRUE.   ! Enable translation
        MOV_ROT   = .FALSE.  ! Disable rotation
        MOVT      = MOVT + 1 ! Increment move counter
      ! Rotation criterion
      ELSE IF( RANDOM_N >= PROB_TRANS_INIT ) THEN
        MOV_ROT   = .TRUE.   ! Enable rotation
        MOV_TRANS = .FALSE.  ! Disable translation
        MOVR      = MOVR + 1 ! Increment move counter
      END IF
    END IF

    ! Pseudorandom number generator (uniform distribution)
    CALL RANF(  )
    ! Random selection of particles of component C
    I = SUM( N_COMPONENT(0:(CI-1)) ) + INT( RANDOM_N * DBLE( N_COMPONENT(CI) ) ) + 1

    ! Assignment of previous configuration (Microstate m)
    RM(:) = RMC(:,I) ! Position
    QM(:) = QMC(:,I) ! Quaternion
    EM(:) = EMC(:,I) ! Orientation

    ! Translation Movement
    IF( MOV_TRANS ) THEN
      ! Random translation along x-axis
      CALL RANF(  )
      RN(1) = RM(1) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
      ! Random translation along y-axis
      CALL RANF(  )
      RN(2) = RM(2) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
      ! Random translation along z-axis
      CALL RANF(  )
      RN(3) = RM(3) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
      ! Minimum Image Convention
      CALL MULTI_MATRIX( BOX_LENGTH_NVT_I, RN, S12 )
      S12 = S12 - ANINT( S12 )
      CALL MULTI_MATRIX( BOX_LENGTH_NVT, S12, RN )
    ! No Translation
    ELSE IF( .NOT. MOV_TRANS ) THEN
      RN(:) = RM(:)
    END IF

    ! Rotation Movement
    IF( MOV_ROT ) THEN
      ! Random Composed Unit Quaternion
      CALL COMPOSED_QUATERNION( QM, QN, ANGMAX )
      ! Active transformation
      CALL ACTIVE_TRANSFORMATION( AXISZ, QN, EN )
    ! No Rotation
    ELSE IF( .NOT. MOV_ROT ) THEN
      QN(:) = QM(:)
      EN(:) = EM(:)
    END IF

    ! Overlap Check
    CALL CHECK_OVERLAP( CI, I, QN, EN, RN, CD, BOX_LENGTH_NVT, BOX_LENGTH_NVT_I, OVERLAP )

    ! Acceptance Criterion
    IF( .NOT. OVERLAP ) THEN
      ! System configuration update
      RMC(:,I) = RN(:) ! Update position
      QMC(:,I) = QN(:) ! Update quaternion
      EMC(:,I) = EN(:) ! Update orientation
      ! Displacement counter update
      IF( MOV_TRANS ) THEN
        NACCT = NACCT + 1  ! Translational move counter
      ELSE IF ( MOV_ROT ) THEN
        NACCR = NACCR + 1  ! Rotational move counter
      END IF
    ELSE
      ! Retrieve old configuration
      RMC(:,I) = RM(:) ! Retrieve position
      QMC(:,I) = QM(:) ! Retrieve quaternion
      EMC(:,I) = EM(:) ! Retrieve orientation
    END IF

  END DO

  ! Adjustment of maximum displacement (Translation and Rotation)
  IF( MOD( (MOVT + MOVR), (N_ADJUST_INIT * N_PARTICLES) ) == 0 ) THEN

    IF( .NOT. DISABLE_TRANS_ADJ ) THEN
      ! Acceptance ratio (translation)
      RATIO = DBLE( NACCT ) / DBLE( MOVT )
      ! Translational adjustment
      IF( RATIO <= R_ACC_T ) THEN
        DRMAX = 0.95D0 * DRMAX
      ELSE
        DRMAX = 1.05D0 * DRMAX
      END IF
    END IF

    ! Limiting displacement (translation)
    IF( DRMAX <= 1.D-2 .AND. .NOT. DISABLE_TRANS_ADJ ) THEN
      IF( MAXVAL(DIAMETER) <= MAXVAL(LENGTH) ) THEN
        DRMAX = 3.5D0 * MAXVAL(LENGTH)
      ELSE
        DRMAX = 3.5D0 * MAXVAL(DIAMETER)
      END IF
      DISABLE_TRANS_ADJ = .TRUE.
    END IF

    IF( .NOT. DISABLE_ROT_ADJ ) THEN
      ! Acceptance ratio (rotation)
      RATIO = DBLE( NACCR ) / DBLE( MOVR )
      ! Rotational adjustment
      IF( RATIO <= R_ACC_R ) THEN
        ANGMAX = 0.95D0 * ANGMAX
      ELSE
        ANGMAX = 1.05D0 * ANGMAX
      END IF
      ! 4π-rotation condition
      IF( ( ANGMAX > 4.D0 * PI ) ) THEN
        ANGMAX = ANGMAX - 2.D0 * PI
      END IF
    END IF

    ! Limiting displacement (rotation)
    IF( ANGMAX <= 1.D-2 .AND. .NOT. DISABLE_ROT_ADJ ) THEN
      ANGMAX = 0.5D0 * PI
      DISABLE_ROT_ADJ = .TRUE.
    END IF

    ! Reset counter
    NACCT = 0
    MOVT  = 0
    NACCR = 0
    MOVR  = 0

  END IF

  ! Iteration
  ATTEMPTS = ATTEMPTS + 1

  ! Overlap check for the proposed initial configuration
  IF( MOD( ATTEMPTS, MAX_ATTEMPTS ) == 0 ) THEN

    ! Progress bar
    CALL PROGRESS_BAR_HITMISS( ATTEMPTS, MAX_ATTEMPTS )

    ! Finalization condition
    IF( ATTEMPTS / MAX_ATTEMPTS >= 9999 ) THEN
      WRITE( *, * ) " "
      WRITE( *, * ) " "
      WRITE( *, "(G0)" ) "Too many attempts! Try decreasing the initial packing fraction of the system. Exiting..."
      CALL EXIT(  )
    END IF

    ! Initial configuration (partial)
    OPEN( UNIT= 55, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/initconf_rnd_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 55, "(I4)" ) N_PARTICLES
    WRITE( 55, * ) " "
    IF( GEOM_SELEC(1) ) THEN
      DO C = 1, COMPONENTS
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 55, * ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
          &              0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
        END DO
      END DO
    ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
      DO C = 1, COMPONENTS
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 55, * ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
          &              0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
        END DO
      END DO
    END IF
    CLOSE( 55 )

    ! Validation loop
    LOOP_VALIDATION_INITIAL_CONF: DO

      ! Initialization
      OVERLAP_VALIDATION = .FALSE.

      ! Anisomorphic molecules (unlike components)
      DO CI = 1, COMPONENTS - 1
        DO CJ = CI + 1, COMPONENTS
          ! First loop represents all particles with indexes i of component Ci
          DO I = SUM( N_COMPONENT(0:(CI-1)) ) + 1, SUM( N_COMPONENT(0:CI) )
            ! Second loop represents all particles with indexes j of component Cj
            DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, SUM( N_COMPONENT(0:CJ) )
              ! Position of particle i
              RI(1)  = RMC(1,I)
              RI(2)  = RMC(2,I)
              RI(3)  = RMC(3,I)
              ! Position of particle j
              RJ(1)  = RMC(1,J)
              RJ(2)  = RMC(2,J)
              RJ(3)  = RMC(3,J)
              ! Orientation of particle i
              EI(1)  = EMC(1,I)
              EI(2)  = EMC(2,I)
              EI(3)  = EMC(3,I)
              ! Orientation of particle j
              EJ(1)  = EMC(1,J)
              EJ(2)  = EMC(2,J)
              EJ(3)  = EMC(3,J)
              ! Quaternion of particle i
              QI(0)  = QMC(0,I)
              QI(1)  = QMC(1,I)
              QI(2)  = QMC(2,I)
              QI(3)  = QMC(3,I)
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
              CALL MULTI_MATRIX( BOX_LENGTH_NVT_I, RIJ, S12 )
              S12 = S12 - ANINT( S12 )
              CALL MULTI_MATRIX( BOX_LENGTH_NVT, S12, RIJ )
              ! Magnitude of the vector distance (squared)
              RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
              ! Cutoff distance (squared)
              CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
              CUTOFF_D = CUTOFF_D * CUTOFF_D
              ! Preliminary test (circumscribing spheres)
              IF( RIJSQ <= CUTOFF_D ) THEN
                ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
                IF( GEOM_SELEC(1) ) THEN
                  CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP_VALIDATION )
                  ! Overlap criterion
                  IF( OVERLAP_VALIDATION ) THEN
                    ! Overlap detected
                    CYCLE HIT_AND_MISS_NVT
                  END IF
                ! Overlap test for spherocylinders (Vega-Lago Method)
                ELSE IF( GEOM_SELEC(2) ) THEN
                  CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_VALIDATION )
                  ! Overlap criterion
                  IF( OVERLAP_VALIDATION ) THEN
                    ! Overlap detected
                    CYCLE HIT_AND_MISS_NVT
                  END IF
                ! Overlap test for cylinders (Lopes et al. Method)
                ELSE IF( GEOM_SELEC(3) ) THEN
                  ! Preliminary test (circumscribing spherocylinders)
                  OVERLAP_PRELIMINAR = .FALSE.
                  CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                  ! Overlap criterion
                  IF( OVERLAP_PRELIMINAR ) THEN
                    ! Retrive position of the particle j after applying the PBC
                    RJ(1) = RI(1) + RIJ(1)
                    RJ(2) = RI(2) + RIJ(2)
                    RJ(3) = RI(3) + RIJ(3)
                    CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP_VALIDATION )
                    ! Overlap criterion
                    IF( OVERLAP_VALIDATION ) THEN
                      ! Overlap detected
                      CYCLE HIT_AND_MISS_NVT
                    END IF
                  END IF
                END IF
              END IF
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
            RI(1)  = RMC(1,I)
            RI(2)  = RMC(2,I)
            RI(3)  = RMC(3,I)
            ! Position of particle j
            RJ(1)  = RMC(1,J)
            RJ(2)  = RMC(2,J)
            RJ(3)  = RMC(3,J)
            ! Orientation of particle i
            EI(1)  = EMC(1,I)
            EI(2)  = EMC(2,I)
            EI(3)  = EMC(3,I)
            ! Orientation of particle j
            EJ(1)  = EMC(1,J)
            EJ(2)  = EMC(2,J)
            EJ(3)  = EMC(3,J)
            ! Quaternion of particle i
            QI(0)  = QMC(0,I)
            QI(1)  = QMC(1,I)
            QI(2)  = QMC(2,I)
            QI(3)  = QMC(3,I)
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
            CALL MULTI_MATRIX( BOX_LENGTH_NVT_I, RIJ, S12 )
            S12 = S12 - ANINT( S12 )
            CALL MULTI_MATRIX( BOX_LENGTH_NVT, S12, RIJ )
            ! Magnitude of the vector distance (squared)
            RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
            ! Cutoff distance (squared)
            CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
            CUTOFF_D = CUTOFF_D * CUTOFF_D
            ! Preliminary test (circumscribing spheres)
            IF( RIJSQ <= CUTOFF_D ) THEN
              ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
              IF( GEOM_SELEC(1) ) THEN
                CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP_VALIDATION )
                ! Overlap criterion
                IF( OVERLAP_VALIDATION ) THEN
                  ! Overlap detected
                  CYCLE HIT_AND_MISS_NVT
                END IF
              ! Overlap test for spherocylinders (Vega-Lago Method)
              ELSE IF( GEOM_SELEC(2) ) THEN
                CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_VALIDATION )
                ! Overlap criterion
                IF( OVERLAP_VALIDATION ) THEN
                  ! Overlap detected
                  CYCLE HIT_AND_MISS_NVT
                END IF
              ! Overlap test for cylinders (Lopes et al. Method)
              ELSE IF( GEOM_SELEC(3) ) THEN
                ! Preliminary test (circumscribing spherocylinders)
                OVERLAP_PRELIMINAR = .FALSE.
                CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                ! Overlap criterion
                IF( OVERLAP_PRELIMINAR ) THEN
                  ! Retrive position of the particle j after applying the PBC
                  RJ(1) = RI(1) + RIJ(1)
                  RJ(2) = RI(2) + RIJ(2)
                  RJ(3) = RI(3) + RIJ(3)
                  CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP_VALIDATION )
                  ! Overlap criterion
                  IF( OVERLAP_VALIDATION ) THEN
                    ! Overlap detected
                    CYCLE HIT_AND_MISS_NVT
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO

      ! Possible initial configuration
      IF( .NOT. OVERLAP_VALIDATION ) THEN
        EXIT HIT_AND_MISS_NVT
      END IF

    END DO LOOP_VALIDATION_INITIAL_CONF

  END IF

END DO HIT_AND_MISS_NVT

! *********************************************************************************************** !
! Reassign positions, rotation quaternions, and orientations                                      !
! *********************************************************************************************** !
Q(:,:) = QMC(:,:) ! Quaternion of particles
R(:,:) = RMC(:,:) ! Position of particles
E(:,:) = EMC(:,:) ! Orientation of particles

! Summary
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( (ATTEMPTS / MAX_ATTEMPTS) == 1 ) THEN
  WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", &
  &                   ATTEMPTS / MAX_ATTEMPTS, " attempt."
ELSE IF( (ATTEMPTS / MAX_ATTEMPTS) > 1 ) THEN
  WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", &
  &                   ATTEMPTS / MAX_ATTEMPTS, " attempts."
END IF
CALL SLEEP( 3 )
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Now running an NPT simulation up to the target packing fraction of ", PACKING_F,"..."
WRITE( *, "(G0)" ) " "
CALL SLEEP( 1 )

! *********************************************************************************************** !
! Monte Carlo parameters (NPT Simulation)                                                         !
! *********************************************************************************************** !
DISABLE_TRANS_ADJ = .FALSE.             ! Translational adjustments             (initial value)
DISABLE_ROT_ADJ   = .FALSE.             ! Rotational adjustments                (initial value)
MOV_TRANS         = .FALSE.             ! Translational move selector           (initial value)
MOV_ROT           = .FALSE.             ! Rotational move selector              (initial value)
MOV_VOL_I         = .FALSE.             ! Volume move selector                  (initial value)
MOV_VOL_A         = .FALSE.             ! Volume move selector                  (initial value)
DRMAX             = DRMAX_INIT          ! Maximum translational displacement    (initial value)
ANGMAX            = ANGMAX_INIT         ! Maximum rotational displacement       (initial value)
DVMAXISO          = DVMAXISO_INIT       ! Maximum isovolumetric displacement    (initial value)
DVMAXANI          = DVMAXANISO_INIT     ! Maximum anisovolumetric displacement  (initial value)
NACCT             = 0                   ! Translational move acceptance counter (initial value)
NACCR             = 0                   ! Rotational move acceptance counter    (initial value)
NACCVI            = 0                   ! Volumetric move acceptance counter    (initial value)
NACCVA            = 0                   ! Volumetric move acceptance counter    (initial value)
MOVT              = 0                   ! Translational move counter            (initial value)
MOVR              = 0                   ! Rotational move counter               (initial value)
MOVVI             = 0                   ! Volume change counter                 (initial value)
MOVVA             = 0                   ! Volume change counter                 (initial value)
QMC(:,:)          = Q(:,:)              ! Quaternion algebra                    (initial value)
RMC(:,:)          = R(:,:)              ! Position of particles                 (initial value)
EMC(:,:)          = E(:,:)              ! Orientation of particles              (initial value)
BOXLMC(:)         = BOX_LENGTH_NVT(:)   ! Box length                            (initial value)
BOXLMC_I(:)       = BOX_LENGTH_NVT_I(:) ! Inverse of box length                 (initial value)
BOXVMC_RND        = BOX_VOLUME_NVT      ! Box volume                            (initial value)
ATTEMPTS          = 0                   ! Number of attempts                    (initial value)
ETA_NPT           = P_FRACTION_NVT      ! Packing fraction                      (initial value)

! Isobaric-Isothermal Monte Carlo Simulation
NPT_SIMULATION: DO

  ! Choose between displacement of molecules or volume change
  CALL RANF(  )
  IF( RANDOM_N < PROB_MOV_INIT ) THEN
    MOV_TRANS = .TRUE.  ! Enable translation
    MOV_ROT   = .TRUE.  ! Enable rotation
    MOV_VOL_I = .FALSE. ! Disable volume change
    MOV_VOL_A = .FALSE. ! Disable volume change
  ELSE IF( RANDOM_N >= PROB_MOV_INIT ) THEN
    MOV_TRANS = .FALSE. ! Disable translation
    MOV_ROT   = .FALSE. ! Disable rotation
    MOV_VOL_I = .TRUE.  ! Enable volume change
    MOV_VOL_A = .TRUE.  ! Enable volume change
  END IF

  ! Movement (translation or rotation)
  IF( MOV_TRANS .OR. MOV_ROT ) THEN

    DO CYCLES = 1, N_PARTICLES

      ! Component index
      IF( COMPONENTS > 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        CALL RANF(  )
        ! Get component index
        CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
      ELSE
        ! Get component index
        CI = 1
      END IF

      ! Avoid components with molar fraction of 0
      DO WHILE( N_COMPONENT(CI) < 1 )
        ! Component index
        IF( COMPONENTS > 1 ) THEN
          ! Pseudorandom number generator (uniform distribution)
          CALL RANF(  )
          ! Get component index
          CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
        ELSE
          ! Get component index
          CI = 1
        END IF
      END DO

      ! Forbid rotation if component is spherical
      IF( GEOM_SELEC(1) .AND. DABS( ASPECT_RATIO(CI) - 1.D0 ) < EPSILON( 1.D0 ) ) THEN
        MOV_TRANS = .TRUE.   ! Enable translation
        MOV_ROT   = .FALSE.  ! Disable rotation
        MOVT      = MOVT + 1 ! Increment move counter
      ELSE IF( GEOM_SELEC(2) .AND. DABS( ASPECT_RATIO(CI) - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
        MOV_TRANS = .TRUE.   ! Enable translation
        MOV_ROT   = .FALSE.  ! Disable rotation
        MOVT      = MOVT + 1 ! Increment move counter
      ! Allow rotation if component is nonspherical
      ELSE
        ! Pseudorandom number generator (uniform distribution)
        CALL RANF(  )
        ! Translation criterion
        IF( RANDOM_N < PROB_TRANS_INIT ) THEN
          MOV_TRANS = .TRUE.   ! Enable translation
          MOV_ROT   = .FALSE.  ! Disable rotation
          MOVT      = MOVT + 1 ! Increment move counter
        ! Rotation criterion
        ELSE IF( RANDOM_N >= PROB_TRANS_INIT ) THEN
          MOV_ROT   = .TRUE.   ! Enable rotation
          MOV_TRANS = .FALSE.  ! Disable translation
          MOVR      = MOVR + 1 ! Increment move counter
        END IF
      END IF

      ! Pseudorandom number generator (uniform distribution)
      CALL RANF(  )
      ! Random selection of particles of component C
      I = SUM( N_COMPONENT(0:(CI-1)) ) + INT( RANDOM_N * DBLE( N_COMPONENT(CI) ) ) + 1

      ! Assignment of previous configuration (Microstate m)
      RM(:) = RMC(:,I) ! Position
      QM(:) = QMC(:,I) ! Quaternion
      EM(:) = EMC(:,I) ! Orientation

      ! Translation Movement
      IF( MOV_TRANS ) THEN
        ! Random translation along x-axis
        CALL RANF(  )
        RN(1) = RM(1) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
        ! Random translation along y-axis
        CALL RANF(  )
        RN(2) = RM(2) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
        ! Random translation along z-axis
        CALL RANF(  )
        RN(3) = RM(3) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
        ! Minimum Image Convention
        CALL MULTI_MATRIX( BOXLMC_I, RN, S12 )
        S12 = S12 - ANINT( S12 )
        CALL MULTI_MATRIX( BOXLMC, S12, RN )
      ! No Translation
      ELSE IF( .NOT. MOV_TRANS ) THEN
        RN(:) = RM(:)
      END IF

      ! Rotation Movement
      IF( MOV_ROT ) THEN
        ! Random Composed Unit Quaternion
        CALL COMPOSED_QUATERNION( QM, QN, ANGMAX )
        ! Active transformation
        CALL ACTIVE_TRANSFORMATION( AXISZ, QN, EN )
      ! No Rotation
      ELSE IF( .NOT. MOV_ROT ) THEN
        QN(:) = QM(:)
        EN(:) = EM(:)
      END IF

      ! Overlap Check
      CALL CHECK_OVERLAP( CI, I, QN, EN, RN, CD, BOXLMC, BOXLMC_I, OVERLAP )

      ! Acceptance Criterion
      IF( .NOT. OVERLAP ) THEN
        ! System configuration update
        RMC(:,I) = RN(:) ! Update position
        QMC(:,I) = QN(:) ! Update quaternion
        EMC(:,I) = EN(:) ! Update orientation
        ! Displacement counter update
        IF( MOV_TRANS ) THEN
          NACCT = NACCT + 1  ! Translational move counter
        ELSE IF ( MOV_ROT ) THEN
          NACCR = NACCR + 1  ! Rotational move counter
        END IF
      ELSE
        ! Retrieve old configuration
        RMC(:,I) = RM(:) ! Retrieve position
        QMC(:,I) = QM(:) ! Retrieve quaternion
        EMC(:,I) = EM(:) ! Retrieve orientation
      END IF

    END DO

  ELSE IF( MOV_VOL_I .OR. MOV_VOL_A ) THEN

    ! Assignment of previous configuration (microstate m)   
    BOXLM(:)   = BOXLMC(:)   ! Box length
    BOXLM_I(:) = BOXLMC_I(:) ! Box length (inverse)
    BOXVM      = BOXVMC_RND  ! Box volume

    ! Expansion/compression type
    CALL RANF(  )

    ! Isotropic volume change
    IF( RANDOM_N < PROB_ISO_INIT ) THEN
      ! Random scaling factor
      CALL RANF(  )
      SCALE_FACTOR = 1.D0 + DVMAXISO * (RANDOM_N - 0.5D0)
      ! Proportional box length
      BOXLN(:) = BOXLM(:) * SCALE_FACTOR
      CALL INVERSE_COF( BOXLN, BOXLN_I, BOXVN )
      ! Movement counter
      MOVVI = MOVVI + 1
      ! Movement type
      MOV_VOL_I = .TRUE.  ! Enable isotropic volume change
      MOV_VOL_A = .FALSE. ! Disable anisotropic volume change
    ! Anisotropic volume change
    ELSE IF( RANDOM_N >= PROB_ISO_INIT ) THEN
      ! Random box component
      CALL RANF(  )
      COMPONENT = INT( RANDOM_N * 6.D0 ) + 1
      IF( COMPONENT == 1 ) THEN
        COMPONENT = 1 ! x vector
      ELSE IF( COMPONENT == 2 ) THEN
        COMPONENT = 4 ! y vector
      ELSE IF( COMPONENT == 3 ) THEN
        COMPONENT = 5 ! y vector
      ELSE IF( COMPONENT == 4 ) THEN
        COMPONENT = 7 ! z vector
      ELSE IF( COMPONENT == 5 ) THEN
        COMPONENT = 8 ! z vector
      ELSE IF( COMPONENT == 6 ) THEN
        COMPONENT = 9 ! z vector
      END IF
      BOXLN(:) = BOXLM(:)
      ! Random factor
      CALL RANF(  )
      BOXLN(COMPONENT) = BOXLM(COMPONENT) + DVMAXANI * (RANDOM_N - 0.5D0)
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL INVERSE_COF( BOXLN, BOXLN_I, BOXVN )
      ! Movement counter
      MOVVA = MOVVA + 1
      ! Movement type
      MOV_VOL_I = .FALSE. ! Disable isotropic volume change
      MOV_VOL_A = .TRUE.  ! Enable anisotropic volume change
    END IF

    ! Reset condition of anisotropic volume change
    IGNORE = .FALSE.

    ! Condition of anisotropic volume change (box distortion)
    IF( MOV_VOL_A ) THEN
      ! Box length
      LBOX(1) = DSQRT( DOT_PRODUCT( BOXLN(1:3), BOXLN(1:3) ) )
      LBOX(2) = DSQRT( DOT_PRODUCT( BOXLN(4:6), BOXLN(4:6) ) )
      LBOX(3) = DSQRT( DOT_PRODUCT( BOXLN(7:9), BOXLN(7:9) ) )
      ! Length ratio
      LBOXR(1) = LBOX(1) / LBOX(2)
      LBOXR(2) = LBOX(1) / LBOX(3)
      LBOXR(3) = LBOX(2) / LBOX(3)
      ! Angle between box vectors
      COSANGLE_VEC(1) = DOT_PRODUCT( BOXLN(1:3), BOXLN(4:6) ) / ( LBOX(1) * LBOX(2) )
      COSANGLE_VEC(2) = DOT_PRODUCT( BOXLN(1:3), BOXLN(7:9) ) / ( LBOX(1) * LBOX(3) )
      COSANGLE_VEC(3) = DOT_PRODUCT( BOXLN(4:6), BOXLN(7:9) ) / ( LBOX(2) * LBOX(3) )
      ! Avoid big distortions of the simulation box
      DO L = 1, 3
        ! Angle distortion
        IF( COSANGLE_VEC(L) < DCOS( (PI / 2.D0) + MAX_ANGLE ) .OR. &
        &   COSANGLE_VEC(L) > DCOS( (PI / 2.D0) - MAX_ANGLE ) ) THEN
          BOXVMC_RND  = BOXVM
          BOXLMC(:)   = BOXLM(:)
          BOXLMC_I(:) = BOXLM_I(:)
          IGNORE = .TRUE.
          EXIT
        END IF
        ! Length distortion
        IF( LBOXR(L) > MAX_LENGTH_RATIO .OR. LBOXR(L) < 1.D0 / MAX_LENGTH_RATIO ) THEN
          BOXVMC_RND  = BOXVM
          BOXLMC(:)   = BOXLM(:)
          BOXLMC_I(:) = BOXLM_I(:)
          IGNORE = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    ! Box not too distorted
    IF( .NOT. IGNORE ) THEN

      ! Enthalpy (weighing function)
      HNM = ( PRESS_RND * ( BOXVN - BOXVM ) ) - ( DBLE( N_PARTICLES ) * DLOG( BOXVN / BOXVM ) )

      ! Random number
      CALL RANF(  )

      ! Enthalpy Criterion
      IF( DEXP( - HNM ) >= RANDOM_N ) THEN

        ! System configuration
        RMCV(:,:) = RMC(:,:) ! Old configuration

        ! Isotropic volume change
        IF( MOV_VOL_I ) THEN
          ! Rescale positions of particles accordingly
          DO K = 1, N_PARTICLES
            RMC(:,K) = RMC(:,K) * SCALE_FACTOR
          END DO
        ! Anisotropic volume change
        ELSE IF( MOV_VOL_A ) THEN
          ! Rescale positions of particles accordingly
          DO K = 1, N_PARTICLES
            ! Scaled coordinates using old dimensions
            CALL MULTI_MATRIX( BOXLM_I, RMC(:,K), S12 )
            ! New real coordinates using new dimensions
            CALL MULTI_MATRIX( BOXLN, S12, RMC(:,K) )
          END DO
        END IF

        ! Overlap check after expansion/compression of the simulation box
        LOOP_OVERLAP_NPT: DO

          ! Initialization
          OVERLAP = .FALSE.

          ! Anisomorphic molecules (unlike components)
          DO CI = 1, COMPONENTS - 1
            DO CJ = CI + 1, COMPONENTS
              ! First loop represents all particles with indexes i of component Ci
              DO I = SUM( N_COMPONENT(0:(CI-1)) ) + 1, SUM( N_COMPONENT(0:CI) )
                ! Second loop represents all particles with indexes j of component Cj
                DO J = SUM( N_COMPONENT(0:(CJ-1)) ) + 1, SUM( N_COMPONENT(0:CJ) )
                  ! Position of particle i
                  RI(1)  = RMC(1,I)
                  RI(2)  = RMC(2,I)
                  RI(3)  = RMC(3,I)
                  ! Position of particle j
                  RJ(1)  = RMC(1,J)
                  RJ(2)  = RMC(2,J)
                  RJ(3)  = RMC(3,J)
                  ! Orientation of particle i
                  EI(1)  = EMC(1,I)
                  EI(2)  = EMC(2,I)
                  EI(3)  = EMC(3,I)
                  ! Orientation of particle j
                  EJ(1)  = EMC(1,J)
                  EJ(2)  = EMC(2,J)
                  EJ(3)  = EMC(3,J)
                  ! Quaternion of particle i
                  QI(0)  = QMC(0,I)
                  QI(1)  = QMC(1,I)
                  QI(2)  = QMC(2,I)
                  QI(3)  = QMC(3,I)
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
                  CALL MULTI_MATRIX( BOXLN_I, RIJ, S12 )
                  S12 = S12 - ANINT( S12 )
                  CALL MULTI_MATRIX( BOXLN, S12, RIJ )
                  ! Magnitude of the vector distance (squared)
                  RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
                  ! Cutoff distance (squared)
                  CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
                  CUTOFF_D = CUTOFF_D * CUTOFF_D
                  ! Preliminary test (circumscribing spheres)
                  IF( RIJSQ <= CUTOFF_D ) THEN
                    ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
                    IF( GEOM_SELEC(1) ) THEN
                      CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
                      ! Overlap criterion
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_OVERLAP_NPT
                      END IF
                    ! Overlap test for spherocylinders (Vega-Lago Method)
                    ELSE IF( GEOM_SELEC(2) ) THEN
                      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                      ! Overlap criterion
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_OVERLAP_NPT
                      END IF
                    ! Overlap test for cylinders (Lopes et al. Method)
                    ELSE IF( GEOM_SELEC(3) ) THEN
                      ! Preliminary test (circumscribing spherocylinders)
                      OVERLAP_PRELIMINAR = .FALSE.
                      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                      ! Overlap criterion
                      IF( OVERLAP_PRELIMINAR ) THEN
                        RJ(1) = RI(1) + RIJ(1)
                        RJ(2) = RI(2) + RIJ(2)
                        RJ(3) = RI(3) + RIJ(3)
                        CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                        ! Overlap criterion
                        IF( OVERLAP ) THEN
                          ! Overlap detected
                          EXIT LOOP_OVERLAP_NPT
                        END IF
                      END IF
                    END IF
                  END IF
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
                RI(1)  = RMC(1,I)
                RI(2)  = RMC(2,I)
                RI(3)  = RMC(3,I)
                ! Position of particle j
                RJ(1)  = RMC(1,J)
                RJ(2)  = RMC(2,J)
                RJ(3)  = RMC(3,J)
                ! Orientation of particle i
                EI(1)  = EMC(1,I)
                EI(2)  = EMC(2,I)
                EI(3)  = EMC(3,I)
                ! Orientation of particle j
                EJ(1)  = EMC(1,J)
                EJ(2)  = EMC(2,J)
                EJ(3)  = EMC(3,J)
                ! Quaternion of particle i
                QI(0)  = QMC(0,I)
                QI(1)  = QMC(1,I)
                QI(2)  = QMC(2,I)
                QI(3)  = QMC(3,I)
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
                CALL MULTI_MATRIX( BOXLN_I, RIJ, S12 )
                S12 = S12 - ANINT( S12 )
                CALL MULTI_MATRIX( BOXLN, S12, RIJ )
                ! Magnitude of the vector distance (squared)
                RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
                ! Cutoff distance (squared)
                CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
                CUTOFF_D = CUTOFF_D * CUTOFF_D
                ! Preliminary test (circumscribing spheres)
                IF( RIJSQ <= CUTOFF_D ) THEN
                  ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
                  IF( GEOM_SELEC(1) ) THEN
                    CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
                    ! Overlap criterion
                    IF( OVERLAP ) THEN
                      ! Overlap detected
                      EXIT LOOP_OVERLAP_NPT
                    END IF
                  ! Overlap test for spherocylinders (Vega-Lago Method)
                  ELSE IF( GEOM_SELEC(2) ) THEN
                    CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                    ! Overlap criterion
                    IF( OVERLAP ) THEN
                      ! Overlap detected
                      EXIT LOOP_OVERLAP_NPT
                    END IF
                  ! Overlap test for cylinders (Lopes et al. Method)
                  ELSE IF( GEOM_SELEC(3) ) THEN
                    ! Preliminary test (circumscribing spherocylinders)
                    OVERLAP_PRELIMINAR = .FALSE.
                    CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                    ! Overlap criterion
                    IF( OVERLAP_PRELIMINAR ) THEN
                      RJ(1) = RI(1) + RIJ(1)
                      RJ(2) = RI(2) + RIJ(2)
                      RJ(3) = RI(3) + RIJ(3)
                      CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                      ! Overlap criterion
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_OVERLAP_NPT
                      END IF
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END DO

          ! No overlaps
          OVERLAP = .FALSE.
          EXIT LOOP_OVERLAP_NPT

        END DO LOOP_OVERLAP_NPT

        ! Acceptance Criterion
        IF( .NOT. OVERLAP ) THEN
          ! Assigns the simulation box properties of a trial volume change to the system configuration.
          BOXVMC_RND  = BOXVN      ! Update volume
          BOXLMC(:)   = BOXLN(:)   ! Update length
          BOXLMC_I(:) = BOXLN_I(:) ! Update length (inverse)
          ! Displacement counter update
          IF( MOV_VOL_I ) THEN
            NACCVI = NACCVI + 1 ! Isotropic move counter
          ELSE IF( MOV_VOL_A ) THEN
            NACCVA = NACCVA + 1 ! Anisotropic move counter
          END IF
          ! Update packing fraction
          ETA_NPT = TOTAL_VP / BOXVN
          ! Re-initialization
          IGNORE = .FALSE.
          ! Lattice reduction
          LATTICER = .FALSE.
          CALL LATTICE_REDUCTION( BOXLMC, DISTORTION, LATTICER )
          IF( LATTICER ) THEN
            ! Calculate the new reciprocal box basis vectors
            CALL INVERSE_COF( BOXLMC, BOXLMC_I, BOXVMC_RND )
            DO K = 1, N_PARTICLES
              ! Minimum Image Convention
              CALL MULTI_MATRIX( BOXLMC_I, RMC(:,K), S12 )
              S12 = S12 - ANINT( S12 )
              CALL MULTI_MATRIX( BOXLMC, S12, RMC(:,K) )
            END DO
            ! Undo box rotation
            IF( DABS( BOXLMC(2) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. DABS( BOXLMC(3) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. &
            &   DABS( BOXLMC(6) - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
              ! Initialization
              BOXLM   = BOXLMC
              BOXLM_I = BOXLMC_I
              ! Box vectors
              V1 = BOXLMC(1:3)
              V2 = BOXLMC(4:6)
              V3 = BOXLMC(7:9)
              ! Angle between x-vector and x-axis
              THETA = DACOS( DOT_PRODUCT( V1, [1.D0,0.D0,0.D0] ) / DSQRT( DOT_PRODUCT( V1, V1 ) ) )
              ! Cross product between x-vector and x-axis (rotation axis)
              RAXIS(1) = 0.D0
              RAXIS(2) = V1(3)
              RAXIS(3) = - V1(2)
              ! Magnitude of rotation axis
              RAXISMAG = DSQRT( DOT_PRODUCT( RAXIS, RAXIS ) )
              ! Avoid null vectors
              IF( DABS( RAXISMAG - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
                RAXIS(:) = 0.D0
              ELSE
                RAXIS(:) = RAXIS(:) / RAXISMAG
              END IF
              ! Rotation quaternion
              QROT(0) = DCOS( THETA * 0.5D0 )            ! Real part
              QROT(1) = DSIN( THETA * 0.5D0 ) * RAXIS(1) ! Imaginary part (Vector)
              QROT(2) = DSIN( THETA * 0.5D0 ) * RAXIS(2) ! Imaginary part (Vector)
              QROT(3) = DSIN( THETA * 0.5D0 ) * RAXIS(3) ! Imaginary part (Vector)
              ! Make box x-vector parallel to x-axis of coordination system
              IF( DABS( RAXISMAG - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V1 / ( DSQRT( DOT_PRODUCT( V1, V1 ) ) ), QROT, AUXV )
                ! New x-vector of the simulation box
                BOXLROT(1:3) = DSQRT( DOT_PRODUCT( BOXLMC(1:3), BOXLMC(1:3) ) ) * AUXV
                BOXLROT(2:3) = 0.D0
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V2 / ( DSQRT( DOT_PRODUCT( V2, V2 ) ) ), QROT, AUXV )
                ! New y-vector of the simulation box
                BOXLROT(4:6) = DSQRT( DOT_PRODUCT( BOXLMC(4:6), BOXLMC(4:6) ) ) * AUXV
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V3 / ( DSQRT( DOT_PRODUCT( V3, V3 ) ) ), QROT, AUXV )
                ! New z-vector of the simulation box
                BOXLROT(7:9) = DSQRT( DOT_PRODUCT( BOXLMC(7:9), BOXLMC(7:9) ) ) * AUXV
                ! Calculate the new reciprocal box basis vectors
                CALL INVERSE_COF( BOXLROT, BOXIROT, BOXVROT )
                ! Rescale positions and orientations of particles accordingly
                DO K = 1, N_PARTICLES
                  ! Scaled coordinates using old dimensions
                  CALL MULTI_MATRIX( BOXLM_I, RMC(:,K), S12 )
                  ! New real coordinates using new dimensions
                  CALL MULTI_MATRIX( BOXLROT, S12, RPROT(:,K) )
                  ! Reorient particles
                  QAUX(:) = QMC(:,K)
                  CALL MULTIPLY_QUATERNIONS( QROT, QAUX, QPROT(:,K) )
                  ! Active transformation (rotation)
                  CALL ACTIVE_TRANSFORMATION( AXISZ, QPROT(:,K), EPROT(:,K) )
                END DO
              ! Box x-vector already parallel to x-axis of coordination system
              ELSE
                ! Retrive old box properties
                BOXLROT(:) = BOXLMC(:)
                BOXIROT(:) = BOXLMC_I(:)
                BOXVROT    = BOXVMC_RND
                ! Retrive old molecular properties
                RPROT(:,:) = RMC(:,:)
                QPROT(:,:) = QMC(:,:)
                EPROT(:,:) = EMC(:,:)
              END IF
              ! Initialization
              BOXLM   = BOXLROT
              BOXLM_I = BOXIROT
              ! Box vectors
              V1 = BOXLROT(1:3)
              V2 = BOXLROT(4:6)
              V3 = BOXLROT(7:9)
              ! Axis of rotation
              RAXIS = V1
              ! Magnitude of rotation axis
              RAXISMAG = DSQRT( DOT_PRODUCT( RAXIS, RAXIS ) )
              ! Avoid null vectors
              IF( DABS( RAXISMAG - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
                RAXIS(:) = 0.D0
              ELSE
                RAXIS(:) = RAXIS(:) / RAXISMAG
              END IF
              ! Projection of the y-vector of the box onto the ZY-plane
              PROJY_XY = V2 - ( DOT_PRODUCT( V2, RAXIS ) ) * RAXIS
              ! Versor of the projection of the y-vector of the box onto the ZY-plane
              PROJY_XY = PROJY_XY / DSQRT( DOT_PRODUCT( PROJY_XY, PROJY_XY ) )
              ! Angle between the projection of the y-vector of the box and the y-axis of the coordination system
              THETA = DACOS( DOT_PRODUCT( PROJY_XY, [0.D0,1.D0,0.D0] ) / DSQRT( DOT_PRODUCT( PROJY_XY, PROJY_XY ) ) )
              ! Direction of rotation
              IF( PROJY_XY(3) < 0.D0 ) THEN
                ! Rotation quaternion (clockwise rotation)
                QROT(0) = DCOS( THETA * 0.5D0 )              ! Real part
                QROT(1) = DSIN( THETA * 0.5D0 ) * RAXIS(1)   ! Imaginary part (Vector)
                QROT(2) = DSIN( THETA * 0.5D0 ) * RAXIS(2)   ! Imaginary part (Vector)
                QROT(3) = DSIN( THETA * 0.5D0 ) * RAXIS(3)   ! Imaginary part (Vector)
              ELSE
                ! Rotation quaternion (counterclockwise rotation)
                QROT(0) = DCOS( - THETA * 0.5D0 )            ! Real part
                QROT(1) = DSIN( - THETA * 0.5D0 ) * RAXIS(1) ! Imaginary part (Vector)
                QROT(2) = DSIN( - THETA * 0.5D0 ) * RAXIS(2) ! Imaginary part (Vector)
                QROT(3) = DSIN( - THETA * 0.5D0 ) * RAXIS(3) ! Imaginary part (Vector)
              END IF
              ! Make box y-vector coplanar with the XY-plane of the coordination system
              IF( DABS( RAXISMAG - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V1 / ( DSQRT( DOT_PRODUCT( V1, V1 ) ) ), QROT, AUXV )
                ! New x-vector of the simulation box
                BOXLROT(1:3) = DSQRT( DOT_PRODUCT( BOXLROT(1:3), BOXLROT(1:3) ) ) * AUXV
                BOXLROT(2:3) = 0.D0
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V2 / ( DSQRT( DOT_PRODUCT( V2, V2 ) ) ), QROT, AUXV )
                ! New y-vector of the simulation box
                BOXLROT(4:6) = DSQRT( DOT_PRODUCT( BOXLROT(4:6), BOXLROT(4:6) ) ) * AUXV
                BOXLROT(6)   = 0.D0
                ! Auxiliary vector
                CALL ACTIVE_TRANSFORMATION( V3 / ( DSQRT( DOT_PRODUCT( V3, V3 ) ) ), QROT, AUXV )
                ! New z-vector of the simulation box
                BOXLROT(7:9) = DSQRT( DOT_PRODUCT( BOXLROT(7:9), BOXLROT(7:9) ) ) * AUXV
                ! Calculate the new reciprocal box basis vectors
                CALL INVERSE_COF( BOXLROT, BOXIROT, BOXVROT )
                ! Rescale positions and orientations of particles accordingly
                DO K = 1, N_PARTICLES
                  ! Scaled coordinates using old dimensions
                  CALL MULTI_MATRIX( BOXLM_I, RPROT(:,K), S12 )
                  ! New real coordinates using new dimensions
                  CALL MULTI_MATRIX( BOXLROT, S12, RPROT(:,K) )
                  ! Reorient particles
                  QAUX(:) = QPROT(:,K)
                  CALL MULTIPLY_QUATERNIONS( QROT, QAUX, QPROT(:,K) )
                  ! Active transformation (rotation)
                  CALL ACTIVE_TRANSFORMATION( AXISZ, QPROT(:,K), EPROT(:,K) )
                END DO
              END IF
              ! Update box properties
              BOXLMC(:)   = BOXLROT(:)
              BOXLMC_I(:) = BOXIROT(:)
              BOXVMC_RND  = BOXVROT
              ! Update molecular properties
              RMC(:,:) = RPROT(:,:)
              QMC(:,:) = QPROT(:,:)
              EMC(:,:) = EPROT(:,:)
            END IF
          END IF
        ! Retrieve old properties of the system configuration and the simulation box
        ELSE IF( OVERLAP ) THEN
          BOXVMC_RND  = BOXVM      ! Retrieve box volume
          BOXLMC(:)   = BOXLM(:)   ! Retrieve box length
          BOXLMC_I(:) = BOXLM_I(:) ! Retrieve box length (inverse)
          RMC(:,:)    = RMCV(:,:)  ! Retrieve position of particles
        END IF

      ! Retrieve old properties of the simulation box
      ELSE

        BOXVMC_RND  = BOXVM      ! Retrieve box volume
        BOXLMC(:)   = BOXLM(:)   ! Retrieve box length
        BOXLMC_I(:) = BOXLM_I(:) ! Retrieve box length (inverse)

      END IF ! Enthalpy criterion

    END IF ! Box distortion criterion

  END IF

  ! Iteration
  ATTEMPTS = ATTEMPTS + 1

  ! Adjustment of maximum displacement (Translation, Rotation, and Volume Change)
  IF( MOD( ATTEMPTS, N_ADJUST_INIT ) == 0 ) THEN

    ! Translational adjustment
    IF( MOVT > 200 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      RATIO = DBLE( NACCT ) / DBLE( MOVT )
      ! Translational adjustment
      IF( RATIO <= R_ACC_T ) THEN
        DRMAX  = 0.95D0 * DRMAX
      ELSE
        DRMAX  = 1.05D0 * DRMAX
      END IF
      ! Reset counter
      NACCT = 0
      MOVT  = 0
    END IF

    ! Avoid multiple turns
    IF( DRMAX >= 2.D0 * MAXVAL( BOXLMC ) ) THEN
      DRMAX = DRMAX - MAXVAL( BOXLMC )
    END IF

    ! Rotational adjustment
    IF( MOVR > 200 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      RATIO = DBLE( NACCR ) / DBLE( MOVR )
      ! Rotation adjustment
      IF( RATIO <= R_ACC_R ) THEN
        ANGMAX = 0.95D0 * ANGMAX
      ELSE
        ANGMAX = 1.05D0 * ANGMAX
      END IF
      ! Reset counter
      NACCR = 0
      MOVR  = 0
    END IF

    ! Avoid multiple turns
    IF( ANGMAX >= 4.D0 * PI ) THEN
      ANGMAX = ANGMAX - 2.D0 * PI
    END IF

    ! Volumetric adjustment (isotropic)
    IF( MOVVI >= 25 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      RATIO = DBLE( NACCVI ) / DBLE( MOVVI )
      ! Volumetric adjustment
      IF( RATIO <= R_ACC_VI ) THEN
        DVMAXISO = 0.95D0 * DVMAXISO
      ELSE
        DVMAXISO = 1.05D0 * DVMAXISO
      END IF
      ! Reset counter
      NACCVI = 0
      MOVVI  = 0
    END IF

    ! Volumetric adjustment (anisotropic)
    IF( MOVVA >= 25 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      RATIO = DBLE( NACCVA ) / DBLE( MOVVA )
      ! Volumetric adjustment
      IF( RATIO <= R_ACC_VA ) THEN
        DVMAXANI = 0.95D0 * DVMAXANI
      ELSE
        DVMAXANI = 1.05D0 * DVMAXANI
      END IF
      ! Reset counter
      NACCVA = 0
      MOVVA  = 0
    END IF

    ! Avoid low volume changes
    IF( DVMAXISO <= DVMIN_INIT ) THEN
      DVMAXISO = DVMIN_INIT
    END IF
    IF( DVMAXANI <= DVMIN_INIT ) THEN
      DVMAXANI = DVMIN_INIT
    END IF

  END IF

  ! Summary
  CALL PROGRESS_BAR_NPT( ATTEMPTS, ETA_NPT, PACKING_F )

  ! Target packing fraction
  IF( ETA_NPT >= PACKING_F .AND. ATTEMPTS == 1 ) THEN
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(3G0,G0.5)" ) "Target packing fraction reached after ", ATTEMPTS, " attempt! Final value: ", ETA_NPT
    WRITE( *, "(G0)" ) " "
    EXIT NPT_SIMULATION
  ELSE IF( ETA_NPT >= PACKING_F .AND. ATTEMPTS > 1 ) THEN
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(3G0,G0.5)" ) "Target packing fraction reached after ", ATTEMPTS, " attempts! Final value: ", ETA_NPT
    WRITE( *, "(G0)" ) " "
    EXIT NPT_SIMULATION
  END IF

  ! Initial configuration (partial)
  OPEN( UNIT= 55, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/initconf_rnd_"// &
  &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  WRITE( 55, "(I4)" ) N_PARTICLES
  WRITE( 55, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 55, * ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
        &              0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 55, * ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
        &              0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 55 )

END DO NPT_SIMULATION

CALL SLEEP( 3 )

! *********************************************************************************************** !
! Monte Carlo parameters (NVT Simulation | Fixing the packing fraction)                           !
! *********************************************************************************************** !
MOV_TRANS         = .FALSE.           ! Translational move selector           (initial value)
MOV_ROT           = .FALSE.           ! Rotational move selector              (initial value)
DRMAX             = DRMAX_INIT        ! Maximum translational displacement    (initial value)
ANGMAX            = ANGMAX_INIT       ! Maximum rotational displacement       (initial value)
NACCT             = 0                 ! Translational move acceptance counter (initial value)
NACCR             = 0                 ! Rotational move acceptance counter    (initial value)
MOVT              = 0                 ! Translational move counter            (initial value)
MOVR              = 0                 ! Rotational move counter               (initial value)
ATTEMPTS          = 0                 ! Number of attempts                    (initial value)

! Scale factor
SCALE_FACTOR = ETA_NPT / PACKING_F
SCALE_FACTOR = SCALE_FACTOR ** ( 1.D0 / 3.D0 )

! Scaled box length
BOX_LENGTH(:) = BOXLMC(:) * SCALE_FACTOR

! Inverse of box length
CALL INVERSE_COF( BOX_LENGTH, BOX_LENGTH_I, BOXVMC )

! Fix packing fraction with a volume expansion
IF( ETA_NPT > PACKING_F ) THEN

  ! Summary
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Attempting to fix the packing fraction of ", ETA_NPT, &
  &                                  " obtained in the NPT simulation to the target value of ", PACKING_F, "..."
  WRITE( *, "(G0)" ) " "
  CALL SLEEP( 2 )

  LOOP_PFRACTION_FIX: DO

    ! Update rotation quaternions and orientations
    Q(:,:) = QMC(:,:)
    E(:,:) = EMC(:,:)

    ! Rescale positions of particles accordingly
    DO K = 1, N_PARTICLES
      ! Scaling coordinates using the old box length
      CALL MULTI_MATRIX( BOXLMC_I, RMC(:,K), S12 )
      ! New real coordinates using the new box length
      CALL MULTI_MATRIX( BOX_LENGTH, S12, R(:,K) )
    END DO

    ! Overlap check after expansion of the simulation box
    LOOP_OVERLAP_NPT_FIX: DO

      ! Initialization
      OVERLAP = .FALSE.

      ! Iteration
      ATTEMPTS = ATTEMPTS + 1

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
              ! Orientation of particle i
              EI(1)  = E(1,I)
              EI(2)  = E(2,I)
              EI(3)  = E(3,I)
              ! Orientation of particle j
              EJ(1)  = E(1,J)
              EJ(2)  = E(2,J)
              EJ(3)  = E(3,J)
              ! Quaternion of particle i
              QI(0)  = Q(0,I)
              QI(1)  = Q(1,I)
              QI(2)  = Q(2,I)
              QI(3)  = Q(3,I)
              ! Quaternion of particle j
              QJ(0)  = Q(0,J)
              QJ(1)  = Q(1,J)
              QJ(2)  = Q(2,J)
              QJ(3)  = Q(3,J)
              ! Vector distance between particles i and j
              RIJ(1) = RJ(1) - RI(1)
              RIJ(2) = RJ(2) - RI(2)
              RIJ(3) = RJ(3) - RI(3)
              ! Minimum Image Convention
              CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
              S12 = S12 - ANINT( S12 )
              CALL MULTI_MATRIX( BOX_LENGTH, S12, RIJ )
              ! Magnitude of the vector distance (squared)
              RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
              ! Cutoff distance (squared)
              CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
              CUTOFF_D = CUTOFF_D * CUTOFF_D
              ! Preliminary test (circumscribing spheres)
              IF( RIJSQ <= CUTOFF_D ) THEN
                ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
                IF( GEOM_SELEC(1) ) THEN
                  CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
                  ! Overlap criterion
                  IF( OVERLAP ) THEN
                    ! Overlap detected
                    EXIT LOOP_OVERLAP_NPT_FIX
                  END IF
                ! Overlap test for spherocylinders (Vega-Lago Method)
                ELSE IF( GEOM_SELEC(2) ) THEN
                  CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                  ! Overlap criterion
                  IF( OVERLAP ) THEN
                    ! Overlap detected
                    EXIT LOOP_OVERLAP_NPT_FIX
                  END IF
                ! Overlap test for cylinders (Lopes et al. Method)
                ELSE IF( GEOM_SELEC(3) ) THEN
                  ! Preliminary test (circumscribing spherocylinders)
                  OVERLAP_PRELIMINAR = .FALSE.
                  CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                  ! Overlap criterion
                  IF( OVERLAP_PRELIMINAR ) THEN
                    RJ(1) = RI(1) + RIJ(1)
                    RJ(2) = RI(2) + RIJ(2)
                    RJ(3) = RI(3) + RIJ(3)
                    CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                    ! Overlap criterion
                    IF( OVERLAP ) THEN
                      ! Overlap detected
                      EXIT LOOP_OVERLAP_NPT_FIX
                    END IF
                  END IF
                END IF
              END IF
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
            ! Orientation of particle i
            EI(1)  = E(1,I)
            EI(2)  = E(2,I)
            EI(3)  = E(3,I)
            ! Orientation of particle j
            EJ(1)  = E(1,J)
            EJ(2)  = E(2,J)
            EJ(3)  = E(3,J)
            ! Quaternion of particle i
            QI(0)  = Q(0,I)
            QI(1)  = Q(1,I)
            QI(2)  = Q(2,I)
            QI(3)  = Q(3,I)
            ! Quaternion of particle j
            QJ(0)  = Q(0,J)
            QJ(1)  = Q(1,J)
            QJ(2)  = Q(2,J)
            QJ(3)  = Q(3,J)
            ! Vector distance between particles i and j
            RIJ(1) = RJ(1) - RI(1)
            RIJ(2) = RJ(2) - RI(2)
            RIJ(3) = RJ(3) - RI(3)
            ! Minimum Image Convention
            CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
            S12 = S12 - ANINT( S12 )
            CALL MULTI_MATRIX( BOX_LENGTH, S12, RIJ )
            ! Magnitude of the vector distance (squared)
            RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
            ! Cutoff distance (squared)
            CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
            CUTOFF_D = CUTOFF_D * CUTOFF_D
            ! Preliminary test (circumscribing spheres)
            IF( RIJSQ <= CUTOFF_D ) THEN
              ! Overlap test for ellipsoids of revolution (Perram-Wertheim Method)
              IF( GEOM_SELEC(1) ) THEN
                CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
                ! Overlap criterion
                IF( OVERLAP ) THEN
                  ! Overlap detected
                  EXIT LOOP_OVERLAP_NPT_FIX
                END IF
              ! Overlap test for spherocylinders (Vega-Lago Method)
              ELSE IF( GEOM_SELEC(2) ) THEN
                CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                ! Overlap criterion
                IF( OVERLAP ) THEN
                  ! Overlap detected
                  EXIT LOOP_OVERLAP_NPT_FIX
                END IF
              ! Overlap test for cylinders (Lopes et al. Method)
              ELSE IF( GEOM_SELEC(3) ) THEN
                ! Preliminary test (circumscribing spherocylinders)
                OVERLAP_PRELIMINAR = .FALSE.
                CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                ! Overlap criterion
                IF( OVERLAP_PRELIMINAR ) THEN
                  RJ(1) = RI(1) + RIJ(1)
                  RJ(2) = RI(2) + RIJ(2)
                  RJ(3) = RI(3) + RIJ(3)
                  CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                  ! Overlap criterion
                  IF( OVERLAP ) THEN
                    ! Overlap detected
                    EXIT LOOP_OVERLAP_NPT_FIX
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO

      ! No overlaps
      OVERLAP = .FALSE.
      EXIT LOOP_OVERLAP_NPT_FIX

    END DO LOOP_OVERLAP_NPT_FIX

    ! Packing fraction fixed
    IF( .NOT. OVERLAP ) THEN

      EXIT LOOP_PFRACTION_FIX

    ! Attempt fixing the packing fraction with a new configuration of particles
    ELSE IF( OVERLAP ) THEN

      ! Displace particles (constant volume)
      DO CYCLES = 1, N_PARTICLES

        ! Component index
        IF( COMPONENTS > 1 ) THEN
          ! Pseudorandom number generator (uniform distribution)
          CALL RANF(  )
          ! Get component index
          CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
        ELSE
          ! Get component index
          CI = 1
        END IF

        ! Avoid components with molar fraction of 0
        DO WHILE( N_COMPONENT(CI) < 1 )
          ! Component index
          IF( COMPONENTS > 1 ) THEN
            ! Pseudorandom number generator (uniform distribution)
            CALL RANF(  )
            ! Get component index
            CI = INT( RANDOM_N * DBLE( COMPONENTS ) ) + 1
          ELSE
            ! Get component index
            CI = 1
          END IF
        END DO

        ! Forbid rotation if component is spherical
        IF( GEOM_SELEC(1) .AND. DABS( ASPECT_RATIO(CI) - 1.D0 ) < EPSILON( 1.D0 ) ) THEN
          MOV_TRANS = .TRUE.   ! Enable translation
          MOV_ROT   = .FALSE.  ! Disable rotation
          MOVT      = MOVT + 1 ! Increment move counter
        ELSE IF( GEOM_SELEC(2) .AND. DABS( ASPECT_RATIO(CI) - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
          MOV_TRANS = .TRUE.   ! Enable translation
          MOV_ROT   = .FALSE.  ! Disable rotation
          MOVT      = MOVT + 1 ! Increment move counter
        ! Allow rotation if component is nonspherical
        ELSE
          ! Pseudorandom number generator (uniform distribution)
          CALL RANF(  )
          ! Translation criterion
          IF( RANDOM_N < PROB_TRANS_INIT ) THEN
            MOV_TRANS = .TRUE.   ! Enable translation
            MOV_ROT   = .FALSE.  ! Disable rotation
            MOVT      = MOVT + 1 ! Increment move counter
          ! Rotation criterion
          ELSE IF( RANDOM_N >= PROB_TRANS_INIT ) THEN
            MOV_ROT   = .TRUE.   ! Enable rotation
            MOV_TRANS = .FALSE.  ! Disable translation
            MOVR      = MOVR + 1 ! Increment move counter
          END IF
        END IF

        ! Pseudorandom number generator (uniform distribution)
        CALL RANF(  )
        ! Random selection of particles of component C
        I = SUM( N_COMPONENT(0:(CI-1)) ) + INT( RANDOM_N * DBLE( N_COMPONENT(CI) ) ) + 1

        ! Assignment of previous configuration (Microstate m)
        RM(:) = RMC(:,I) ! Position
        QM(:) = QMC(:,I) ! Quaternion
        EM(:) = EMC(:,I) ! Orientation

        ! Translation Movement
        IF( MOV_TRANS ) THEN
          ! Random translation along x-axis
          CALL RANF(  )
          RN(1) = RM(1) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
          ! Random translation along y-axis
          CALL RANF(  )
          RN(2) = RM(2) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
          ! Random translation along z-axis
          CALL RANF(  )
          RN(3) = RM(3) + ( ( 2.D0 * RANDOM_N ) - 1.D0 ) * DRMAX  ! Range [-drmax,drmax]
          ! Minimum Image Convention
          CALL MULTI_MATRIX( BOXLMC_I, RN, S12 )
          S12 = S12 - ANINT( S12 )
          CALL MULTI_MATRIX( BOXLMC, S12, RN )
        ! No Translation
        ELSE IF( .NOT. MOV_TRANS ) THEN
          RN(:) = RM(:)
        END IF

        ! Rotation Movement
        IF( MOV_ROT ) THEN
          ! Random Composed Unit Quaternion
          CALL COMPOSED_QUATERNION( QM, QN, ANGMAX )
          ! Active transformation
          CALL ACTIVE_TRANSFORMATION( AXISZ, QN, EN )
        ! No Rotation
        ELSE IF( .NOT. MOV_ROT ) THEN
          QN(:) = QM(:)
          EN(:) = EM(:)
        END IF

        ! Overlap Check
        CALL CHECK_OVERLAP( CI, I, QN, EN, RN, CD, BOXLMC, BOXLMC_I, OVERLAP )

        ! Acceptance Criterion
        IF( .NOT. OVERLAP ) THEN
          ! System configuration update
          RMC(:,I) = RN(:) ! Update position
          QMC(:,I) = QN(:) ! Update quaternion
          EMC(:,I) = EN(:) ! Update orientation
          ! Displacement counter update
          IF( MOV_TRANS ) THEN
            NACCT = NACCT + 1  ! Translational move counter
          ELSE IF ( MOV_ROT ) THEN
            NACCR = NACCR + 1  ! Rotational move counter
          END IF
        ELSE
          ! Retrieve old configuration
          RMC(:,I) = RM(:) ! Retrieve position
          QMC(:,I) = QM(:) ! Retrieve quaternion
          EMC(:,I) = EM(:) ! Retrieve orientation
        END IF

      END DO

      ! Adjustment of maximum displacement (Translation and Rotation)
      IF( MOD( (MOVT + MOVR), (N_ADJUST_INIT * N_PARTICLES) ) == 0 ) THEN

        ! Acceptance ratio (translation)
        RATIO = DBLE( NACCT ) / DBLE( MOVT )
        ! Translational adjustment
        IF( RATIO <= R_ACC_T ) THEN
          DRMAX = 0.95D0 * DRMAX
        ELSE
          DRMAX = 1.05D0 * DRMAX
        END IF

        ! Acceptance ratio (rotation)
        RATIO = DBLE( NACCR ) / DBLE( MOVR )
        ! Rotational adjustment
        IF( RATIO <= R_ACC_R ) THEN
          ANGMAX = 0.95D0 * ANGMAX
        ELSE
          ANGMAX = 1.05D0 * ANGMAX
        END IF

        ! 4π-rotation condition
        IF( ( ANGMAX > 4.D0 * PI ) ) THEN
          ANGMAX = ANGMAX - 2.D0 * PI
        END IF

        ! Reset counter
        NACCT = 0
        MOVT  = 0
        NACCR = 0
        MOVR  = 0

      END IF

      ! Summary
      CALL PROGRESS_BAR_RND( ATTEMPTS )

      CYCLE LOOP_PFRACTION_FIX

    END IF

  END DO LOOP_PFRACTION_FIX

END IF

! Deallocation
DEALLOCATE( RMCV, QPROT, RPROT, EPROT )

! Summary
IF( ATTEMPTS > 1 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
END IF
IF( ETA_NPT > PACKING_F .AND. ATTEMPTS == 1 ) THEN
  WRITE( *, "(G0,G0,G0)" ) "Packing fraction fixed after ", ATTEMPTS, " attempt."
  WRITE( *, "(G0)" ) " "
ELSE IF( ETA_NPT > PACKING_F .AND. ATTEMPTS > 1 ) THEN
  WRITE( *, "(G0,G0,G0)" ) "Packing fraction fixed after ", ATTEMPTS, " attempts."
  WRITE( *, "(G0)" ) " "
END IF
CALL SLEEP( 2 )

RETURN

END SUBROUTINE CONFIG_RND

! *********************************************************************************************** !
!           This subroutine allocates particles according to a packed-box configuration           !
! *********************************************************************************************** !
SUBROUTINE CONFIG_PB(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: C ! Component index

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( COMPONENTS ) :: CUTOFF ! Cutoff diameter

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

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AXIS_SELEC(1) ) THEN
  AXISN(:) = AXISX(:)
ELSE IF( AXIS_SELEC(2) ) THEN
  AXISN(:) = AXISY(:)
ELSE IF( AXIS_SELEC(3) ) THEN
  AXISN(:) = AXISZ(:)
END IF

! *********************************************************************************************** !
! Convert degrees to radians                                                                      !
! *********************************************************************************************** !
QUATERNION_ANGLE = QUATERNION_ANGLE * PI / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
Q(0,:) = DCOS( QUATERNION_ANGLE * 0.5D0 )             ! Real part
Q(1,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(1)  ! Imaginary part (Vector)
Q(2,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(2)  ! Imaginary part (Vector)
Q(3,:) = DSIN( QUATERNION_ANGLE * 0.5D0 ) * AXISN(3)  ! Imaginary part (Vector)

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 15 )//"PACKED-BOX CONFIGURATION"//REPEAT( " ", 16 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Under construction. Exiting..."

CALL EXIT(  )

RETURN

END SUBROUTINE CONFIG_PB

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!          This subroutine generates a progress bar for the HIT-AND-MISS NVT algorithm.           !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR_HITMISS( I, J )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K ! Counters

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 29 ) :: BAR ! Progress bar

! Attempts
K = I / J

! *********************************************************************************************** !
! Progress bar (FORMAT)                                                                           !
! *********************************************************************************************** !
BAR = "Cycle #?????? | Attempt #????"

! *********************************************************************************************** !
! Progress bar (replace character positions)                                                      !
! *********************************************************************************************** !
WRITE( UNIT= BAR(8:13), FMT= "(I0.6)" ) I
WRITE( UNIT= BAR(26:29), FMT= "(I0.4)" ) K

! *********************************************************************************************** !
! Print progress bar                                                                              !
! *********************************************************************************************** !
IF( K >= 100 ) THEN
  WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A29,A65)", ADVANCE= "NO" ) CHAR(13), BAR, " (Too many attempts! Try decreasing the"// &
  &                                                                             " initial packing fraction)"
ELSE
  WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A29)", ADVANCE= "NO" ) CHAR(13), BAR
END IF

! *********************************************************************************************** !
! Flush standard output unit                                                                      !
! *********************************************************************************************** !
FLUSH( UNIT= OUTPUT_UNIT )

RETURN

END SUBROUTINE PROGRESS_BAR_HITMISS

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!          This subroutine generates a progress bar for the HIT-AND-MISS NPT algorithm.           !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR_NPT( I, ETA_NPT, ETA )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ) :: ETA_NPT, ETA ! Packing fraction

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 67 ) :: BAR ! Progress bar

! *********************************************************************************************** !
! Progress bar (FORMAT)                                                                           !
! *********************************************************************************************** !
IF( I < 10 ) THEN
  BAR = "Attempts: ? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 100 ) THEN
  BAR = "Attempts: ?? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 1000 ) THEN
  BAR = "Attempts: ??? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 10000 ) THEN
  BAR = "Attempts: ???? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 100000 ) THEN
  BAR = "Attempts: ????? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 1000000 ) THEN
  BAR = "Attempts: ?????? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 10000000 ) THEN
  BAR = "Attempts: ??????? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 100000000 ) THEN
  BAR = "Attempts: ???????? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I < 1000000000 ) THEN
  BAR = "Attempts: ????????? | Packing fraction: ??????? (TARGET = ???????)"
ELSE IF( I >= 1000000000 ) THEN
  BAR = "Attempts: ?????????? | Packing fraction: ??????? (TARGET = ???????)"
END IF

! *********************************************************************************************** !
! Progress bar (replace character positions)                                                      !
! *********************************************************************************************** !
IF( I < 10 ) THEN
  WRITE( UNIT= BAR(11:11), FMT= "(I1)" ) I
  WRITE( UNIT= BAR(33:39), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(51:57), FMT= "(F7.5)" ) ETA
ELSE IF( I < 100 ) THEN
  WRITE( UNIT= BAR(11:12), FMT= "(I2)" ) I
  WRITE( UNIT= BAR(34:40), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(52:58), FMT= "(F7.5)" ) ETA
ELSE IF( I < 1000 ) THEN
  WRITE( UNIT= BAR(11:13), FMT= "(I3)" ) I
  WRITE( UNIT= BAR(35:41), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(53:59), FMT= "(F7.5)" ) ETA
ELSE IF( I < 10000 ) THEN
  WRITE( UNIT= BAR(11:14), FMT= "(I4)" ) I
  WRITE( UNIT= BAR(36:42), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(54:60), FMT= "(F7.5)" ) ETA
ELSE IF( I < 100000 ) THEN
  WRITE( UNIT= BAR(11:15), FMT= "(I5)" ) I
  WRITE( UNIT= BAR(37:43), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(55:61), FMT= "(F7.5)" ) ETA
ELSE IF( I < 1000000 ) THEN
  WRITE( UNIT= BAR(11:16), FMT= "(I6)" ) I
  WRITE( UNIT= BAR(38:44), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(56:62), FMT= "(F7.5)" ) ETA
ELSE IF( I < 10000000 ) THEN
  WRITE( UNIT= BAR(11:17), FMT= "(I7)" ) I
  WRITE( UNIT= BAR(39:45), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(57:63), FMT= "(F7.5)" ) ETA
ELSE IF( I < 100000000 ) THEN
  WRITE( UNIT= BAR(11:18), FMT= "(I8)" ) I
  WRITE( UNIT= BAR(40:46), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(58:64), FMT= "(F7.5)" ) ETA
ELSE IF( I < 1000000000 ) THEN
  WRITE( UNIT= BAR(11:19), FMT= "(I9)" ) I
  WRITE( UNIT= BAR(41:47), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(59:65), FMT= "(F7.5)" ) ETA
ELSE IF( I >= 1000000000 ) THEN
  WRITE( UNIT= BAR(11:20), FMT= "(I10)" ) I
  WRITE( UNIT= BAR(42:48), FMT= "(F7.5)" ) ETA_NPT
  WRITE( UNIT= BAR(60:66), FMT= "(F7.5)" ) ETA
END IF

! *********************************************************************************************** !
! Print progress bar                                                                              !
! *********************************************************************************************** !
WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A67)", ADVANCE= "NO" ) CHAR(13), BAR

! *********************************************************************************************** !
! Flush standard output unit                                                                      !
! *********************************************************************************************** !
FLUSH( UNIT= OUTPUT_UNIT )

RETURN

END SUBROUTINE PROGRESS_BAR_NPT

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!                 This subroutine generates a progress bar for the RND algorithm.                 !
! *********************************************************************************************** !
SUBROUTINE PROGRESS_BAR_RND( J )

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: J ! Counter

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 16 ) :: BAR ! Progress bar

! *********************************************************************************************** !
! Progress bar (FORMAT)                                                                           !
! *********************************************************************************************** !
BAR = "Attempts: ??????"

! *********************************************************************************************** !
! Progress bar (replace character positions)                                                      !
! *********************************************************************************************** !
WRITE( UNIT= BAR(11:16), FMT= "(I6.6)" ) J

! *********************************************************************************************** !
! Print progress bar                                                                              !
! *********************************************************************************************** !
WRITE( UNIT= OUTPUT_UNIT, FMT= "(A1,A16)", ADVANCE= "NO" ) CHAR(13), BAR

! *********************************************************************************************** !
! Flush standard output unit                                                                      !
! *********************************************************************************************** !
FLUSH( UNIT= OUTPUT_UNIT )

RETURN

END SUBROUTINE PROGRESS_BAR_RND

! *********************************************************************************************** !
!      This subroutine creates a file containing all properties of the initial configuration      !
! *********************************************************************************************** !
SUBROUTINE CONFIG_OUT(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, C ! Counter

! *********************************************************************************************** !
! Simple cubic structure                                                                          !
! *********************************************************************************************** !
IF( CONFIG_SELEC(1) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_sc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_sc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) N_PARTICLES
  WRITE( 10, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (SC) | All Properties
  IF( .NOT. INIT_CONF ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/"//TRIM( DESCRIPTOR_DATE )//""//TRIM( DESCRIPTOR_HOUR )//"_initconf_sc_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GEOM_INQ
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", CONFIG_INQ
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PACKING_F
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", COMPONENTS
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", C, ": ", DIAMETER(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", C, ": ", LENGTH(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Elongation_of_Component_#", C, ": ", ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", C, ": ", MOLAR_F(C)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", N_PARTICLES
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TOTAL_VP, " Å³"
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", C, ": ", N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", C, ": ", PARTICLE_VOL(C), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BOX_VOLUME, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BOX_LENGTH(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BOX_LENGTH(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BOX_LENGTH(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BOX_LENGTH_I(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BOX_LENGTH_I(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", C, ": ", RHO_PARTICLE(C), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TOTAL_RHO, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", TEMP, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", PRESS
    WRITE( 10, * ) " "
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        FORMAT_SELF = "(G0,1X,7(G0.15,1X))"
        WRITE( 10, FORMAT_SELF ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! *********************************************************************************************** !
! Body-centered cubic structure                                                                   !
! *********************************************************************************************** !
ELSE IF( CONFIG_SELEC(2) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_bcc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_bcc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) N_PARTICLES
  WRITE( 10, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (BCC) | All Properties
  IF( .NOT. INIT_CONF ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/"//TRIM( DESCRIPTOR_DATE )//""//TRIM( DESCRIPTOR_HOUR )//"_initconf_bcc_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GEOM_INQ
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", CONFIG_INQ
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PACKING_F
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", COMPONENTS
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", C, ": ", DIAMETER(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", C, ": ", LENGTH(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Elongation_of_Component_#", C, ": ", ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", C, ": ", MOLAR_F(C)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", N_PARTICLES
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TOTAL_VP, " Å³"
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", C, ": ", N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", C, ": ", PARTICLE_VOL(C), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BOX_VOLUME, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BOX_LENGTH(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BOX_LENGTH(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BOX_LENGTH(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BOX_LENGTH_I(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BOX_LENGTH_I(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", C, ": ", RHO_PARTICLE(C), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TOTAL_RHO, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", TEMP, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", PRESS
    WRITE( 10, * ) " "
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        FORMAT_SELF = "(G0,1X,7(G0.15,1X))"
        WRITE( 10, FORMAT_SELF ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! *********************************************************************************************** !
! Face-centered cubic structure                                                                   !
! *********************************************************************************************** !
ELSE IF( CONFIG_SELEC(3) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_fcc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_fcc_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) N_PARTICLES
  WRITE( 10, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (FCC) | All Properties
  IF( .NOT. INIT_CONF ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/"//TRIM( DESCRIPTOR_DATE )//""//TRIM( DESCRIPTOR_HOUR )//"_initconf_fcc_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GEOM_INQ
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", CONFIG_INQ
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PACKING_F
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", COMPONENTS
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", C, ": ", DIAMETER(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", C, ": ", LENGTH(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Elongation_of_Component_#", C, ": ", ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", C, ": ", MOLAR_F(C)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", N_PARTICLES
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TOTAL_VP, " Å³"
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", C, ": ", N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", C, ": ", PARTICLE_VOL(C), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BOX_VOLUME, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BOX_LENGTH(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BOX_LENGTH(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BOX_LENGTH(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BOX_LENGTH_I(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BOX_LENGTH_I(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", C, ": ", RHO_PARTICLE(C), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TOTAL_RHO, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", TEMP, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", PRESS
    WRITE( 10, * ) " "
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        FORMAT_SELF = "(G0,1X,7(G0.15,1X))"
        WRITE( 10, FORMAT_SELF ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! *********************************************************************************************** !
! Random cubic structure                                                                          !
! *********************************************************************************************** !
ELSE IF( CONFIG_SELEC(4) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_rnd_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_rnd_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) N_PARTICLES
  WRITE( 10, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (RND) | All Properties
  IF( .NOT. INIT_CONF ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/"//TRIM( DESCRIPTOR_DATE )//""//TRIM( DESCRIPTOR_HOUR )//"_initconf_rnd_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GEOM_INQ
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", CONFIG_INQ
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PACKING_F
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", COMPONENTS
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", C, ": ", DIAMETER(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", C, ": ", LENGTH(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Elongation_of_Component_#", C, ": ", ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", C, ": ", MOLAR_F(C)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", N_PARTICLES
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TOTAL_VP, " Å³"
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", C, ": ", N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", C, ": ", PARTICLE_VOL(C), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BOX_VOLUME, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BOX_LENGTH(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BOX_LENGTH(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BOX_LENGTH(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BOX_LENGTH_I(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BOX_LENGTH_I(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", C, ": ", RHO_PARTICLE(C), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TOTAL_RHO, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", TEMP, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", PRESS
    WRITE( 10, * ) " "
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        FORMAT_SELF = "(G0,1X,7(G0.15,1X))"
        WRITE( 10, FORMAT_SELF ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! *********************************************************************************************** !
! Packed cubic box structure                                                                      !
! *********************************************************************************************** !
ELSE IF( CONFIG_SELEC(5) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_pb_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_initconf_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_pb_"// &
    &                     TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) N_PARTICLES
  WRITE( 10, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
      END DO
    END DO
  ELSE IF( GEOM_SELEC(2) .OR. GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        WRITE( 10, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), 0.5D0 * DIAMETER(C), &
        &                          0.5D0 * DIAMETER(C), LENGTH(C)
      END DO
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (PB) | All Properties
  IF( .NOT. INIT_CONF ) THEN
    OPEN( UNIT= 10, FILE= "Initial_Configuration/"//TRIM( DESCRIPTOR_DATE )//""//TRIM( DESCRIPTOR_HOUR )//"_initconf_pb_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GEOM_INQ
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", CONFIG_INQ
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PACKING_F
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", COMPONENTS
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", C, ": ", DIAMETER(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", C, ": ", LENGTH(C), " Å"
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Elongation_of_Component_#", C, ": ", ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", C, ": ", MOLAR_F(C)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", N_PARTICLES
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TOTAL_VP, " Å³"
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", C, ": ", N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", C, ": ", PARTICLE_VOL(C), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BOX_VOLUME, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BOX_LENGTH(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BOX_LENGTH(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BOX_LENGTH(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BOX_LENGTH_I(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BOX_LENGTH_I(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", C, ": ", RHO_PARTICLE(C), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TOTAL_RHO, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", TEMP, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", PRESS
    WRITE( 10, * ) " "
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        FORMAT_SELF = "(G0,1X,7(G0.15,1X))"
        WRITE( 10, FORMAT_SELF ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    CLOSE( 10 )
  END IF

END IF

RETURN

END SUBROUTINE CONFIG_OUT

END MODULE INITCONFIG
