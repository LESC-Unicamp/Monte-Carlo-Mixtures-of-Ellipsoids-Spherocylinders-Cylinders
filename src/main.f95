! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!     The Perram-Wertheim Method (J. W. Perram & M. S. Wertheim, 1985) is used to search for      !
!             molecular overlaps between ellipsoids of revolution after a trial move.             !
!              The Vega-Lago Method (C. Vega & S. Lago, 1993) is used to search for               !
!                 molecular overlaps between spherocylinders after a trial move.                  !
!  The algorithm of Lopes et al. (J. Lopes, F. Romano, E. Grelet, L. Franco, A. Giacometti, 2021) !
!         is used to search for molecular overlaps between cylinders after a trial move.          !
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
! Main References:                 J. W. Perram, M. S. Wertheim                                   !
!                               J. Comput. Phys 15, 409-416 (1985)                                !
!                                DOI: 10.1016/0021-9991(85)90171-8                                !
!                             --------------------------------------                              !
!                                   J. W. Perram, J. Rasmussen                                    !
!                                  Phys. Rev. E 54, 6565 (1996)                                   !
!                                 DOI: 10.1103/PhysRevE.54.6565                                   !
!                             --------------------------------------                              !
!                                        C. Vega, S. Lago                                         !
!                                 Computers Chem. 18, 55-59 (1993)                                !
!                                DOI: 10.1016/0097-8485(94)80023-5                                !
!                             --------------------------------------                              !
!                    J. Lopes, F. Romano, E. Grelet, L. Franco, A. Giacometti                     !
!                                 Chem. Phys. 154, 104902 (2021)                                  !
!                                     DOI: 10.1063/5.0040942                                      !
!                             --------------------------------------                              !
!                                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                      N. Metropolis et al.                                       !
!                                 J. Chem. Phys. 21, 1087 (1953)                                  !
!                                     DOI: 10.1063/1.1699114                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
PROGRAM MAIN

! USES FOUR MODULES: global variables, variable initialization, initial configuration, and directory creator
USE GLOBALVAR
USE INITVAR
USE INITCONFIG
USE FOLDERS

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SIZE_SEED   ! Seed array size
INTEGER :: SEED_COMP   ! Seed array component

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: RANS

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K, L ! Counters
INTEGER( KIND= INT64 ) :: CI, CJ, C  ! Counters (component)
INTEGER( KIND= INT64 ) :: C_LAMB     ! Counter
INTEGER( KIND= INT64 ) :: COMPONENT  ! Box matrix component
INTEGER( KIND= INT64 ) :: CYCLES     ! Cycles
INTEGER( KIND= INT64 ) :: NACCT      ! Move acceptance counter: Translation
INTEGER( KIND= INT64 ) :: NACCR      ! Move acceptance counter: Rotation
INTEGER( KIND= INT64 ) :: NACCVI     ! Move acceptance counter: Volume change
INTEGER( KIND= INT64 ) :: NACCVA     ! Move acceptance counter: Volume change
INTEGER( KIND= INT64 ) :: MOVT       ! Move counter (Translation)
INTEGER( KIND= INT64 ) :: MOVR       ! Move counter (Rotation)
INTEGER( KIND= INT64 ) :: MOVVI      ! Move counter (Volume change)
INTEGER( KIND= INT64 ) :: MOVVA      ! Move counter (Volume change)
INTEGER( KIND= INT64 ) :: REMAINDER  ! Remainder of division

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                      :: RIJSQ            ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                      :: CD               ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( KIND= REAL64 )                      :: CUTOFF_D         ! Cutoff distance
REAL( KIND= REAL64 )                      :: RATIO            ! Acceptance ratio (simulation)
REAL( KIND= REAL64 )                      :: S                ! Nematic order parameter
REAL( KIND= REAL64 )                      :: BOXVM, BOXVN     ! Volume of simulation box (before/after a trial move)
REAL( KIND= REAL64 )                      :: HNM              ! Enthalpy criterion (reduced)
REAL( KIND= REAL64 )                      :: SCALE_FACTOR     ! Scaling factor
REAL( KIND= REAL64 )                      :: DISTORTION       ! Box distortion
REAL( KIND= REAL64 )                      :: EXEC_TIME        ! Execution time
REAL( KIND= REAL64 )                      :: BOXVROT          ! Volume of simulation box (after undoing box rotation)
REAL( KIND= REAL64 )                      :: THETA            ! Angle between box vector and coordination system
REAL( KIND= REAL64 )                      :: RAXISMAG         ! Magnitude of rotation axis
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: PROJY_XY         ! Projection of the y-vector of the box onto the ZY-plane and the unit vector of the y-axis
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: S12              ! Position (unit box)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: RIJ              ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: EM, EN           ! Orientation (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: EI, EJ           ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: RM, RN           ! Position (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: RI, RJ           ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 9 )      :: BOXLM, BOXLN     ! Length of simulation box (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 9 )      :: BOXLM_I, BOXLN_I ! Length of simulation box (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 9 )      :: BOXLMC           ! Length (x,y,z) of the simulation box
REAL( KIND= REAL64 ), DIMENSION( 9 )      :: BOXLMC_I         ! Length (x,y,z) of simulation box (inverse)
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: COSANGLE_VEC     ! Cossine of angle between box vectors
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: LBOX             ! Length of box edges
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: LBOXR            ! Length ratio of box edges
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: RAXIS            ! Rotation axis
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: AUXV             ! Auxiliary vector
REAL( KIND= REAL64 ), DIMENSION( 3 )      :: V1, V2, V3       ! Box vectors
REAL( KIND= REAL64 ), DIMENSION( 9 )      :: BOXLROT, BOXIROT ! Length of simulation box (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )    :: QROT             ! Rotation quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 )    :: QAUX             ! Auxiliary quaternion
REAL( KIND= REAL64 ), DIMENSION( 0:3 )    :: QM, QN           ! Quaternion (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )    :: QI, QJ           ! Quaternion of particles i and j

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: CUTOFF     ! Cutoff distance
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: VM, VN     ! Perturbed potential energy (before/after a trial move)
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: DV         ! Perturbed potential energy difference between two microstates m and n
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: A1, A2     ! First- and second-order TPT coefficients
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: APERT      ! Full Helmholtz free energy of the perturbed system
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: DPA1, DPA2 ! First- and second-order TPT coefficients (standard deviation)
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: DPAPERT    ! Full Helmholtz free energy of the perturbed system (standard deviation)
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: STRSH      ! String size
REAL( KIND= REAL64 ), DIMENSION( : ), ALLOCATABLE    :: STRSS      ! String size
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: STRST      ! String size
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: RMCV       ! Old position of particles
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: RPROT      ! Position of the center of mass (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: QPROT      ! Quaternion of the center of mass (after undoing box rotation)
REAL( KIND= REAL64 ), DIMENSION( :, : ), ALLOCATABLE :: EPROT      ! Orientation of the center of mass (after undoing box rotation)

! *********************************************************************************************** !
! CHARACTER STRINGS (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
CHARACTER( LEN= 68 ), DIMENSION( : ), ALLOCATABLE     :: LOG_STRINGS_H ! Strings used in the simulation log file
CHARACTER( LEN= 68 ), DIMENSION( : ), ALLOCATABLE     :: LOG_STRINGS_S ! Strings used in the simulation log file
CHARACTER( LEN= 68 ), DIMENSION( :, : ), ALLOCATABLE  :: LOG_STRINGS_T ! Strings used in the simulation log file
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: LOG_STRINGS_L ! Strings used in the simulation log file
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: CHAR_LABEL    ! Simulation log
CHARACTER( LEN= 100 ), DIMENSION( :, : ), ALLOCATABLE :: CHAR_LABEL_L  ! Simulation log

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP            ! Detects overlap between two particles (PW) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAP_PRELIMINAR ! Detects overlap between two particles : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PARALLEL           ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: MOV_ROT            ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_TRANS          ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_VOL_I          ! Volumetric movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MOV_VOL_A          ! Volumetric movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: IGNORE_VOL_ATTEMPT ! Detects if a box deformation is valid or not : TRUE = ignore box deformation; FALSE = consider box deformation
LOGICAL :: LATTICER           ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: FLAG               ! Generic TRUE/FALSE flag

! Title
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 6 )//"NVT/NPT-MONTE CARLO SIMULATION OF MIXTURES"//REPEAT( " ", 7 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR

! Molecular geometry selection (see 'INITCONFIG' module)
CALL GEOM_SELECTION(  )

! Molecular configuration selection (see 'INITCONFIG' module)
CALL CONFIG_SELECTION(  )

! Initialization of Monte Carlo parameters (see 'INITVAR' module)
CALL MONTECARLO_VAR(  )

! Initialization of Inquiry/Control variables (see 'INITVAR' module)
CALL INQUERY_VAR(  )

! Initialization of common variables (see 'INITVAR' module)
CALL COMMON_VAR(  )

! Initialization of potential variables (see 'INITVAR' module)
CALL POTENTIAL_VAR(  )

! Pseudorandom number generator seed
IF( FSEED ) THEN
  SEED = 123456789
ELSE IF( .NOT. FSEED ) THEN
  CALL RANDOM_SEED( SIZE= SIZE_SEED )
  ALLOCATE( RANS(SIZE_SEED) )
  CALL RANDOM_SEED( GET= RANS )
  CALL RANDOM_NUMBER( RANDOM_N )
  SEED_COMP = INT( RANDOM_N * DBLE( SIZE_SEED ) ) + 1
  SEED = ABS( RANS(SEED_COMP) )
END IF

! Allocation
ALLOCATE( Q(0:3,N_PARTICLES), QMC(0:3,N_PARTICLES) )
ALLOCATE( R(3,N_PARTICLES), RMC(3,N_PARTICLES) )
ALLOCATE( E(3,N_PARTICLES), EMC(3,N_PARTICLES) )
ALLOCATE( RMCV(3,N_PARTICLES) )
ALLOCATE( CUTOFF(COMPONENTS) )
ALLOCATE( A1(N_LAMBDA), A2(N_LAMBDA), APERT(N_LAMBDA) )
ALLOCATE( DPA1(N_LAMBDA), DPA2(N_LAMBDA), DPAPERT(N_LAMBDA) )
ALLOCATE( V(N_LAMBDA), VMC(N_LAMBDA) )
ALLOCATE( VN(N_LAMBDA), VM(N_LAMBDA), DV(N_LAMBDA) )
ALLOCATE( SWRANGE(COMPONENTS,N_LAMBDA) )
ALLOCATE( QPROT(0:3,N_PARTICLES) )
ALLOCATE( RPROT(3,N_PARTICLES) )
ALLOCATE( EPROT(3,N_PARTICLES) )

! Atom ID (Required in some visualization and analysis software)
ALLOCATE( INDEX_P(COMPONENTS) )
DO C = 1, COMPONENTS
  INDEX_P(C) = C
END DO

! Unrotated reference orientation (Allen and Tildesley, 2nd Edition, pages 106-111)
AXISX(1) = 1.D0 ! X-axis
AXISX(2) = 0.D0 ! X-axis
AXISX(3) = 0.D0 ! X-axis
AXISY(1) = 0.D0 ! Y-axis
AXISY(2) = 1.D0 ! Y-axis
AXISY(3) = 0.D0 ! Y-axis
AXISZ(1) = 0.D0 ! Z-axis
AXISZ(2) = 0.D0 ! Z-axis
AXISZ(3) = 1.D0 ! Z-axis

! CPU Clock
CALL DATE_AND_TIME( VALUES= DATE_TIME )

! Date format (YYYY/MM/DD)
FORMAT_DATE = "(I4,2I2.2)"
WRITE( DESCRIPTOR_DATE, FORMAT_DATE ) DATE_TIME(1), DATE_TIME(2), DATE_TIME(3)

! Time format (HH:MM:SS)
FORMAT_HOUR = "(3I2.2)"
! Hour descriptor
WRITE( DESCRIPTOR_HOUR, FORMAT_HOUR ) DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)

! Output file descriptors (p. fraction/pressure [1], # of components [2], and geometry [3])
FORMAT_FILE1 = "(F0.5)"
IF( MC_ENSEMBLE == "NVT" ) THEN
  WRITE( DESCRIPTOR_FILE1, FORMAT_FILE1 ) PACKING_F
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  WRITE( DESCRIPTOR_FILE1, FORMAT_FILE1 ) PRESS
END IF
FORMAT_FILE2 = "(I0.3)"
WRITE( DESCRIPTOR_FILE2, FORMAT_FILE2 ) COMPONENTS
FORMAT_FILE3 = "(A3)"
WRITE( DESCRIPTOR_FILE3, FORMAT_FILE3 ) GEO_ACRONYM

! Initial configuration (see 'INITCONFIG' module)
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 17 )//"INITIAL CONFIGURATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Setting up initial configuration folder..."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Initial configuration folder (see 'FOLDERS' module)
CALL INITFOLDER(  )

! Calls 'CONFIG_SC' subroutine if the user chooses a simple cubic structure
IF( CONFIG_SELEC(1) ) THEN
  IF( .NOT. INIT_CONF ) THEN
    DO C = 1, COMPONENTS - 1
      IF( DABS( DIAMETER(C) - DIAMETER(C+1) ) >= EPSILON( 1.D0 ) .OR. DABS( LENGTH(C) - LENGTH(C+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The simple cubic structure might create overlaping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
      IF( ( SPHERCOMP(C) .AND. .NOT. SPHERCOMP(C+1) ) .OR. ( .NOT. SPHERCOMP(C) .AND. SPHERCOMP(C+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The simple cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
    END DO
    CALL CONFIG_SC(  )
    CALL CONFIG_OUT(  )
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be overwritten.", &
    &                   " Continue? [Y/N]"
    READ( *, * ) DUMMY
    CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
    IF( DUMMY /= "Y" ) THEN
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) CODE_DESCRIPTOR
    INQUIRE( FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_sc_"//TRIM( DESCRIPTOR_FILE3 )//".xyz", &
    &        EXIST= FILE_EXIST )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FILE_EXIST ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL SLEEP( 1 )
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL SLEEP( 1 )
    OPEN( UNIT= 10, ACTION= "READ", FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_sc_" &
    &                                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    READ( 10, * ) DUMMY, GEOM_INQ
    ! Extended configuration name
    GEOM_SELEC(:) = .FALSE.
    IF( GEOM_INQ == "EOR" ) THEN
      GEOMETRY    = "Ellipsoids of revolution"
      GEO_ACRONYM = "eor"
      GEOM_SELEC(1) = .TRUE.
    ELSE IF( GEOM_INQ == "SPC" ) THEN
      GEOMETRY    = "Spherocylinders"
      GEO_ACRONYM = "spc"
      GEOM_SELEC(2) = .TRUE.
    ELSE IF( GEOM_INQ == "CYL" ) THEN
      GEOMETRY    = "Cylinders"
      GEO_ACRONYM = "cyl"
      GEOM_SELEC(3) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, CONFIG_INQ
    ! Extended configuration name
    CONFIG_SELEC(:) = .FALSE.
    IF( CONFIG_INQ == "SC" ) THEN
      CONFIGURATION   = "Simple Cube"
      CONFIG_SELEC(1) = .TRUE.
    ELSE IF( CONFIG_INQ == "BCC" ) THEN
      CONFIGURATION   = "Body-Centered Cube"
      CONFIG_SELEC(2) = .TRUE.
    ELSE IF( CONFIG_INQ == "FCC" ) THEN
      CONFIGURATION   = "Face-Centered Cube"
      CONFIG_SELEC(3) = .TRUE.
    ELSE IF( CONFIG_INQ == "RND" ) THEN
      CONFIGURATION   = "Random"
      CONFIG_SELEC(4) = .TRUE.
    ELSE IF( CONFIG_INQ == "PB" ) THEN
      CONFIGURATION   = "Packed Box"
      CONFIG_SELEC(5) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, PACKING_F
    READ( 10, * ) DUMMY, COMPONENTS
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, SPHCOMP_INQ(C)
      IF( SPHCOMP_INQ(C) == "T" ) THEN
        SPHERCOMP(C) = .TRUE.
      ELSE
        SPHERCOMP(C) = .FALSE.
      END IF
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, DIAMETER(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, LENGTH(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, MOLAR_F(C)
    END DO
    READ( 10, * ) DUMMY, N_PARTICLES
    READ( 10, * ) DUMMY, TOTAL_VP, DUMMY
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, PARTICLE_VOL(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, BOX_VOLUME, DUMMY
    READ( 10, * ) DUMMY, BOX_LENGTH(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH(7:9)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, RHO_PARTICLE(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, TOTAL_RHO, DUMMY
    READ( 10, * ) DUMMY, TEMP, DUMMY
    READ( 10, * ) DUMMY, PRESS
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        READ( 10, * ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
    CALL SLEEP( 1 )
    WRITE( *, "(G0)" ) " "
    CALL CONFIG_OUT(  )
  END IF
! Calls 'CONFIG_BCC' subroutine if the user chooses a body-centered cubic structure
ELSE IF( CONFIG_SELEC(2) ) THEN
  IF( .NOT. INIT_CONF ) THEN
    DO C = 1, COMPONENTS - 1
      IF( DABS( DIAMETER(C) - DIAMETER(C+1) ) >= EPSILON( 1.D0 ) .OR. DABS( LENGTH(C) - LENGTH(C+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The body-centered cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
      IF( ( SPHERCOMP(C) .AND. .NOT. SPHERCOMP(C+1) ) .OR. ( .NOT. SPHERCOMP(C) .AND. SPHERCOMP(C+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The body-centered cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
    END DO
    CALL CONFIG_BCC(  )
    CALL CONFIG_OUT(  )
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be overwritten.", &
    &                   " Continue? [Y/N]"
    READ( *, * ) DUMMY
    CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
    IF( DUMMY /= "Y" ) THEN
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) CODE_DESCRIPTOR
    INQUIRE( FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_bcc_"//TRIM( DESCRIPTOR_FILE3 )//".xyz", &
    &        EXIST= FILE_EXIST )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FILE_EXIST ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL SLEEP( 1 )
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL SLEEP( 1 )
    OPEN( UNIT= 10, ACTION= "READ", FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_bcc_" &
    &                                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    READ( 10, * ) DUMMY, GEOM_INQ
    ! Extended configuration name
    GEOM_SELEC(:) = .FALSE.
    IF( GEOM_INQ == "EOR" ) THEN
      GEOMETRY    = "Ellipsoids of revolution"
      GEO_ACRONYM = "eor"
      GEOM_SELEC(1) = .TRUE.
    ELSE IF( GEOM_INQ == "SPC" ) THEN
      GEOMETRY    = "Spherocylinders"
      GEO_ACRONYM = "spc"
      GEOM_SELEC(2) = .TRUE.
    ELSE IF( GEOM_INQ == "CYL" ) THEN
      GEOMETRY    = "Cylinders"
      GEO_ACRONYM = "cyl"
      GEOM_SELEC(3) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, CONFIG_INQ
    ! Extended configuration name
    CONFIG_SELEC(:) = .FALSE.
    IF( CONFIG_INQ == "SC" ) THEN
      CONFIGURATION   = "Simple Cube"
      CONFIG_SELEC(1) = .TRUE.
    ELSE IF( CONFIG_INQ == "BCC" ) THEN
      CONFIGURATION   = "Body-Centered Cube"
      CONFIG_SELEC(2) = .TRUE.
    ELSE IF( CONFIG_INQ == "FCC" ) THEN
      CONFIGURATION   = "Face-Centered Cube"
      CONFIG_SELEC(3) = .TRUE.
    ELSE IF( CONFIG_INQ == "RND" ) THEN
      CONFIGURATION   = "Random"
      CONFIG_SELEC(4) = .TRUE.
    ELSE IF( CONFIG_INQ == "PB" ) THEN
      CONFIGURATION   = "Packed Box"
      CONFIG_SELEC(5) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, PACKING_F
    READ( 10, * ) DUMMY, COMPONENTS
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, SPHCOMP_INQ(C)
      IF( SPHCOMP_INQ(C) == "T" ) THEN
        SPHERCOMP(C) = .TRUE.
      ELSE
        SPHERCOMP(C) = .FALSE.
      END IF
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, DIAMETER(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, LENGTH(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, MOLAR_F(C)
    END DO
    READ( 10, * ) DUMMY, N_PARTICLES
    READ( 10, * ) DUMMY, TOTAL_VP, DUMMY
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, PARTICLE_VOL(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, BOX_VOLUME, DUMMY
    READ( 10, * ) DUMMY, BOX_LENGTH(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH(7:9)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, RHO_PARTICLE(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, TOTAL_RHO, DUMMY
    READ( 10, * ) DUMMY, TEMP, DUMMY
    READ( 10, * ) DUMMY, PRESS
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        READ( 10, * ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
    CALL SLEEP( 1 )
    WRITE( *, "(G0)" ) " "
    CALL CONFIG_OUT(  )
  END IF
! Calls 'CONFIG_FCC' subroutine if the user chooses a face-centered cubic structure
ELSE IF( CONFIG_SELEC(3) ) THEN
  IF( .NOT. INIT_CONF ) THEN
    DO C = 1, COMPONENTS - 1
      IF( DABS( DIAMETER(C) - DIAMETER(C+1) ) >= EPSILON( 1.D0 ) .OR. DABS( LENGTH(C) - LENGTH(C+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The face-centered cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
      IF( ( SPHERCOMP(C) .AND. .NOT. SPHERCOMP(C+1) ) .OR. ( .NOT. SPHERCOMP(C) .AND. SPHERCOMP(C+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", C, " are not isomorphic with molecules of component ", C + 1, "! ", &
        &                         "The face-centered cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) DUMMY
        CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
        IF( DUMMY /= "Y" ) THEN
          CALL SLEEP( 1 )
          CALL EXIT(  )
        END IF
        WRITE( *, "(G0)" ) " "
      END IF
    END DO
    CALL CONFIG_FCC(  )
    CALL CONFIG_OUT(  )
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be overwritten.", &
    &                   " Continue? [Y/N]"
    READ( *, * ) DUMMY
    CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
    IF( DUMMY /= "Y" ) THEN
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) CODE_DESCRIPTOR
    INQUIRE( FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_fcc_"//TRIM( DESCRIPTOR_FILE3 )//".xyz", &
    &        EXIST= FILE_EXIST )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FILE_EXIST ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL SLEEP( 1 )
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL SLEEP( 1 )
    OPEN( UNIT= 10, ACTION= "READ", FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_fcc_" &
    &                                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    READ( 10, * ) DUMMY, GEOM_INQ
    ! Extended configuration name
    GEOM_SELEC(:) = .FALSE.
    IF( GEOM_INQ == "EOR" ) THEN
      GEOMETRY    = "Ellipsoids of revolution"
      GEO_ACRONYM = "eor"
      GEOM_SELEC(1) = .TRUE.
    ELSE IF( GEOM_INQ == "SPC" ) THEN
      GEOMETRY    = "Spherocylinders"
      GEO_ACRONYM = "spc"
      GEOM_SELEC(2) = .TRUE.
    ELSE IF( GEOM_INQ == "CYL" ) THEN
      GEOMETRY    = "Cylinders"
      GEO_ACRONYM = "cyl"
      GEOM_SELEC(3) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, CONFIG_INQ
    ! Extended configuration name
    CONFIG_SELEC(:) = .FALSE.
    IF( CONFIG_INQ == "SC" ) THEN
      CONFIGURATION   = "Simple Cube"
      CONFIG_SELEC(1) = .TRUE.
    ELSE IF( CONFIG_INQ == "BCC" ) THEN
      CONFIGURATION   = "Body-Centered Cube"
      CONFIG_SELEC(2) = .TRUE.
    ELSE IF( CONFIG_INQ == "FCC" ) THEN
      CONFIGURATION   = "Face-Centered Cube"
      CONFIG_SELEC(3) = .TRUE.
    ELSE IF( CONFIG_INQ == "RND" ) THEN
      CONFIGURATION   = "Random"
      CONFIG_SELEC(4) = .TRUE.
    ELSE IF( CONFIG_INQ == "PB" ) THEN
      CONFIGURATION   = "Packed Box"
      CONFIG_SELEC(5) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, PACKING_F
    READ( 10, * ) DUMMY, COMPONENTS
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, SPHCOMP_INQ(C)
      IF( SPHCOMP_INQ(C) == "T" ) THEN
        SPHERCOMP(C) = .TRUE.
      ELSE
        SPHERCOMP(C) = .FALSE.
      END IF
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, DIAMETER(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, LENGTH(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, MOLAR_F(C)
    END DO
    READ( 10, * ) DUMMY, N_PARTICLES
    READ( 10, * ) DUMMY, TOTAL_VP, DUMMY
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, PARTICLE_VOL(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, BOX_VOLUME, DUMMY
    READ( 10, * ) DUMMY, BOX_LENGTH(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH(7:9)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, RHO_PARTICLE(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, TOTAL_RHO, DUMMY
    READ( 10, * ) DUMMY, TEMP, DUMMY
    READ( 10, * ) DUMMY, PRESS
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        READ( 10, * ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
    CALL SLEEP( 1 )
    WRITE( *, "(G0)" ) " "
    CALL CONFIG_OUT(  )
  END IF
! Calls 'CONFIG_RND' subroutine if the user chooses a random structure
ELSE IF( CONFIG_SELEC(4) ) THEN
  IF( .NOT. INIT_CONF ) THEN
    CALL CONFIG_RND(  )
    CALL CONFIG_OUT(  )
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be overwritten.", &
    &                   " Continue? [Y/N]"
    READ( *, * ) DUMMY
    CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
    IF( DUMMY /= "Y" ) THEN
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) CODE_DESCRIPTOR
    INQUIRE( FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_rnd_"//TRIM( DESCRIPTOR_FILE3 )//".xyz", &
    &        EXIST= FILE_EXIST )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FILE_EXIST ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL SLEEP( 1 )
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL SLEEP( 1 )
    OPEN( UNIT= 10, ACTION= "READ", FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_rnd_" &
    &                                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    READ( 10, * ) DUMMY, GEOM_INQ
    ! Extended configuration name
    GEOM_SELEC(:) = .FALSE.
    IF( GEOM_INQ == "EOR" ) THEN
      GEOMETRY    = "Ellipsoids of revolution"
      GEO_ACRONYM = "eor"
      GEOM_SELEC(1) = .TRUE.
    ELSE IF( GEOM_INQ == "SPC" ) THEN
      GEOMETRY    = "Spherocylinders"
      GEO_ACRONYM = "spc"
      GEOM_SELEC(2) = .TRUE.
    ELSE IF( GEOM_INQ == "CYL" ) THEN
      GEOMETRY    = "Cylinders"
      GEO_ACRONYM = "cyl"
      GEOM_SELEC(3) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, CONFIG_INQ
    ! Extended configuration name
    CONFIG_SELEC(:) = .FALSE.
    IF( CONFIG_INQ == "SC" ) THEN
      CONFIGURATION   = "Simple Cube"
      CONFIG_SELEC(1) = .TRUE.
    ELSE IF( CONFIG_INQ == "BCC" ) THEN
      CONFIGURATION   = "Body-Centered Cube"
      CONFIG_SELEC(2) = .TRUE.
    ELSE IF( CONFIG_INQ == "FCC" ) THEN
      CONFIGURATION   = "Face-Centered Cube"
      CONFIG_SELEC(3) = .TRUE.
    ELSE IF( CONFIG_INQ == "RND" ) THEN
      CONFIGURATION   = "Random"
      CONFIG_SELEC(4) = .TRUE.
    ELSE IF( CONFIG_INQ == "PB" ) THEN
      CONFIGURATION   = "Packed Box"
      CONFIG_SELEC(5) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, PACKING_F
    READ( 10, * ) DUMMY, COMPONENTS
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, SPHCOMP_INQ(C)
      IF( SPHCOMP_INQ(C) == "T" ) THEN
        SPHERCOMP(C) = .TRUE.
      ELSE
        SPHERCOMP(C) = .FALSE.
      END IF
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, DIAMETER(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, LENGTH(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, MOLAR_F(C)
    END DO
    READ( 10, * ) DUMMY, N_PARTICLES
    READ( 10, * ) DUMMY, TOTAL_VP, DUMMY
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, PARTICLE_VOL(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, BOX_VOLUME, DUMMY
    READ( 10, * ) DUMMY, BOX_LENGTH(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH(7:9)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, RHO_PARTICLE(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, TOTAL_RHO, DUMMY
    READ( 10, * ) DUMMY, TEMP, DUMMY
    READ( 10, * ) DUMMY, PRESS
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        READ( 10, * ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
    CALL SLEEP( 1 )
    WRITE( *, "(G0)" ) " "
    CALL CONFIG_OUT(  )
  END IF
! Calls 'CONFIG_PB' subroutine if the user chooses a packed-box structure
ELSE IF( CONFIG_SELEC(5) ) THEN
  IF( .NOT. INIT_CONF ) THEN
    CALL CONFIG_PB(  )
    CALL CONFIG_OUT(  )
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be overwritten.", &
    &                   " Continue? [Y/N]"
    READ( *, * ) DUMMY
    CALL TO_UPPER( DUMMY, LEN_TRIM( DUMMY ), DUMMY )
    IF( DUMMY /= "Y" ) THEN
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) CODE_DESCRIPTOR
    INQUIRE( FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_pb_"//TRIM( DESCRIPTOR_FILE3 )//".xyz", &
    &        EXIST= FILE_EXIST )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FILE_EXIST ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL SLEEP( 1 )
      CALL EXIT(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL SLEEP( 1 )
    OPEN( UNIT= 10, ACTION= "READ", FILE= "Initial_Configuration/"//TRIM( CODE_DESCRIPTOR )//"_initconf_pb_" &
    &                                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
    READ( 10, * ) DUMMY, GEOM_INQ
    ! Extended configuration name
    GEOM_SELEC(:) = .FALSE.
    IF( GEOM_INQ == "EOR" ) THEN
      GEOMETRY    = "Ellipsoids of revolution"
      GEO_ACRONYM = "eor"
      GEOM_SELEC(1) = .TRUE.
    ELSE IF( GEOM_INQ == "SPC" ) THEN
      GEOMETRY    = "Spherocylinders"
      GEO_ACRONYM = "spc"
      GEOM_SELEC(2) = .TRUE.
    ELSE IF( GEOM_INQ == "CYL" ) THEN
      GEOMETRY    = "Cylinders"
      GEO_ACRONYM = "cyl"
      GEOM_SELEC(3) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, CONFIG_INQ
    ! Extended configuration name
    CONFIG_SELEC(:) = .FALSE.
    IF( CONFIG_INQ == "SC" ) THEN
      CONFIGURATION   = "Simple Cube"
      CONFIG_SELEC(1) = .TRUE.
    ELSE IF( CONFIG_INQ == "BCC" ) THEN
      CONFIGURATION   = "Body-Centered Cube"
      CONFIG_SELEC(2) = .TRUE.
    ELSE IF( CONFIG_INQ == "FCC" ) THEN
      CONFIGURATION   = "Face-Centered Cube"
      CONFIG_SELEC(3) = .TRUE.
    ELSE IF( CONFIG_INQ == "RND" ) THEN
      CONFIGURATION   = "Random"
      CONFIG_SELEC(4) = .TRUE.
    ELSE IF( CONFIG_INQ == "PB" ) THEN
      CONFIGURATION   = "Packed Box"
      CONFIG_SELEC(5) = .TRUE.
    END IF
    READ( 10, * ) DUMMY, PACKING_F
    READ( 10, * ) DUMMY, COMPONENTS
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, SPHCOMP_INQ(C)
      IF( SPHCOMP_INQ(C) == "T" ) THEN
        SPHERCOMP(C) = .TRUE.
      ELSE
        SPHERCOMP(C) = .FALSE.
      END IF
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, DIAMETER(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, LENGTH(C), DUMMY
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, ASPECT_RATIO(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, MOLAR_F(C)
    END DO
    READ( 10, * ) DUMMY, N_PARTICLES
    READ( 10, * ) DUMMY, TOTAL_VP, DUMMY
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, N_COMPONENT(C)
    END DO
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, PARTICLE_VOL(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, BOX_VOLUME, DUMMY
    READ( 10, * ) DUMMY, BOX_LENGTH(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH(7:9)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(1:3)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(4:6)
    READ( 10, * ) DUMMY, BOX_LENGTH_I(7:9)
    DO C = 1, COMPONENTS
      READ( 10, * ) DUMMY, RHO_PARTICLE(C), DUMMY
    END DO
    READ( 10, * ) DUMMY, TOTAL_RHO, DUMMY
    READ( 10, * ) DUMMY, TEMP, DUMMY
    READ( 10, * ) DUMMY, PRESS
    DO C = 1, COMPONENTS
      DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
        READ( 10, * ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I)
      END DO
    END DO
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
    CALL SLEEP( 1 )
    WRITE( *, "(G0)" ) " "
    CALL CONFIG_OUT(  )
  END IF
END IF

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"FOLDER ORGANIZER"//REPEAT( " ", 20 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Setting up folders..."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Create output directories (see 'FOLDERS' module)
CALL INQUIRE_FOLDERS(  )

! Create date subfolders (see 'FOLDERS' module)
CALL DATE_FOLDERS(  )

! Initialization of the attractive range subfolders (see 'FOLDERS' module)
IF( POTENTIAL_SELEC(2) ) THEN
  CALL LAMBDA_FOLDERS(  )
END IF

! Hard-core volumetric relation (EOR/SPC/HC and SPHERES)
IF( GEOM_SELEC(1) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      SIGSPHERE(C) = DIAMETER(C) * ( ( ASPECT_RATIO(C) ) ** ( 1.D0 / 3.D0 ) ) ! Ellipsoids of revolution
    ELSE
      SIGSPHERE(C) = DIAMETER(C) ! Spheres
    END IF
  END DO
ELSE IF( GEOM_SELEC(2) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      SIGSPHERE(C) = DIAMETER(C) * ( ( 1.D0 + ( 1.5D0 * ASPECT_RATIO(C) ) ) ** ( 1.D0 / 3.D0 ) ) ! Spherocylinders
    ELSE
      SIGSPHERE(C) = DIAMETER(C) ! Spheres
    END IF
  END DO
ELSE IF( GEOM_SELEC(3) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      SIGSPHERE(C) = DIAMETER(C) * ( ( 1.5D0 * ASPECT_RATIO(C) ) ** ( 1.D0 / 3.D0 ) ) ! Cylinders
    ELSE
      SIGSPHERE(C) = DIAMETER(C) ! Spheres
    END IF
  END DO
END IF

! Effective range of attraction
IF( POTENTIAL_SELEC(2) ) THEN
  DO C = 1, COMPONENTS
    SWRANGE(C,:) = L_RANGE(:) * SIGSPHERE(C)
  END DO
END IF

! Active transformation (orientation of particles)
DO I = 1, N_PARTICLES
  CALL ACTIVE_TRANSFORMATION( AXISZ, Q(:,I), E(:,I) )
END DO

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 6 )//"OVERLAP CHECK FOR THE INITIAL CONFIGURATION"//REPEAT( " ", 6 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)" ) "Verifying initial configuration..."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Diameter of circumscribing sphere
IF( GEOM_SELEC(1) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      IF( ASPECT_RATIO(C) > 0.D0 .AND. ASPECT_RATIO(C) < 1.D0 ) THEN
        CUTOFF(C) = DIAMETER(C)
      ELSE IF( ASPECT_RATIO(C) > 1.D0 ) THEN
        CUTOFF(C) = LENGTH(C)
      END IF
    ELSE
      CUTOFF(C) = DIAMETER(C)
    END IF
  END DO
ELSE IF( GEOM_SELEC(2) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      CUTOFF(C) = DIAMETER(C) + LENGTH(C)
    ELSE
      CUTOFF(C) = DIAMETER(C)
    END IF
  END DO
ELSE IF( GEOM_SELEC(3) ) THEN
  DO C = 1, COMPONENTS
    IF( .NOT. SPHERCOMP(C) ) THEN
      CUTOFF(C) = DIAMETER(C) + LENGTH(C)
    ELSE
      CUTOFF(C) = DIAMETER(C)
    END IF
  END DO
END IF

! Overlap check (initial configuration)
OVERLAP = .FALSE.
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
        ! Minimum image convention
        CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
        S12 = S12 - ANINT(S12)
        CALL MULTI_MATRIX( BOX_LENGTH, S12, RIJ )
        ! Magnitude of the vector distance (squared)
        RIJSQ = ( RIJ(1) * RIJ(1) ) + ( RIJ(2) * RIJ(2) ) + ( RIJ(3) * RIJ(3) )
        ! Cutoff distance (squared)
        CUTOFF_D = 0.5D0 * ( CUTOFF(CI) + CUTOFF(CJ) )
        CUTOFF_D = CUTOFF_D * CUTOFF_D
        ! Preliminary test (circumscribing spheres)
        IF( RIJSQ <= CUTOFF_D ) THEN
          ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
          IF( GEOM_SELEC(1) ) THEN
            CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
            ! Overlap criterion
            IF( OVERLAP ) THEN
              ! Overlap detected
              IF( COMPONENTS > 1 ) THEN
                WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
              ELSE IF( COMPONENTS == 1 ) THEN
                WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
              END IF
              CALL SLEEP( 1 )
              CALL EXIT(  )
            END IF
          ! Overlap test for spherocylinders (Vega-Lago method)
          ELSE IF( GEOM_SELEC(2) ) THEN
            CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
            ! Overlap criterion
            IF( OVERLAP ) THEN
              ! Overlap detected
              IF( COMPONENTS > 1 ) THEN
                WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
              ELSE IF( COMPONENTS == 1 ) THEN
                WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
              END IF
              CALL SLEEP( 1 )
              CALL EXIT(  )
            END IF
          ! Overlap test for cylinders and/or spheres
          ELSE IF( GEOM_SELEC(3) ) THEN
            IF( .NOT. SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
              ! Initialization
              OVERLAP_PRELIMINAR = .FALSE.
              ! Preliminary test (circumscribing spherocylinders)
              CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
              ! Overlap criterion
              IF( OVERLAP_PRELIMINAR ) THEN
                ! Retrive position of the particle j after applying the PBC
                RJ(1) = RI(1) + RIJ(1)
                RJ(2) = RI(2) + RIJ(2)
                RJ(3) = RI(3) + RIJ(3)
                ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                ! Overlap criterion
                IF( OVERLAP ) THEN
                  ! Overlap detected
                  IF( COMPONENTS > 1 ) THEN
                    WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                    &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
                  ELSE IF( COMPONENTS == 1 ) THEN
                    WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, &
                    &                   ". Exiting..."
                  END IF
                  CALL SLEEP( 1 )
                  CALL EXIT(  )
                END IF
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( .NOT. SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
              ! Retrive position of the particle j after applying the PBC
              RJ(1) = RI(1) + RIJ(1)
              RJ(2) = RI(2) + RIJ(2)
              RJ(3) = RI(3) + RIJ(3)
              CALL CYLINDERSPHERE_OVERLAP( CI, CJ, QI, RI, RJ, OVERLAP )
              ! Overlap criterion
              IF( OVERLAP ) THEN
                ! Overlap detected
                IF( COMPONENTS > 1 ) THEN
                  WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                  &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
                ELSE IF( COMPONENTS == 1 ) THEN
                  WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, &
                  &                   ". Exiting..."
                END IF
                CALL SLEEP( 1 )
                CALL EXIT(  )
              END IF
            ! Overlap test for cylinders and spheres
            ELSE IF( SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
              ! Retrive position of the particle j after applying the PBC
              RJ(1) = RI(1) + RIJ(1)
              RJ(2) = RI(2) + RIJ(2)
              RJ(3) = RI(3) + RIJ(3)
              CALL CYLINDERSPHERE_OVERLAP( CJ, CI, QJ, RJ, RI, OVERLAP )
              ! Overlap criterion
              IF( OVERLAP ) THEN
                ! Overlap detected
                IF( COMPONENTS > 1 ) THEN
                  WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                  &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
                ELSE IF( COMPONENTS == 1 ) THEN
                  WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, &
                  &                   ". Exiting..."
                END IF
                CALL SLEEP( 1 )
                CALL EXIT(  )
              END IF
            ! Overlap test for spheres
            ELSE IF( SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
              ! Overlap detected
              IF( COMPONENTS > 1 ) THEN
                WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
              ELSE IF( COMPONENTS == 1 ) THEN
                WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
              END IF
              CALL SLEEP( 1 )
              CALL EXIT(  )
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
      ! Minimum image convention
      CALL MULTI_MATRIX( BOX_LENGTH_I, RIJ, S12 )
      S12 = S12 - ANINT(S12)
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
            IF( COMPONENTS > 1 ) THEN
              WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
              &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
            ELSE IF( COMPONENTS == 1 ) THEN
              WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
            END IF
            CALL SLEEP( 1 )
            CALL EXIT(  )
          END IF
        ! Overlap test for spherocylinders (Vega-Lago Method)
        ELSE IF( GEOM_SELEC(2) ) THEN
          CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
          ! Overlap criterion
          IF( OVERLAP ) THEN
            ! Overlap detected
            IF( COMPONENTS > 1 ) THEN
              WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
              &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
            ELSE IF( COMPONENTS == 1 ) THEN
              WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
            END IF
            CALL SLEEP( 1 )
            CALL EXIT(  )
          END IF
        ! Overlap test for cylinders and/or spheres
        ELSE IF( GEOM_SELEC(3) ) THEN
          IF( .NOT. SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
            ! Initialization
            OVERLAP_PRELIMINAR = .FALSE.
            ! Preliminary test (circumscribing spherocylinders)
            CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
            ! Overlap criterion
            IF( OVERLAP_PRELIMINAR ) THEN
              ! Retrive position of the particle j after applying the PBC
              RJ(1) = RI(1) + RIJ(1)
              RJ(2) = RI(2) + RIJ(2)
              RJ(3) = RI(3) + RIJ(3)
              ! Overlap test for cylinders (modified algorithm of Lopes et al.)
              CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
              ! Overlap criterion
              IF( OVERLAP ) THEN
                ! Overlap detected
                IF( COMPONENTS > 1 ) THEN
                  WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                  &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
                ELSE IF( COMPONENTS == 1 ) THEN
                  WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
                END IF
                CALL SLEEP( 1 )
                CALL EXIT(  )
              END IF
            END IF
            ! Overlap test for cylinders and spheres
          ELSE IF( .NOT. SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
            ! Retrive position of the particle j after applying the PBC
            RJ(1) = RI(1) + RIJ(1)
            RJ(2) = RI(2) + RIJ(2)
            RJ(3) = RI(3) + RIJ(3)
            CALL CYLINDERSPHERE_OVERLAP( CI, CJ, QI, RI, RJ, OVERLAP )
            ! Overlap criterion
            IF( OVERLAP ) THEN
              ! Overlap detected
              IF( COMPONENTS > 1 ) THEN
                WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
              ELSE IF( COMPONENTS == 1 ) THEN
                WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, &
                &                   ". Exiting..."
              END IF
              CALL SLEEP( 1 )
              CALL EXIT(  )
            END IF
          ! Overlap test for cylinders and spheres
          ELSE IF( SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
            ! Retrive position of the particle j after applying the PBC
            RJ(1) = RI(1) + RIJ(1)
            RJ(2) = RI(2) + RIJ(2)
            RJ(3) = RI(3) + RIJ(3)
            CALL CYLINDERSPHERE_OVERLAP( CJ, CI, QJ, RJ, RI, OVERLAP )
            ! Overlap criterion
            IF( OVERLAP ) THEN
              ! Overlap detected
              IF( COMPONENTS > 1 ) THEN
                WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
                &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
              ELSE IF( COMPONENTS == 1 ) THEN
                WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, &
                &                   ". Exiting..."
              END IF
              CALL SLEEP( 1 )
              CALL EXIT(  )
            END IF
          ! Overlap test for spheres
          ELSE IF( SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
            ! Overlap detected
            IF( COMPONENTS > 1 ) THEN
              WRITE( *, "(9G0)" ) "Overlap found in the initial configuration between particles ", I, " of component ", &
              &                   CI, " and ", J, " of component ", CJ, ". Exiting..."
            ELSE IF( COMPONENTS == 1 ) THEN
              WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", I, " and ", J, ". Exiting..."
            END IF
            CALL SLEEP( 1 )
            CALL EXIT(  )
          END IF
        END IF
      END IF
    END DO
  END DO
END DO

! Status
WRITE( *, "(G0)" ) "No overlaps found in the initial configuration. Resuming..."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Start simulation timer
CALL CPU_TIME( START_TIMER )

! Active transformation
DO I = 1, N_PARTICLES
  CALL ACTIVE_TRANSFORMATION( AXISZ, Q(:,I), E(:,I) )
END DO

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 20 )//"FILE ORGANIZER"//REPEAT( " ", 21 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating files..."

! *********************************************************************************************** !
! Output file units                                                                               !
! *********************************************************************************************** !

! Trajectory file (depends on user's choice)
IF( TRAJ_CHECK ) THEN
  IF( MC_ENSEMBLE == "NVT" ) THEN
    OPEN( UNIT= 20, FILE= "Trajectories/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_traj_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
    OPEN( UNIT= 20, FILE= "Trajectories/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
    &                     "_traj_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
    &                     //TRIM( DESCRIPTOR_FILE3 )//".xyz" )
  END IF
  WRITE( 20, "(G0)" ) N_PARTICLES
  WRITE( 20, * ) " "
  IF( GEOM_SELEC(1) ) THEN
    DO C = 1, COMPONENTS
      IF( .NOT. SPHERCOMP(C) ) THEN
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
        END DO
      ELSE
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
        END DO
      END IF
    END DO
  ELSE IF( GEOM_SELEC(2) ) THEN
    DO C = 1, COMPONENTS
      IF( .NOT. SPHERCOMP(C) ) THEN
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
        END DO
      ELSE
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
        END DO
      END IF
    END DO
  ELSE IF( GEOM_SELEC(3) ) THEN
    DO C = 1, COMPONENTS
      IF( .NOT. SPHERCOMP(C) ) THEN
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
        END DO
      ELSE
        DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
          WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), R(1,I), R(2,I), R(3,I), Q(0,I), Q(1,I), Q(2,I), Q(3,I), &
          &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
        END DO
      END IF
    END DO
  END IF
  FLUSH( 20 )
END IF

! Ratio file (translation)
IF( MC_ENSEMBLE == "NVT" ) THEN
  OPEN( UNIT= 30, FILE= "Ratio/Translation/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 30, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å]"//'"', ",", &
  &                    '"'//"Acceptance Ratio Threshold"//'"'
  FLUSH( 30 )
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 30, FILE= "Ratio/Translation/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 30, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å]"//'"', ",", &
  &                    '"'//"Acceptance Ratio Threshold"//'"'
  FLUSH( 30 )
END IF

! Ratio file (rotation)
IF( MC_ENSEMBLE == "NVT" ) THEN
  OPEN( UNIT= 40, FILE= "Ratio/Rotation/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 40, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [rad]"//'"', ",", &
  &                    '"'//"Acceptance Ratio Threshold"//'"'
  FLUSH( 40 )
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 40, FILE= "Ratio/Rotation/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 40, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [rad]"//'"', ",", &
  &                    '"'//"Acceptance Ratio Threshold"//'"'
  FLUSH( 40 )
END IF

! Ratio file (volume change)
IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 50, FILE= "Ratio/Volume/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 50, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å³]"//'"', ",", &
  &                    '"'//"Acceptance Ratio Threshold"//'"', ",", '"'//"Type of Volume Change"//'"'
  FLUSH( 50 )
END IF

! Ratio file (box properties)
IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 55, FILE= "Ratio/Box/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_ratio_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 55, "(23G0)" ) '"'//"Cycles"//'"', ",", '"'//"Box Distortion"//'"', ",", &
  &                     '"'//"Box Volume [Å³]"//'"', ",", '"'//"Box Length 1 [Å]"//'"', ",", &
  &                     '"'//"Box Length 2 [Å]"//'"', ",", '"'//"Box Length 3 [Å]"//'"', ",", &
  &                     '"'//"Box Length 4 [Å]"//'"', ",", '"'//"Box Length 5 [Å]"//'"', ",", &
  &                     '"'//"Box Length 6 [Å]"//'"', ",", '"'//"Box Length 7 [Å]"//'"', ",", &
  &                     '"'//"Box Length 8 [Å]"//'"', ",", '"'//"Box Length 9 [Å]"//'"'
  FLUSH( 55 )
END IF

! Order parameter file
IF( MC_ENSEMBLE == "NVT" ) THEN
  OPEN( UNIT= 60, FILE= "Order_Parameter/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_order_η"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 60, "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Nematic Order Parameter"//'"'
  FLUSH( 60 )
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 60, FILE= "Order_Parameter/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_order_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 60, "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Nematic Order Parameter"//'"'
  FLUSH( 60 )
END IF

! Results file
IF( MC_ENSEMBLE == "NPT" ) THEN
  OPEN( UNIT= 70, FILE= "Results/"//TRIM( DESCRIPTOR_DATE )//"/"//TRIM( DESCRIPTOR_HOUR )// &
  &                     "_results_P"//TRIM( DESCRIPTOR_FILE1 )//"_C"//TRIM( DESCRIPTOR_FILE2 )//"_" &
  &                     //TRIM( DESCRIPTOR_FILE3 )//".dat" )
  WRITE( 70, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Packing Fraction"//'"', ",", '"'//"Number Density [Å⁻³]"//'"', ",", &
  &                    '"'//"Box Volume"//'"', ",", '"'//"Reduced Pressure"//'"'
  FLUSH( 70 )
END IF

! Potential file
IF( POTENTIAL_SELEC(2) ) THEN
  DO C_LAMB = 1, N_LAMBDA
    WRITE( DESCRIPTOR_LAMB, FORMAT_LAMB ) L_RANGE(C_LAMB)
    ! Potential files
    OPEN( UNIT= (80 + C_LAMB), FILE= "Potential/"//TRIM( DESCRIPTOR_DATE )//"/Lambda_"//TRIM( DESCRIPTOR_LAMB )//"/"&
    &                                //TRIM( DESCRIPTOR_HOUR )//"_thermo_η"//TRIM( DESCRIPTOR_FILE1 )//"_C" &
    &                                //TRIM( DESCRIPTOR_FILE2 )//"_"//TRIM( DESCRIPTOR_FILE3 )//".dat" )
    WRITE( (80 + C_LAMB), "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Potential Energy"//'"'
    FLUSH( (80 + C_LAMB) )
  END DO
END IF

! Status
CALL SLEEP( 1 )
WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
WRITE( *, "(G0)" ) " "

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 16 )//"MONTE CARLO SIMULATION"//REPEAT( " ", 17 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
IF( MC_ENSEMBLE == "NVT" ) THEN
  WRITE( *, "(G0)" ) "Starting up the NVT-Monte Carlo simulation..."
ELSE IF( MC_ENSEMBLE == "NPT" ) THEN
  WRITE( *, "(G0)" ) "Starting up the NPT-Monte Carlo simulation..."
END IF
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Computation of total potential energy (initial configuration)
IF( POTENTIAL_SELEC(2) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Computing total potential energy of the initial configuration..."
  CALL COMPUTE_TOTAL_ENERGY(  )
  ! Status
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Pseudorandom number generator seed
IF( FSEED ) THEN
  SEED = 123456789
ELSE IF( .NOT. FSEED ) THEN
  CALL RANDOM_NUMBER( RANDOM_N )
  SEED_COMP = INT( RANDOM_N * DBLE( SIZE_SEED ) ) + 1
  SEED = ABS( RANS(SEED_COMP) )
END IF

! Initialize box distortion parameter
CALL LATTICE_REDUCTION( BOX_LENGTH, DISTORTION, LATTICER )

! *********************************************************************************************** !
! Monte Carlo parameters                                                                          !
! *********************************************************************************************** !
MOV_TRANS    = .FALSE.         ! Translational move selector                       (initial value)
MOV_ROT      = .FALSE.         ! Rotational move selector                          (initial value)
MOV_VOL_I    = .FALSE.         ! Volume move selector (Isotropic)                  (initial value)
MOV_VOL_A    = .FALSE.         ! Volume move selector (Anisotropic)                (initial value)
DRMAX        = MAX_TRANS       ! Maximum translational displacement                (initial value)
ANGMAX       = MAX_ROT         ! Maximum rotational displacement                   (initial value)
DVMAXISO     = MAX_VOLI        ! Maximum isovolumetric displacement                (initial value)
DVMAXANI     = MAX_VOLA        ! Maximum anisovolumetric displacement              (initial value)
NACCT        = 0               ! Translational move acceptance counter             (initial value)
NACCR        = 0               ! Rotational move acceptance counter                (initial value)
NACCVI       = 0               ! Volumetric move acceptance counter (Isotropic)    (initial value)
NACCVA       = 0               ! Volumetric move acceptance counter (Anisotropic)  (initial value)
MOVT         = 0               ! Translational move counter                        (initial value)
MOVR         = 0               ! Rotational move counter                           (initial value)
MOVVI        = 0               ! Volume change counter (Isotropic)                 (initial value)
MOVVA        = 0               ! Volume change counter (Anisotropic)               (initial value)
QMC(:,:)     = Q(:,:)          ! Rotation quaternions                              (initial value)
RMC(:,:)     = R(:,:)          ! Position of particles                             (initial value)
EMC(:,:)     = E(:,:)          ! Orientation of particles                          (initial value)
VMC(:)       = V(:)            ! Total potential energy                            (initial value)
BOXLMC(:)    = BOX_LENGTH(:)   ! Box length                                        (initial value)
BOXLMC_I(:)  = BOX_LENGTH_I(:) ! Box length (inverse)                              (initial value)
BOXVMC       = BOX_VOLUME      ! Box volume                                        (initial value)
SCALE_FACTOR = 1.D0            ! Scaling factor                                    (initial value)
LATTICER     = .FALSE.         ! Lattice reduction                                 (initial value)

! Metropolis Algorithm - Importance Sampling
WRITE( *, "(G0)" ) "Running Metropolis algorithm..."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! *********************************************************************************************** !
! Simulation cycles                                                                               !
! *********************************************************************************************** !
!  A 'cycle' is characterized by N trial moves (rotation or translation) of a random particle or  !
!  by a single change of the simulation box volume.                                               !
! *********************************************************************************************** !
DO CYCLES = 1, MAX_CYCLES

  ! Simulation progress
  CALL PROGRESS_BAR_MC( CYCLES, MAX_CYCLES, MC_ENSEMBLE )

  ! Save all potential data (equilibration and production)
  IF( POTENTIAL_SELEC(2) ) THEN
    IF( .NOT. POTENTIAL_CHECK ) THEN
      IF( MOD( CYCLES, N_SAVE ) == 0 ) THEN
        DO C_LAMB = 1, N_LAMBDA
          WRITE( (80 + C_LAMB), "(3G0)" ) CYCLES, ",", VMC(C_LAMB)
          FLUSH( (80 + C_LAMB) )
        END DO
      END IF
    ! Save potential data (production-related only)
    ELSE IF( POTENTIAL_CHECK ) THEN
      IF( CYCLES > N_EQUIL .AND. MOD( CYCLES, N_SAVE ) == 0 ) THEN
        DO C_LAMB = 1, N_LAMBDA
          WRITE( (80 + C_LAMB), "(3G0)" ) CYCLES, ",", VMC(C_LAMB)
          FLUSH( (80 + C_LAMB) )
        END DO
      END IF
    END IF
  END IF

  ! Generates a random number for the NPT-simulation
  IF( MC_ENSEMBLE == "NPT" ) THEN
    CALL RANF(  )
  END IF

  ! Movement (Translation or Rotation)
  IF( RANDOM_N <= PROB_MOV .OR. MC_ENSEMBLE == "NVT" ) THEN

    ! Disable volume change
    MOV_VOL_I = .FALSE.
    MOV_VOL_A = .FALSE.

    ! Particle loop
    DO K = 1, N_PARTICLES

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
      IF( SPHERCOMP(CI) ) THEN
        MOV_TRANS = .TRUE.   ! Enable translation
        MOV_ROT   = .FALSE.  ! Disable rotation
        MOVT      = MOVT + 1 ! Increment move counter
      ! Allow rotation if component is nonspherical
      ELSE
        ! Pseudorandom number generator (uniform distribution)
        CALL RANF(  )
        ! Translation criterion
        IF( RANDOM_N < PROB_TRANS ) THEN
          MOV_TRANS = .TRUE.   ! Enable translation
          MOV_ROT   = .FALSE.  ! Disable rotation
          MOVT      = MOVT + 1 ! Increment move counter
        ! Rotation criterion
        ELSE IF( RANDOM_N >= PROB_TRANS ) THEN
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

      ! Translational movement
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
      ! No translation
      ELSE IF( .NOT. MOV_TRANS ) THEN
        RN(:) = RM(:)
      END IF

      ! Rotational movement
      IF( MOV_ROT ) THEN
        ! Random Composed Unit Quaternion
        CALL COMPOSED_QUATERNION( QM, QN, ANGMAX )
        ! Active transformation
        CALL ACTIVE_TRANSFORMATION( AXISZ, QN, EN )
      ! No rotation
      ELSE IF( .NOT. MOV_ROT ) THEN
        QN(:) = QM(:)
        EN(:) = EM(:)
      END IF

      ! Overlap check
      CALL CHECK_OVERLAP( CI, I, QN, EN, RN, CD, BOXLMC, BOXLMC_I, OVERLAP )

      ! Acceptance criterion
      IF( .NOT. OVERLAP ) THEN
        ! System configuration update
        RMC(:,I) = RN(:) ! Update position
        QMC(:,I) = QN(:) ! Update quaternion
        EMC(:,I) = EN(:) ! Update orientation
        ! Update total potential energy
        IF( POTENTIAL_SELEC(2) ) THEN
          ! Computation of potential energy of particle i (microstate m)
          CALL COMPUTE_PARTICLE_ENERGY( CI, I, RM, VM, BOXLMC, BOXLMC_I )
          ! Computation of potential energy of particle i (microstate n)
          CALL COMPUTE_PARTICLE_ENERGY( CI, I, RN, VN, BOXLMC, BOXLMC_I )
          ! Computation of energy difference of microstates n and m
          DV(:)  = VN(:) - VM(:)
          ! System energy update
          VMC(:) = VMC(:) + DV(:)
        END IF
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

  ! Volume change
  ELSE IF( RANDOM_N > PROB_MOV .AND. MC_ENSEMBLE == "NPT" ) THEN

    ! Disable translation and rotation
    MOV_TRANS = .FALSE.
    MOV_ROT   = .FALSE.

    ! Assignment of previous configuration
    BOXLM(:)   = BOXLMC(:)   ! Box length
    BOXLM_I(:) = BOXLMC_I(:) ! Box length (inverse)
    BOXVM      = BOXVMC      ! Box volume

    ! Expansion/compression type
    CALL RANF(  )

    ! Isotropic volume change
    IF( RANDOM_N <= PROB_VOL_ISO ) THEN
      ! Random scaling factor
      CALL RANF(  )
      SCALE_FACTOR = 1.0D0 + DVMAXISO * (RANDOM_N - 0.5D0)
      ! Proportional box length
      BOXLN(:) = BOXLM(:) * SCALE_FACTOR
      VOL_TYPE = "IS"
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL INVERSE_COF( BOXLN, BOXLN_I, BOXVN )
      ! Movement counter
      MOVVI = MOVVI + 1
      ! Movement type
      MOV_VOL_I = .TRUE.  ! Enable isotropic volume change
      MOV_VOL_A = .FALSE. ! Disable anisotropic volume change
    ! Anisotropic volume change
    ELSE IF( RANDOM_N > PROB_VOL_ISO ) THEN
      ! Random box component
      CALL RANF(  )
      COMPONENT = INT( RANDOM_N * 6.D0 ) + 1
      IF( COMPONENT == 1 ) THEN
        COMPONENT = 1 ! xx component
      ELSE IF( COMPONENT == 2 ) THEN
        COMPONENT = 4 ! xy component
      ELSE IF( COMPONENT == 3 ) THEN
        COMPONENT = 5 ! yy component
      ELSE IF( COMPONENT == 4 ) THEN
        COMPONENT = 7 ! xz component
      ELSE IF( COMPONENT == 5 ) THEN
        COMPONENT = 8 ! yz component
      ELSE IF( COMPONENT == 6 ) THEN
        COMPONENT = 9 ! zz component
      END IF
      BOXLN(:) = BOXLM(:)
      ! Random stretching/shortening of the box length component
      CALL RANF(  )
      BOXLN(COMPONENT) = BOXLM(COMPONENT) + DVMAXANI * (RANDOM_N - 0.5D0)
      VOL_TYPE = "AN"
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL INVERSE_COF( BOXLN, BOXLN_I, BOXVN )
      ! Movement counter
      MOVVA = MOVVA + 1
      ! Movement type
      MOV_VOL_I = .FALSE. ! Disable isotropic volume change
      MOV_VOL_A = .TRUE.  ! Enable anisotropic volume change
    END IF

    ! Reset condition of anisotropic volume change
    IGNORE_VOL_ATTEMPT = .FALSE.

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
        IF( COSANGLE_VEC(L) < DCOS( (PI / 2.D0) + MAX_ANGLE ) .OR. COSANGLE_VEC(L) > DCOS( (PI / 2.D0) - MAX_ANGLE ) ) THEN
          BOXVMC      = BOXVM
          BOXLMC(:)   = BOXLM(:)
          BOXLMC_I(:) = BOXLM_I(:)
          IGNORE_VOL_ATTEMPT = .TRUE.
          EXIT
        END IF
        ! Length distortion
        IF( LBOXR(L) > MAX_LENGTH_RATIO .OR. LBOXR(L) < 1.D0 / MAX_LENGTH_RATIO ) THEN
          BOXVMC      = BOXVM
          BOXLMC(:)   = BOXLM(:)
          BOXLMC_I(:) = BOXLM_I(:)
          IGNORE_VOL_ATTEMPT = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    ! Box not too distorted
    IF( .NOT. IGNORE_VOL_ATTEMPT ) THEN

      ! Enthalpy (weighing function)
      HNM = ( PRESS * ( BOXVN - BOXVM ) ) - ( DBLE( N_PARTICLES ) * DLOG( BOXVN / BOXVM ) )

      ! Random number
      CALL RANF(  )

      ! Enthalpy Criterion
      IF( DEXP( - HNM ) >= RANDOM_N ) THEN

        ! System configuration update
        RMCV(:,:) = RMC(:,:) ! Old configuration

        ! Rescale positions of all particles accordingly
        IF( MOV_VOL_I ) THEN
          ! Isotropic volume change
          DO K = 1, N_PARTICLES
            RMC(:,K) = RMC(:,K) * SCALE_FACTOR
          END DO
        ELSE IF( MOV_VOL_A ) THEN
          ! Anisotropic volume change
          DO K = 1, N_PARTICLES
            ! Scaling coordinates using the old box length
            CALL MULTI_MATRIX( BOXLM_I, RMC(:,K), S12 )
            ! New real coordinates using the new box length
            CALL MULTI_MATRIX( BOXLN, S12, RMC(:,K) )
          END DO
        END IF

        ! Overlap check after expansion/compression of the simulation box
        LOOP_ALLOVERLAP_NPT: DO

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
                    ! Overlap test for ellipsoids of revolution (Perram-Wertheim method)
                    IF( GEOM_SELEC(1) ) THEN
                      CALL ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CD, OVERLAP )
                      ! Overlap criterion
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_ALLOVERLAP_NPT
                      END IF
                    ! Overlap test for spherocylinders (Vega-Lago method)
                    ELSE IF( GEOM_SELEC(2) ) THEN
                      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                      ! Overlap criterion
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_ALLOVERLAP_NPT
                      END IF
                    ! Overlap test for cylinders and/or spheres
                    ELSE IF( GEOM_SELEC(3) ) THEN
                      IF( .NOT. SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
                        ! Initialization
                        OVERLAP_PRELIMINAR = .FALSE.
                        ! Preliminary test (circumscribing spherocylinders)
                        CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                        ! Overlap criterion
                        IF( OVERLAP_PRELIMINAR ) THEN
                          ! Retrive position of the particle j after applying the PBC
                          RJ(1) = RI(1) + RIJ(1)
                          RJ(2) = RI(2) + RIJ(2)
                          RJ(3) = RI(3) + RIJ(3)
                          ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                          CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                          ! Overlap criterion
                          IF( OVERLAP ) THEN
                            ! Overlap detected
                            EXIT LOOP_ALLOVERLAP_NPT
                          END IF
                        END IF
                      ! Overlap test for cylinders and spheres
                      ELSE IF( .NOT. SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
                        ! Retrive position of the particle j after applying the PBC
                        RJ(1) = RI(1) + RIJ(1)
                        RJ(2) = RI(2) + RIJ(2)
                        RJ(3) = RI(3) + RIJ(3)
                        CALL CYLINDERSPHERE_OVERLAP( CI, CJ, QI, RI, RJ, OVERLAP )
                        IF( OVERLAP ) THEN
                          ! Overlap detected
                          EXIT LOOP_ALLOVERLAP_NPT
                        END IF
                      ! Overlap test for cylinders and spheres
                      ELSE IF( SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
                        ! Retrive position of the particle j after applying the PBC
                        RJ(1) = RI(1) + RIJ(1)
                        RJ(2) = RI(2) + RIJ(2)
                        RJ(3) = RI(3) + RIJ(3)
                        CALL CYLINDERSPHERE_OVERLAP( CJ, CI, QJ, RJ, RI, OVERLAP )
                        IF( OVERLAP ) THEN
                          ! Overlap detected
                          EXIT LOOP_ALLOVERLAP_NPT
                        END IF
                      ! Overlap test for spheres
                      ELSE IF( SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
                        ! Overlap detected
                        OVERLAP = .TRUE.
                        EXIT LOOP_ALLOVERLAP_NPT
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
                      EXIT LOOP_ALLOVERLAP_NPT
                    END IF
                  ! Overlap test for spherocylinders (Vega-Lago Method)
                  ELSE IF( GEOM_SELEC(2) ) THEN
                    CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP )
                    ! Overlap criterion
                    IF( OVERLAP ) THEN
                      ! Overlap detected
                      EXIT LOOP_ALLOVERLAP_NPT
                    END IF
                  ! Overlap test for cylinders (Lopes et al. Method)
                  ELSE IF( GEOM_SELEC(3) ) THEN
                    IF( .NOT. SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
                      ! Initialization
                      OVERLAP_PRELIMINAR = .FALSE.
                      ! Preliminary test (circumscribing spherocylinders)
                      CALL SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, CD, PARALLEL, OVERLAP_PRELIMINAR )
                      ! Overlap criterion
                      IF( OVERLAP_PRELIMINAR ) THEN
                        ! Retrive position of the particle j after applying the PBC
                        RJ(1) = RI(1) + RIJ(1)
                        RJ(2) = RI(2) + RIJ(2)
                        RJ(3) = RI(3) + RIJ(3)
                        ! Overlap test for cylinders (modified algorithm of Lopes et al.)
                        CALL CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, PARALLEL, OVERLAP )
                        ! Overlap criterion
                        IF( OVERLAP ) THEN
                          ! Overlap detected
                          EXIT LOOP_ALLOVERLAP_NPT
                        END IF
                      END IF
                    ! Overlap test for cylinders and spheres
                    ELSE IF( .NOT. SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
                      ! Retrive position of the particle j after applying the PBC
                      RJ(1) = RI(1) + RIJ(1)
                      RJ(2) = RI(2) + RIJ(2)
                      RJ(3) = RI(3) + RIJ(3)
                      CALL CYLINDERSPHERE_OVERLAP( CI, CJ, QI, RI, RJ, OVERLAP )
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_ALLOVERLAP_NPT
                      END IF
                    ! Overlap test for cylinders and spheres
                    ELSE IF( SPHERCOMP(CI) .AND. .NOT. SPHERCOMP(CJ) ) THEN
                      ! Retrive position of the particle j after applying the PBC
                      RJ(1) = RI(1) + RIJ(1)
                      RJ(2) = RI(2) + RIJ(2)
                      RJ(3) = RI(3) + RIJ(3)
                      CALL CYLINDERSPHERE_OVERLAP( CJ, CI, QJ, RJ, RI, OVERLAP )
                      IF( OVERLAP ) THEN
                        ! Overlap detected
                        EXIT LOOP_ALLOVERLAP_NPT
                      END IF
                    ! Overlap test for spheres
                    ELSE IF( SPHERCOMP(CI) .AND. SPHERCOMP(CJ) ) THEN
                      ! Overlap detected
                      OVERLAP = .TRUE.
                      EXIT LOOP_ALLOVERLAP_NPT
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END DO

          ! No overlaps
          EXIT LOOP_ALLOVERLAP_NPT

        END DO LOOP_ALLOVERLAP_NPT

        ! Acceptance criterion
        IF( .NOT. OVERLAP ) THEN
          ! System configuration update
          BOXVMC      = BOXVN      ! Update volume
          BOXLMC(:)   = BOXLN(:)   ! Update length
          BOXLMC_I(:) = BOXLN_I(:) ! Update length (inverse)
          ! Displacement counter update
          IF( MOV_VOL_I ) THEN
            NACCVI = NACCVI + 1 ! Isotropic move counter
          ELSE IF( MOV_VOL_A ) THEN
            NACCVA = NACCVA + 1 ! Anisotropic move counter
          END IF
          ! Update packing fraction and reduced number density
          PACKING_F = ( TOTAL_VP / BOXVN )
          TOTAL_RHO = ( DBLE( N_PARTICLES ) / BOXVN )
          ! Re-initialization
          IGNORE_VOL_ATTEMPT = .FALSE.
          ! Lattice reduction
          LATTICER = .FALSE.
          CALL LATTICE_REDUCTION( BOXLMC, DISTORTION, LATTICER )
          IF( LATTICER ) THEN
            ! Calculate the new reciprocal box basis vectors
            CALL INVERSE_COF( BOXLMC, BOXLMC_I, BOXVMC )
            DO K = 1, N_PARTICLES
              ! Minimum image convention
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
                BOXVROT    = BOXVMC
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
              BOXVMC      = BOXVROT
              ! Update molecular properties
              RMC(:,:) = RPROT(:,:)
              QMC(:,:) = QPROT(:,:)
              EMC(:,:) = EPROT(:,:)
            END IF
          END IF
        ! Retrieve old properties of the system configuration
        ELSE IF( OVERLAP ) THEN
          BOXVMC      = BOXVM      ! Retrieve box volume
          BOXLMC(:)   = BOXLM(:)   ! Retrieve box length
          BOXLMC_I(:) = BOXLM_I(:) ! Retrieve box length (inverse)
          RMC(:,:)    = RMCV(:,:)  ! Retrieve position of particles
        END IF

      ! Retrieve old properties of the simulation box
      ELSE

        BOXVMC      = BOXVM      ! Retrieve box volume
        BOXLMC(:)   = BOXLM(:)   ! Retrieve box length
        BOXLMC_I(:) = BOXLM_I(:) ! Retrieve box length (inverse)

      END IF ! Enthalpy criterion

    END IF ! Box distortion criterion

  END IF

  ! Adjustment of maximum displacements
  IF( CYCLES <= N_EQUIL ) THEN  ! During equilibration only

    ! Adjustment of maximum displacement (translation)
    IF( DBLE( MOVT ) >= DBLE( N_ADJUST * N_PARTICLES ) * R_ACC_T ) THEN

      ! Acceptance ratio (translation)
      IF( MOVT > 0 ) THEN
        RATIO = DBLE( NACCT ) / DBLE( MOVT )
        ! Translational adjustment
        IF( RATIO <= R_ACC_T ) THEN
          DRMAX = 0.95D0 * DRMAX
        ELSE
          DRMAX = 1.05D0 * DRMAX
        END IF
        ! Ratio data (translation)
        WRITE( 30, "(7G0)" ) CYCLES, ",", RATIO, ",", DRMAX, ",", R_ACC_T
        FLUSH( 30 )
        ! Reset counter
        NACCT = 0
        MOVT  = 0
      END IF

      ! Avoid large translations
      IF( DRMAX > 5.D0 * MAXVAL( BOXLMC ) ) THEN
        DRMAX = DRMAX - 2.5D0 * MAXVAL( BOXLMC )
      END IF

    END IF

    ! Adjustment of maximum displacement (rotation)
    IF( DBLE( MOVR ) >= DBLE( N_ADJUST * N_PARTICLES ) * R_ACC_R ) THEN

      ! Acceptance ratio (rotation)
      IF( MOVR > 0 ) THEN
        RATIO = DBLE( NACCR ) / DBLE( MOVR )
        ! Rotational adjustment
        IF( RATIO <= R_ACC_R ) THEN
          ANGMAX = 0.95D0 * ANGMAX
        ELSE
          ANGMAX = 1.05D0 * ANGMAX
        END IF
        ! Ratio data (rotation)
        WRITE( 40, "(7G0)" ) CYCLES, ",", RATIO, ",", ANGMAX, ",", R_ACC_R
        FLUSH( 40 )
        ! Reset counter
        NACCR = 0
        MOVR  = 0
      END IF

      ! Avoid 4π-rotations
      IF( ANGMAX > 4.D0 * PI ) THEN
        ANGMAX = ANGMAX - 2.D0 * PI
      END IF

    END IF

    ! Adjustment of maximum displacement (isotropic volume change)
    IF( DBLE( MOVVI ) >= DBLE( N_ADJUST ) * R_ACC_VI ) THEN

      IF( MOVVI > 0 ) THEN
        ! Acceptance ratio (isotropic volume change)
        RATIO = DBLE( NACCVI ) / DBLE( MOVVI )
        ! Volumetric adjustment
        IF( RATIO <= R_ACC_VI ) THEN
          DVMAXISO = 0.95D0 * DVMAXISO
        ELSE
          DVMAXISO = 1.05D0 * DVMAXISO
        END IF
        ! Ratio data (volume change)
        IF( MC_ENSEMBLE == "NPT" ) THEN
          WRITE( 50, "(9G0)" ) CYCLES, ",", RATIO, ",", DVMAXISO, ",", R_ACC_VI, ",", VOL_TYPE
          FLUSH( 50 )
        END IF
        ! Reset counter
        NACCVI = 0
        MOVVI  = 0
      END IF

      ! Ratio data (box properties)
      IF( MC_ENSEMBLE == "NPT" ) THEN
        WRITE( 55, "(23G0)" ) CYCLES, ",", DISTORTION, ",", BOXVMC, ",", BOXLMC(1), ",", BOXLMC(2), ",", BOXLMC(3), ",", &
        &                     BOXLMC(4), ",", BOXLMC(5), ",", BOXLMC(6), ",", BOXLMC(7), ",", BOXLMC(8), ",", BOXLMC(9)
        FLUSH( 55 )
      END IF

    END IF

    ! Adjustment of maximum displacement (anisotropic volume change)
    IF( DBLE( MOVVA ) >= DBLE( N_ADJUST ) * R_ACC_VA ) THEN

      IF( MOVVA > 0 ) THEN
        ! Acceptance ratio (isotropic volume change)
        RATIO = DBLE( NACCVA ) / DBLE( MOVVA )
        ! Volumetric adjustment
        IF( RATIO <= R_ACC_VA ) THEN
          DVMAXANI = 0.95D0 * DVMAXANI
        ELSE
          DVMAXANI = 1.05D0 * DVMAXANI
        END IF
        ! Ratio data (volume change)
        IF( MC_ENSEMBLE == "NPT" ) THEN
          WRITE( 50, "(9G0)" ) CYCLES, ",", RATIO, ",", DVMAXANI, ",", R_ACC_VA, ",", VOL_TYPE
          FLUSH( 50 )
        END IF
        ! Reset counter
        NACCVA = 0
        MOVVA  = 0
      END IF

      ! Ratio data (box properties)
      IF( MC_ENSEMBLE == "NPT" ) THEN
        WRITE( 55, "(23G0)" ) CYCLES, ",", DISTORTION, ",", BOXVMC, ",", BOXLMC(1), ",", BOXLMC(2), ",", BOXLMC(3), ",", &
        &                     BOXLMC(4), ",", BOXLMC(5), ",", BOXLMC(6), ",", BOXLMC(7), ",", BOXLMC(8), ",", BOXLMC(9)
        FLUSH( 55 )
      END IF

    END IF

  END IF

  ! Order parameter data
  IF( MOD( CYCLES, N_SAVE ) == 0 ) THEN
    ! Nematic order parameter (Q-tensor method)
    CALL NEMATIC_ORDER_PARAMETER( S, EMC )
    ! Only production-related data
    WRITE( 60, "(3G0)" ) CYCLES, ",", S
    FLUSH( 60 )
  END IF

  ! Trajectory data
  IF( TRAJ_CHECK ) THEN
    IF( MOD ( CYCLES, N_SAVE ) == 0 ) THEN
      WRITE( 20, "(G0)" ) N_PARTICLES
      WRITE( 20, * ) " "
      IF( GEOM_SELEC(1) ) THEN
        DO C = 1, COMPONENTS
          IF( .NOT. SPHERCOMP(C) ) THEN
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * LENGTH(C)
            END DO
          ELSE
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
            END DO
          END IF
        END DO
      ELSE IF( GEOM_SELEC(2) ) THEN
        DO C = 1, COMPONENTS
          IF( .NOT. SPHERCOMP(C) ) THEN
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
            END DO
          ELSE
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
            END DO
          END IF
        END DO
      ELSE IF( GEOM_SELEC(3) ) THEN
        DO C = 1, COMPONENTS
          IF( .NOT. SPHERCOMP(C) ) THEN
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), LENGTH(C)
            END DO
          ELSE
            DO I = SUM( N_COMPONENT(0:(C-1)) ) + 1, SUM( N_COMPONENT(0:C) )
              WRITE( 20, "(11(G0,1X))" ) INDEX_P(C), RMC(1,I), RMC(2,I), RMC(3,I), QMC(0,I), QMC(1,I), QMC(2,I), QMC(3,I), &
              &                          0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C), 0.5D0 * DIAMETER(C)
            END DO
          END IF
        END DO
      END IF
      FLUSH( 20 )
    END IF
  END IF

  ! Results data (packing fraction, number density, box volume, and pressure)
  IF( MOD( CYCLES, N_SAVE ) == 0 .AND. MC_ENSEMBLE == "NPT" ) THEN
    WRITE( 70, "(9G0)" ) CYCLES, ",", PACKING_F, ",", TOTAL_RHO, ",", BOXVMC, ",", PRESS
    FLUSH( 70 )
  END IF

END DO

! End of Metropolis algorithm
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Monte Carlo simulation finished successfully! See directories for results."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Output unit
IF( TRAJ_CHECK ) THEN
  CLOSE( 20 )
END IF
CLOSE( 30 )
CLOSE( 40 )
CLOSE( 60 )
IF( MC_ENSEMBLE == "NPT" ) THEN
  CLOSE( 50 )
  CLOSE( 55 )
  CLOSE( 70 )
END IF
IF( POTENTIAL_SELEC(2) ) THEN
  DO C_LAMB = 1, N_LAMBDA
    CLOSE ( 80 + C_LAMB )
  END DO
END IF

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"SIMULATION LENGTH"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR

! End simulation timer
CALL CPU_TIME( STOP_TIMER )
WRITE( *, "(G0,G0.5,G0)" ) "Elapsed Time: ", (STOP_TIMER - START_TIMER), "s."
CALL SLEEP( 1 )
WRITE( *, "(G0)" ) " "

! Deallocation of arrays
DEALLOCATE( QMC, RMC, EMC, RMCV )
DEALLOCATE( V, VMC, VN, VM, DV )
DEALLOCATE( SWRANGE )
DEALLOCATE( Q, R, E )
DEALLOCATE( QPROT, RPROT, EPROT )

! Calculation of the first- and second-order TPT coefficients
IF( COEF_CHECK .AND. MC_ENSEMBLE == "NVT" .AND. POTENTIAL_SELEC(2) ) THEN
  ! Status
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
  WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"TPT COEFFICIENTS"//REPEAT( " ", 20 )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR
  WRITE( *, "(G0)" ) "Initializing calculation of the TPT coefficients..."
  WRITE( *, "(G0)" ) " "
  ! Calculation of the TPT coefficients for every potential range
  BLOCK_AV_LOOP: DO C_LAMB = 1, N_LAMBDA
    ! Condition
    IF( ( MIN_BLOCKS > MAX_BLOCKS ) ) THEN
      WRITE( *, "(5G0)" ) "Minimum number of blocks (", MIN_BLOCKS, ") greater than the maximum number of blocks (", MAX_BLOCKS, &
      &                   "). Skipping..."
      WRITE( *, "(G0)" ) " "
      ! Terminate calculation if min. block greater than max. blocks
      EXIT BLOCK_AV_LOOP
    END IF
    ! Status
    WRITE( *, "(G0,G0.5,G0)", ADVANCE= "NO" ) "Current Potential Range: ", L_RANGE(C_LAMB), ". Calculating..."
    WRITE( DESCRIPTOR_LAMB, FORMAT_LAMB ) L_RANGE(C_LAMB)
    ! Initialization
    FLAG = .FALSE.
    ! Block-averaging method
    CALL BLOCK_AVERAGING( FLAG, A1(C_LAMB), A2(C_LAMB), APERT(C_LAMB), DPA1(C_LAMB), DPA2(C_LAMB), DPAPERT(C_LAMB), EXEC_TIME )
    ! Stop criterion
    IF( FLAG ) THEN
      WRITE( *, "(G0)", ADVANCE= "YES" ) " "
      WRITE( *, "(5G0)" ) "Number of blocks (", MAX_BLOCKS - MIN_BLOCKS, ") exceed the maximum number of blocks (", 10000, &
      &                   "). Skipping..."
      WRITE( *, "(G0)" ) " "
      ! Terminate calculation if min. block greater than max. blocks
      EXIT BLOCK_AV_LOOP
    END IF
    ! Status
    WRITE( *, "(G0,G0.5,G0)", ADVANCE= "YES" ) " Done in ", EXEC_TIME, "s."
  END DO BLOCK_AV_LOOP
  ! Subroutine returns no errors
  IF( .NOT. FLAG .AND. ( MIN_BLOCKS <= MAX_BLOCKS ) ) THEN
    DO C_LAMB = 1, N_LAMBDA
      WRITE ( DESCRIPTOR_LAMB, FORMAT_LAMB ) L_RANGE(C_LAMB)
      OPEN( UNIT= 150, FILE= "Perturbed_Coefficient/"//TRIM( DESCRIPTOR_DATE )//"/Lambda_"//TRIM( DESCRIPTOR_LAMB )//"/" &
      &                      //TRIM( DESCRIPTOR_HOUR )//"_TPT_coefficients_η"//TRIM( DESCRIPTOR_FILE1 )//"_C" &
      &                      //TRIM( DESCRIPTOR_FILE2 )//"_"//TRIM( DESCRIPTOR_FILE3 )//".csv")
      WRITE( 150, "(G0,G0)" ) "# of Particles: ", N_PARTICLES
      WRITE( 150, "(G0,G0.5)" ) "Reduced Temperature: ", RED_TEMP
      WRITE( 150, "(G0,G0.5)" ) "Reduced Number Density: ", TOTAL_RHO
      WRITE( 150, "(G0,G0.5)" ) "Aspect Ratio (L/D): ", ASPECT_RATIO
      WRITE( 150, "(G0,G0.5)" ) "Attractive Range: ", L_RANGE(C_LAMB)
      WRITE( 150, "(G0)" ) " "
      WRITE( 150, "(11G0)" ) "'1st_Order_TPTCoefficient'", ",", "'1st_Order_TPTCoefficient_STD'", ",", &
      &                      "'2nd_Order_TPTCoefficient'", ",", "'2nd_Order_TPTCoefficient_STD'", ",", &
      &                      "'Perturbed_Helmholtz_FEnergy'", ",", "'Perturbed_Helmholtz_FEnergy_STD'"
      WRITE( 150, "(G0.7,5(G0,G0.7))" ) A1(C_LAMB), ",", DPA1(C_LAMB), ",", A2(C_LAMB), ",", DPA2(C_LAMB), ",", APERT(C_LAMB), &
      &                                 ",", DPAPERT(C_LAMB)
      CLOSE( 150 )
    END DO
    WRITE( *, "(G0)", ADVANCE= "YES" ) " "
  END IF
END IF

! Write down results
WRITE( *, "(G0)", ADVANCE= "NO" ) "Writing simulation log..."
CALL SLEEP( 1 )

! Allocation
ALLOCATE( CHAR_LABEL(72,COMPONENTS) )
ALLOCATE( CHAR_LABEL_L(1,N_LAMBDA ) )

! Simulation log descriptors
WRITE( CHAR_LABEL(1,1), "(G0)"    ) COMPONENTS
WRITE( CHAR_LABEL(2,1), "(G0.5)"  ) PACKING_F
DO C = 1, COMPONENTS
  WRITE( CHAR_LABEL(3,C), "(G0.5)"  ) DIAMETER(C)
  WRITE( CHAR_LABEL(4,C), "(G0.5)"  ) LENGTH(C)
  WRITE( CHAR_LABEL(5,C), "(G0.5)"  ) ASPECT_RATIO(C)
  WRITE( CHAR_LABEL(6,C), "(G0.5)"  ) MOLAR_F(C)
  WRITE( CHAR_LABEL(7,C), "(G0)"    ) N_COMPONENT(C)
  WRITE( CHAR_LABEL(8,C), "(G0.5)"  ) PARTICLE_VOL(C)
  WRITE( CHAR_LABEL(9,C), "(G0.5)"  ) RHO_PARTICLE(C)
END DO
WRITE( CHAR_LABEL(10,1), "(G0)"   ) N_PARTICLES
WRITE( CHAR_LABEL(11,1), "(G0.5)" ) TOTAL_VP
WRITE( CHAR_LABEL(12,1), "(G0.5)" ) TOTAL_RHO
WRITE( CHAR_LABEL(13,1), "(G0.5)" ) TEMP
WRITE( CHAR_LABEL(14,1), "(G0.5)" ) PRESS
WRITE( CHAR_LABEL(15,1), "(G0.5)" ) BOXVMC
WRITE( CHAR_LABEL(16,1), "(2(E11.4,G0),E11.4)" ) BOXLMC(1), ", ", BOXLMC(2), ", ", BOXLMC(3)
WRITE( CHAR_LABEL(17,1), "(2(E11.4,G0),E11.4)" ) BOXLMC(4), ", ", BOXLMC(5), ", ", BOXLMC(6)
WRITE( CHAR_LABEL(18,1), "(2(E11.4,G0),E11.4)" ) BOXLMC(7), ", ", BOXLMC(8), ", ", BOXLMC(9)
WRITE( CHAR_LABEL(19,1), "(G0)"   ) MAX_CYCLES
WRITE( CHAR_LABEL(20,1), "(G0)"   ) N_EQUIL
WRITE( CHAR_LABEL(21,1), "(G0)"   ) MAX_CYCLES - N_EQUIL
WRITE( CHAR_LABEL(22,1), "(G0)"   ) N_SAVE
WRITE( CHAR_LABEL(23,1), "(G0)"   ) N_ADJUST
WRITE( CHAR_LABEL(24,1), "(G0.5)" ) R_ACC_T
WRITE( CHAR_LABEL(25,1), "(G0.5)" ) R_ACC_R
WRITE( CHAR_LABEL(26,1), "(G0.5)" ) R_ACC_VI
WRITE( CHAR_LABEL(54,1), "(G0.5)" ) R_ACC_VA
WRITE( CHAR_LABEL(60,1), "(G0.5)" ) PROB_MOV
WRITE( CHAR_LABEL(27,1), "(G0.5)" ) PROB_TRANS
WRITE( CHAR_LABEL(28,1), "(G0.5)" ) PROB_ROT
WRITE( CHAR_LABEL(29,1), "(G0.5)" ) PROB_VOL
WRITE( CHAR_LABEL(30,1), "(G0.5)" ) PROB_VOL_ISO
WRITE( CHAR_LABEL(55,1), "(G0.5)" ) PROB_VOL_ANISO
WRITE( CHAR_LABEL(31,1), "(G0.5)" ) MAX_TRANS
WRITE( CHAR_LABEL(32,1), "(G0.5)" ) MAX_ROT
WRITE( CHAR_LABEL(33,1), "(G0.5)" ) MAX_VOLI
WRITE( CHAR_LABEL(56,1), "(G0.5)" ) MAX_VOLA
WRITE( CHAR_LABEL(34,1), "(G0)"   ) UNROT_AXIS
WRITE( CHAR_LABEL(35,1), "(G0.5)" ) QUATERNION_ANGLE
WRITE( CHAR_LABEL(60,1), "(G0.5)" ) BOX_DIST
WRITE( CHAR_LABEL(61,1), "(G0.5)" ) MAX_LENGTH_RATIO
WRITE( CHAR_LABEL(62,1), "(G0.5)" ) MAX_ANGLE * 180.D0 / PI
IF( LRED_SELEC(1) ) THEN
  WRITE( CHAR_LABEL(63,1), "(G0)" ) "Gottwald"
ELSE IF( LRED_SELEC(2) ) THEN
  WRITE( CHAR_LABEL(63,1), "(G0)" ) "Lenstra-Lenstra-Lovász"
END IF
WRITE( CHAR_LABEL(36,1), "(G0)"   ) CONFIGURATION
IF( ( STOP_TIMER - START_TIMER ) < 60.D0 ) THEN
  WRITE( CHAR_LABEL(37,1), "(2G0)" ) FLOOR( STOP_TIMER - START_TIMER ), "s"
ELSE IF( ( STOP_TIMER - START_TIMER ) < 3600.D0 ) THEN
  WRITE( CHAR_LABEL(37,1), "(4G0)" ) FLOOR( ( STOP_TIMER - START_TIMER ) / 60.D0 ), "m:", &
  &                                  FLOOR( ( STOP_TIMER - START_TIMER ) - 60.D0 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / &
  &                                  60.D0 ) ) ), "s"
ELSE IF( ( STOP_TIMER - START_TIMER ) < 86400.D0 ) THEN
  WRITE( CHAR_LABEL(37,1), "(6G0)" ) FLOOR( ( STOP_TIMER - START_TIMER ) / 3.6d3 ), "h:", &
  &                                  FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 60.D0 ) ) - 60.D0 * &
  &                                  DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 3.6d3 ) ) ), "m:", &
  &                                  FLOOR( ( STOP_TIMER - START_TIMER ) - 3600.D0 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / &
  &                                  3.6d3 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 60.D0 ) ) - &
  &                                  60.D0 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 3.6d3 ) ) ) ) ), "s"
ELSE IF( ( STOP_TIMER - START_TIMER ) >= 86400.D0 ) THEN
  WRITE( CHAR_LABEL(37,1), "(8G0)" ) FLOOR( ( STOP_TIMER - START_TIMER ) / 864.d2 ), "d:", &
  &                                  FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 3600.D0 ) ) - 24.D0 * &
  &                                  DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 864.d2 ) ) ), "h:", &
  &                                  FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 60.D0 ) ) - 1440.D0 * &
  &                                  DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 864.d2 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( &
  &                                  STOP_TIMER - START_TIMER ) / 3600.D0 ) ) - 24.D0 * DBLE( FLOOR( ( STOP_TIMER - &
  &                                  START_TIMER ) / 864.d2 ) ) ) ) ), "m:", &
  &                                  FLOOR( ( STOP_TIMER - START_TIMER ) - 864.D2 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / &
  &                                  864.d2 ) ) - 3600.D0 * DBLE( FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 3600.D0 ) ) - &
  &                                  24.D0 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 864.d2 ) ) ) ) - 60.D0 * &
  &                                  DBLE( FLOOR( ( ( STOP_TIMER - START_TIMER ) / 60.D0 ) - 1440.D0 * DBLE( FLOOR( ( STOP_TIMER - &
  &                                  START_TIMER ) / 864.d2 ) ) - 60.D0 * DBLE( FLOOR( DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / &
  &                                  3600.D0 ) ) - 24.D0 * DBLE( FLOOR( ( STOP_TIMER - START_TIMER ) / 864.d2 ) ) ) ) ) ) ), "s"
END IF
FORMAT_SELF = "(I4,G0,I2.2,G0,I2.2,G0,I2.2,G0,I2.2,G0,I2.2)"
WRITE( CHAR_LABEL(38,1), FORMAT_SELF ) DATE_TIME(1), "/", DATE_TIME(2), "/", DATE_TIME(3), " ", &
&                                      DATE_TIME(5), ":", DATE_TIME(6), ":", DATE_TIME(7)
DO C = 1, COMPONENTS
  WRITE( CHAR_LABEL(39,C), "(G0)" ) C
END DO
IF( CONFIG_SELEC(4) ) THEN
  WRITE( CHAR_LABEL(40,1), "(G0.5)" ) ETA_INI
  WRITE( CHAR_LABEL(42,1), "(G0.5)" ) PRESS_RND
  WRITE( CHAR_LABEL(43,1), "(G0)"   ) N_ADJUST_INIT
  WRITE( CHAR_LABEL(44,1), "(G0.5)" ) DRMAX_INIT
  WRITE( CHAR_LABEL(45,1), "(G0.5)" ) ANGMAX_INIT
  WRITE( CHAR_LABEL(46,1), "(G0.5)" ) DVMAXISO_INIT
  WRITE( CHAR_LABEL(57,1), "(G0.5)" ) DVMAXANISO_INIT
  WRITE( CHAR_LABEL(47,1), "(G0.5)" ) DVMIN_INIT
  WRITE( CHAR_LABEL(48,1), "(G0.5)" ) PROB_MOV_INIT
  WRITE( CHAR_LABEL(49,1), "(G0.5)" ) PROB_VOL_INIT
  WRITE( CHAR_LABEL(50,1), "(G0.5)" ) PROB_TRANS_INIT
  WRITE( CHAR_LABEL(51,1), "(G0.5)" ) PROB_ROT_INIT
  WRITE( CHAR_LABEL(58,1), "(G0.5)" ) PROB_ISO_INIT
  WRITE( CHAR_LABEL(59,1), "(G0.5)" ) PROB_ANISO_INIT
END IF
IF( INIT_CONF ) THEN
  WRITE( CHAR_LABEL(52,1), "(G0)" ) "[YES]"
ELSE IF( .NOT. INIT_CONF ) THEN
  WRITE( CHAR_LABEL(52,1), "(G0)" ) "[NO]"
END IF
IF( FSEED ) THEN
  WRITE( CHAR_LABEL(53,1), "(G0)" ) "[YES]"
ELSE IF( .NOT. FSEED ) THEN
  WRITE( CHAR_LABEL(53,1), "(G0)" ) "[NO]"
END IF
IF( POTENTIAL_SELEC(1) ) THEN
  WRITE( CHAR_LABEL(64,1), "(G0)" ) "Purely Repulsive Hard-Core"
ELSE IF( POTENTIAL_SELEC(2) ) THEN
  WRITE( CHAR_LABEL(64,1), "(G0)" ) "Spherical Square-Well"
END IF
IF( POTENTIAL_SELEC(2) ) THEN
  WRITE( CHAR_LABEL(65,1), "(G0)"   ) N_LAMBDA
  DO I = 1, N_LAMBDA
    WRITE( CHAR_LABEL_L(1,I), "(G0.5)" ) L_RANGE(I)
  END DO
  WRITE( CHAR_LABEL(67,1), "(G0.5)" ) RED_TEMP
  IF( POTENTIAL_CHECK ) THEN
    WRITE( CHAR_LABEL(68,1), "(G0)" ) "Only Production Cycles"
  ELSE IF( .NOT. POTENTIAL_CHECK ) THEN
    WRITE( CHAR_LABEL(68,1), "(G0)" ) "Equilibration and Production Cycles"
  END IF
  WRITE( CHAR_LABEL(69,1), "(G0)"   ) MIN_BLOCKS
  WRITE( CHAR_LABEL(70,1), "(G0)"   ) MAX_BLOCKS
END IF
DO C = 1, COMPONENTS
  WRITE( CHAR_LABEL(71,C), "(G0)"   ) SPHERCOMP(C)
END DO
WRITE( CHAR_LABEL(72,C), "(G0)"   ) SEED

! Log strings
ALLOCATE( LOG_STRINGS_H(6), LOG_STRINGS_T(78,COMPONENTS), LOG_STRINGS_S(6) )
ALLOCATE( LOG_STRINGS_L( 1,( INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1 ) ) )
ALLOCATE( STRSH(6), STRST(78,COMPONENTS), STRSS(6) )

! String name (header)
LOG_STRINGS_H(1)  = "MONTE CARLO SIMULATION LOG"
LOG_STRINGS_H(2)  = "NVT/NPT-Monte Carlo algorithm for cylindrycally-symmetric molecules"
LOG_STRINGS_H(3)  = "Nathan Barros de Souza"
LOG_STRINGS_H(4)  = "University of Campinas"
LOG_STRINGS_H(5)  = "School of Chemical Engineering"
LOG_STRINGS_H(6)  = "Supervisor: Luis Fernando Mercier Franco"

! String name (text)
LOG_STRINGS_T(1,1) = "Execution Date: "//TRIM( CHAR_LABEL(38,1) )
LOG_STRINGS_T(2,1) = "Ensemble: "//TRIM( MC_ENSEMBLE )
LOG_STRINGS_T(3,1) = "Molecular Shape: "//TRIM( GEOMETRY )
LOG_STRINGS_T(4,1) = "Number of Components: "//TRIM( CHAR_LABEL(1,1) )
LOG_STRINGS_T(5,1) = "Packing Fraction: "//TRIM( CHAR_LABEL(2,1) )
DO C = 1, COMPONENTS
  LOG_STRINGS_T(6,C)  = "● COMPONENT #"//TRIM( CHAR_LABEL(39,C) )
  IF( SPHERCOMP(C) ) THEN
    LOG_STRINGS_T(77,C) = "Spherical Component"
  ELSE
    LOG_STRINGS_T(77,C) = "Non-Spherical Component"
  END IF
  LOG_STRINGS_T(7,C)  = "Diameter: "//TRIM( CHAR_LABEL(3,C) )//"Å"
  LOG_STRINGS_T(8,C)  = "Length: "//TRIM( CHAR_LABEL(4,C) )//"Å"
  LOG_STRINGS_T(9,C)  = "Elongation: "//TRIM( CHAR_LABEL(5,C) )
  LOG_STRINGS_T(10,C) = "Molar Fraction: "//TRIM( CHAR_LABEL(6,C) )
  LOG_STRINGS_T(11,C) = "Number of Particles: "//TRIM( CHAR_LABEL(7,C) )
  LOG_STRINGS_T(12,C) = "Particle Volume: "//TRIM( CHAR_LABEL(8,C) )//"Å³"
  LOG_STRINGS_T(13,C) = "Number Density: "//TRIM( CHAR_LABEL(9,C) )//"Å⁻³"
END DO
LOG_STRINGS_T(14,1) = "Total Number of Particles: "//TRIM( CHAR_LABEL(10,1) )
LOG_STRINGS_T(15,1) = "Total Particle Volume: "//TRIM( CHAR_LABEL(11,1) )//"Å³"
LOG_STRINGS_T(16,1) = "Total Number Density: "//TRIM( CHAR_LABEL(12,1) )//"Å⁻³"
LOG_STRINGS_T(17,1) = "Absolute Temperature: "//TRIM( CHAR_LABEL(13,1) )//"K"
LOG_STRINGS_T(18,1) = "Reduced Pressure: "//TRIM( CHAR_LABEL(14,1) )
LOG_STRINGS_T(19,1) = "Box Volume: "//TRIM( CHAR_LABEL(15,1) )//"Å³"
LOG_STRINGS_T(20,1) = "Box Length (Å):"
LOG_STRINGS_T(21,1) = TRIM( CHAR_LABEL(16,1) )
LOG_STRINGS_T(22,1) = TRIM( CHAR_LABEL(17,1) )
LOG_STRINGS_T(23,1) = TRIM( CHAR_LABEL(18,1) )
LOG_STRINGS_T(24,1) = "Total Number of Cycles: "//TRIM( CHAR_LABEL(19,1) )
LOG_STRINGS_T(25,1) = "Equilibration Cycles: "//TRIM( CHAR_LABEL(20,1) )
LOG_STRINGS_T(26,1) = "Production Cycles: "//TRIM( CHAR_LABEL(21,1) )
LOG_STRINGS_T(27,1) = "Saving Frequency (Cycles): "//TRIM( CHAR_LABEL(22,1) )
LOG_STRINGS_T(28,1) = "Adjustment Frequency (Cycles): "//TRIM( CHAR_LABEL(23,1) )
LOG_STRINGS_T(29,1) = "Ratio Threshold (Translation): "//TRIM( CHAR_LABEL(24,1) )
LOG_STRINGS_T(30,1) = "Ratio Threshold (Rotation): "//TRIM( CHAR_LABEL(25,1) )
LOG_STRINGS_T(31,1) = "Ratio Threshold (Isotropic): "//TRIM( CHAR_LABEL(26,1) )
LOG_STRINGS_T(58,1) = "Ratio Threshold (Anisotropic): "//TRIM( CHAR_LABEL(54,1) )
LOG_STRINGS_T(59,1) = "Probability (Movement): "//TRIM( CHAR_LABEL(60,1) )
LOG_STRINGS_T(32,1) = "Probability (Translation): "//TRIM( CHAR_LABEL(27,1) )
LOG_STRINGS_T(33,1) = "Probability (Rotation): "//TRIM( CHAR_LABEL(28,1) )
LOG_STRINGS_T(34,1) = "Probability (Volume): "//TRIM( CHAR_LABEL(29,1) )
LOG_STRINGS_T(35,1) = "Probability (Isotropic): "//TRIM( CHAR_LABEL(30,1) )
LOG_STRINGS_T(60,1) = "Probability (Anisotropic): "//TRIM( CHAR_LABEL(55,1) )
LOG_STRINGS_T(36,1) = "Initial Maximum Displacement (Translation): "//TRIM( CHAR_LABEL(31,1) )
LOG_STRINGS_T(37,1) = "Initial Maximum Displacement (Rotation): "//TRIM( CHAR_LABEL(32,1) )
LOG_STRINGS_T(38,1) = "Initial Maximum Displacement (Isotropic): "//TRIM( CHAR_LABEL(33,1) )
LOG_STRINGS_T(61,1) = "Initial Maximum Displacement (Anisotropic): "//TRIM( CHAR_LABEL(56,1) )
LOG_STRINGS_T(65,1) = "Maximum Box Distortion: "//TRIM( CHAR_LABEL(60,1) )
LOG_STRINGS_T(66,1) = "Maximum Box Length Distortion: "//TRIM( CHAR_LABEL(61,1) )
LOG_STRINGS_T(67,1) = "Maximum Box Angle Distortion: "//TRIM( CHAR_LABEL(62,1) )//"°"
LOG_STRINGS_T(68,1) = "Lattice Reduction Algorithm: "//TRIM( CHAR_LABEL(63,1) )
LOG_STRINGS_T(39,1) = "Unrotated Axis (Initial Configuration): "//TRIM( CHAR_LABEL(34,1) )
LOG_STRINGS_T(40,1) = "Quaternion Angle (Initial Configuration): "//TRIM( CHAR_LABEL(35,1) )
LOG_STRINGS_T(41,1) = "● Configuration: "//TRIM( CHAR_LABEL(36,1) )//" Structure"
LOG_STRINGS_T(42,1) = "Initial Packing Fraction (Random Configuration): "//TRIM( CHAR_LABEL(40,1) )
LOG_STRINGS_T(44,1) = "Target Reduced Pressure (Random Configuration): "//TRIM( CHAR_LABEL(42,1) )
LOG_STRINGS_T(45,1) = "Adjustment Frequency (Random Configuration): Every "//TRIM( CHAR_LABEL(43,1) )//" Cycle(s)"
LOG_STRINGS_T(46,1) = "Maximum Rotation (Random Configuration): "//TRIM( CHAR_LABEL(44,1) )
LOG_STRINGS_T(47,1) = "Maximum Translation (Random Configuration): "//TRIM( CHAR_LABEL(45,1) )
LOG_STRINGS_T(48,1) = "Maximum Isotropic Volume Change (Random Configuration): "//TRIM( CHAR_LABEL(46,1) )
LOG_STRINGS_T(62,1) = "Maximum Anisotropic Volume Change (Random Configuration): "//TRIM( CHAR_LABEL(57,1) )
LOG_STRINGS_T(49,1) = "Minimum Volume Change (Random Configuration): "//TRIM( CHAR_LABEL(47,1) )
LOG_STRINGS_T(50,1) = "Movement Probability (Random Configuration): "//TRIM( CHAR_LABEL(48,1) )
LOG_STRINGS_T(51,1) = "Volume Change Probability (Random Configuration): "//TRIM( CHAR_LABEL(49,1) )
LOG_STRINGS_T(52,1) = "Translation Probability (Random Configuration): "//TRIM( CHAR_LABEL(50,1) )
LOG_STRINGS_T(53,1) = "Rotation Probability (Random Configuration): "//TRIM( CHAR_LABEL(51,1) )
LOG_STRINGS_T(63,1) = "Isotropic Probability (Random Configuration): "//TRIM( CHAR_LABEL(58,1) )
LOG_STRINGS_T(64,1) = "Anisotropic Probability (Random Configuration): "//TRIM( CHAR_LABEL(59,1) )
LOG_STRINGS_T(54,1) = "Preset Initial Configuration: "//TRIM( CHAR_LABEL(52,1) )
IF( TRAJ_CHECK ) THEN
  LOG_STRINGS_T(55,1) = "Trajectory of particles computed."
ELSE IF( .NOT. TRAJ_CHECK ) THEN
  LOG_STRINGS_T(55,1) = "Trajectory of particles not computed."
END IF
LOG_STRINGS_T(56,1) = "Fixed seed: "//TRIM( CHAR_LABEL(53,1) )
LOG_STRINGS_T(57,1) = "Simulation length: "//TRIM( CHAR_LABEL(37,1) )
LOG_STRINGS_T(69,1) = "● Potential Type: "//TRIM( CHAR_LABEL(64,1) )//" Potential"
LOG_STRINGS_T(70,1) = "Number of Attractive Range Points: "//TRIM( CHAR_LABEL(65,1) )
IF( N_LAMBDA >= 4 ) THEN
  REMAINDER = MOD( N_LAMBDA, 4 )
ELSE
  REMAINDER = N_LAMBDA
END IF
IF( N_LAMBDA >= 4 .AND. REMAINDER == 0 ) THEN
  DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 ) - 1
    LOG_STRINGS_L(1,I) = TRIM( CHAR_LABEL_L(1,4*I-3) )//", "//TRIM( CHAR_LABEL_L(1,4*I-2) )//", "// &
    &                    TRIM( CHAR_LABEL_L(1,4*I-1) )//", "//TRIM( CHAR_LABEL_L(1,4*I-0) )//", "
  END DO
  I = INT( DBLE( N_LAMBDA ) / 4.D0 )
  LOG_STRINGS_L(1,I) = TRIM( CHAR_LABEL_L(1,4*I-3) )//", "//TRIM( CHAR_LABEL_L(1,4*I-2) )//", "// &
  &                    TRIM( CHAR_LABEL_L(1,4*I-1) )//", "//TRIM( CHAR_LABEL_L(1,4*I-0) )
ELSE IF( N_LAMBDA >= 4 .AND. REMAINDER /= 0 ) THEN
  DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 )
    LOG_STRINGS_L(1,I) = TRIM( CHAR_LABEL_L(1,4*I-3) )//", "//TRIM( CHAR_LABEL_L(1,4*I-2) )//", "// &
    &                    TRIM( CHAR_LABEL_L(1,4*I-1) )//", "//TRIM( CHAR_LABEL_L(1,4*I-0) )//", "
  END DO
END IF
IF( N_LAMBDA < 4 .OR. (N_LAMBDA >= 4 .AND. REMAINDER /= 0) ) THEN
  IF( REMAINDER == 1 ) THEN
    LOG_STRINGS_L(1,INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1) = TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+1) )
  ELSE IF( REMAINDER == 2 ) THEN
    LOG_STRINGS_L(1,INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1) = TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+1) )//", "// &
    &                                                     TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+2) )
  ELSE IF( REMAINDER == 3 ) THEN
    LOG_STRINGS_L(1,INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1) = TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+1) )//", "// &
    &                                                     TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+2) )//", "// &
    &                                                     TRIM( CHAR_LABEL_L(1,N_LAMBDA-REMAINDER+3) )
  END IF
END IF
LOG_STRINGS_T(71,1) = "Attractive Range Points:"
LOG_STRINGS_T(72,1) = "Reduced Temperature: "//TRIM( CHAR_LABEL(67,1) )
LOG_STRINGS_T(73,1) = "Potential Files: "//TRIM( CHAR_LABEL(68,1) )
IF( COEF_CHECK ) THEN
  LOG_STRINGS_T(74,1) = "Perturbation coefficients computed."
ELSE IF( .NOT. COEF_CHECK ) THEN
  LOG_STRINGS_T(74,1) = "Perturbation coefficients not computed."
END IF
LOG_STRINGS_T(75,1) = "Minimum Number of Blocks (Block Averaging): "//TRIM( CHAR_LABEL(69,1) )
LOG_STRINGS_T(76,1) = "Maximum Number of Blocks (Block Averaging): "//TRIM( CHAR_LABEL(70,1) )
LOG_STRINGS_T(78,1) = "Seed Value: "//TRIM( CHAR_LABEL(72,1) )

! String name (subtitle)
IF( COMPONENTS > 1 ) THEN
  LOG_STRINGS_S(1) = "PARAMETERS OF THE COMPONENTS"
ELSE IF( COMPONENTS == 1 ) THEN
  LOG_STRINGS_S(1) = "PARAMETERS OF THE COMPONENT"
END IF
LOG_STRINGS_S(2) = "PARAMETERS OF THE SYSTEM"
LOG_STRINGS_S(3) = "PARAMETERS OF THE SIMULATION"
LOG_STRINGS_S(4) = "PARAMETERS OF THE INITIAL CONFIGURATION"
LOG_STRINGS_S(5) = "OTHER PARAMETERS"
LOG_STRINGS_S(6) = "PARAMETERS OF THE POTENTIAL"

! String size
DO I = 1, 6
  STRSH(I) = ( 70.D0 - DBLE( LEN( TRIM( LOG_STRINGS_H(I) ) ) ) ) * 0.5D0
END DO
DO I = 1, 78
  DO C = 1, COMPONENTS
    STRST(I,C) = ( 69 - LEN( TRIM( LOG_STRINGS_T(I,C) ) ) )
  END DO
END DO
DO I = 1, 6
  STRSS(I) = ( 70.D0 - DBLE( LEN( TRIM( LOG_STRINGS_S(I) ) ) ) ) * 0.5D0
END DO

! *********************************************************************************************** !
! File inquiry                                                                                    !
! *********************************************************************************************** !
INQUIRE( FILE= "Simulation_Log.txt", EXIST= FILE_EXIST )

! *********************************************************************************************** !
! Simulation log                                                                                  !
! *********************************************************************************************** !
IF( .NOT. FILE_EXIST ) THEN
  OPEN( UNIT= 95, FILE= "Simulation_Log.txt" )
  WRITE( 95, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(1) ) )//TRIM( LOG_STRINGS_H(1) )//REPEAT( " ", CEILING( STRSH(1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(2) ) )//TRIM( LOG_STRINGS_H(2) )//REPEAT( " ", CEILING( STRSH(2) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(3) ) )//TRIM( LOG_STRINGS_H(3) )//REPEAT( " ", CEILING( STRSH(3) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(4) ) )//TRIM( LOG_STRINGS_H(4) )//REPEAT( " ", CEILING( STRSH(4) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(5) ) )//TRIM( LOG_STRINGS_H(5) )//REPEAT( " ", CEILING( STRSH(5) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSH(6) ) )//TRIM( LOG_STRINGS_H(6) )//REPEAT( " ", CEILING( STRSH(6) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VL//REPEAT( CH_HS, 70 )//CH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(1,1) )//REPEAT( " ", NINT( STRST(1,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(2,1) )//REPEAT( " ", NINT( STRST(2,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(3,1) )//REPEAT( " ", NINT( STRST(3,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(4,1) )//REPEAT( " ", NINT( STRST(4,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(5,1) )//REPEAT( " ", NINT( STRST(5,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
  IF( COMPONENTS > 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(1) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
    &                   SS_BR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    DO C = 1, COMPONENTS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(6,C) )//REPEAT( " ", NINT( STRST(6,C) ) + 2 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(77,C) )//REPEAT( " ", NINT( STRST(77,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(7,C) )//REPEAT( " ", NINT( STRST(7,C) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(8,C) )//REPEAT( " ", NINT( STRST(8,C) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(9,C) )//REPEAT( " ", NINT( STRST(9,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(10,C) )//REPEAT( " ", NINT( STRST(10,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(11,C) )//REPEAT( " ", NINT( STRST(11,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(12,C) )//REPEAT( " ", NINT( STRST(12,C) ) + 2 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(13,C) )//REPEAT( " ", NINT( STRST(13,C) ) + 4 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    END DO
  ELSE IF( COMPONENTS == 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(1) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
    &                   SS_BR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(77,C) )//REPEAT( " ", NINT( STRST(77,C) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(7,1) )//REPEAT( " ", NINT( STRST(7,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(8,1) )//REPEAT( " ", NINT( STRST(8,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(9,1) )//REPEAT( " ", NINT( STRST(9,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(10,1) )//REPEAT( " ", NINT( STRST(10,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(11,1) )//REPEAT( " ", NINT( STRST(11,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(12,1) )//REPEAT( " ", NINT( STRST(12,1) ) + 2 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(13,1) )//REPEAT( " ", NINT( STRST(13,1) ) + 4 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(2) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(2) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(2) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(2) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(2) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(2) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(2) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(14,1) )//REPEAT( " ", NINT( STRST(14,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(15,1) )//REPEAT( " ", NINT( STRST(15,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(16,1) )//REPEAT( " ", NINT( STRST(16,1) ) + 4 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(17,1) )//REPEAT( " ", NINT( STRST(17,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(18,1) )//REPEAT( " ", NINT( STRST(18,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(19,1) )//REPEAT( " ", NINT( STRST(19,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )// &
  &                   SS_UL//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(21,1) ) ) )// &
  &                   SS_UR//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )//SS_VS//TRIM( LOG_STRINGS_T(21,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(20,1) )//REPEAT( " ", 1 )//SS_VS//TRIM( LOG_STRINGS_T(22,1) )// &
  &                   SS_VS//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )//SS_VS//TRIM( LOG_STRINGS_T(23,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )// &
  &                   SS_BL//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(21,1) ) ) )// &
  &                   SS_BR//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(3) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(3) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(3) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(3) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(3) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(3) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(3) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(24,1) )//REPEAT( " ", NINT( STRST(24,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(25,1) )//REPEAT( " ", NINT( STRST(25,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(26,1) )//REPEAT( " ", NINT( STRST(26,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(27,1) )//REPEAT( " ", NINT( STRST(27,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(28,1) )//REPEAT( " ", NINT( STRST(28,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(29,1) )//REPEAT( " ", NINT( STRST(29,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(30,1) )//REPEAT( " ", NINT( STRST(30,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(31,1) )//REPEAT( " ", NINT( STRST(31,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(58,1) )//REPEAT( " ", NINT( STRST(58,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(59,1) )//REPEAT( " ", NINT( STRST(59,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(32,1) )//REPEAT( " ", NINT( STRST(32,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(33,1) )//REPEAT( " ", NINT( STRST(33,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(34,1) )//REPEAT( " ", NINT( STRST(34,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(35,1) )//REPEAT( " ", NINT( STRST(35,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(60,1) )//REPEAT( " ", NINT( STRST(60,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(36,1) )//REPEAT( " ", NINT( STRST(36,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(37,1) )//REPEAT( " ", NINT( STRST(37,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(38,1) )//REPEAT( " ", NINT( STRST(38,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(61,1) )//REPEAT( " ", NINT( STRST(61,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(65,1) )//REPEAT( " ", NINT( STRST(65,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(66,1) )//REPEAT( " ", NINT( STRST(66,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(67,1) )//REPEAT( " ", NINT( STRST(67,1) ) + 1 )//CH_VS
  IF( LRED_SELEC(1) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(68,1) )//REPEAT( " ", NINT( STRST(68,1) ) )//CH_VS
  ELSE IF( LRED_SELEC(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(68,1) )//REPEAT( " ", NINT( STRST(68,1) ) + 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(4) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(4) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(4) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(4) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(4) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(4) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(4) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(39,1) )//REPEAT( " ", NINT( STRST(39,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(40,1) )//REPEAT( " ", NINT( STRST(40,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(41,1) )//REPEAT( " ", NINT( STRST(41,1) ) + 2 )//CH_VS
  IF( CONFIG_SELEC(4) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(42,1) )//REPEAT( " ", NINT( STRST(42,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(44,1) )//REPEAT( " ", NINT( STRST(44,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(45,1) )//REPEAT( " ", NINT( STRST(45,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(46,1) )//REPEAT( " ", NINT( STRST(46,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(47,1) )//REPEAT( " ", NINT( STRST(47,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(48,1) )//REPEAT( " ", NINT( STRST(48,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(62,1) )//REPEAT( " ", NINT( STRST(62,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(49,1) )//REPEAT( " ", NINT( STRST(49,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(50,1) )//REPEAT( " ", NINT( STRST(50,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(51,1) )//REPEAT( " ", NINT( STRST(51,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(52,1) )//REPEAT( " ", NINT( STRST(52,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(53,1) )//REPEAT( " ", NINT( STRST(53,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(63,1) )//REPEAT( " ", NINT( STRST(63,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(64,1) )//REPEAT( " ", NINT( STRST(64,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(54,1) )//REPEAT( " ", NINT( STRST(54,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(6) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(6) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(6) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(6) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(6) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(6) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(6) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(69,1) )//REPEAT( " ", NINT( STRST(69,1) ) + 2 )//CH_VS
  IF( POTENTIAL_SELEC(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(70,1) )//REPEAT( " ", NINT( STRST(70,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(71,1) )//REPEAT( " ", NINT( STRST(71,1) ) - 1 )//CH_VS
    IF( N_LAMBDA == 4 ) THEN
      DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 )
        WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 )//CH_VS
      END DO
    ELSE IF( N_LAMBDA /= 4 ) THEN
      DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 )
        IF( I == INT( DBLE( N_LAMBDA ) / 4.D0 ) .AND. REMAINDER == 0 ) THEN
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 )//CH_VS
        ELSE
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 35 )//CH_VS
        END IF
      END DO
    END IF
    IF( REMAINDER /= 0 ) THEN
      I = INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1
      REMAINDER = 4 - REMAINDER
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 + (8 * REMAINDER) )//CH_VS
    END IF
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(72,1) )//REPEAT( " ", NINT( STRST(72,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(73,1) )//REPEAT( " ", NINT( STRST(73,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(74,1) )//REPEAT( " ", NINT( STRST(74,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(75,1) )//REPEAT( " ", NINT( STRST(75,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(76,1) )//REPEAT( " ", NINT( STRST(76,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(5) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(5) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(5) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(5) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(5) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(5) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(5) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(55,1) )//REPEAT( " ", NINT( STRST(55,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(56,1) )//REPEAT( " ", NINT( STRST(56,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(78,1) )//REPEAT( " ", NINT( STRST(78,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(57,1) )//REPEAT( " ", NINT( STRST(57,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  CLOSE( 95 )
! *********************************************************************************************** !
! Simulation log (appending)                                                                      !
! *********************************************************************************************** !
ELSE IF( FILE_EXIST ) THEN
  OPEN( UNIT= 95, FILE= "Simulation_Log.txt", POSITION= "APPEND" )
  WRITE( 95, "(G0)" ) REPEAT( " ", 1 )//REPEAT( C_FUL, 70 )//REPEAT( " ", 1 )
  WRITE( 95, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(1,1) )//REPEAT( " ", NINT( STRST(1,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(2,1) )//REPEAT( " ", NINT( STRST(2,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(3,1) )//REPEAT( " ", NINT( STRST(3,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(4,1) )//REPEAT( " ", NINT( STRST(4,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(5,1) )//REPEAT( " ", NINT( STRST(5,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
  IF( COMPONENTS > 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(1) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
    &                   SS_BR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    DO C = 1, COMPONENTS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(6,C) )//REPEAT( " ", NINT( STRST(6,C) ) + 2 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(77,C) )//REPEAT( " ", NINT( STRST(77,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(7,C) )//REPEAT( " ", NINT( STRST(7,C) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(8,C) )//REPEAT( " ", NINT( STRST(8,C) ) + 1 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(9,C) )//REPEAT( " ", NINT( STRST(9,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(10,C) )//REPEAT( " ", NINT( STRST(10,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(11,C) )//REPEAT( " ", NINT( STRST(11,C) ) - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(12,C) )//REPEAT( " ", NINT( STRST(12,C) ) + 2 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(13,C) )//REPEAT( " ", NINT( STRST(13,C) ) + 4 - 1 )//CH_VS
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    END DO
  ELSE IF( COMPONENTS == 1 ) THEN
    WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(1) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(1) )// &
    &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(1) - 1 ) )//SH_VR
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(1) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(1) ) ) )// &
    &                   SS_BR//REPEAT( " ", CEILING( STRSS(1) - 1 ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(77,C) )//REPEAT( " ", NINT( STRST(77,C) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(7,1) )//REPEAT( " ", NINT( STRST(7,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(8,1) )//REPEAT( " ", NINT( STRST(8,1) ) + 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(9,1) )//REPEAT( " ", NINT( STRST(9,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(10,1) )//REPEAT( " ", NINT( STRST(10,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(11,1) )//REPEAT( " ", NINT( STRST(11,1) ) )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(12,1) )//REPEAT( " ", NINT( STRST(12,1) ) + 2 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(13,1) )//REPEAT( " ", NINT( STRST(13,1) ) + 4 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(2) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(2) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(2) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(2) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(2) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(2) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(2) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(2) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(14,1) )//REPEAT( " ", NINT( STRST(14,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(15,1) )//REPEAT( " ", NINT( STRST(15,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(16,1) )//REPEAT( " ", NINT( STRST(16,1) ) + 4 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(17,1) )//REPEAT( " ", NINT( STRST(17,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(18,1) )//REPEAT( " ", NINT( STRST(18,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(19,1) )//REPEAT( " ", NINT( STRST(19,1) ) + 2 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )// &
  &                   SS_UL//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(21,1) ) ) )// &
  &                   SS_UR//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )//SS_VS//TRIM( LOG_STRINGS_T(21,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(20,1) )//REPEAT( " ", 1 )//SS_VS//TRIM( LOG_STRINGS_T(22,1) )// &
  &                   SS_VS//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )//SS_VS//TRIM( LOG_STRINGS_T(23,1) )//SS_VS// &
  &                   REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + 1 )// &
  &                   SS_BL//REPEAT( " ", LEN( TRIM( LOG_STRINGS_T(21,1) ) ) )// &
  &                   SS_BR//REPEAT( " ", 67 - ( LEN( TRIM( LOG_STRINGS_T(20,1) ) ) + LEN( TRIM( LOG_STRINGS_T(21,1) ) ) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(3) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(3) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(3) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(3) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(3) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(3) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(3) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(3) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(24,1) )//REPEAT( " ", NINT( STRST(24,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(25,1) )//REPEAT( " ", NINT( STRST(25,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(26,1) )//REPEAT( " ", NINT( STRST(26,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(27,1) )//REPEAT( " ", NINT( STRST(27,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(28,1) )//REPEAT( " ", NINT( STRST(28,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(29,1) )//REPEAT( " ", NINT( STRST(29,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(30,1) )//REPEAT( " ", NINT( STRST(30,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(31,1) )//REPEAT( " ", NINT( STRST(31,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(58,1) )//REPEAT( " ", NINT( STRST(58,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(59,1) )//REPEAT( " ", NINT( STRST(59,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(32,1) )//REPEAT( " ", NINT( STRST(32,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(33,1) )//REPEAT( " ", NINT( STRST(33,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(34,1) )//REPEAT( " ", NINT( STRST(34,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(35,1) )//REPEAT( " ", NINT( STRST(35,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(60,1) )//REPEAT( " ", NINT( STRST(60,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(36,1) )//REPEAT( " ", NINT( STRST(36,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(37,1) )//REPEAT( " ", NINT( STRST(37,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(38,1) )//REPEAT( " ", NINT( STRST(38,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(61,1) )//REPEAT( " ", NINT( STRST(61,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(65,1) )//REPEAT( " ", NINT( STRST(65,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(66,1) )//REPEAT( " ", NINT( STRST(66,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(67,1) )//REPEAT( " ", NINT( STRST(67,1) ) + 1 )//CH_VS
  IF( LRED_SELEC(1) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(68,1) )//REPEAT( " ", NINT( STRST(68,1) ) )//CH_VS
  ELSE IF( LRED_SELEC(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(68,1) )//REPEAT( " ", NINT( STRST(68,1) ) + 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(4) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(4) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(4) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(4) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(4) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(4) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(4) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(4) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(39,1) )//REPEAT( " ", NINT( STRST(39,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(40,1) )//REPEAT( " ", NINT( STRST(40,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(41,1) )//REPEAT( " ", NINT( STRST(41,1) ) + 2 )//CH_VS
  IF( CONFIG_SELEC(4) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(42,1) )//REPEAT( " ", NINT( STRST(42,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(43,1) )//REPEAT( " ", NINT( STRST(43,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(44,1) )//REPEAT( " ", NINT( STRST(44,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(45,1) )//REPEAT( " ", NINT( STRST(45,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(46,1) )//REPEAT( " ", NINT( STRST(46,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(47,1) )//REPEAT( " ", NINT( STRST(47,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(48,1) )//REPEAT( " ", NINT( STRST(48,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(62,1) )//REPEAT( " ", NINT( STRST(62,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(49,1) )//REPEAT( " ", NINT( STRST(49,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(50,1) )//REPEAT( " ", NINT( STRST(50,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(51,1) )//REPEAT( " ", NINT( STRST(51,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(52,1) )//REPEAT( " ", NINT( STRST(52,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(53,1) )//REPEAT( " ", NINT( STRST(53,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(63,1) )//REPEAT( " ", NINT( STRST(63,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(64,1) )//REPEAT( " ", NINT( STRST(64,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(54,1) )//REPEAT( " ", NINT( STRST(54,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(6) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(6) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(6) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(6) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(6) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(6) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(6) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(6) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(69,1) )//REPEAT( " ", NINT( STRST(69,1) ) + 2 )//CH_VS
  IF( POTENTIAL_SELEC(2) ) THEN
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(70,1) )//REPEAT( " ", NINT( STRST(70,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(71,1) )//REPEAT( " ", NINT( STRST(71,1) ) - 1 )//CH_VS
    IF( N_LAMBDA == 4 ) THEN
      DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 )
        WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 )//CH_VS
      END DO
    ELSE IF( N_LAMBDA /= 4 ) THEN
      DO I = 1, INT( DBLE( N_LAMBDA ) / 4.D0 )
        IF( I == INT( DBLE( N_LAMBDA ) / 4.D0 ) .AND. REMAINDER == 0 ) THEN
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 )//CH_VS
        ELSE
          WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 35 )//CH_VS
        END IF
      END DO
    END IF
    IF( REMAINDER /= 0 ) THEN
      I = INT( DBLE( N_LAMBDA ) / 4.D0 ) + 1
      REMAINDER = 4 - REMAINDER
      WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 4 )//TRIM( LOG_STRINGS_L(1,I) )//REPEAT( " ", 36 + (8 * REMAINDER) )//CH_VS
    END IF
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(72,1) )//REPEAT( " ", NINT( STRST(72,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(73,1) )//REPEAT( " ", NINT( STRST(73,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(74,1) )//REPEAT( " ", NINT( STRST(74,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(75,1) )//REPEAT( " ", NINT( STRST(75,1) ) - 1 )//CH_VS
    WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 2 )//TRIM( LOG_STRINGS_T(76,1) )//REPEAT( " ", NINT( STRST(76,1) ) - 1 )//CH_VS
  END IF
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(5) - 1 ) )//SS_UL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(5) ) ) )// &
  &                   SS_UR//REPEAT( " ", CEILING( STRSS(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) SH_VL//REPEAT( SS_HS, FLOOR( STRSS(5) - 1 ) )//SS_VL//TRIM( LOG_STRINGS_S(5) )// &
  &                   SS_VR//REPEAT( SS_HS, CEILING( STRSS(5) - 1 ) )//SH_VR
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", FLOOR( STRSS(5) - 1 ) )//SS_BL//REPEAT( SS_HS, LEN( TRIM( LOG_STRINGS_S(5) ) ) )// &
  &                   SS_BR//REPEAT( " ", CEILING( STRSS(5) - 1 ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(55,1) )//REPEAT( " ", NINT( STRST(55,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(56,1) )//REPEAT( " ", NINT( STRST(56,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(78,1) )//REPEAT( " ", NINT( STRST(78,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 70 )//CH_VS
  WRITE( 95, "(G0)" ) CH_VS//REPEAT( " ", 1 )//TRIM( LOG_STRINGS_T(57,1) )//REPEAT( " ", NINT( STRST(57,1) ) )//CH_VS
  WRITE( 95, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  CLOSE( 95 )
END IF

! Status
WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"

END PROGRAM MAIN
