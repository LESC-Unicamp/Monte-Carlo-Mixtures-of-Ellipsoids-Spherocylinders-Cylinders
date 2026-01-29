! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!          This module allows the user to choose one of the molecular configurations for          !
!              isomorphic molecules - Simple Cube (SC), Body-Centered Cube (BCC), or              !
!         Face-Centered Cube (FCC) -, and for anisomorphic molecules - Random Box (RND).          !
!        Molecules will be then allocated in accordance to the selected crystal structure.        !
!                       See Macpherson et al. (2007) for some information.                        !
!     This module also writes out a file containing all particles' positions and quaternions.     !
!                                                                                                 !
! Version number: 2.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       January 28th, 2026                                        !
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
MODULE InitialSystemConfiguration

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BodyFixedAxis ! Body-fixed axis of rotation (for the initial configuration only)


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
! *********************************************************************************************** !

CONTAINS

! *********************************************************************************************** !
!             This subroutine allows the user to choose the geometry of the molecules             !
! *********************************************************************************************** !
SUBROUTINE GeometrySelection(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = EOR                                                                                      !
!  (2) = SPC                                                                                      !
!  (3) = CYL                                                                                      !
! *********************************************************************************************** !
GeometryType = .FALSE.

! Molecular geometry file
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )
READ( 100, * ) Dummy, GeometryInquiry
CALL ToUpper( GeometryInquiry, LEN_TRIM( GeometryInquiry ), GeometryInquiry )
CLOSE( 100 )

! Extended molecular geometry name
IF( GeometryInquiry == "EOR" ) THEN
  MolecularGeometry = "Ellipsoids of revolution"
  GeometryAcronym   = "eor"
ELSE IF( GeometryInquiry == "SPC" ) THEN
  MolecularGeometry = "Spherocylinders"
  GeometryAcronym   = "spc"
ELSE IF( GeometryInquiry == "CYL" ) THEN
  MolecularGeometry = "Cylinders"
  GeometryAcronym   = "cyl"
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined ", TRIM( GeometryInquiry ), " is not an available molecular geometry. Exiting... "
  CALL Exit(  )
END IF

! Molecular geometry inquiry
WRITE( *, "(G0)" ) "The molecules are: "//TRIM( MolecularGeometry )//". "
IF( GeometryInquiry == "EOR" ) THEN
  GeometryType(1) = .TRUE.
ELSE IF( GeometryInquiry == "SPC" ) THEN
  GeometryType(2) = .TRUE.
ELSE IF( GeometryInquiry == "CYL" ) THEN
  GeometryType(3) = .TRUE.
END IF
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE GeometrySelection

! *********************************************************************************************** !
!        This subroutine allows the user to choose the initial configuration of the system        !
! *********************************************************************************************** !
SUBROUTINE InitialConfigurationSelection(  )

IMPLICIT NONE

! *********************************************************************************************** !
! Initialization of logical array                                                                 !
!  (1) = SC                                                                                       !
!  (2) = BCC                                                                                      !
!  (3) = FCC                                                                                      !
!  (4) = RND                                                                                      !
! *********************************************************************************************** !
ConfigurationSelection = .FALSE.

! Initial configuration file
OPEN( Unit= 100, File= "ini_config.ini", Action= "READ" )
READ( 100, * ) Dummy, Dummy
READ( 100, * ) Dummy, ConfigurationInquiry
CALL ToUpper( ConfigurationInquiry, LEN_TRIM( ConfigurationInquiry ), ConfigurationInquiry )
CLOSE( 100 )

! Extended configuration name
IF( ConfigurationInquiry == "SC" ) THEN
  InitialConfiguration = "Simple Cube"
ELSE IF( ConfigurationInquiry == "BCC" ) THEN
  InitialConfiguration = "Body-Centered Cube"
ELSE IF( ConfigurationInquiry == "FCC" ) THEN
  InitialConfiguration = "Face-Centered Cube"
ELSE IF( ConfigurationInquiry == "RND" ) THEN
  InitialConfiguration = "Random"
ELSE ! Stop condition
  WRITE( *, "(3G0)" ) "The user-defined ", TRIM( ConfigurationInquiry ), " is not an available initial configuration. Exiting... "
  CALL Exit(  )
END IF

! Initial configuration inquiry
WRITE( *, "(G0)") "Initial configuration is: "//TRIM( InitialConfiguration )//". "
IF( ConfigurationInquiry == "SC" ) THEN
  ConfigurationSelection(1) = .TRUE.
ELSE IF( ConfigurationInquiry == "BCC" ) THEN
  ConfigurationSelection(2) = .TRUE.
ELSE IF( ConfigurationInquiry == "FCC" ) THEN
  ConfigurationSelection(3) = .TRUE.
ELSE IF( ConfigurationInquiry == "RND" ) THEN
  ConfigurationSelection(4) = .TRUE.
END IF
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE InitialConfigurationSelection

! *********************************************************************************************** !
!              This subroutine initializes the molecular configuration of the system              !
! *********************************************************************************************** !
SUBROUTINE InitialConfigurationStructure(  )

! Uses one module: folders
USE Folders, ONLY: InitialConfigurationFolders

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: cComponent ! Counter (component)
INTEGER( Kind= Int64 ) :: pParticle  ! Counter (particle)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FileExist ! Checks whether a file exists or not

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 14 ) :: DescriptorInitialConfig ! Descriptor for output file (initial configuration)

! Initial configuration folder (see 'Folders' module)
CALL InitialConfigurationFolders(  )

! Calls 'SimpleCubicConfiguration' subroutine if the user chooses a simple cubic structure
IF( ConfigurationSelection(1) ) THEN
  IF( .NOT. PresetInitialConfiguration ) THEN
    DO cComponent = 1, nComponents - 1
      IF( DABS( cDiameter(cComponent) - cDiameter(cComponent+1) ) >= EPSILON( 1.D0 ) .OR. &
      &   DABS( cLength(cComponent) - cLength(cComponent+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The simple cubic structure might create overlaping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
      IF( ( SphericalComponentLogical(cComponent) .AND. .NOT. SphericalComponentLogical(cComponent+1) ) .OR. &
      &   ( .NOT. SphericalComponentLogical(cComponent) .AND. SphericalComponentLogical(cComponent+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The simple cubic structure might create overlapping configurations. ", &
        &                         "Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
    END DO
    CALL SimpleCubicConfiguration(  )
    CALL ConfigurationOutput(  )
    RETURN
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be ", &
    &                   "overwritten. Continue? [Y/N]"
    READ( *, * ) Dummy
    CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
    IF( Dummy /= "Y" ) THEN
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) DescriptorInitialConfig
    INQUIRE( File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_sc_"//TRIM( DescriptorFileGeometry )// &
    &              ".xyz", EXIST= FileExist )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the selected ", &
      &                     "molecular geometry. Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL Sleep( 1 )
    OPEN( Unit= 10, Action= "READ", File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_sc_" &
    &                                     //TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
! Calls 'BodyCentredCubicConfiguration' subroutine if the user chooses a body-centered cubic structure
ELSE IF( ConfigurationSelection(2) ) THEN
  IF( .NOT. PresetInitialConfiguration ) THEN
    DO cComponent = 1, nComponents - 1
      IF( DABS( cDiameter(cComponent) - cDiameter(cComponent+1) ) >= EPSILON( 1.D0 ) .OR. &
      &   DABS( cLength(cComponent) - cLength(cComponent+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The body-centered cubic structure might create overlapping ", &
        &                         "configurations. Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
      IF( ( SphericalComponentLogical(cComponent) .AND. .NOT. SphericalComponentLogical(cComponent+1) ) .OR. &
      &   ( .NOT. SphericalComponentLogical(cComponent) .AND. SphericalComponentLogical(cComponent+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The body-centered cubic structure might create overlapping ", &
        &                         "configurations. Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
    END DO
    CALL BodyCentredCubicConfiguration(  )
    CALL ConfigurationOutput(  )
    RETURN
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be ", &
    &                   "overwritten. Continue? [Y/N]"
    READ( *, * ) Dummy
    CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
    IF( Dummy /= "Y" ) THEN
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) DescriptorInitialConfig
    INQUIRE( File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_bcc_"//TRIM( DescriptorFileGeometry )// &
    &              ".xyz", EXIST= FileExist )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL Sleep( 1 )
    OPEN( Unit= 10, Action= "READ", File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_bcc_" &
    &                                     //TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
! Calls 'FaceCentredCubicConfiguration' subroutine if the user chooses a face-centered cubic structure
ELSE IF( ConfigurationSelection(3) ) THEN
  IF( .NOT. PresetInitialConfiguration ) THEN
    DO cComponent = 1, nComponents - 1
      IF( DABS( cDiameter(cComponent) - cDiameter(cComponent+1) ) >= EPSILON( 1.D0 ) .OR. &
      &   DABS( cLength(cComponent) - cLength(cComponent+1) ) >= EPSILON( 1.D0 ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The face-centered cubic structure might create overlapping ", &
        &                         "configurations. Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
      IF( ( SphericalComponentLogical(cComponent) .AND. .NOT. SphericalComponentLogical(cComponent+1) ) .OR. &
      &   ( .NOT. SphericalComponentLogical(cComponent) .AND. SphericalComponentLogical(cComponent+1) ) ) THEN
        WRITE( *, "(5G0,/,2G0)" ) "Molecules of component ", cComponent, " are not isomorphic with molecules of component ", &
        &                         cComponent + 1, "! ", "The face-centered cubic structure might create overlapping ", &
        &                         "configurations. Do you wish to continue?"
        READ( *, * ) Dummy
        CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
        IF( Dummy /= "Y" ) THEN
          CALL Sleep( 1 )
          CALL Exit(  )
        END IF
        WRITE( *, "(G0)" ) " "
        EXIT
      END IF
    END DO
    CALL FaceCentredCubicConfiguration(  )
    CALL ConfigurationOutput(  )
    RETURN
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be ", &
    &                   "overwritten. Continue? [Y/N]"
    READ( *, * ) Dummy
    CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
    IF( Dummy /= "Y" ) THEN
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) DescriptorInitialConfig
    INQUIRE( File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_fcc_"//TRIM( DescriptorFileGeometry )// &
    &              ".xyz", EXIST= FileExist )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL Sleep( 1 )
    OPEN( Unit= 10, Action= "READ", File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_fcc_" &
    &                                     //TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
! Calls 'RandomConfiguration' subroutine if the user chooses a random structure
ELSE IF( ConfigurationSelection(4) ) THEN
  IF( .NOT. PresetInitialConfiguration ) THEN
    CALL RandomConfiguration(  )
    CALL ConfigurationOutput(  )
    RETURN
  ELSE
    WRITE( *, "(2G0)" ) "Preset initial configuration selected! Simulation parameters previously calculated will be ", &
    &                   "overwritten. Continue? [Y/N]"
    READ( *, * ) Dummy
    CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
    IF( Dummy /= "Y" ) THEN
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Please enter the 14-character descriptor code of the file: "
    READ( *, * ) DescriptorInitialConfig
    INQUIRE( File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_rnd_"//TRIM( DescriptorFileGeometry )// &
    &              ".xyz", EXIST= FileExist )
    WRITE( *, "(G0)" ) " "
    IF( .NOT. FileExist ) THEN
      WRITE( *, "(G0,G0)" ) "Preset initial configuration not found for the specified configuration or the ", &
      &                     "selected molecular geometry. Exiting..."
      CALL Sleep( 1 )
      CALL Exit(  )
    END IF
    WRITE( *, "(G0)" ) "Preset initial configuration found for the specified configuration and molecular geometry!"
    CALL Sleep( 1 )
    OPEN( Unit= 10, Action= "READ", File= "Initial_Configuration/"//TRIM( DescriptorInitialConfig )//"_initconf_rnd_" &
    &                                     //TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
END IF

! Read variables from preset initial configuration file
READ( 10, * ) Dummy, GeometryInquiry
GeometryType = .FALSE.
IF( GeometryInquiry == "EOR" ) THEN
  MolecularGeometry = "Ellipsoids of revolution"
  GeometryAcronym   = "eor"
  GeometryType(1)   = .TRUE.
ELSE IF( GeometryInquiry == "SPC" ) THEN
  MolecularGeometry = "Spherocylinders"
  GeometryAcronym   = "spc"
  GeometryType(2)   = .TRUE.
ELSE IF( GeometryInquiry == "CYL" ) THEN
  MolecularGeometry = "Cylinders"
  GeometryAcronym   = "cyl"
  GeometryType(3)   = .TRUE.
END IF
READ( 10, * ) Dummy, ConfigurationInquiry
ConfigurationSelection = .FALSE.
IF( ConfigurationInquiry == "SC" ) THEN
  InitialConfiguration      = "Simple Cube"
  ConfigurationSelection(1) = .TRUE.
ELSE IF( ConfigurationInquiry == "BCC" ) THEN
  InitialConfiguration      = "Body-Centered Cube"
  ConfigurationSelection(2) = .TRUE.
ELSE IF( ConfigurationInquiry == "FCC" ) THEN
  InitialConfiguration      = "Face-Centered Cube"
  ConfigurationSelection(3) = .TRUE.
ELSE IF( ConfigurationInquiry == "RND" ) THEN
  InitialConfiguration      = "Random"
  ConfigurationSelection(4) = .TRUE.
END IF
READ( 10, * ) Dummy, PackingFraction
READ( 10, * ) Dummy, nComponents
! Deallocate
IF( ALLOCATED( SphericalComponentInquiry ) ) DEALLOCATE( SphericalComponentInquiry )
IF( ALLOCATED( SphericalComponentLogical ) ) DEALLOCATE( SphericalComponentLogical )
IF( ALLOCATED( cDiameter ) ) DEALLOCATE( cDiameter )
IF( ALLOCATED( cLength ) ) DEALLOCATE( cLength )
IF( ALLOCATED( cAspectRatio ) ) DEALLOCATE( cAspectRatio )
IF( ALLOCATED( cMolarFraction ) ) DEALLOCATE( cMolarFraction )
IF( ALLOCATED( cParticles ) ) DEALLOCATE( cParticles )
IF( ALLOCATED( cMolecularVolume ) ) DEALLOCATE( cMolecularVolume )
IF( ALLOCATED( cNumberDensity ) ) DEALLOCATE( cNumberDensity )
IF( ALLOCATED( cIndex ) ) DEALLOCATE( cIndex )
IF( ALLOCATED( cPotentialRange ) ) DEALLOCATE( cPotentialRange )
IF( ALLOCATED( cCircumscribingPotentialRange ) ) DEALLOCATE( cCircumscribingPotentialRange )
IF( ALLOCATED( cCircumscribingSphereDiameter ) ) DEALLOCATE( cCircumscribingSphereDiameter )
IF( ALLOCATED( cDiameterEquivalentSphere ) ) DEALLOCATE( cDiameterEquivalentSphere )
! Reallocate
ALLOCATE( cDiameter(nComponents), cLength(nComponents), cAspectRatio(nComponents), cMolecularVolume(nComponents) )
ALLOCATE( cParticles(0:nComponents), cMolarFraction(nComponents), cNumberDensity(nComponents) )
ALLOCATE( cCircumscribingSphereDiameter(nComponents) )
ALLOCATE( cDiameterEquivalentSphere(nComponents) )
ALLOCATE( SphericalComponentInquiry(nComponents) )
ALLOCATE( SphericalComponentLogical(nComponents) )
ALLOCATE( cIndex(nComponents) )
ALLOCATE( cPotentialRange(nComponents,nRange) )
ALLOCATE( cCircumscribingPotentialRange(nComponents,nRange) )
! Resume
cParticles = 0
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, SphericalComponentInquiry(cComponent)
  IF( SphericalComponentInquiry(cComponent) == "T" ) THEN
    SphericalComponentLogical(cComponent) = .TRUE.
  ELSE
    SphericalComponentLogical(cComponent) = .FALSE.
  END IF
END DO
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cDiameter(cComponent), Dummy
END DO
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cLength(cComponent), Dummy
END DO
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cAspectRatio(cComponent)
END DO
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cMolarFraction(cComponent)
END DO
READ( 10, * ) Dummy, nParticles
! Deallocate
IF( ALLOCATED( pComponents ) ) DEALLOCATE( pComponents )
IF( ALLOCATED( pQuaternion ) ) DEALLOCATE( pQuaternion )
IF( ALLOCATED( pQuaternionMC ) ) DEALLOCATE( pQuaternionMC )
IF( ALLOCATED( pPosition ) ) DEALLOCATE( pPosition )
IF( ALLOCATED( pPositionMC ) ) DEALLOCATE( pPositionMC )
IF( ALLOCATED( pOrientation ) ) DEALLOCATE( pOrientation )
IF( ALLOCATED( pOrientationMC ) ) DEALLOCATE( pOrientationMC )
! Reallocate
ALLOCATE( pComponents(nParticles) )
ALLOCATE( pQuaternion(0:3,nParticles), pQuaternionMC(0:3,nParticles) )
ALLOCATE( pPosition(3,nParticles), pPositionMC(3,nParticles) )
ALLOCATE( pOrientation(3,nParticles), pOrientationMC(3,nParticles) )
! Resume
READ( 10, * ) Dummy, TotalParticleVolume, Dummy
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cParticles(cComponent)
END DO
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cMolecularVolume(cComponent), Dummy
END DO
READ( 10, * ) Dummy, BoxVolume, Dummy
READ( 10, * ) Dummy, BoxLength(1:3)
READ( 10, * ) Dummy, BoxLength(4:6)
READ( 10, * ) Dummy, BoxLength(7:9)
READ( 10, * ) Dummy, BoxLengthInverse(1:3)
READ( 10, * ) Dummy, BoxLengthInverse(4:6)
READ( 10, * ) Dummy, BoxLengthInverse(7:9)
DO cComponent = 1, nComponents
  READ( 10, * ) Dummy, cNumberDensity(cComponent), Dummy
END DO
READ( 10, * ) Dummy, TotalNumberDensity, Dummy
READ( 10, * ) Dummy, Dummy, Dummy
READ( 10, * ) Dummy, Dummy
DO cComponent = 1, nComponents
  cIndex(cComponent) = cComponent
END DO
DO cComponent = 1, nComponents
  DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
    READ( 10, * ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
    &             pQuaternion(0,pParticle), pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle)
  END DO
END DO
CLOSE( 10 )

! Diameter of circumscribing sphere
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      IF( cAspectRatio(cComponent) > 0.D0 .AND. cAspectRatio(cComponent) < 1.D0 ) THEN
        cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
      ELSE IF( cAspectRatio(cComponent) > 1.D0 ) THEN
        cCircumscribingSphereDiameter(cComponent) = cLength(cComponent)
      END IF
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent) + cLength(cComponent)
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cCircumscribingSphereDiameter(cComponent) = DSQRT( cDiameter(cComponent) * cDiameter(cComponent) + &
      &                                                  cLength(cComponent) * cLength(cComponent) )
    ELSE
      cCircumscribingSphereDiameter(cComponent) = cDiameter(cComponent)
    END IF
  END DO
END IF

! Hard-core volumetric relation (EOR/SPC/HC and SPHERES)
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ) ! Ellipsoids of revolution
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( 1.D0 + ( 1.5D0 * cAspectRatio(cComponent) ) ) ** &
      &                                       ( 1.D0 / 3.D0 ) ) ! Spherocylinders
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) * ( ( 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ) ! Cylinders
    ELSE
      cDiameterEquivalentSphere(cComponent) = cDiameter(cComponent) ! Spheres
    END IF
  END DO
END IF

! Effective range of attraction
IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN ! Spherical square-well potential
  DO cComponent = 1, nComponents
    cPotentialRange(cComponent,:) = PotentialRange(:) * cDiameterEquivalentSphere(cComponent)
  END DO
ELSE IF( PerturbedPotentialTypeLogical(3) .OR. FullPotentialTypeLogical(3) ) THEN ! Anisotropic square-well potential
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      cPotentialRange(cComponent,:) = PotentialRange(:) * cDiameter(cComponent)
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        IF( cAspectRatio(cComponent) > 0.D0 .AND. cAspectRatio(cComponent) < 1.D0 ) THEN
          cCircumscribingPotentialRange(cComponent,:) = cDiameter(cComponent) + PotentialRange(:) * cDiameter(cComponent)
        ELSE IF( cAspectRatio(cComponent) > 1.D0 ) THEN
          cCircumscribingPotentialRange(cComponent,:) = cLength(cComponent) + PotentialRange(:) * cDiameter(cComponent)
        END IF
      ELSE
        cCircumscribingPotentialRange(cComponent,:) = cDiameter(cComponent) + PotentialRange(:) * cDiameter(cComponent)
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      cPotentialRange(cComponent,:) = PotentialRange(:) * cDiameter(cComponent)
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        cCircumscribingPotentialRange(cComponent,:) = cDiameter(cComponent) + cLength(cComponent) + PotentialRange(:) * &
        &                                             cDiameter(cComponent)
      ELSE
        cCircumscribingPotentialRange(cComponent,:) = cDiameter(cComponent) + PotentialRange(:) * cDiameter(cComponent)
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      cPotentialRange(cComponent,:) = PotentialRange(:) * cDiameter(cComponent)
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        cCircumscribingPotentialRange(cComponent,:) = DSQRT( cLength(cComponent) * cLength(cComponent) + cDiameter(cComponent) * &
        &                                             cDiameter(cComponent) * ( 1.D0 + 2.D0 * PotentialRange(:) + 2.D0 * &
        &                                             PotentialRange(:) * PotentialRange(:) ) + 2.D0 * cDiameter(cComponent) * &
        &                                             cLength(cComponent) * PotentialRange(:) )
      ELSE
        cCircumscribingPotentialRange(cComponent,:) = cDiameter(cComponent) + PotentialRange(:) * cDiameter(cComponent)
      END IF
    END DO
  END IF
END IF

! Component index of a particle
DO cComponent = 1, nComponents
  DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
    pComponents(pParticle) = cComponent
  END DO
END DO

! Summary (for preset initial configuration)
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "File sucessfully read! Resuming..."
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "

! Print out initial configuration file
CALL ConfigurationOutput(  )

RETURN

END SUBROUTINE InitialConfigurationStructure

! *********************************************************************************************** !
!          This subroutine allocates particles according to a simple cubic configuration          !
! *********************************************************************************************** !
SUBROUTINE SimpleCubicConfiguration(  )

! Uses one module: linked lists
USE LinkedLists, ONLY: MakeList

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nCells              ! Number of unit cells
INTEGER( Kind= Int64 ) :: iCell, jCell, kCell ! Counter of cells
INTEGER( Kind= Int64 ) :: iParticle           ! Counter of particles
INTEGER( Kind= Int64 ) :: Counter             ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: CellLength             ! Length of unit cell (cubic structure)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: ScalingDistanceUnitBox ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BoxCutoff              ! Box cutoff (x-, y-, and z-directions)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AxisSelection(1) ) THEN
  BodyFixedAxis = xAxis
ELSE IF( AxisSelection(2) ) THEN
  BodyFixedAxis = yAxis
ELSE IF( AxisSelection(3) ) THEN
  BodyFixedAxis = zAxis
END IF

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
pQuaternion(0,:) = DCOS( QuaternionAngle * 0.5D0 )                    ! Real part
pQuaternion(1,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(3) ! Imaginary part (Vector)

! Number of unit cells per axis (simple cube)
nCells = NINT( DBLE( nParticles ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (simple cube)
CellLength = ( 1.D0 / TotalNumberDensity ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BoxLength    = 0.D0
BoxLength(1) = CellLength * DBLE( nCells )
BoxLength(5) = CellLength * DBLE( nCells )
BoxLength(9) = CellLength * DBLE( nCells )

! Simulation box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Position of particles (centers of mass)
Counter = 1
DO iCell = 1, nCells
  DO jCell = 1, nCells
    DO kCell = 1, nCells
      ! Particles on the right vertex of unit cell
      pPosition(1,Counter) = DBLE( iCell - 1 ) * CellLength
      pPosition(2,Counter) = DBLE( jCell - 1 ) * CellLength
      pPosition(3,Counter) = DBLE( kCell - 1 ) * CellLength
      Counter = Counter + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO iParticle = 1, nParticles
  ! Spatial transformation
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,iParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,iParticle) )
END DO

! Cell list
IF( CellListLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLength(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLength(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLength(9)
  CALL MakeList( BoxCutoff, pPosition, BoxLengthInverse )
END IF

RETURN

END SUBROUTINE SimpleCubicConfiguration

! *********************************************************************************************** !
!      This subroutine allocates particles according to the body-centred cubic configuration      !
! *********************************************************************************************** !
SUBROUTINE BodyCentredCubicConfiguration(  )

! Uses one module: linked lists
USE LinkedLists, ONLY: MakeList

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nCells              ! Number of unit cells
INTEGER( Kind= Int64 ) :: iCell, jCell, kCell ! Counter of cells
INTEGER( Kind= Int64 ) :: iParticle           ! Counter of particles
INTEGER( Kind= Int64 ) :: Counter             ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: CellLength             ! Length of unit cell (cubic structure)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: ScalingDistanceUnitBox ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BoxCutoff              ! Box cutoff (x-, y-, and z-directions)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AxisSelection(1) ) THEN
  BodyFixedAxis = xAxis
ELSE IF( AxisSelection(2) ) THEN
  BodyFixedAxis = yAxis
ELSE IF( AxisSelection(3) ) THEN
  BodyFixedAxis = zAxis
END IF

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
pQuaternion(0,:) = DCOS( QuaternionAngle * 0.5D0 )                    ! Real part
pQuaternion(1,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(3) ! Imaginary part (Vector)

! Number of unit cells per axis (body-centered cube)
nCells = NINT( ( 0.5D0 * DBLE( nParticles ) ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (body-centered cube)
CellLength = ( 2.D0 / TotalNumberDensity ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BoxLength    = 0.D0
BoxLength(1) = CellLength * DBLE( nCells )
BoxLength(5) = CellLength * DBLE( nCells )
BoxLength(9) = CellLength * DBLE( nCells )

! Simulation box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Positioning of particles (centers of mass)
Counter = 1
DO iCell = 1, nCells
  DO jCell = 1, nCells
    DO kCell = 1, nCells
      ! Particles on the right vertex of unit cell
      pPosition(1,Counter) = DBLE( iCell - 1 ) * CellLength
      pPosition(2,Counter) = DBLE( jCell - 1 ) * CellLength
      pPosition(3,Counter) = DBLE( kCell - 1 ) * CellLength
      Counter = Counter + 1
      ! Particles on the center of unit cell
      pPosition(1,Counter) = ( DBLE( iCell ) - 0.5D0 ) * CellLength
      pPosition(2,Counter) = ( DBLE( jCell ) - 0.5D0 ) * CellLength
      pPosition(3,Counter) = ( DBLE( kCell ) - 0.5D0 ) * CellLength
      Counter = Counter + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO iParticle = 1, nParticles
  ! Spatial transformation
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,iParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,iParticle) )
END DO

! Cell list
IF( CellListLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLength(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLength(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLength(9)
  CALL MakeList( BoxCutoff, pPosition, BoxLengthInverse )
END IF

RETURN

END SUBROUTINE BodyCentredCubicConfiguration

! *********************************************************************************************** !
!      This subroutine allocates particles according to the face-centred cubic configuration      !
! *********************************************************************************************** !
SUBROUTINE FaceCentredCubicConfiguration(  )

! Uses one module: linked lists
USE LinkedLists, ONLY: MakeList

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: nCells              ! Number of unit cells
INTEGER( Kind= Int64 ) :: iCell, jCell, kCell ! Counter of cells
INTEGER( Kind= Int64 ) :: iParticle           ! Counter of particles
INTEGER( Kind= Int64 ) :: Counter             ! Counter

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: CellLength             ! Length of unit cell (cubic structure)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: ScalingDistanceUnitBox ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BoxCutoff              ! Box cutoff (x-, y-, and z-directions)

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AxisSelection(1) ) THEN
  BodyFixedAxis = xAxis
ELSE IF( AxisSelection(2) ) THEN
  BodyFixedAxis = yAxis
ELSE IF( AxisSelection(3) ) THEN
  BodyFixedAxis = zAxis
END IF

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
pQuaternion(0,:) = DCOS( QuaternionAngle * 0.5D0 )                    ! Real part
pQuaternion(1,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(3) ! Imaginary part (Vector)

! Number of unit cells per axis (face-centered cube)
nCells = NINT( ( 0.25D0 * DBLE( nParticles ) ) ** ( 1.D0 / 3.D0 ) )

! Unit cell length (face-centered cube)
CellLength = ( 4.D0 / TotalNumberDensity ) ** ( 1.D0 / 3.D0 )

! Simulation box length
BoxLength    = 0.D0
BoxLength(1) = CellLength * DBLE( nCells )
BoxLength(5) = CellLength * DBLE( nCells )
BoxLength(9) = CellLength * DBLE( nCells )

! Simulation box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Positioning of particles (centers of mass)
Counter = 1
DO iCell = 1, nCells
  DO jCell = 1, nCells
    DO kCell = 1, nCells
      ! Particles on the right vertex of unit cell
      pPosition(1,Counter) = DBLE( iCell - 1 ) * CellLength
      pPosition(2,Counter) = DBLE( jCell - 1 ) * CellLength
      pPosition(3,Counter) = DBLE( kCell - 1 ) * CellLength
      Counter = Counter + 1
      ! Particles on the front face of unit cell
      pPosition(1,Counter) = DBLE( iCell - 1 ) * CellLength
      pPosition(2,Counter) = ( DBLE( jCell ) - 0.5D0 ) * CellLength
      pPosition(3,Counter) = ( DBLE( kCell ) - 0.5D0 ) * CellLength
      Counter = Counter + 1
      ! Particles on the left face of unit cell
      pPosition(1,Counter) = ( DBLE( iCell ) - 0.5D0 ) * CellLength
      pPosition(2,Counter) = DBLE( jCell - 1 ) * CellLength
      pPosition(3,Counter) = ( DBLE( kCell ) - 0.5D0 ) * CellLength
      Counter = Counter + 1
      ! Particles on the lower face of unit cell
      pPosition(1,Counter) = ( DBLE( iCell ) - 0.5D0 ) * CellLength
      pPosition(2,Counter) = ( DBLE( jCell ) - 0.5D0 ) * CellLength
      pPosition(3,Counter) = DBLE( kCell - 1 ) * CellLength
      Counter = Counter + 1
    END DO
  END DO
END DO

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO iParticle = 1, nParticles
  ! Spatial transformation
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,iParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,iParticle) )
END DO

! Cell list
IF( CellListLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLength(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLength(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLength(9)
  CALL MakeList( BoxCutoff, pPosition, BoxLengthInverse )
END IF

RETURN

END SUBROUTINE FaceCentredCubicConfiguration

! *********************************************************************************************** !
! This subroutine allocates particles according to a random orthorhombic/triclinic configuration  !
! *********************************************************************************************** !
SUBROUTINE RandomConfiguration(  )

! Uses two modules: linked lists and overlap check
USE LinkedLists
USE OverlapCheck

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: bEdge                  ! Counter (box edges)
INTEGER( Kind= Int64 ) :: iParticle, pParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iComponent, cComponent ! Counters (component)
INTEGER( Kind= Int64 ) :: BoxMatrixComponent     ! Box matrix component
INTEGER( Kind= Int64 ) :: OverlappingParticles   ! Counter of overlapping particles

! *********************************************************************************************** !
! INTEGER VARIABLES (PARAMETER)                                                                   !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER :: FixPFraction = 1 ! Control variable for the algorithm that fixes the packing fraction of the system

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pCycle                             ! Counter of cycles
INTEGER( Kind= Int64 ) :: nAttempts                          ! Counter
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation             ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation                ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nAcceptanceIsotropicVolumeChange   ! Move acceptance counter: Isotropic volume scaling
INTEGER( Kind= Int64 ) :: nAcceptanceAnisotropicVolumeChange ! Move acceptance counter: Anistropic volume scaling
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter        ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter           ! Move counter (Rotation)
INTEGER( Kind= Int64 ) :: nMovementIsoVolumeChangeCounter    ! Move counter (Isotropic volume scaling)
INTEGER( Kind= Int64 ) :: nMovementAnisoVolumeChangeCounter  ! Move counter (Anisotropic volume scaling)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: ContactDistance                  ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 )                   :: BoxVolumeNVT                     ! Box volume (NVT simulation)
REAL( Kind= Real64 )                   :: PackingFractionNVT               ! Packing fraction (NVT simulation)
REAL( Kind= Real64 )                   :: BoxVolumeRandomConfiguration     ! Box volume (NPT simulation)
REAL( Kind= Real64 )                   :: OldBoxVolume, NewBoxVolume       ! Box volume (before/after a trial move)
REAL( Kind= Real64 )                   :: VolumeScalingFactor              ! Scaling factor of the volume of the simulation box
REAL( Kind= Real64 )                   :: EnthalpyChange                   ! Enthalpy change between microstates (reduced)
REAL( Kind= Real64 )                   :: PackingFractionNPT               ! Packing fraction (NPT simulation)
REAL( Kind= Real64 )                   :: Ratio                            ! Acceptance ratio (simulation)
REAL( Kind= Real64 )                   :: CurrentBoxDistortion             ! Box distortion
REAL( Kind= Real64 )                   :: MaxTranslationalDisplacement     ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                   :: MaxAngularDisplacement           ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                   :: MaxIsoVolumetricDisplacement     ! Maximum displacement [+/-] (Isotropic volume scaling)
REAL( Kind= Real64 )                   :: MaxAnisoVolumetricDisplacement   ! Maximum displacement [+/-] (Anisotropic volume scaling)
REAL( Kind= Real64 )                   :: ScalingBoxVolume                 ! Scaled volume of the simulation box
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxCutoff                        ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthNVT                     ! Box length (NVT simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthInverseNVT              ! Inverse of box length (NVT simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthMC                      ! Box length (NPT simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: BoxLengthInverseMC               ! Inverse of box length (NPT simulation)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: OldBoxLength, NewBoxLength       ! Box length (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: OldBoxLengthInverse              ! Inverse of box length (before a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: NewBoxLengthInverse              ! Inverse of box length (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxVectorAngle                   ! Cossine of angle between box vectors
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxEdgeLength                    ! Length of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxEdgeRatio                     ! Length ratio of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox           ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOldPosition, iNewPosition       ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOldOrientation, iNewOrientation ! Orientation of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iOldQuaternion, iNewQuaternion   ! Quaternion of particle (before/after a trial move)

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: PositionSaveMC ! Old position of all particles

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap                          ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CheckBoxDistortion               ! Detects if a box shape deformation is valid or not : TRUE = ignore box shape deformation; FALSE = consider box shape deformation
LOGICAL :: LatticeReductionLogical          ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: MovementRotationLogical          ! Rotational move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical       ! Translational movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementIsoVolumeChangeLogical   ! Isotropic volume scaling selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementAnisoVolumeChangeLogical ! Anisotropic volume scaling selection : TRUE = movement selected; FALSE = movement not selected

! *********************************************************************************************** !
! LOGICAL VARIABLES (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
LOGICAL, DIMENSION( : ), ALLOCATABLE :: OverlapCounterLogical ! Checks how many particles are overlapping each other during the hit-and-miss algorithm

! Allocation
ALLOCATE( PositionSaveMC(3,nParticles) )
ALLOCATE( OverlapCounterLogical(nParticles) )

! Initialization
VolumeScalingFactor   = 1.D0
OverlapCounterLogical = .TRUE.
OverlappingParticles  = COUNT( OverlapCounterLogical, Dim= 1 ) - 1

! Chosen unrotated reference (x-, y-, or z-axis)
IF( AxisSelection(1) ) THEN
  BodyFixedAxis = xAxis
ELSE IF( AxisSelection(2) ) THEN
  BodyFixedAxis = yAxis
ELSE IF( AxisSelection(3) ) THEN
  BodyFixedAxis = zAxis
END IF

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! *********************************************************************************************** !
! Quaternion Algebra                                                                              !
! *********************************************************************************************** !
!  See 'Quaternion algebras (2021)' book by John Voight.                                          !
!  Available at: <https://math.dartmouth.edu/~jvoight/quat-book.pdf>.                             !
! *********************************************************************************************** !
pQuaternion(0,:) = DCOS( QuaternionAngle * 0.5D0 )                    ! Real part
pQuaternion(1,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:) = DSIN( QuaternionAngle * 0.5D0 ) * BodyFixedAxis(3) ! Imaginary part (Vector)

! Box length (cube)
BoxLength    = 0.D0
BoxLength(1) = BoxVolume ** (1.D0 / 3.D0)
BoxLength(5) = BoxVolume ** (1.D0 / 3.D0)
BoxLength(9) = BoxVolume ** (1.D0 / 3.D0)

! Simulation box length (inverse)
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Active transformation (orientation)
DO iParticle = 1, nParticles
  CALL VectorRotation( zAxis, pQuaternion(:,iParticle), pOrientation(:,iParticle) )
END DO

! Packing fraction (NVT simulation)
PackingFractionNVT = PackingFractionInitialConfiguration

! Box length (NVT simulation)
BoxVolumeNVT = TotalParticleVolume / PackingFractionNVT

! Box length (NVT simulation)
BoxLengthNVT    = 0.D0
BoxLengthNVT(1) = (BoxVolumeNVT) ** (1.D0 / 3.D0)
BoxLengthNVT(5) = (BoxVolumeNVT) ** (1.D0 / 3.D0)
BoxLengthNVT(9) = (BoxVolumeNVT) ** (1.D0 / 3.D0)

! Inverse of box length (NVT simulation)
CALL InverseMatrixCofactorVec( BoxLengthNVT, BoxLengthInverseNVT, BoxVolumeNVT )

! Positioning of particles (centers of mass)
pPosition = 0.D0

! *********************************************************************************************** !
! Monte Carlo parameters (NVT simulation)                                                         !
! *********************************************************************************************** !
MovementTranslationLogical  = .FALSE.      ! Translational move selector           (initial value)
MovementRotationLogical     = .FALSE.      ! Rotational move selector              (initial value)
nAcceptanceTranslation      = 0            ! Translational move acceptance counter (initial value)
nAcceptanceRotation         = 0            ! Rotational move acceptance counter    (initial value)
nMovementTranslationCounter = 0            ! Translational move counter            (initial value)
nMovementRotationCounter    = 0            ! Rotational move counter               (initial value)
nAttempts                   = 0            ! Number of attempts                    (initial value)
MaxAngularDisplacement      = 0.1D0        ! Maximum rotational displacement       (initial value)
pQuaternionMC               = pQuaternion  ! Quaternion algebra                    (initial value)
pPositionMC                 = pPosition    ! Position of particles                 (initial value)
pOrientationMC              = pOrientation ! Orientation of particles              (initial value)

! Maximum translational displacement (initial value)
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  IF( MAXVAL( cDiameter ) <= MAXVAL( cLength ) ) THEN
    MaxTranslationalDisplacement = 1.05D0 * MAXVAL( cLength )   
  ELSE
    MaxTranslationalDisplacement = 1.05D0 * MAXVAL( cDiameter )
  END IF
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  MaxTranslationalDisplacement = 1.05D0 * ( MAXVAL( cLength ) + MAXVAL( cDiameter ) )
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  MaxTranslationalDisplacement = 1.05D0 * ( MAXVAL( cLength ) + MAXVAL( cDiameter ) )
END IF

! Summary
WRITE( *, "(G0)" ) "Attempting to randomly distribute particles inside a cubic box. It may take a while..."
WRITE( *, "(G0)" ) " "
CALL Sleep( 1 )

! Open initial configuration file
OPEN( Unit= 55, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/initconf_rnd_"// &
&                     TRIM( DescriptorFileGeometry )//".xyz", Status= "REPLACE" )

! *********************************************************************************************** !
! Hit-and-miss (NVT simulation)                                                                   !
! *********************************************************************************************** !
HitAndMissNVT: DO

  ! Loop over particles
  DO pCycle = 1, nParticles

    ! Rewind output unit
    REWIND( 55 )

    ! Component index
    IF( nComponents > 1 ) THEN
      ! Pseudorandom number generator (uniform distribution)
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      ! Get component index
      iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
    ELSE
      ! Get component index
      iComponent = 1
    END IF

    ! Avoid components with molar fraction of 0
    DO WHILE( cParticles(iComponent) < 1 )
      ! Component index
      IF( nComponents > 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        ! Get component index
        iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
      ELSE
        ! Get component index
        iComponent = 1
      END IF
    END DO

    ! Forbid rotation if component is spherical
    IF( SphericalComponentLogical(iComponent) ) THEN
      MovementTranslationLogical  = .TRUE.  ! Enable translation
      MovementRotationLogical     = .FALSE. ! Disable rotation
      nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
    ! Allow rotation if component is nonspherical
    ELSE
      ! Pseudorandom number generator (uniform distribution)
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      ! Translation criterion
      IF( RandomNumber < TranslationalProbabilityRandomConfig ) THEN
        MovementTranslationLogical  = .TRUE.  ! Enable translation
        MovementRotationLogical     = .FALSE. ! Disable rotation
        nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
      ! Rotation criterion
      ELSE IF( RandomNumber >= TranslationalProbabilityRandomConfig ) THEN
        MovementRotationLogical    = .TRUE.  ! Enable rotation
        MovementTranslationLogical = .FALSE. ! Disable translation
        nMovementRotationCounter   = nMovementRotationCounter + 1 ! Increment move counter
      END IF
    END IF

    ! Pseudorandom number generator (uniform distribution)
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
    ! Random selection of particles of component i
    iParticle = SUM( cParticles(0:(iComponent-1)) ) + INT( RandomNumber * DBLE( cParticles(iComponent) ) ) + 1

    ! Assignment of previous configuration (microstate m)
    iOldPosition(:)    = pPositionMC(:,iParticle)    ! Old position
    iOldQuaternion(:)  = pQuaternionMC(:,iParticle)  ! Old quaternion
    iOldOrientation(:) = pOrientationMC(:,iParticle) ! Old orientation

    ! Translational movement
    IF( MovementTranslationLogical ) THEN
      ! Random translation along x-axis
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
      ! Random translation along y-axis
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
      ! Random translation along z-axis
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
      ! Minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverseNVT, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLengthNVT, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
    ! Disable translation
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition = iOldPosition
    END IF

    ! Rotational movement
    IF( MovementRotationLogical ) THEN
      ! Random composed unit quaternion
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
      ! Active transformation
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ! Disable rotation
    ELSE IF( .NOT. MovementRotationLogical ) THEN
      iNewQuaternion  = iOldQuaternion
      iNewOrientation = iOldOrientation
    END IF

    ! Overlap check
    CALL ParticleOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
    &                          BoxLengthNVT, BoxLengthInverseNVT, Overlap )

    ! Acceptance criterion
    IF( .NOT. Overlap ) THEN
      ! System configuration update
      pPositionMC(:,iParticle)    = iNewPosition(:)    ! Update position
      pQuaternionMC(:,iParticle)  = iNewQuaternion(:)  ! Update quaternion
      pOrientationMC(:,iParticle) = iNewOrientation(:) ! Update orientation
      ! Update counter of overlapping configurations
      IF( OverlapCounterLogical(iParticle) ) THEN
        OverlapCounterLogical(iParticle) = .FALSE.
        OverlappingParticles = OverlappingParticles - 1
        IF( OverlappingParticles < 0 ) OverlappingParticles = 0
      END IF
      ! Displacement counter update
      IF( MovementTranslationLogical ) THEN
        nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
      ELSE IF ( MovementRotationLogical ) THEN
        nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
      END IF
    ELSE
      ! Retrieve old configuration
      pPositionMC(:,iParticle)    = iOldPosition(:)    ! Retrieve old position
      pQuaternionMC(:,iParticle)  = iOldQuaternion(:)  ! Retrieve old quaternion
      pOrientationMC(:,iParticle) = iOldOrientation(:) ! Retrieve old orientation
    END IF

  END DO

  ! Adjustment of maximum displacement (translation and rotation)
  IF( MOD( (nMovementTranslationCounter + nMovementRotationCounter), (nAdjustmentMovementRandomConfig * nParticles) ) == 0 ) THEN

    ! Acceptance ratio (translation)
    IF( nMovementTranslationCounter > 0 ) THEN
      Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
      ! Translational adjustment
      IF( Ratio <= AcceptanceRatioTranslation ) THEN
        MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
      ELSE
        MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
      END IF
      ! Set minimum translational displacement (arbitrary)
      IF( MAXVAL( cDiameter ) <= MAXVAL( cLength ) ) THEN
        IF( MaxTranslationalDisplacement <= MAXVAL( cLength ) ) MaxTranslationalDisplacement = MAXVAL( cLength )
      ELSE
        IF( MaxTranslationalDisplacement <= MAXVAL( cDiameter ) ) MaxTranslationalDisplacement = MAXVAL( cDiameter )
      END IF
    END IF

    ! Acceptance ratio (rotation)
    IF( nMovementRotationCounter > 0 ) THEN
      Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
      ! Rotational adjustment
      IF( Ratio <= AcceptanceRatioRotation ) THEN
        MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
      ELSE
        MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
      END IF
      ! Avoid multiple turns (arbitrary)
      IF( ( MaxAngularDisplacement > 4.D0 * cPi ) ) THEN
        MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
      END IF
    END IF

    ! Reset counters
    nAcceptanceTranslation      = 0
    nMovementTranslationCounter = 0
    nAcceptanceRotation         = 0
    nMovementRotationCounter    = 0

  END IF

  ! Iteration
  nAttempts = nAttempts + 1

  ! Initial configuration (partial)
  WRITE( 55, "(G0)" ) nParticles
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 55, DescriptorString ) 'Lattice="', BoxLengthNVT(1:9), '" Origin="', -0.5D0 * ( BoxLengthNVT(1) + BoxLengthNVT(4) + &
  &                             BoxLengthNVT(7) ), -0.5D0 * ( BoxLengthNVT(2) + BoxLengthNVT(5) + BoxLengthNVT(8) ), -0.5D0 * &
  &                             ( BoxLengthNVT(3) + BoxLengthNVT(6) + BoxLengthNVT(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  FLUSH( 55 )

  ! Progress bar
  CALL ProgressBarHitAndMiss( nAttempts, OverlappingParticles )

  ! Check number of overlapping particles
  IF( OverlappingParticles == 0 ) THEN
    ! Possible initial configuration
    EXIT HitAndMissNVT
  ELSE
    ! Attempt another configuration
    CYCLE HitAndMissNVT
  END IF

END DO HitAndMissNVT

! Initial configuration (partial)
REWIND( 55 )
WRITE( 55, "(G0)" ) nParticles
DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
WRITE( 55, DescriptorString ) 'Lattice="', BoxLengthNVT(1:9), '" Origin="', -0.5D0 * ( BoxLengthNVT(1) + BoxLengthNVT(4) + &
&                             BoxLengthNVT(7) ), -0.5D0 * ( BoxLengthNVT(2) + BoxLengthNVT(5) + BoxLengthNVT(8) ), -0.5D0 * &
&                             ( BoxLengthNVT(3) + BoxLengthNVT(6) + BoxLengthNVT(9) ), '" ', &
&                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
END IF

! Close output unit
CLOSE( 55 )

! Reassign positions, rotation quaternions, and orientations of all particles
pQuaternion  = pQuaternionMC  ! Quaternion of particles
pPosition    = pPositionMC    ! Position of particles
pOrientation = pOrientationMC ! Orientation of particles

! Deallocation
DEALLOCATE( OverlapCounterLogical )

! Summary
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( nAttempts == 1 ) THEN
  WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", nAttempts, " attempt."
ELSE IF( nAttempts > 1 ) THEN
  WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", nAttempts, " attempts."
END IF
CALL Sleep( 1 )
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Now running an NPT simulation up to the target packing fraction of ", PackingFraction, "..."
WRITE( *, "(G0)" ) " "
CALL Sleep( 1 )

! *********************************************************************************************** !
! Monte Carlo parameters (NPT simulation)                                                         !
! *********************************************************************************************** !
MovementTranslationLogical         = .FALSE.                                    ! Translational move selector           (initial value)
MovementRotationLogical            = .FALSE.                                    ! Rotational move selector              (initial value)
MovementIsoVolumeChangeLogical     = .FALSE.                                    ! Volume move selector                  (initial value)
MovementAnisoVolumeChangeLogical   = .FALSE.                                    ! Volume move selector                  (initial value)
nAcceptanceTranslation             = 0                                          ! Translational move acceptance counter (initial value)
nAcceptanceRotation                = 0                                          ! Rotational move acceptance counter    (initial value)
nAcceptanceIsotropicVolumeChange   = 0                                          ! Volumetric move acceptance counter    (initial value)
nAcceptanceAnisotropicVolumeChange = 0                                          ! Volumetric move acceptance counter    (initial value)
nMovementTranslationCounter        = 0                                          ! Translational move counter            (initial value)
nMovementRotationCounter           = 0                                          ! Rotational move counter               (initial value)
nMovementIsoVolumeChangeCounter    = 0                                          ! Volume scaling counter                (initial value)
nMovementAnisoVolumeChangeCounter  = 0                                          ! Volume scaling counter                (initial value)
nAttempts                          = 0                                          ! Number of attempts                    (initial value)
pQuaternionMC                      = pQuaternion                                ! Quaternion algebra                    (initial value)
pPositionMC                        = pPosition                                  ! Position of particles                 (initial value)
pOrientationMC                     = pOrientation                               ! Orientation of particles              (initial value)
BoxLengthMC                        = BoxLengthNVT                               ! Box length                            (initial value)
BoxLengthInverseMC                 = BoxLengthInverseNVT                        ! Inverse of box length                 (initial value)
BoxVolumeRandomConfiguration       = BoxVolumeNVT                               ! Box volume                            (initial value)
PackingFractionNPT                 = PackingFractionNVT                         ! Packing fraction                      (initial value)
MaxTranslationalDisplacement       = MaxTranslationalDisplacementRandomConfig   ! Maximum translational displacement    (initial value)
MaxAngularDisplacement             = MaxAngularDisplacementRandomConfig         ! Maximum rotational displacement       (initial value)
MaxIsoVolumetricDisplacement       = MaxIsoVolumetricDisplacementRandomConfig   ! Maximum isovolumetric displacement    (initial value)
MaxAnisoVolumetricDisplacement     = MaxAnisoVolumetricDisplacementRandomConfig ! Maximum anisovolumetric displacement  (initial value)

! Initialize cell list
IF( CellListLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLengthMC(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLengthMC(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLengthMC(9)
  OldBoxLength = BoxLengthMC
  CALL MakeList( BoxCutoff, pPositionMC, BoxLengthInverseMC )
END IF

! Open initial configuration file
OPEN( Unit= 55, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/initconf_rnd_"// &
&                     TRIM( DescriptorFileGeometry )//".xyz", Status= "REPLACE" )

! Isothermal-isobaric Monte Carlo simulation
SimulationNPT: DO

  ! Rewind output unit
  REWIND( 55 )

  ! Prepare simulation box for single-particle moves
  IF( CellListLogical .AND. (MovementAnisoVolumeChangeLogical .OR. MovementIsoVolumeChangeLogical) ) THEN ! Check cells only after a volume change
    CALL BoxCheckNPT( pPositionMC, OldBoxLength, BoxLengthMC, BoxLengthInverseMC, .FALSE. )
  END IF

  ! Choose between displacement of molecules or volume scaling
  IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
  IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
  IF( RandomNumber < MovementProbabilityRandomConfig ) THEN
    MovementTranslationLogical       = .TRUE.  ! Enable translation
    MovementRotationLogical          = .TRUE.  ! Enable rotation
    MovementIsoVolumeChangeLogical   = .FALSE. ! Disable volume scaling
    MovementAnisoVolumeChangeLogical = .FALSE. ! Disable volume scaling
  ELSE IF( RandomNumber >= MovementProbabilityRandomConfig ) THEN
    MovementTranslationLogical       = .FALSE. ! Disable translation
    MovementRotationLogical          = .FALSE. ! Disable rotation
    MovementIsoVolumeChangeLogical   = .TRUE.  ! Enable volume scaling
    MovementAnisoVolumeChangeLogical = .TRUE.  ! Enable volume scaling
  END IF

  ! Movement (translation or rotation)
  IF( MovementTranslationLogical .OR. MovementRotationLogical ) THEN

    ! Loop over particles
    DO pCycle = 1, nParticles

      ! Component index
      IF( nComponents > 1 ) THEN
        ! Pseudorandom number generator (uniform distribution)
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        ! Get component index
        iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
      ELSE
        ! Get component index
        iComponent = 1
      END IF

      ! Avoid components with molar fraction of 0
      DO WHILE( cParticles(iComponent) < 1 )
        ! Component index
        IF( nComponents > 1 ) THEN
          ! Pseudorandom number generator (uniform distribution)
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          ! Get component index
          iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
        ELSE
          ! Get component index
          iComponent = 1
        END IF
      END DO

      ! Forbid rotation if component is spherical
      IF( SphericalComponentLogical(iComponent) ) THEN
        MovementTranslationLogical  = .TRUE.  ! Enable translation
        MovementRotationLogical     = .FALSE. ! Disable rotation
        nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
      ! Allow rotation if component is nonspherical
      ELSE
        ! Pseudorandom number generator (uniform distribution)
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        ! Translation criterion
        IF( RandomNumber < TranslationalProbabilityRandomConfig ) THEN
          MovementTranslationLogical  = .TRUE.  ! Enable translation
          MovementRotationLogical     = .FALSE. ! Disable rotation
          nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
        ! Rotation criterion
        ELSE IF( RandomNumber >= TranslationalProbabilityRandomConfig ) THEN
          MovementRotationLogical    = .TRUE.  ! Enable rotation
          MovementTranslationLogical = .FALSE. ! Disable translation
          nMovementRotationCounter   = nMovementRotationCounter + 1 ! Increment move counter
        END IF
      END IF

      ! Pseudorandom number generator (uniform distribution)
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      ! Random selection of particles of component cComponent
      iParticle = SUM( cParticles(0:(iComponent-1)) ) + INT( RandomNumber * DBLE( cParticles(iComponent) ) ) + 1

      ! Assignment of previous configuration (microstate m)
      iOldPosition(:)    = pPositionMC(:,iParticle)    ! Old position
      iOldQuaternion(:)  = pQuaternionMC(:,iParticle)  ! Old quaternion
      iOldOrientation(:) = pOrientationMC(:,iParticle) ! Old orientation

      ! Translational movement
      IF( MovementTranslationLogical ) THEN
        ! Random translation along x-axis
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along y-axis
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along z-axis
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
      ! Disable translation
      ELSE IF( .NOT. MovementTranslationLogical ) THEN
        iNewPosition = iOldPosition
      END IF

      ! Rotational movement
      IF( MovementRotationLogical ) THEN
        ! Random composed unit quaternion
        CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
        ! Active transformation
        CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
      ! Disable rotation
      ELSE IF( .NOT. MovementRotationLogical ) THEN
        iNewQuaternion  = iOldQuaternion
        iNewOrientation = iOldOrientation
      END IF

      ! Overlap check after displacement of a particle
      IF( .NOT. CellListControl ) THEN
        ! Whole system
        CALL ParticleOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
        &                          BoxLengthMC, BoxLengthInverseMC, Overlap )
      ELSE
        ! Linked lists
        CALL ListOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
        &                      BoxLengthMC, BoxLengthInverseMC, Overlap, .FALSE. )
      END IF

      ! Acceptance criterion
      IF( .NOT. Overlap ) THEN
        ! System configuration update
        pPositionMC(:,iParticle)    = iNewPosition(:)    ! Update position
        pQuaternionMC(:,iParticle)  = iNewQuaternion(:)  ! Update quaternion
        pOrientationMC(:,iParticle) = iNewOrientation(:) ! Update orientation
        ! Displacement counter update
        IF( MovementTranslationLogical ) THEN
          IF( CellListControl ) CALL ParticleTranslationNVT( iParticle, ScalingDistanceUnitBox, .FALSE. ) ! Update cell
          nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
        ELSE IF ( MovementRotationLogical ) THEN
          nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
        END IF
      ELSE
        ! Retrieve old configuration
        pPositionMC(:,iParticle)    = iOldPosition(:)    ! Retrieve old position
        pQuaternionMC(:,iParticle)  = iOldQuaternion(:)  ! Retrieve old quaternion
        pOrientationMC(:,iParticle) = iOldOrientation(:) ! Retrieve old orientation
      END IF

    END DO

  ! Volume scaling (isotropic or anisotropic)
  ELSE IF( MovementIsoVolumeChangeLogical .OR. MovementAnisoVolumeChangeLogical ) THEN

    ! Assignment of previous configuration (microstate m)
    OldBoxLength        = BoxLengthMC                  ! Box length
    OldBoxLengthInverse = BoxLengthInverseMC           ! Box length (inverse)
    OldBoxVolume        = BoxVolumeRandomConfiguration ! Box volume

    ! Expansion/compression type
    IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
    IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

    ! Isotropic volume scaling
    IF( RandomNumber < IsoVolumetricProbabilityRandomConfig ) THEN
      ! Random scaling factor
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      ! Random walk on the logarithm of the volume
      VolumeScalingFactor = DLOG( OldBoxVolume ) + (RandomNumber - 0.5D0) * MaxIsoVolumetricDisplacement
      VolumeScalingFactor = DEXP( VolumeScalingFactor )
      VolumeScalingFactor = (VolumeScalingFactor / OldBoxVolume) ** (1.D0 / 3.D0)
      ! Proportional box length
      NewBoxLength = OldBoxLength * VolumeScalingFactor
      CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
      ! Movement counter
      nMovementIsoVolumeChangeCounter = nMovementIsoVolumeChangeCounter + 1
      ! Movement type
      MovementIsoVolumeChangeLogical   = .TRUE.  ! Enable isotropic volume scaling
      MovementAnisoVolumeChangeLogical = .FALSE. ! Disable anisotropic volume scaling
    ! Anisotropic volume scaling
    ELSE IF( RandomNumber >= IsoVolumetricProbabilityRandomConfig ) THEN
      ! Random box component
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      BoxMatrixComponent = INT( RandomNumber * 6.D0 ) + 1
      IF( BoxMatrixComponent == 1 ) THEN
        BoxMatrixComponent = 1 ! XX component
      ELSE IF( BoxMatrixComponent == 2 ) THEN
        BoxMatrixComponent = 4 ! YX component
      ELSE IF( BoxMatrixComponent == 3 ) THEN
        BoxMatrixComponent = 5 ! YY component
      ELSE IF( BoxMatrixComponent == 4 ) THEN
        BoxMatrixComponent = 7 ! ZX component
      ELSE IF( BoxMatrixComponent == 5 ) THEN
        BoxMatrixComponent = 8 ! ZY component
      ELSE IF( BoxMatrixComponent == 6 ) THEN
        BoxMatrixComponent = 9 ! ZZ component
      END IF
      NewBoxLength = OldBoxLength
      ! Random factor
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
      NewBoxLength(BoxMatrixComponent) = OldBoxLength(BoxMatrixComponent) + MaxAnisoVolumetricDisplacement * (RandomNumber - 0.5D0)
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
      ! Movement counter
      nMovementAnisoVolumeChangeCounter = nMovementAnisoVolumeChangeCounter + 1
      ! Movement type
      MovementIsoVolumeChangeLogical   = .FALSE. ! Disable isotropic volume scaling
      MovementAnisoVolumeChangeLogical = .TRUE.  ! Enable anisotropic volume scaling
    END IF

    ! Reset condition of anisotropic volume scaling (we must be careful not to induce a bias in the system)
    CheckBoxDistortion = .FALSE.

    ! Condition of anisotropic volume scaling (box distortion)
    IF( MovementAnisoVolumeChangeLogical ) THEN
      ! Box length
      BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(1:3) ) )
      BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(4:6) ) )
      BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( NewBoxLength(7:9), NewBoxLength(7:9) ) )
      ! Length ratio
      BoxEdgeRatio(1) = BoxEdgeLength(1) / BoxEdgeLength(2)
      BoxEdgeRatio(2) = BoxEdgeLength(1) / BoxEdgeLength(3)
      BoxEdgeRatio(3) = BoxEdgeLength(2) / BoxEdgeLength(3)
      ! Angle between box vectors
      BoxVectorAngle(1) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(4:6) ) / ( BoxEdgeLength(1) * BoxEdgeLength(2) )
      BoxVectorAngle(2) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(7:9) ) / ( BoxEdgeLength(1) * BoxEdgeLength(3) )
      BoxVectorAngle(3) = DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(7:9) ) / ( BoxEdgeLength(2) * BoxEdgeLength(3) )
      ! Avoid big distortions of the simulation box
      DO bEdge = 1, 3
        ! Angle distortion
        IF( BoxVectorAngle(bEdge) < DCOS( (cPi / 2.D0) + BoxVectorMaxAngle ) .OR. &
        &   BoxVectorAngle(bEdge) > DCOS( (cPi / 2.D0) - BoxVectorMaxAngle ) ) THEN
          BoxLengthMC                  = OldBoxLength
          BoxLengthInverseMC           = OldBoxLengthInverse
          BoxVolumeRandomConfiguration = OldBoxVolume
          CheckBoxDistortion           = .TRUE.
          EXIT
        END IF
        ! Length distortion
        IF( BoxEdgeRatio(bEdge) > BoxEdgeMaxRatio .OR. BoxEdgeRatio(bEdge) < 1.D0 / BoxEdgeMaxRatio ) THEN
          BoxLengthMC                  = OldBoxLength
          BoxLengthInverseMC           = OldBoxLengthInverse
          BoxVolumeRandomConfiguration = OldBoxVolume
          CheckBoxDistortion           = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    ! Box not too distorted
    IF( .NOT. CheckBoxDistortion ) THEN

      ! Enthalpy change (weighing function)
      EnthalpyChange = ( PressureRandomConfig * ( NewBoxVolume - OldBoxVolume ) ) - ( DBLE( nParticles + 1 ) * &
      &                DLOG( NewBoxVolume / OldBoxVolume ) )

      ! Random number
      IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
      IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )

      ! Enthalpy change criterion
      IF( DEXP( - EnthalpyChange ) >= RandomNumber ) THEN

        ! System configuration
        PositionSaveMC = pPositionMC ! Old configuration

        ! Isotropic volume scaling
        IF( MovementIsoVolumeChangeLogical ) THEN
          ! Rescale positions of particles accordingly
          DO pParticle = 1, nParticles
            pPositionMC(:,pParticle) = pPositionMC(:,pParticle) * VolumeScalingFactor
          END DO
        ! Anisotropic volume scaling
        ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
          ! Rescale positions of particles accordingly
          DO pParticle = 1, nParticles
            ! Transform spatial coordinates using old box dimensions
            CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
            ! New spatial coordinates using new box dimensions
            CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceUnitBox, pPositionMC(:,pParticle) )
          END DO
        END IF

        ! Overlap check after expansion/compression of the simulation box
        IF( .NOT. CellListControl ) THEN
          ! Whole system
          CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        ELSE
          ! Linked lists
          CALL FullListOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap, HalfNeighboursControl )
          ! In case the number of cells in one direction (x, y, or z) becomes less than 3
          IF( .NOT. CellListControl ) CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        END IF

        ! Acceptance criterion
        IF( .NOT. Overlap ) THEN
          ! Assigns the simulation box properties of a trial volume scaling to the system configuration
          BoxVolumeRandomConfiguration = NewBoxVolume        ! Update box volume
          BoxLengthMC                  = NewBoxLength        ! Update box length
          BoxLengthInverseMC           = NewBoxLengthInverse ! Update box length (inverse)
          ! Displacement counter update
          IF( MovementIsoVolumeChangeLogical ) THEN
            nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Isotropic move counter
          ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
            nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Anisotropic move counter
          END IF
          ! Update packing fraction
          PackingFractionNPT = TotalParticleVolume / NewBoxVolume
          ! Re-initialization
          CheckBoxDistortion = .FALSE.
          ! Lattice reduction
          LatticeReductionLogical = .FALSE.
          CALL LatticeReduction( BoxLengthMC, CurrentBoxDistortion, LatticeReductionLogical )
          IF( LatticeReductionLogical ) THEN
            ! Calculate the new reciprocal box basis vectors
            CALL InverseMatrixCofactorVec( BoxLengthMC, BoxLengthInverseMC, BoxVolumeRandomConfiguration )
            DO pParticle = 1, nParticles
              ! Minimum image convention (the spatial distribution of particles remains unchanged)
              CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pParticle), ScalingDistanceUnitBox ) ! Spatial transformation
              ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
              CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, pPositionMC(:,pParticle) ) ! Spatial transformation
            END DO
            ! Check orientation of the box (eliminate box rotations)
            IF( DABS( BoxLengthMC(2) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. DABS( BoxLengthMC(3) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. &
            &   DABS( BoxLengthMC(6) - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
              ! Undo box rotation
              CALL UndoBoxRotation( BoxLengthMC, BoxLengthInverseMC, BoxVolumeRandomConfiguration )
            END IF
          END IF
        ! Retrieve old properties of the system configuration and the simulation box
        ELSE
          BoxVolumeRandomConfiguration = OldBoxVolume        ! Retrieve old box volume
          BoxLengthMC                  = OldBoxLength        ! Retrieve old box length
          BoxLengthInverseMC           = OldBoxLengthInverse ! Retrieve old box length (inverse)
          pPositionMC                  = PositionSaveMC      ! Retrieve old position of particles
        END IF

      ! Retrieve old properties of the simulation box
      ELSE

        BoxVolumeRandomConfiguration = OldBoxVolume        ! Retrieve old box volume
        BoxLengthMC                  = OldBoxLength        ! Retrieve old box length
        BoxLengthInverseMC           = OldBoxLengthInverse ! Retrieve old box length (inverse)

      END IF ! Enthalpy criterion

    END IF ! Box distortion criterion

  END IF

  ! Iteration
  nAttempts = nAttempts + 1

  ! Adjustment of maximum displacement (translation and rotation)
  IF( MOD( nAttempts, nAdjustmentMovementRandomConfig ) == 0 ) THEN

    ! Translational adjustment
    IF( nMovementTranslationCounter > 0 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
      ! Translational adjustment
      IF( Ratio <= AcceptanceRatioTranslation ) THEN
        MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
      ELSE
        MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
      END IF
      ! Reset counter
      nAcceptanceTranslation      = 0
      nMovementTranslationCounter = 0
    END IF

    ! Avoid multiple turns (arbitrary)
    BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( BoxLengthMC(1:3), BoxLengthMC(1:3) ) )
    BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( BoxLengthMC(4:6), BoxLengthMC(4:6) ) )
    BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( BoxLengthMC(7:9), BoxLengthMC(7:9) ) )
    IF( MaxTranslationalDisplacement >= 2.D0 * MAXVAL( BoxEdgeLength ) ) THEN
      MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxEdgeLength )
    END IF

    ! Rotational adjustment
    IF( nMovementRotationCounter > 0 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
      ! Rotational adjustment
      IF( Ratio <= AcceptanceRatioRotation ) THEN
        MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
      ELSE
        MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
      END IF
      ! Reset counter
      nAcceptanceRotation      = 0
      nMovementRotationCounter = 0
    END IF

    ! Avoid multiple turns (arbitrary)
    IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
      MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
    END IF

  END IF

  ! Adjustment of maximum displacement (volume scaling)
  IF( MOD( nAttempts, nAdjustmentVolumeRandomConfig ) == 0 ) THEN

    ! Volumetric adjustment (isotropic)
    IF( nMovementIsoVolumeChangeCounter > 0 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      Ratio = DBLE( nAcceptanceIsotropicVolumeChange ) / DBLE( nMovementIsoVolumeChangeCounter )
      ! Volumetric adjustment
      IF( Ratio <= AcceptanceRatioIsoVolumeChange ) THEN
        MaxIsoVolumetricDisplacement = 0.95D0 * MaxIsoVolumetricDisplacement
      ELSE
        MaxIsoVolumetricDisplacement = 1.05D0 * MaxIsoVolumetricDisplacement
      END IF
      ! Reset counter
      nAcceptanceIsotropicVolumeChange = 0
      nMovementIsoVolumeChangeCounter  = 0
    END IF

    ! Volumetric adjustment (anisotropic)
    IF( nMovementAnisoVolumeChangeCounter > 0 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      Ratio = DBLE( nAcceptanceAnisotropicVolumeChange ) / DBLE( nMovementAnisoVolumeChangeCounter )
      ! Volumetric adjustment
      IF( Ratio <= AcceptanceRatioAnisoVolumeChange ) THEN
        MaxAnisoVolumetricDisplacement = 0.95D0 * MaxAnisoVolumetricDisplacement
      ELSE
        MaxAnisoVolumetricDisplacement = 1.05D0 * MaxAnisoVolumetricDisplacement
      END IF
      ! Reset counter
      nAcceptanceAnisotropicVolumeChange = 0
      nMovementAnisoVolumeChangeCounter  = 0
    END IF

    ! Avoid low volume changes (arbitrary)
    IF( MaxIsoVolumetricDisplacement <= MinVolumetricDisplacementRandomConfig ) THEN
      MaxIsoVolumetricDisplacement = MinVolumetricDisplacementRandomConfig
    END IF
    IF( MaxAnisoVolumetricDisplacement <= MinVolumetricDisplacementRandomConfig ) THEN
      MaxAnisoVolumetricDisplacement = MinVolumetricDisplacementRandomConfig
    END IF

  END IF

  ! Summary
  CALL ProgressBarRandomConfigNPT( nAttempts, PackingFractionNPT, PackingFraction )

  ! Target packing fraction
  IF( PackingFractionNPT >= PackingFraction .AND. nAttempts == 1 ) THEN
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(3G0,G0.7)" ) "Target packing fraction reached after ", nAttempts, " attempt! Final value: ", PackingFractionNPT
    WRITE( *, "(G0)" ) " "
    EXIT SimulationNPT
  ELSE IF( PackingFractionNPT >= PackingFraction .AND. nAttempts > 1 ) THEN
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(3G0,G0.7)" ) "Target packing fraction reached after ", nAttempts, " attempts! Final value: ", PackingFractionNPT
    WRITE( *, "(G0)" ) " "
    EXIT SimulationNPT
  END IF

  ! Initial configuration (partial)
  WRITE( 55, "(G0)" ) nParticles
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 55, DescriptorString ) 'Lattice="', BoxLengthMC(1:9), '" Origin="', -0.5D0 * ( BoxLengthMC(1) + BoxLengthMC(4) + &
  &                             BoxLengthMC(7) ), -0.5D0 * ( BoxLengthMC(2) + BoxLengthMC(5) + BoxLengthMC(8) ), -0.5D0 * &
  &                             ( BoxLengthMC(3) + BoxLengthMC(6) + BoxLengthMC(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              cLength(cComponent)
        END DO
      ELSE
        DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
          &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &              0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  FLUSH( 55 )

END DO SimulationNPT

! Initial configuration (partial)
REWIND( 55 )
WRITE( 55, "(G0)" ) nParticles
DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
WRITE( 55, DescriptorString ) 'Lattice="', BoxLengthMC(1:9), '" Origin="', -0.5D0 * ( BoxLengthMC(1) + BoxLengthMC(4) + &
&                             BoxLengthMC(7) ), -0.5D0 * ( BoxLengthMC(2) + BoxLengthMC(5) + BoxLengthMC(8) ), -0.5D0 * &
&                             ( BoxLengthMC(3) + BoxLengthMC(6) + BoxLengthMC(9) ), '" ', &
&                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
ELSE IF( GeometryType(3) ) THEN ! Cylinders
  DO cComponent = 1, nComponents
    IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              cLength(cComponent)
      END DO
    ELSE
      DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
        &              pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
        &              pQuaternionMC(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
        &              0.5D0 * cDiameter(cComponent)
      END DO
    END IF
  END DO
END IF

! Close output unit
CLOSE( 55 )

! Deallocation
DEALLOCATE( PositionSaveMC )

! *********************************************************************************************** !
! Monte Carlo parameters (NVT Simulation | Fixing the packing fraction)                           !
! *********************************************************************************************** !
nAcceptanceTranslation       = 0                                        ! Translational move acceptance counter (initial value)
nAcceptanceRotation          = 0                                        ! Rotational move acceptance counter    (initial value)
nMovementTranslationCounter  = 0                                        ! Translational move counter            (initial value)
nMovementRotationCounter     = 0                                        ! Rotational move counter               (initial value)
nAttempts                    = 0                                        ! Number of attempts                    (initial value)
MovementTranslationLogical   = .FALSE.                                  ! Translational move selector           (initial value)
MovementRotationLogical      = .FALSE.                                  ! Rotational move selector              (initial value)
MaxTranslationalDisplacement = MaxTranslationalDisplacementRandomConfig ! Maximum translational displacement    (initial value)
MaxAngularDisplacement       = MaxAngularDisplacementRandomConfig       ! Maximum rotational displacement       (initial value)

! Scaled box volume
ScalingBoxVolume = PackingFractionNPT / PackingFraction
ScalingBoxVolume = ScalingBoxVolume ** ( 1.D0 / 3.D0 )

! Scaled box length
BoxLength = BoxLengthMC * ScalingBoxVolume

! Inverse of box length
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Fix packing fraction with a volume expansion
IF( PackingFractionNPT >= PackingFraction ) THEN

  ! Summary
  CALL Sleep( 1 )
  WRITE( *, "(G0,G0.7,G0,G0.7,G0)" ) "Attempting to fix the packing fraction of ", PackingFractionNPT, &
  &                                  " obtained in the NPT simulation to the target value of ", PackingFraction, "..."
  WRITE( *, "(G0)" ) " "
  CALL Sleep( 1 )

  ! Open initial configuration file
  OPEN( Unit= 55, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/initconf_rnd_"// &
  &                     TRIM( DescriptorFileGeometry )//".xyz", Status= "REPLACE" )

  ! Correction of the packing fraction
  PackingFractionCorrection: DO

    ! Rewind output unit
    REWIND( 55 )

    ! Update rotation quaternions and orientations
    pQuaternion  = pQuaternionMC
    pOrientation = pOrientationMC

    ! Update positions of particles accordingly
    DO pParticle = 1, nParticles
      ! Transform spatial coordinates using the box dimensions from the NPT simulation
      CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
      ! New spatial coordinates using the actual box dimensions
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,pParticle) )
    END DO

    ! Initialization
    Overlap = .FALSE.

    ! Iteration
    nAttempts = nAttempts + 1

    ! Summary
    CALL ProgressBarRandomConfigPackingFractionCorrection( nAttempts )

    ! Overlap check after expansion/compression of the simulation box
    CALL FullOverlapCheckFixPFraction( ContactDistance, BoxLength, BoxLengthInverse, Overlap )

    ! Initial configuration (partial)
    WRITE( 55, "(G0)" ) nParticles
    DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
    WRITE( 55, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
    &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
    &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
    &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPosition(1,iParticle), pPosition(2,iParticle), pPosition(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              0.5D0 * cLength(cComponent)
          END DO
        ELSE
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPosition(1,iParticle), pPosition(2,iParticle), pPosition(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              cLength(cComponent)
          END DO
        ELSE
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    ELSE IF( GeometryType(3) ) THEN ! Cylinders
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              cLength(cComponent)
          END DO
        ELSE
          DO iParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 55, * ) cIndex(cComponent), pPositionMC(1,iParticle), pPositionMC(2,iParticle), pPositionMC(3,iParticle), &
            &              pQuaternion(1,iParticle), pQuaternion(2,iParticle), pQuaternion(3,iParticle), &
            &              pQuaternion(0,iParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
            &              0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    END IF
    FLUSH( 55 )

    ! Packing fraction fixed
    IF( .NOT. Overlap ) THEN

      EXIT PackingFractionCorrection

    ! Attempt to fix the packing fraction with a new configuration of particles
    ELSE

      ! Displace particles (constant volume)
      DO pCycle = 1, nParticles

        ! Component index
        IF( nComponents > 1 ) THEN
          ! Pseudorandom number generator (uniform distribution)
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          ! Get component index
          iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
        ELSE
          ! Get component index
          iComponent = 1
        END IF

        ! Avoid components with molar fraction of 0
        DO WHILE( cParticles(iComponent) < 1 )
          ! Component index
          IF( nComponents > 1 ) THEN
            ! Pseudorandom number generator (uniform distribution)
            IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
            IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
            ! Get component index
            iComponent = INT( RandomNumber * DBLE( nComponents ) ) + 1
          ELSE
            ! Get component index
            iComponent = 1
          END IF
        END DO

        ! Forbid rotation if component is spherical
        IF( SphericalComponentLogical(iComponent) ) THEN
          MovementTranslationLogical  = .TRUE.  ! Enable translation
          MovementRotationLogical     = .FALSE. ! Disable rotation
          nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
        ! Allow rotation if component is nonspherical
        ELSE
          ! Pseudorandom number generator (uniform distribution)
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          ! Translation criterion
          IF( RandomNumber < TranslationalProbabilityRandomConfig ) THEN
            MovementTranslationLogical  = .TRUE.  ! Enable translation
            MovementRotationLogical     = .FALSE. ! Disable rotation
            nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
          ! Rotation criterion
          ELSE IF( RandomNumber >= TranslationalProbabilityRandomConfig ) THEN
            MovementRotationLogical    = .TRUE.  ! Enable rotation
            MovementTranslationLogical = .FALSE. ! Disable translation
            nMovementRotationCounter   = nMovementRotationCounter + 1 ! Increment move counter
          END IF
        END IF

        ! Pseudorandom number generator (uniform distribution)
        IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
        IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
        ! Random selection of particles of component i
        iParticle = SUM( cParticles(0:(iComponent-1)) ) + INT( RandomNumber * DBLE( cParticles(iComponent) ) ) + 1

        ! Assignment of previous configuration (microstate m)
        iOldPosition(:)    = pPositionMC(:,iParticle)    ! Old position
        iOldQuaternion(:)  = pQuaternionMC(:,iParticle)  ! Old quaternion
        iOldOrientation(:) = pOrientationMC(:,iParticle) ! Old orientation

        ! Translational movement
        IF( MovementTranslationLogical ) THEN
          ! Random translation along x-axis
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Random translation along y-axis
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Random translation along z-axis
          IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
          IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
          iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Minimum image convention
          CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
        ! Disable translation
        ELSE IF( .NOT. MovementTranslationLogical ) THEN
          iNewPosition = iOldPosition
        END IF

        ! Rotational movement
        IF( MovementRotationLogical ) THEN
          ! Random Composed Unit Quaternion
          CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
          ! Active transformation
          CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
        ! Disable rotation
        ELSE IF( .NOT. MovementRotationLogical ) THEN
          iNewQuaternion  = iOldQuaternion
          iNewOrientation = iOldOrientation
        END IF

        ! Overlap check after displacement of a particle
        CALL ParticleOverlapCheck( iComponent, iParticle, iNewQuaternion, iNewOrientation, iNewPosition, ContactDistance, &
        &                          BoxLengthMC, BoxLengthInverseMC, Overlap )

        ! Acceptance criterion
        IF( .NOT. Overlap ) THEN
          ! System configuration update
          pPositionMC(:,iParticle)    = iNewPosition(:)    ! Update position
          pQuaternionMC(:,iParticle)  = iNewQuaternion(:)  ! Update quaternion
          pOrientationMC(:,iParticle) = iNewOrientation(:) ! Update orientation
          ! Displacement counter update
          IF( MovementTranslationLogical ) THEN
            nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
          ELSE IF ( MovementRotationLogical ) THEN
            nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
          END IF
        ELSE
          ! Retrieve old configuration
          pPositionMC(:,iParticle)    = iOldPosition(:)    ! Retrieve old position
          pQuaternionMC(:,iParticle)  = iOldQuaternion(:)  ! Retrieve old quaternion
          pOrientationMC(:,iParticle) = iOldOrientation(:) ! Retrieve old orientation
        END IF

      END DO

      ! Adjustment of maximum displacement (translation and rotation)
      IF( MOD( (nMovementTranslationCounter + nMovementRotationCounter), &
      &        (nAdjustmentMovementRandomConfig * nParticles) ) == 0 ) THEN

        ! Translational adjustment
        IF( nMovementTranslationCounter > 0 ) THEN
          ! Acceptance ratio (translation)
          Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
          ! Translational adjustment
          IF( Ratio <= AcceptanceRatioTranslation ) THEN
            MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
          ELSE
            MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
          END IF
        END IF

        ! Avoid multiple turns (arbitrary)
        BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( BoxLengthMC(1:3), BoxLengthMC(1:3) ) )
        BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( BoxLengthMC(4:6), BoxLengthMC(4:6) ) )
        BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( BoxLengthMC(7:9), BoxLengthMC(7:9) ) )
        IF( MaxTranslationalDisplacement >= 2.D0 * MAXVAL( BoxEdgeLength ) ) THEN
          MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxEdgeLength )
        END IF

        ! Rotational adjustment
        IF( nMovementRotationCounter > 0 ) THEN
          ! Acceptance ratio (rotation)
          Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
          ! Rotational adjustment
          IF( Ratio <= AcceptanceRatioRotation ) THEN
            MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
          ELSE
            MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
          END IF
        END IF

        ! Avoid multiple turns (arbitrary)
        IF( ( MaxAngularDisplacement > 4.D0 * cPi ) ) THEN
          MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
        END IF

        ! Reset counter
        nAcceptanceTranslation      = 0
        nMovementTranslationCounter = 0
        nAcceptanceRotation         = 0
        nMovementRotationCounter    = 0

      END IF

      ! Check new configuration
      CYCLE PackingFractionCorrection

    END IF

  END DO PackingFractionCorrection

  ! Close output unit
  CLOSE( 55 )

END IF

! Summary
IF( PackingFractionNPT >= PackingFraction .AND. nAttempts == 1 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0,G0)" ) "Packing fraction fixed after ", nAttempts, " attempt."
  WRITE( *, "(G0)" ) " "
ELSE IF( PackingFractionNPT >= PackingFraction .AND. nAttempts > 1 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0,G0)" ) "Packing fraction fixed after ", nAttempts, " attempts."
  WRITE( *, "(G0)" ) " "
END IF
CALL Sleep( 1 )

! Deallocate list arrays
CALL FinalizeList(  )

RETURN

END SUBROUTINE RandomConfiguration

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!          This subroutine generates a progress bar for the hit-and-miss NVT algorithm.           !
! *********************************************************************************************** !
SUBROUTINE ProgressBarHitAndMiss( iParticle, OverlappingParticles )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle            ! Counter (particle)
INTEGER( Kind= Int64 ) :: OverlappingParticles ! Counter of overlapping particles
INTEGER( Kind= Int64 ) :: AuxiliarInt1         ! Auxiliar
INTEGER( Kind= Int64 ) :: AuxiliarInt2         ! Auxiliar

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 59 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0
AuxiliarInt2 = 0

! Progress bar (FORMAT)
IF( iParticle < 10 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I1)" ) iParticle
ELSE IF( iParticle < 100 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I2)" ) iParticle
ELSE IF( iParticle < 1000 ) THEN
  AuxiliarInt1 = 2
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I3)" ) iParticle
ELSE IF( iParticle < 10000 ) THEN
  AuxiliarInt1 = 3
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I4)" ) iParticle
ELSE IF( iParticle < 100000 ) THEN
  AuxiliarInt1 = 4
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I5)" ) iParticle
ELSE IF( iParticle < 1000000 ) THEN
  AuxiliarInt1 = 5
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I6)" ) iParticle
ELSE IF( iParticle < 10000000 ) THEN
  AuxiliarInt1 = 6
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I7)" ) iParticle
ELSE IF( iParticle < 100000000 ) THEN
  AuxiliarInt1 = 7
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I8)" ) iParticle
ELSE IF( iParticle < 1000000000 ) THEN
  AuxiliarInt1 = 8
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I9)" ) iParticle
ELSE IF( iParticle >= 1000000000 ) THEN
  AuxiliarInt1 = 10
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: > 1 billion "
END IF
IF( OverlappingParticles < 10 ) THEN
  AuxiliarInt2 = 0
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ?"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I1)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 100 ) THEN
  AuxiliarInt2 = 1
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ??"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I2)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 1000 ) THEN
  AuxiliarInt2 = 2
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ???"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I3)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 10000 ) THEN
  AuxiliarInt2 = 3
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I4)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 100000 ) THEN
  AuxiliarInt2 = 4
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ?????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I5)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 1000000 ) THEN
  AuxiliarInt2 = 5
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ??????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I6)" ) OverlappingParticles
ELSE IF( OverlappingParticles >= 1000000 ) THEN
  AuxiliarInt2 = 10
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: > 1 million"
END IF
ProgressBar((39+AuxiliarInt1+AuxiliarInt2):59) = REPEAT( " ", ( (20 - AuxiliarInt1 - AuxiliarInt2) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 38 + AuxiliarInt1 + AuxiliarInt2 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(38+AuxiliarInt1+AuxiliarInt2+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarHitAndMiss

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!       This subroutine generates a progress bar for the algorithm that compress/expand the       !
!              volume of the simulation box to the desirable target packing fraction              !
! *********************************************************************************************** !
SUBROUTINE ProgressBarRandomConfigNPT( iParticle, PackingFractionNPT, TargetPackingFraction )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle    ! Counter
INTEGER( Kind= Int64 ) :: AuxiliarInt1 ! Auxiliar

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: PackingFractionNPT, TargetPackingFraction ! Packing fraction

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 78 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0

! Progress bar (FORMAT)
IF( iParticle < 10 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I1)" ) iParticle
ELSE IF( iParticle < 100 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I2)" ) iParticle
ELSE IF( iParticle < 1000 ) THEN
  AuxiliarInt1 = 2
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I3)" ) iParticle
ELSE IF( iParticle < 10000 ) THEN
  AuxiliarInt1 = 3
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I4)" ) iParticle
ELSE IF( iParticle < 100000 ) THEN
  AuxiliarInt1 = 4
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I5)" ) iParticle
ELSE IF( iParticle < 1000000 ) THEN
  AuxiliarInt1 = 5
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I6)" ) iParticle
ELSE IF( iParticle < 10000000 ) THEN
  AuxiliarInt1 = 6
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I7)" ) iParticle
ELSE IF( iParticle < 100000000 ) THEN
  AuxiliarInt1 = 7
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I8)" ) iParticle
ELSE IF( iParticle < 1000000000 ) THEN
  AuxiliarInt1 = 8
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I9)" ) iParticle
ELSE IF( iParticle >= 1000000000 ) THEN
  AuxiliarInt1 = 10
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: > 1 billion "
END IF
ProgressBar((13+AuxiliarInt1):(45+AuxiliarInt1)) = "| Packing fraction: ???????????? "
WRITE( Unit= ProgressBar((33+AuxiliarInt1):(44+AuxiliarInt1)), Fmt= "(E12.6)" ) PackingFractionNPT
ProgressBar((46+AuxiliarInt1):(67+AuxiliarInt1)) = "(TARGET: ????????????)"
WRITE( Unit= ProgressBar((55+AuxiliarInt1):(66+AuxiliarInt1)), Fmt= "(E12.6)" ) TargetPackingFraction
ProgressBar((68+AuxiliarInt1):78) = REPEAT( " ", ( (10 - AuxiliarInt1) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 67 + AuxiliarInt1 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(67+AuxiliarInt1+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarRandomConfigNPT

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!            This subroutine generates a progress bar for the algorithm that fixes the            !
!                                 packing fraction of the system                                  !
! *********************************************************************************************** !
SUBROUTINE ProgressBarRandomConfigPackingFractionCorrection( iParticle )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle   ! Counter
INTEGER( Kind= Int64 ) :: AuxiliarInt ! Auxiliar

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 22 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Progress bar (FORMAT)
IF( iParticle < 10 ) THEN
  AuxiliarInt = 0
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ?"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I1)" ) iParticle
ELSE IF( iParticle < 100 ) THEN
  AuxiliarInt = 1
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ??"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I2)" ) iParticle
ELSE IF( iParticle < 1000 ) THEN
  AuxiliarInt = 2
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ???"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I3)" ) iParticle
ELSE IF( iParticle < 10000 ) THEN
  AuxiliarInt = 3
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I4)" ) iParticle
ELSE IF( iParticle < 100000 ) THEN
  AuxiliarInt = 4
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ?????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I5)" ) iParticle
ELSE IF( iParticle < 1000000 ) THEN
  AuxiliarInt = 5
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ??????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I6)" ) iParticle
ELSE IF( iParticle < 10000000 ) THEN
  AuxiliarInt = 6
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ???????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I7)" ) iParticle
ELSE IF( iParticle < 100000000 ) THEN
  AuxiliarInt = 7
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ????????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I8)" ) iParticle
ELSE IF( iParticle < 1000000000 ) THEN
  AuxiliarInt = 8
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: ?????????"
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt)), Fmt= "(I9)" ) iParticle
ELSE IF( iParticle >= 1000000000 ) THEN
  AuxiliarInt = 10
  ProgressBar(1:(11+AuxiliarInt)) = "Attempts: > 1 billion"
END IF
ProgressBar((12+AuxiliarInt):22) = REPEAT( " ", ( (10 - AuxiliarInt) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 11 + AuxiliarInt + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(11+AuxiliarInt+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarRandomConfigPackingFractionCorrection

! *********************************************************************************************** !
!      This subroutine creates a file containing all properties of the initial configuration      !
! *********************************************************************************************** !
SUBROUTINE ConfigurationOutput(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle  ! Counter (particles)
INTEGER( Kind= Int64 ) :: cComponent ! Counter (component)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 64 )  :: FormatPosition   ! String format (general)
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Format (position)
FormatPosition = "(G0,1X,7(G0.15,1X))"

! Simple cubic structure
IF( ConfigurationSelection(1) ) THEN

  ! Initial configuration (OVITO) | Positions, orientations, and geometric details
  IF( EnsembleMC == "NVT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_initconf_η"// &
    &                      TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_sc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  ELSE IF( EnsembleMC == "NPT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_initconf_P"// &
    &                      TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_sc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) nParticles
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (SC) | All Properties
  IF( .NOT. PresetInitialConfiguration ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_initconf_sc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GeometryInquiry
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", ConfigurationInquiry
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PackingFraction
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", nComponents
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,3G0)" ) "Spherical_Component_#", cComponent, ": ", SphericalComponentInquiry(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", cComponent, ": ", cDiameter(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", cComponent, ": ", cLength(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Aspect_Ratio_of_Component_#", cComponent, ": ", cAspectRatio(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", cComponent, ": ", cMolarFraction(cComponent)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", nParticles
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TotalParticleVolume, " Å³"
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", cComponent, ": ", cParticles(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", cComponent, ": ", cMolecularVolume(cComponent), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BoxVolume, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BoxLength(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BoxLength(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BoxLength(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BoxLengthInverse(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BoxLengthInverse(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BoxLengthInverse(7:9)
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", cComponent, ": ", cNumberDensity(cComponent), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TotalNumberDensity, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", AbsoluteTemperature, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", ReducedPressure
    WRITE( 10, * ) " "
    DO cComponent = 1, nComponents
      DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 10, FormatPosition ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
        &                           pQuaternion(0,pParticle), pQuaternion(1,pParticle), pQuaternion(2,pParticle), &
        &                           pQuaternion(3,pParticle)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! Body-centered cubic structure
ELSE IF( ConfigurationSelection(2) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( EnsembleMC == "NVT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_bcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  ELSE IF( EnsembleMC == "NPT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_bcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) nParticles
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (BCC) | All Properties
  IF( .NOT. PresetInitialConfiguration ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_initconf_bcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GeometryInquiry
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", ConfigurationInquiry
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PackingFraction
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", nComponents
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,3G0)" ) "Spherical_Component_#", cComponent, ": ", SphericalComponentInquiry(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", cComponent, ": ", cDiameter(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", cComponent, ": ", cLength(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Aspect_Ratio_of_Component_#", cComponent, ": ", cAspectRatio(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", cComponent, ": ", cMolarFraction(cComponent)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", nParticles
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TotalParticleVolume, " Å³"
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", cComponent, ": ", cParticles(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", cComponent, ": ", cMolecularVolume(cComponent), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BoxVolume, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BoxLength(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BoxLength(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BoxLength(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BoxLengthInverse(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BoxLengthInverse(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BoxLengthInverse(7:9)
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", cComponent, ": ", cNumberDensity(cComponent), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TotalNumberDensity, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", AbsoluteTemperature, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", ReducedPressure
    WRITE( 10, * ) " "
    DO cComponent = 1, nComponents
      DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 10, FormatPosition ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
        &                           pQuaternion(0,pParticle), pQuaternion(1,pParticle), pQuaternion(2,pParticle), &
        &                           pQuaternion(3,pParticle)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! Face-centered cubic structure
ELSE IF( ConfigurationSelection(3) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( EnsembleMC == "NVT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_fcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  ELSE IF( EnsembleMC == "NPT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_fcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) nParticles
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (FCC) | All Properties
  IF( .NOT. PresetInitialConfiguration ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_initconf_fcc_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GeometryInquiry
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", ConfigurationInquiry
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PackingFraction
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", nComponents
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,3G0)" ) "Spherical_Component_#", cComponent, ": ", SphericalComponentInquiry(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", cComponent, ": ", cDiameter(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", cComponent, ": ", cLength(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Aspect_Ratio_of_Component_#", cComponent, ": ", cAspectRatio(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", cComponent, ": ", cMolarFraction(cComponent)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", nParticles
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TotalParticleVolume, " Å³"
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", cComponent, ": ", cParticles(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", cComponent, ": ", cMolecularVolume(cComponent), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BoxVolume, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BoxLength(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BoxLength(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BoxLength(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BoxLengthInverse(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BoxLengthInverse(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BoxLengthInverse(7:9)
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", cComponent, ": ", cNumberDensity(cComponent), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TotalNumberDensity, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", AbsoluteTemperature, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", ReducedPressure
    WRITE( 10, * ) " "
    DO cComponent = 1, nComponents
      DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 10, FormatPosition ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
        &                           pQuaternion(0,pParticle), pQuaternion(1,pParticle), pQuaternion(2,pParticle), &
        &                           pQuaternion(3,pParticle)
      END DO
    END DO
    CLOSE( 10 )
  END IF

! Random cubic structure
ELSE IF( ConfigurationSelection(4) ) THEN

  ! Initial configuration (OVITO) | Positions, Orientations, and Geometric Details
  IF( EnsembleMC == "NVT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_rnd_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  ELSE IF( EnsembleMC == "NPT" ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                      "_initconf_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_rnd_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
  END IF
  WRITE( 10, "(G0)" ) nParticles
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  ELSE IF( GeometryType(3) ) THEN ! Cylinders
    DO cComponent = 1, nComponents
      IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          cLength(cComponent)
        END DO
      ELSE
        DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
          WRITE( 10, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
          &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
          &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent), &
          &                          0.5D0 * cDiameter(cComponent)
        END DO
      END IF
    END DO
  END IF
  CLOSE( 10 )

  ! Initial configuration (RND) | All Properties
  IF( .NOT. PresetInitialConfiguration ) THEN
    OPEN( Unit= 10, File= "Initial_Configuration/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_initconf_rnd_"// &
    &                      TRIM( DescriptorFileGeometry )//".xyz" )
    WRITE( 10, "(G0,G0)" ) "Molecular_geometry: ", GeometryInquiry
    WRITE( 10, "(G0,G0)" ) "Initial_configuration: ", ConfigurationInquiry
    WRITE( 10, "(G0,G0.5)" ) "Packing_fraction: ", PackingFraction
    WRITE( 10, "(G0,G0)" ) "Number_of_Components: ", nComponents
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,3G0)" ) "Spherical_Component_#", cComponent, ": ", SphericalComponentInquiry(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Diameter_of_Component_#", cComponent, ": ", cDiameter(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Length_of_Component_#", cComponent, ": ", cLength(cComponent), " Å"
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Aspect_Ratio_of_Component_#", cComponent, ": ", cAspectRatio(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5)" ) "Molar_Fraction_of_Component_#", cComponent, ": ", cMolarFraction(cComponent)
    END DO
    WRITE( 10, "(G0,G0)" ) "Number_of_Particles: ", nParticles
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Molecular_Volume: ", TotalParticleVolume, " Å³"
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0)" ) "Number_of_Particles_of_Component_#", cComponent, ": ", cParticles(cComponent)
    END DO
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Molecular_Volume_of_Component_#", cComponent, ": ", cMolecularVolume(cComponent), " Å³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Volume_of_the_Simulation_Box_(Cubic): ", BoxVolume, " Å³"
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_1)_in_Å: ", BoxLength(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_2)_in_Å: ", BoxLength(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Cubic_3)_in_Å: ", BoxLength(7:9)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_1)_in_Å⁻¹: ", BoxLengthInverse(1:3)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_2)_in_Å⁻¹: ", BoxLengthInverse(4:6)
    WRITE( 10, "(G0,3(G0.15,1X))" ) "Length_of_the_Simulation_Box_(Inverse_3)_in_Å⁻¹: ", BoxLengthInverse(7:9)
    DO cComponent = 1, nComponents
      WRITE( 10, "(G0,G0,G0,G0.5,G0)" ) "Number_Density_of_Component_#", cComponent, ": ", cNumberDensity(cComponent), " Å⁻³"
    END DO
    WRITE( 10, "(G0,G0.5,G0)" ) "Total_Number_Density: ", TotalNumberDensity, " Å⁻³"
    WRITE( 10, "(G0,G0.5,G0)" ) "Absolute_Temperature: ", AbsoluteTemperature, " K"
    WRITE( 10, "(G0,G0.5)" ) "Reduced_Pressure: ", ReducedPressure
    WRITE( 10, * ) " "
    DO cComponent = 1, nComponents
      DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
        WRITE( 10, FormatPosition ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), pPosition(3,pParticle), &
        &                           pQuaternion(0,pParticle), pQuaternion(1,pParticle), pQuaternion(2,pParticle), &
        &                           pQuaternion(3,pParticle)
      END DO
    END DO
    CLOSE( 10 )
  END IF

END IF

RETURN

END SUBROUTINE ConfigurationOutput

END MODULE InitialSystemConfiguration
