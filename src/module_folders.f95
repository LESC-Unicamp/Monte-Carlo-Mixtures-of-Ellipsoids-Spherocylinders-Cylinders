! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!           This module creates folders and subfolders to organize simulation results.            !
! Directories are created by executing a shell command via an intrinsic function called 'system'. !
!               Please note that which shell is used to invoke the command line is                !
!                           system-dependent and environment-dependent.                           !
!         See <https://gcc.gnu.org/onlinedocs/gfortran/SYSTEM.html> for more information.         !
!                  The code below is meant for Ubuntu 20.04 operational systems.                  !
!                   We have not provided an alternative code for Windows users.                   !
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
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

MODULE FOLDERS

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL, DIMENSION( 5 ) :: FEXIST  ! Checks whether folder exists or not
LOGICAL, DIMENSION( 7 ) :: DFEXIST ! Checks whether date folders exist or not
LOGICAL, DIMENSION( 5 ) :: SFEXIST ! Checks whether subfolder exists or not

CONTAINS

! *********************************************************************************************** !
!                             Initialization of parental directories                              !
! *********************************************************************************************** !
SUBROUTINE INITFOLDER(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( FILE= "Initial_Configuration", EXIST= FEXIST(1) )

! *********************************************************************************************** !
! Initial configuration folder (holds information on the initial molecular structure)             !
! *********************************************************************************************** !
IF( .NOT. FEXIST(1) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating initial configuration folder..."
  CALL SYSTEM( "mkdir Initial_Configuration" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! *********************************************************************************************** !
! Inquires whether a subfolder exists and stores the inquiry result in a logical variable         !
! *********************************************************************************************** !
!  The initial molecular structure at 'OVITO' subfolder is properly formatted to be analyzed      !
!  by that software.                                                                              !
! *********************************************************************************************** !
INQUIRE( FILE= "Initial_Configuration/OVITO/", EXIST= SFEXIST(1) )

! *********************************************************************************************** !
! Initial configuration subfolder (OVITO)                                                         !
! *********************************************************************************************** !
IF( .NOT. SFEXIST(1) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating initial configuration subfolder..."
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( FILE= "Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/", EXIST= DFEXIST(1) )

! *********************************************************************************************** !
! Date subfolder                                                                                  !
! *********************************************************************************************** !
IF( .NOT. DFEXIST(1) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating date subfolder..."
  CALL SYSTEM( "mkdir Initial_Configuration/OVITO/"//TRIM( DESCRIPTOR_DATE )//"/" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE INITFOLDER

! *********************************************************************************************** !
!                             Initialization of parental directories                              !
! *********************************************************************************************** !
SUBROUTINE INQUIRE_FOLDERS(  )

IMPLICIT NONE

! Inquires whether a folder exists and stores the inquiry result in a logical variable
INQUIRE( FILE= "Trajectories", EXIST= FEXIST(2) )
INQUIRE( FILE= "Ratio", EXIST= FEXIST(3) )
INQUIRE( FILE= "Order_Parameter", EXIST= FEXIST(4) )
INQUIRE( FILE= "Results", EXIST= FEXIST(5) )

! *********************************************************************************************** !
! Trajectory folder (holds information on orientation and position of particles)                  !
! *********************************************************************************************** !
IF( .NOT. FEXIST(2) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating trajectory folder..."
  CALL SYSTEM( "mkdir Trajectories" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! *********************************************************************************************** !
! Ratio folder (holds information on the equilibration cycles)                                    !
! *********************************************************************************************** !
IF( .NOT. FEXIST(3) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating ratio folder..."
  CALL SYSTEM( "mkdir Ratio" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! Inquires whether a subfolder exists and stores the inquiry result in a logical variable
INQUIRE( FILE= "Ratio/Translation/", EXIST= SFEXIST(2) )
INQUIRE( FILE= "Ratio/Rotation/", EXIST= SFEXIST(3) )
INQUIRE( FILE= "Ratio/Volume/", EXIST= SFEXIST(4) )
INQUIRE( FILE= "Ratio/Box/", EXIST= SFEXIST(5) )

! *********************************************************************************************** !
! Ratio subfolders                                                                                !
! *********************************************************************************************** !
IF( .NOT. SFEXIST(2) .OR. .NOT. SFEXIST(3) .OR. .NOT. SFEXIST(4) .OR. .NOT. SFEXIST(5) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating ratio subfolders..."
END IF
IF( .NOT. SFEXIST(2) ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/" )
END IF
IF( .NOT. SFEXIST(3) ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/" )
END IF
IF( .NOT. SFEXIST(4) ) THEN
  CALL SYSTEM( "mkdir Ratio/Volume/" )
END IF
IF( .NOT. SFEXIST(5) ) THEN
  CALL SYSTEM( "mkdir Ratio/Box/" )
END IF
IF( .NOT. SFEXIST(2) .OR. .NOT. SFEXIST(3) .OR. .NOT. SFEXIST(4) .OR. .NOT. SFEXIST(5) ) THEN
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! *********************************************************************************************** !
! Order parameter folder (holds information on the nematic order parameter)                       !
! *********************************************************************************************** !
IF( .NOT. FEXIST(4) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating order parameter folder..."
  CALL SYSTEM( "mkdir Order_Parameter" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

! *********************************************************************************************** !
! Results folder (holds information on the packing fraction and the simulation box)               !
! *********************************************************************************************** !
IF( .NOT. FEXIST(5) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating results folder..."
  CALL SYSTEM( "mkdir Results" )
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE INQUIRE_FOLDERS

! *********************************************************************************************** !
!                                Initialization of date subfolders                                !
! *********************************************************************************************** !
SUBROUTINE DATE_FOLDERS(  )

IMPLICIT NONE

! Inquires whether the date subfolder exists and stores the inquiry result in a logical variable
INQUIRE( FILE= "Trajectories/"//TRIM( DESCRIPTOR_DATE )//"/", EXIST= DFEXIST(1) )
INQUIRE( FILE= "Ratio/Translation/"//TRIM( DESCRIPTOR_DATE )//"/", EXIST= DFEXIST(2) )
INQUIRE( FILE= "Ratio/Rotation/"//TRIM( DESCRIPTOR_DATE )//"/", EXIST= DFEXIST(3) )
INQUIRE( FILE= "Ratio/Volume/"//TRIM(DESCRIPTOR_DATE)//"/", EXIST= DFEXIST(4) )
INQUIRE( FILE= "Ratio/Box/"//TRIM(DESCRIPTOR_DATE)//"/", EXIST= DFEXIST(5) )
INQUIRE( FILE= "Order_Parameter/"//TRIM( DESCRIPTOR_DATE )//"/", EXIST= DFEXIST(6) )
INQUIRE( FILE= "Results/"//TRIM(DESCRIPTOR_DATE)//"/", EXIST= DFEXIST(7) )

! *********************************************************************************************** !
! Date subfolders                                                                                 !
! *********************************************************************************************** !
IF( .NOT. DFEXIST(1) .OR. .NOT. DFEXIST(2) .OR. .NOT. DFEXIST(3) .OR. .NOT. DFEXIST(4) .OR. .NOT. DFEXIST(5) .OR. &
&   .NOT. DFEXIST(6) .OR. .NOT. DFEXIST(7) ) THEN
  WRITE( *, "(G0)", ADVANCE= "NO" ) "Creating date subfolders..."
END IF
IF( .NOT. DFEXIST(1) ) THEN
  CALL SYSTEM( "mkdir Trajectories/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(2) ) THEN
  CALL SYSTEM( "mkdir Ratio/Translation/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(3) ) THEN
  CALL SYSTEM( "mkdir Ratio/Rotation/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(4) ) THEN
  CALL SYSTEM( "mkdir Ratio/Volume/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(5) ) THEN
  CALL SYSTEM( "mkdir Ratio/Box/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(6) ) THEN
  CALL SYSTEM( "mkdir Order_Parameter/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(7) ) THEN
  CALL SYSTEM( "mkdir Results/"//TRIM( DESCRIPTOR_DATE )//"/" )
END IF
IF( .NOT. DFEXIST(1) .OR. .NOT. DFEXIST(2) .OR. .NOT. DFEXIST(3) .OR. .NOT. DFEXIST(4) .OR. .NOT. DFEXIST(5) .OR. &
&   .NOT. DFEXIST(6) .OR. .NOT. DFEXIST(7) ) THEN
  CALL SLEEP( 1 )
  WRITE( *, "(G0)", ADVANCE= "YES" ) " Done!"
  WRITE( *, "(G0)" ) " "
END IF

RETURN

END SUBROUTINE DATE_FOLDERS

END MODULE FOLDERS