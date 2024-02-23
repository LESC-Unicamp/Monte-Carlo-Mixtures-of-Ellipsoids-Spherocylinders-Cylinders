! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                    This code contains a subroutine used in the main program                     !
!                      that opens all file units and creates their headers.                       !
!                                                                                                 !
! Version number: 1.3.1                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       February 9th, 2024                                        !
! ############################################################################################### !
! Main References:                          A. Stukowski                                          !
!                      Modelling Simul. Mater. Sci. Eng. 18, 015012O (2010)                       !
!                               DOI: 10.1088/0965-0393/18/1/015012                                !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!             This subroutine opens all file units to be written in the main program              !
! *********************************************************************************************** !
SUBROUTINE FileHandler( LastLine )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: FileStatus ! Status of file

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: pParticle    ! Counter (particles)
INTEGER( Kind= Int64 )                 :: rRange       ! Counter (potential range)
INTEGER( Kind= Int64 )                 :: cComponent   ! Counter (component)
INTEGER( Kind= Int64 )                 :: LastLineTemp ! Last line of a file (temporary)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: LastLine     ! Last line of a file

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FileExist ! Checks whether a file exists or not

! Initialize last line of equilibration files
LastLine  = 0
FileExist = .FALSE.

! *********************************************************************************************** !
! OVITO Visualization Tool                                                                        !
! *********************************************************************************************** !
! The trajectory files are pre-formatted to be visualized in OVITO.                               !
! Please visit <https://www.ovito.org/> for more details.                                         !
! *********************************************************************************************** !

! Trajectory file (depends on user's choice)
IF( TrajectoryLogical ) THEN
  IF( EnsembleMC == "NVT" ) THEN
    INQUIRE( File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_traj_η"// &
    &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
    &              TRIM( DescriptorFileGeometry )//".xyz", Exist= FileExist )
    IF( .NOT. FileExist ) THEN
      OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
      &                     "_traj_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
      &                     //TRIM( DescriptorFileGeometry )//".xyz" )
    ELSE
      OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
      &                     "_traj_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
      &                     //TRIM( DescriptorFileGeometry )//".xyz", Position= "APPEND" )
    END IF
  ELSE IF( EnsembleMC == "NPT" ) THEN
    INQUIRE( File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_traj_P"// &
    &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
    &              TRIM( DescriptorFileGeometry )//".xyz", Exist= FileExist )
    IF( .NOT. FileExist ) THEN
      OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
      &                     "_traj_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
      &                     //TRIM( DescriptorFileGeometry )//".xyz" )
    ELSE
      OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
      &                     "_traj_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
      &                     //TRIM( DescriptorFileGeometry )//".xyz", Position= "APPEND" )
    END IF
  END IF
  IF( .NOT. FileExist ) THEN
    WRITE( 20, "(G0)" ) nParticles
    WRITE( 20, * ) " "
    IF( GeometryType(1) ) THEN ! Ellipsoids-of-revolution
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cLength(cComponent)
          END DO
        ELSE
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    ELSE IF( GeometryType(2) ) THEN ! Spherocylinders
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), cLength(cComponent)
          END DO
        ELSE
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    ELSE IF( GeometryType(3) ) THEN ! Cylinders
      DO cComponent = 1, nComponents
        IF( .NOT. SphericalComponentLogical(cComponent) ) THEN
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), cLength(cComponent)
          END DO
        ELSE
          DO pParticle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
            WRITE( 20, "(11(G0,1X))" ) cIndex(cComponent), pPosition(1,pParticle), pPosition(2,pParticle), &
            &                          pPosition(3,pParticle), pQuaternion(0,pParticle), pQuaternion(1,pParticle), &
            &                          pQuaternion(2,pParticle), pQuaternion(3,pParticle), 0.5D0 * cDiameter(cComponent), &
            &                          0.5D0 * cDiameter(cComponent), 0.5D0 * cDiameter(cComponent)
          END DO
        END IF
      END DO
    END IF
    FLUSH( 20 )
  END IF
END IF

! Ratio file (translation)
IF( EnsembleMC == "NVT" ) THEN
  INQUIRE( File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_η"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 30, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å]"//'"', ",", &
    &                    '"'//"Acceptance Ratio Threshold"//'"'
    FLUSH( 30 )
  ELSE
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
    ! Skip header
    READ( 30, * )
    DO
      READ( 30, *, IOStat= FileStatus ) LastLineTemp
      IF( FileStatus /= 0 ) EXIT
      LastLine(1) = LastLineTemp
    END DO
    CLOSE( 30 )
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
ELSE IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 30, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å]"//'"', ",", &
    &                    '"'//"Acceptance Ratio Threshold"//'"'
    FLUSH( 30 )
  ELSE
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
    ! Skip header
    READ( 30, * )
    DO
      READ( 30, *, IOStat= FileStatus ) LastLineTemp
      IF( FileStatus /= 0 ) EXIT
      LastLine(1) = LastLineTemp
    END DO
    CLOSE( 30 )
    OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Ratio file (rotation)
IF( EnsembleMC == "NVT" ) THEN
  INQUIRE( File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_η"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 40, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [rad]"//'"', ",", &
    &                    '"'//"Acceptance Ratio Threshold"//'"'
    FLUSH( 40 )
  ELSE
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
    ! Skip header
    READ( 40, * )
    DO
      READ( 40, *, IOStat= FileStatus ) LastLineTemp
      IF( FileStatus /= 0 ) EXIT
      LastLine(2) = LastLineTemp
    END DO
    CLOSE( 40 )
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
ELSE IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 40, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [rad]"//'"', ",", &
    &                    '"'//"Acceptance Ratio Threshold"//'"'
    FLUSH( 40 )
  ELSE
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
    ! Skip header
    READ( 40, * )
    DO
      READ( 40, *, IOStat= FileStatus ) LastLineTemp
      IF( FileStatus /= 0 ) EXIT
      LastLine(2) = LastLineTemp
    END DO
    CLOSE( 40 )
    OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Ratio file (volume scaling)
IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 50, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å³]"//'"', ",", &
    &                    '"'//"Acceptance Ratio Threshold"//'"', ",", '"'//"Type of Volume Scaling"//'"'
    FLUSH( 50 )
  ELSE
    OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
    ! Skip header
    READ( 50, * )
    DO
      READ( 50, *, IOStat= FileStatus ) LastLineTemp
      IF( FileStatus /= 0 ) EXIT
      LastLine(3) = LastLineTemp
    END DO
    CLOSE( 50 )
    OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_ratio_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Box file
IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Box/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"//TRIM( DescriptorFileThermoVariable )// &
  &              "_C"//TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 55, File= "Box/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"// &
    &                     TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
    &                     TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 55, "(23G0)" ) '"'//"Cycles"//'"', ",", '"'//"Box Distortion"//'"', ",", &
    &                     '"'//"Box Volume [Å³]"//'"', ",", '"'//"Box Length 1 [Å]"//'"', ",", &
    &                     '"'//"Box Length 2 [Å]"//'"', ",", '"'//"Box Length 3 [Å]"//'"', ",", &
    &                     '"'//"Box Length 4 [Å]"//'"', ",", '"'//"Box Length 5 [Å]"//'"', ",", &
    &                     '"'//"Box Length 6 [Å]"//'"', ",", '"'//"Box Length 7 [Å]"//'"', ",", &
    &                     '"'//"Box Length 8 [Å]"//'"', ",", '"'//"Box Length 9 [Å]"//'"'
    FLUSH( 55 )
  ELSE
    OPEN( Unit= 55, File= "Box/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_P"// &
    &                     TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
    &                     TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Order parameter file
IF( EnsembleMC == "NVT" ) THEN
  INQUIRE( File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_order_η"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 60, File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_order_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 60, "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Nematic Order Parameter"//'"'
    FLUSH( 60 )
  ELSE
    OPEN( Unit= 60, File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_order_η"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
ELSE IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_order_P"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 60, File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_order_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 60, "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Nematic Order Parameter"//'"'
    FLUSH( 60 )
  ELSE
    OPEN( Unit= 60, File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_order_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Results file
IF( EnsembleMC == "NPT" ) THEN
  INQUIRE( File= "Results/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_results_P"// &
  &              TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_"// &
  &              TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    OPEN( Unit= 70, File= "Results/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_results_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat" )
    WRITE( 70, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Packing Fraction"//'"', ",", '"'//"Number Density [Å⁻³]"//'"', ",", &
    &                    '"'//"Box Volume"//'"', ",", '"'//"Reduced Pressure"//'"'
    FLUSH( 70 )
  ELSE
    OPEN( Unit= 70, File= "Results/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
    &                     "_results_P"//TRIM( DescriptorFileThermoVariable )//"_C"//TRIM( DescriptorFileComponents )//"_" &
    &                     //TRIM( DescriptorFileGeometry )//".dat", Position= "APPEND" )
  END IF
END IF

! Potential file
IF( PerturbedPotentialTypeLogical(2) .OR. FullPotentialTypeLogical(2) ) THEN
  IF( EnsembleMC == "NVT" ) THEN
    DO rRange = 1, nRange
      WRITE( DescriptorRange, FormatRange ) PotentialRange(rRange)
      ! Potential files
      INQUIRE( File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
      &              TRIM( DescriptorHour )//"_thermo_η"//TRIM( DescriptorFileThermoVariable )//"_C"// &
      &              TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
      IF( .NOT. FileExist ) THEN
        OPEN( Unit= (80 + rRange), File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
        &                                 TRIM( DescriptorHour )//"_thermo_η"//TRIM( DescriptorFileThermoVariable )//"_C"// &
        &                                 TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat" )
        WRITE( (80 + rRange), "(5G0)" ) '"'//"Cycles"//'"', ",", '"'//"Potential Energy"//'"', ",", '"'//"Temperature"//'"'
        FLUSH( (80 + rRange) )
      ELSE
        OPEN( Unit= (80 + rRange), File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
        &                                 TRIM( DescriptorHour )//"_thermo_η"//TRIM( DescriptorFileThermoVariable )//"_C"// &
        &                                 TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", &
        &                                 Position= "APPEND" )
      END IF
    END DO
  ELSE IF( EnsembleMC == "NPT" ) THEN
    DO rRange = 1, nRange
      WRITE( DescriptorRange, FormatRange ) PotentialRange(rRange)
      ! Potential files
      INQUIRE( File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
      &              TRIM( DescriptorHour )//"_thermo_P"//TRIM( DescriptorFileThermoVariable )//"_C"// &
      &              TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", Exist= FileExist )
      IF( .NOT. FileExist ) THEN
        OPEN( Unit= (80 + rRange), File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
        &                                 TRIM( DescriptorHour )//"_thermo_P"//TRIM( DescriptorFileThermoVariable )//"_C"// &
        &                                 TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat" )
        WRITE( (80 + rRange), "(5G0)" ) '"'//"Cycles"//'"', ",", '"'//"Potential Energy"//'"', ",", '"'//"Temperature"//'"'
        FLUSH( (80 + rRange) )
      ELSE
        OPEN( Unit= (80 + rRange), File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
        &                                 TRIM( DescriptorHour )//"_thermo_P"//TRIM( DescriptorFileThermoVariable )//"_C"// &
        &                                 TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", &
        &                                 Position= "APPEND" )
      END IF
    END DO
  END IF
END IF

RETURN

END SUBROUTINE FileHandler
