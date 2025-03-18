! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!                    This code contains a subroutine used in the main program                     !
!                      that opens all file units and creates their headers.                       !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       February 15th, 2024                                       !
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
SUBROUTINE FileHandler(  )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle ! Counter (particles)
INTEGER( Kind= Int64 ) :: cCylinder ! Counter (cylinders)
INTEGER( Kind= Int64 ) :: pImage    ! Counter (images)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: TrajectoryPosition ! Position array (cylinders)

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Trajectory file (depends on user's choice)
IF( TrajectoryLogical ) THEN
  OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_traj_N" &
  &                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".xyz", Action= "WRITE" )
  WRITE( 20, "(G0)" ) nParticles * 4 + nParticles * nImages * 4
  ! Descriptor string
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 20, DescriptorString ) 'Lattice="', (1.D0 + 2.D0 * DBLE( nLayers) ) * BoxLength(1:9), '" Origin="', -0.5D0 * &
  &                             (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLength(1) + BoxLength(4) + BoxLength(7) ), &
  &                             -0.5D0 * (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLength(2) + BoxLength(5) + &
  &                             BoxLength(8) ), -0.5D0 * (1.D0 + 2.D0 * DBLE( nLayers) ) * ( BoxLength(3) + &
  &                             BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO pParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of cylinders
      TrajectoryPosition(1) = cPosition(1,cCylinder,pParticle)
      TrajectoryPosition(2) = cPosition(2,cCylinder,pParticle)
      TrajectoryPosition(3) = cPosition(3,cCylinder,pParticle)
      WRITE( 20, "(11(G0,1X))" ) Atom(1), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), &
      &                          pQuaternion(1,pParticle), pQuaternion(2,pParticle), pQuaternion(3,pParticle), &
      &                          pQuaternion(0,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      DO pImage = 1, nImages
        TrajectoryPosition(1) = imcPosition(1,cCylinder,pImage,pParticle)
        TrajectoryPosition(2) = imcPosition(2,cCylinder,pImage,pParticle)
        TrajectoryPosition(3) = imcPosition(3,cCylinder,pImage,pParticle)
        WRITE( 20, "(11(G0,1X))" ) Atom(2), TrajectoryPosition(1), TrajectoryPosition(2), TrajectoryPosition(3), &
        &              imQuaternion(1,pImage,pParticle), imQuaternion(2,pImage,pParticle), imQuaternion(3,pImage,pParticle), &
        &              imQuaternion(0,pImage,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      END DO
    END DO
  END DO
  FLUSH( 20 )
END IF

! Ratio file (translation)
OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_N" &
&                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )
WRITE( 30, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å]"//'"', ",", &
&                    '"'//"Acceptance Ratio Threshold"//'"'
FLUSH( 30 )

! Ratio file (rotation)
OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_N" &
&                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )
WRITE( 40, "(7G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [rad]"//'"', ",", &
&                    '"'//"Acceptance Ratio Threshold"//'"'
FLUSH( 40 )

! Ratio file (volume)
OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_N" &
&                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )
WRITE( 50, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Acceptance Ratio"//'"', ",", '"'//"Maximum Displacement [Å³]"//'"', ",", &
&                    '"'//"Acceptance Ratio Threshold"//'"', ",", '"'//"Type of Volume Scaling"//'"'
FLUSH( 50 )

! Results file
OPEN( Unit= 60, File= "Results/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_results_N" &
&                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )
WRITE( 60, "(9G0)" ) '"'//"Cycles"//'"', ",", '"'//"Packing Fraction"//'"', ",", '"'//"Number Density [Å⁻³]"//'"', ",", &
&                    '"'//"Box Volume"//'"', ",", '"'//"Reduced Pressure"//'"'
FLUSH( 60 )

! Floppy box file
OPEN( Unit= 70, File= "Box/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_box_N" &
&                     //TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )

! Order parameter file
OPEN( Unit= 80, File= "Order_Parameter/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )// &
&                     "_order_N"//TRIM( DescriptorFileParticles )//"_AR"//TRIM( DescriptorFileGeometry )//".dat", Action= "WRITE" )
WRITE( 80, "(3G0)" ) '"'//"Cycles"//'"', ",", '"'//"Nematic Order Parameter"//'"'
FLUSH( 80 )

RETURN

END SUBROUTINE FileHandler