! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains the subroutines used in the main program                    !
!                         accounting for quaternion algebraic operations.                         !
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
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!     This subroutine takes a body-fixed orientation/position and a rotation quaternion and       !
!            generates a space-fixed orientation/position using a 3D-rotation matrix.             !
!        See Allen and Tildesley, 2nd Edition (2017), pages 106-111 for more information.         !
! *********************************************************************************************** !
SUBROUTINE ActiveTransformation( BodyFixedFrame, RotationQuaternion, SpaceFixedFrame )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iRow, jCol ! Counters

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BodyFixedFrame     ! Body-fixed orientation/position
REAL( Kind= Real64 ), DIMENSION( 3 )    :: SpaceFixedFrame    ! Space-fixed orientation/position
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: RotationQuaternion ! Rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 3, 3 ) :: RotationMatrix     ! Rotation matrix
REAL( Kind= Real64 ), DIMENSION( 3, 3 ) :: tRotationMatrix    ! Transpose of rotation matrix

! *********************************************************************************************** !
! Rotation Matrix - Allen and Tildesley, 2nd edition, page 110                                    !
! *********************************************************************************************** !
! First row
RotationMatrix(1,1) = ( RotationQuaternion(0) * RotationQuaternion(0) ) + ( RotationQuaternion(1) * RotationQuaternion(1) ) - &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) - ( RotationQuaternion(3) * RotationQuaternion(3) )
RotationMatrix(1,2) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(2) + RotationQuaternion(0) * RotationQuaternion(3) )
RotationMatrix(1,3) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(3) - RotationQuaternion(0) * RotationQuaternion(2) )
! Second row
RotationMatrix(2,1) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(2) - RotationQuaternion(0) * RotationQuaternion(3) )
RotationMatrix(2,2) = ( RotationQuaternion(0) * RotationQuaternion(0) ) - ( RotationQuaternion(1) * RotationQuaternion(1) ) + &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) - ( RotationQuaternion(3) * RotationQuaternion(3) )
RotationMatrix(2,3) = 2.D0 * ( RotationQuaternion(2) * RotationQuaternion(3) + RotationQuaternion(0) * RotationQuaternion(1) )
! Third row
RotationMatrix(3,1) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(3) + RotationQuaternion(0) * RotationQuaternion(2) )
RotationMatrix(3,2) = 2.D0 * ( RotationQuaternion(2) * RotationQuaternion(3) - RotationQuaternion(0) * RotationQuaternion(1) )
RotationMatrix(3,3) = ( RotationQuaternion(0) * RotationQuaternion(0) ) - ( RotationQuaternion(1) * RotationQuaternion(1) ) - &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) + ( RotationQuaternion(3) * RotationQuaternion(3) )

! Transpose of rotation matrix
DO iRow = 1, 3
  DO jCol = 1, 3
    tRotationMatrix(iRow,jCol) = RotationMatrix(jCol,iRow)
  END DO
END DO

! Active tranformation (dot product of body-fixed vector and transpose of rotation matrix)
SpaceFixedFrame(1) = ( BodyFixedFrame(1) * tRotationMatrix(1,1) ) + ( BodyFixedFrame(2) * tRotationMatrix(1,2) ) + &
&                    ( BodyFixedFrame(3) * tRotationMatrix(1,3) )
SpaceFixedFrame(2) = ( BodyFixedFrame(1) * tRotationMatrix(2,1) ) + ( BodyFixedFrame(2) * tRotationMatrix(2,2) ) + &
&                    ( BodyFixedFrame(3) * tRotationMatrix(2,3) )
SpaceFixedFrame(3) = ( BodyFixedFrame(1) * tRotationMatrix(3,1) ) + ( BodyFixedFrame(2) * tRotationMatrix(3,2) ) + &
&                    ( BodyFixedFrame(3) * tRotationMatrix(3,3) )

RETURN

END SUBROUTINE ActiveTransformation

! *********************************************************************************************** !
!              This subroutine rotates a vector using Hamilton's quaternion product.              !
! *********************************************************************************************** !
SUBROUTINE VectorRotation( PointVector, RotationQuaternion, RotatedVector )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 )   :: PointVector         ! Point vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 3 )   :: RotatedVector       ! Rotated vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RotationQuaternion  ! Rotation quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ConjugateQuaternion ! Conjugate quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RotatedQuaternion   ! Rotated quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: TempQuaternion      ! Temporary quaternion [WXYZ]

! Multiplication of quaternions to find the rotated quaternion
CALL QuaternionMultiplication( RotationQuaternion, [ 0.D0, PointVector ], TempQuaternion )

! Conjugate quaternion
ConjugateQuaternion = [ RotationQuaternion(0), -RotationQuaternion(1), -RotationQuaternion(2), -RotationQuaternion(3) ]

! Multiplication of the rotated quaternion and the conjugate of the rotation quaternion to find the rotated vector
CALL QuaternionMultiplication( TempQuaternion, ConjugateQuaternion, RotatedQuaternion )

! Rotated vector
RotatedVector(1:3) = RotatedQuaternion(1:3)

RETURN

END SUBROUTINE VectorRotation

! *********************************************************************************************** !
!        This subroutine takes a rotation quaternion (qm) and combines it with a randomly         !
!        generated quaternion (qr) through quaternion multiplication, creating a randomly         !
!                                composed rotation quaternion (qn)                                !
! *********************************************************************************************** !
SUBROUTINE QuaternionCombination( ReferenceUnitQuaternion, ComposedUnitQuaternion, Angle )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: Angle                   ! Maximum angular displacement
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ReferenceUnitQuaternion ! Reference rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion    ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ComposedUnitQuaternion  ! Composed rotation quaternion

! Random quaternion generator
CALL RandomQuaternionGenerator( RandomUnitQuaternion, Angle )

! Quaternion multiplication (composed rotation)
CALL QuaternionMultiplication( RandomUnitQuaternion, ReferenceUnitQuaternion, ComposedUnitQuaternion )

RETURN

END SUBROUTINE QuaternionCombination

! *********************************************************************************************** !
!            This subroutine creates a composed rotation quaternion by multiplying the            !
!                         reference quaternion with a random quaternion                           !
! *********************************************************************************************** !
SUBROUTINE QuaternionMultiplication( RandomUnitQuaternion, ReferenceUnitQuaternion, ComposedUnitQuaternion )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion    ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ReferenceUnitQuaternion ! Reference quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ComposedUnitQuaternion  ! Composed quaternion

! Cross product of quaternions (qr × qm)
ComposedUnitQuaternion(0) = ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(0) ) - &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(1) ) - &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(2) ) - &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(1) = ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(0) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(1) ) - &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(2) ) + &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(2) = ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(0) ) + &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(1) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(2) ) - &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(3) = ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(0) ) - &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(1) ) + &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(2) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(3) )

RETURN

END SUBROUTINE QuaternionMultiplication
