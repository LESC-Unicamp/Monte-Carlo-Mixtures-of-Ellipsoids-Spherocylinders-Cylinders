! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains all subroutines used in the main program                    !
!                         accounting for quaternion algebraic operations.                         !
!                                                                                                 !
! Version number: 1.2.1                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       October 31st, 2023                                        !
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
!        This subroutine generates a random quaternion from a random angle and random axis        !
! *********************************************************************************************** !
SUBROUTINE RandomQuaternionGenerator( RandomUnitQuaternion, Angle )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: Angle                ! Maximum angular displacement
REAL( Kind= Real64 )                   :: RandomAngle          ! Random angle
REAL( Kind= Real64 ), DIMENSION( 3 )   :: RandomUnitVector     ! Random vector
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion ! Random unit quaternion

! Random angle
CALL RandomNumberGenLCG(  )
RandomAngle = ( (2.D0 * RandomNumber) - 1.D0 ) * Angle ! Range [-angmax,angmax]

! Random unit vector
CALL RandomVectorGenerator( RandomUnitVector )

! Quaternion algebra
RandomUnitQuaternion(0) = DCOS( RandomAngle * 0.5D0 )                       ! Real part
RandomUnitQuaternion(1) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(1) ! Imaginary part (Vector)
RandomUnitQuaternion(2) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(2) ! Imaginary part (Vector)
RandomUnitQuaternion(3) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(3) ! Imaginary part (Vector)

RETURN

END SUBROUTINE RandomQuaternionGenerator

! *********************************************************************************************** !
!            This subroutine generates a random vector on the surface of a unit sphere            !
!           See Allen and Tildesley, 2nd Edition (2017), page 514 for more information.           !
!        (Routine 'maths_module.f90' of Code A.1. in Marsaglia, Ann. Math. Statist., 1972)        !
! *********************************************************************************************** !   
SUBROUTINE RandomVectorGenerator( RandomUnitVector )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: Alpha, Beta, Zeta ! Random numbers
REAL( Kind= Real64 ), DIMENSION( 3 ) :: RandomUnitVector  ! Random vector

! Marsaglia's routine
UnitVectorCriterion: DO
  ! Uniform random number, α
  CALL RandomNumberGenLCG(  )
  Alpha = (2.D0 * RandomNumber) - 1.D0
  ! Uniform random number, β
  CALL RandomNumberGenLCG(  )
  Beta = (2.D0 * RandomNumber) - 1.D0
  ! Sum of squares, ζ
  Zeta  = (Alpha * Alpha) + (Beta * Beta)
  ! Marseglia's criterion
  IF( Zeta < 1.D0 ) THEN
    EXIT UnitVectorCriterion
  END IF
END DO UnitVectorCriterion

! Random vector
RandomUnitVector(1) = 2.D0 * Alpha * DSQRT( 1.D0 - Zeta )
RandomUnitVector(2) = 2.D0 * Beta * DSQRT( 1.D0 - Zeta )
RandomUnitVector(3) = 1.D0 - (2.D0 * Zeta)

RETURN

END SUBROUTINE RandomVectorGenerator

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