! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains the subroutines used in the main program                    !
!                       to generate random quaternions and random vectors.                        !
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
IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
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
  IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
  IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
  Alpha = (2.D0 * RandomNumber) - 1.D0
  ! Uniform random number, β
  IF( RNGeneratorLogical(1) ) CALL Random_Number( RandomNumber )
  IF( RNGeneratorLogical(2) ) CALL RandomNumberGenBitwise(  )
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