! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains the subroutines used in the main program                    !
!                          to compute the order parameter of the system.                          !
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
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!    This subroutine calculates the order parameter of a nematic phase via the Q-tensor method    !
! *********************************************************************************************** !
!                               Programmed by: Joyce Tavares Lopes                                !
!                               Modified by: Nathan Barros de Souza                               !
!                     University of Campinas, School of Chemical Engineering                      !
! *********************************************************************************************** !
!        See O. K. Smith, Communications of the ACM, 4(4), 168 (1961) for more information        !
! *********************************************************************************************** !
SUBROUTINE UniaxialNematicOrderParameter( OrderParameterS2, pCurrentOrientation )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: cComponent            ! Counter (component)
INTEGER( Kind= Int64 ) :: Particle              ! Counter (particle)
INTEGER( Kind= Int64 ) :: pNonSphericalParticle ! Counter (non-spherical particles)
INTEGER( Kind= Int64 ) :: Alpha, Beta           ! Unit vector specifiers (î, ĵ, k̂)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                             :: OrderParameterS2    ! Nematic order parameter
REAL( Kind= Real64 )                             :: KroneckerDelta      ! Kronecker delta
REAL( Kind= Real64 )                             :: M                   ! One-third of the trace (matrix)
REAL( Kind= Real64 )                             :: HQ                  ! One-half of the determinant (matrix)
REAL( Kind= Real64 )                             :: P                   ! Sum of squares of elements of a matrix
REAL( Kind= Real64 )                             :: Phi                 ! Angle
REAL( Kind= Real64 )                             :: TestCondition       ! Condition for real eigenvalues
REAL( Kind= Real64 ), DIMENSION( 3 )             :: Eigenvector         ! Eigenvector of the order tensor Q
REAL( Kind= Real64 ), DIMENSION( 3, 3 )          :: pTensorQ            ! Order tensor Q of particle i (3 x 3 Matrix)
REAL( Kind= Real64 ), DIMENSION( 3, 3 )          :: TensorQ             ! Order tensor Q of particle i (Average)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles ) :: pCurrentOrientation ! Order tensor Q of particle i (Average)

! Initialization
pTensorQ = 0.D0
pNonSphericalParticle = 0

DO cComponent = 1, nComponents
  ! Skip if component is spherical
  IF( SphericalComponentLogical(cComponent) ) THEN
    CYCLE
  END IF
  ! Q-tensor construction
  DO Particle = SUM( cParticles(0:(cComponent-1)) ) + 1, SUM( cParticles(0:cComponent) )
    DO Alpha = 1, 3
      DO Beta = 1, 3
        ! Dyadic product (Kronecker delta)
        IF ( Alpha == Beta ) THEN
          KroneckerDelta = 1.D0
        ELSE
          KroneckerDelta = 0.D0
        END IF
        ! Second-order Legendre polynomial (P₂)
        pTensorQ(Alpha,Beta) = pTensorQ(Alpha,Beta) + ( 1.5D0 * pCurrentOrientation(Alpha,Particle) * &
        &                      pCurrentOrientation(Beta,Particle) ) - (0.5D0 * KroneckerDelta)
      END DO
    END DO
    ! Increment counter of non-spherical particles
    pNonSphericalParticle = pNonSphericalParticle + 1
  END DO
END DO

! All components are spherical
IF( pNonSphericalParticle == 0 ) THEN
  OrderParameterS2 = 0.D0
  RETURN
END IF

! Averaged Q-tensor
TensorQ = pTensorQ / DBLE( pNonSphericalParticle )

! One-third of the trace of symmetric Q matrix
M = TensorQ(1,1) + TensorQ(2,2) + TensorQ(3,3)
M = M / 3.D0

! Half of the determinant of (Q - mI), where I is the identity matrix
HQ = ( TensorQ(1,1) - M ) * ( TensorQ(2,2) - M ) * ( TensorQ(3,3) - M ) + ( TensorQ(1,2) * TensorQ(2,3) * TensorQ(3,1) ) + &
&    ( TensorQ(1,3) * TensorQ(3,2) * TensorQ(2,1) ) - ( TensorQ(1,3) * ( TensorQ(2,2) - M ) * TensorQ(3,1) ) - &
&    ( TensorQ(2,3) * TensorQ(3,2) * ( TensorQ(1,1) - M ) ) - ( ( TensorQ(3,3) - M ) * TensorQ(2,1) * TensorQ(1,2) )
HQ = 0.5D0 * HQ

! One-sixth of the sum of squares of elements of (Q - mI)
P = ( TensorQ(1,1) - M ) * ( TensorQ(1,1) - M ) + ( TensorQ(1,2) * TensorQ(1,2) ) + ( TensorQ(1,3) * TensorQ(1,3) ) + &
&   ( TensorQ(2,2) - M ) * ( TensorQ(2,2) - M ) + ( TensorQ(2,1) * TensorQ(2,1) ) + ( TensorQ(2,3) * TensorQ(2,3) ) + &
&   ( TensorQ(3,3) - M ) * ( TensorQ(3,3) - M ) + ( TensorQ(3,1) * TensorQ(3,1) ) + ( TensorQ(3,2) * TensorQ(3,2) )
P = P / 6.D0

! Test condition
TestCondition = ( P * P * P ) - ( HQ * HQ )

! Real eigenvalues condition (P³ ≥ HQ²)
IF ( TestCondition >= 0.D0 ) THEN
  Phi = DATAN( DSQRT( TestCondition ) / HQ ) / 3.D0 ! 0 ≤ ϕ ≤ π
ELSE
  Phi = 0.D0
END IF

! Eigenvector of Q (phase director)
Eigenvector(1) = M + ( 2.D0 * DSQRT( P ) * DCOS( Phi ) )
Eigenvector(2) = M - DSQRT( P ) * ( DCOS( Phi ) + ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )
Eigenvector(3) = M - DSQRT( P ) * ( DCOS( Phi ) - ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )

! Nematic order parameter
OrderParameterS2 = MAXVAL( Eigenvector ) ! Largest eigenvalue

RETURN

END SUBROUTINE UniaxialNematicOrderParameter