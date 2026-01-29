! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                  This code contains the subroutines used in the main program                    !
!                          to compute the order parameter of the system.                          !
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
SUBROUTINE UniaxialNematicOrderParameter( OrderParameterS2, pCurrentOrientation, PhaseDirector )

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
REAL( Kind= Real64 ), DIMENSION( 3 )             :: PhaseDirector       ! Eigenvector of the order tensor Q associated with the largest eigenvalue
REAL( Kind= Real64 ), DIMENSION( 3 )             :: Eigenvalues         ! Eigenvalues of the order tensor Q
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

! Eigenvalues of Q
Eigenvalues(1) = M + ( 2.D0 * DSQRT( P ) * DCOS( Phi ) )
Eigenvalues(2) = M - DSQRT( P ) * ( DCOS( Phi ) + ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )
Eigenvalues(3) = M - DSQRT( P ) * ( DCOS( Phi ) - ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )

! Nematic order parameter
OrderParameterS2 = MAXVAL( Eigenvalues ) ! Largest eigenvalue

! Phase director
CALL NematicPhaseDirector( TensorQ, OrderParameterS2, PhaseDirector )

RETURN

END SUBROUTINE UniaxialNematicOrderParameter

! *********************************************************************************************** !
!         This subroutine solves a system of linear equations using Cramer's rule method.         !
! *********************************************************************************************** !
! OBS.: This calculation does not differ between antiparallel nematic phase directors, n and -n., !
!       which are usually distinguished for ferroelectric liquid crystals and related to the      !
!       dipole moment.                                                                            !
! *********************************************************************************************** !
SUBROUTINE NematicPhaseDirector( Matrix, OrderParameter, Eigenvector )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: i, j, k   ! Counters
INTEGER( Kind= Int64 ) :: ArraySize ! Size of the matrix
INTEGER( Kind= Int64 ) :: row, col  ! Rows and columns

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Determinant, DeterminantCramer ! Determinants
REAL( Kind= Real64 ) :: OrderParameter                 ! Nematic order parameter

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: Eigenvector ! Eigenvector associated with the largest eigenvalue
REAL( Kind= Real64 ), DIMENSION( 3, 3 ) :: Matrix      ! Input tensor Q matrix

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: AuxMatrix     ! Auxiliary matrix
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ReducedMatrix ! Reduced matrix
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: CramerMatrix  ! Cramer's matrix
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: ResultVector  ! Result vector

! Check if the matrix is square
IF( SIZE( Matrix, Dim= 1 ) /= SIZE( Matrix, Dim= 2 ) ) THEN
  WRITE( *, '(A)' ) 'The matrix is not square!'
  STOP
END IF

! Check if the matrix is 3x3
IF( SIZE( Matrix, Dim= 1 ) /= 3 ) THEN
  WRITE( *, '(A)' ) 'The matrix is not 3x3!'
  STOP
END IF

! Get the size of the matrix
ArraySize = SIZE( Matrix, Dim= 1 )

! Allocation
IF( ALLOCATED( ReducedMatrix ) ) DEALLOCATE( ReducedMatrix )
IF( ALLOCATED( AuxMatrix ) ) DEALLOCATE( AuxMatrix )
IF( ALLOCATED( ResultVector ) ) DEALLOCATE( ResultVector )
IF( ALLOCATED( CramerMatrix ) ) DEALLOCATE( CramerMatrix )
ALLOCATE( ReducedMatrix( ArraySize - 1, ArraySize - 1 ) )
ReducedMatrix = 0.D0
ALLOCATE( CramerMatrix( ArraySize - 1, ArraySize - 1 ) )
CramerMatrix = 0.D0
ALLOCATE( AuxMatrix( ArraySize, ArraySize ) )
AuxMatrix = 0.D0
ALLOCATE( ResultVector( ArraySize - 1 ) )
ResultVector = 0.D0

! Initialize the eigenvector
Eigenvector = 1.D0

! Rewrite input matrix by subtracting the diagonal matrix of the eigenvalues (A - λI)
DO i = 1, ArraySize
  DO j = 1, ArraySize
    AuxMatrix(i,j) = Matrix(i,j) - OrderParameter * IdentityMatrix(i,j)
  END DO
END DO

! Reduced matrix
k = 1
DO
  ! Stop if all reduced matrices are singular
  IF( k > 3 ) THEN
    WRITE( *, '(A)' ) 'Cannot find a the eigenvector of the nematic phase!'
    STOP
  END IF
  ! Trim the k-th row and k-th column of the auxiliary matrix
  row = 1
  DO i = 1, ArraySize
    IF( i == k ) CYCLE
    col = 1
    DO j = 1, ArraySize
      IF( j == k ) CYCLE
      ReducedMatrix(row,col) = AuxMatrix(i,j)
      col = col + 1
    END DO
    row = row + 1
  END DO
  ! Determinant of the reduced 3x3 matrix
  Determinant = ReducedMatrix(1,1) * ReducedMatrix(2,2) - ReducedMatrix(1,2) * ReducedMatrix(2,1)
  IF( DABS( Determinant ) <= 1.D-12 ) THEN ! Singular matrix
    k = k + 1 ! Trim another row and column
    CYCLE
  END IF
  EXIT
END DO

! Result vector (last column of the matrix)
j = 1
DO i = 1, ArraySize
  IF( i == k ) CYCLE
  ResultVector(j) = - Matrix(i,k) * Eigenvector(k) ! Xk is assumed to be 1 (could be any value)
  j = j + 1
END DO

! Cramer's matrix (substitution of the result vector in ith column)
j = 1
DO i = 1, SIZE( Eigenvector )
  IF( i == k ) CYCLE
  CramerMatrix = ReducedMatrix
  CramerMatrix(:,j) = ResultVector(:)
  ! Determinant of the Cramer's matrix
  DeterminantCramer = CramerMatrix(1,1) * CramerMatrix(2,2) - CramerMatrix(1,2) * CramerMatrix(2,1)
  ! Eigenvector elements (solution of the system of linear equations)
  Eigenvector(i) = DeterminantCramer / Determinant
  ! Next column of the Cramer's matrix
  j = j + 1
END DO

! Normalize the eigenvector
Eigenvector = Eigenvector / DSQRT( DOT_PRODUCT( Eigenvector, Eigenvector ) )

RETURN

END SUBROUTINE NematicPhaseDirector
