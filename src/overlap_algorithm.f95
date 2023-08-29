! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!         This code contains the overlap search subroutine for ellipsoids of revolution,          !
!                                 spherocylinders, and cylinders.                                 !
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
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!    This subroutine takes the rotation quaternions of two molecular ellipsoids of revolution     !
!     i and j and the distance vector joining their centers of mass and determinates whether      !
!                                 they overlap each other or not.                                 !
!          See Perram and Rasmussen, Phys. Rev. E 54, 6565 (1996), for more information.          !
! *********************************************************************************************** !
SUBROUTINE ELLIPSOID_OVERLAP( QI, QJ, RIJ, RIJSQ, CI, CJ, CONTACT_D, OVERLAP_HER )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: CI, CJ ! Counters

! *********************************************************************************************** !
! REAL VARIABLES (PARAMETER)                                                                      !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), PARAMETER :: TOLERANCE = 1.D-10  ! Numerical tolerance

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                    :: LAMBDAA     ! Interpolation parameter (extremum) [a]
REAL( KIND= REAL64 )                    :: LAMBDAB     ! Interpolation parameter (extremum) [b]
REAL( KIND= REAL64 )                    :: LAMBDAC     ! Interpolation parameter (midpoint) [c]
REAL( KIND= REAL64 )                    :: FA, FB      ! Derivatives of the interpolating function at the extrema of the interval: λ ∈ [0,1]
REAL( KIND= REAL64 )                    :: FC          ! Derivative of the interpolating function at the midpoint of the considered interval
REAL( KIND= REAL64 )                    :: AUX         ! Auxiliary variable
REAL( KIND= REAL64 )                    :: INTFUNC     ! Interpolating function (Perram-Wertheim method)
REAL( KIND= REAL64 )                    :: DFUNC_DLAMB ! First derivative of the interpolating function with respect to λ (Perram-Wertheim method)
REAL( KIND= REAL64 )                    :: RIJSQ       ! Magnitude of the vector distance between particles i and j (squared)
REAL( KIND= REAL64 )                    :: CONTACT_D   ! Contact distance (Perram-Wertheim method)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RIJ         ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EXI, EXJ    ! Orientation of particles i and j along x-direction
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EYI, EYJ    ! Orientation of particles i and j along y-direction
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EZI, EZJ    ! Orientation of particles i and j along z-direction
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: SEMIAI      ! Semiaxes of particle i
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: SEMIBJ      ! Semiaxes of particle j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DAUX1_DLAMB ! Auxiliary variable (derivative)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DAUX2_DLAMB ! Auxiliary variable (derivative)
REAL( KIND= REAL64 ), DIMENSION( 0:3 )  :: QI, QJ      ! Rotation quaternions of particle i and j
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EXIEXI      ! Outer product of the orientation of particle i along x-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EYIEYI      ! Outer product of the orientation of particle i along y-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EZIEZI      ! Outer product of the orientation of particle i along z-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EXJEXJ      ! Outer product of the orientation of particle j along x-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EYJEYJ      ! Outer product of the orientation of particle j along y-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EZJEZJ      ! Outer product of the orientation of particle j along z-direction
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: AI, BJ      ! Matrices of the semiaxes of the ellipsoids
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: AIINV       ! Inverse matrix of the semiaxis of ellipsoid i
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: BJINV       ! Inverse matrix of the semiaxis of ellipsoid j
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: G, GINV     ! Auxiliary matrix G and inverse matrix G
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: DG_DLAMB    ! Derivative of matrix G with respect to λ
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: IDMATRIX    ! Identity matrix

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP_HER ! Detects overlap between two ellipsoids of revolution: TRUE = overlap detected; FALSE = overlap not detected

! Ellipsoid centers of mass coincide (FA = FB = 0 | Bissection method cannot be used)
IF( DABS( RIJSQ - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
  OVERLAP_HER = .TRUE.
  INTFUNC     = 0.D0
  CONTACT_D   = 0.D0
  RETURN
END IF

! Initialization
OVERLAP_HER = .FALSE.

! Identity matrix
IDMATRIX(:,:) = 0.D0
IDMATRIX(1,1) = 1.D0
IDMATRIX(2,2) = 1.D0
IDMATRIX(3,3) = 1.D0

! Orientation of particle i along x-, y-, and z-directions
CALL ACTIVE_TRANSFORMATION( AXISX, QI, EXI )
CALL ACTIVE_TRANSFORMATION( AXISY, QI, EYI )
CALL ACTIVE_TRANSFORMATION( AXISZ, QI, EZI )

! Orientation of particle j along x-, y-, and z-direction
CALL ACTIVE_TRANSFORMATION( AXISX, QJ, EXJ )
CALL ACTIVE_TRANSFORMATION( AXISY, QJ, EYJ )
CALL ACTIVE_TRANSFORMATION( AXISZ, QJ, EZJ )

! Outer product of the orientation of particle i along x-direction
CALL OUTER_PRODUCT( EXI, EXI, EXIEXI )
! Outer product of the orientation of particle i along y-direction
CALL OUTER_PRODUCT( EYI, EYI, EYIEYI )
! Outer product of the orientation of particle i along z-direction
CALL OUTER_PRODUCT( EZI, EZI, EZIEZI )
! Outer product of the orientation of particle j along x-direction
CALL OUTER_PRODUCT( EXJ, EXJ, EXJEXJ )
! Outer product of the orientation of particle j along y-direction
CALL OUTER_PRODUCT( EYJ, EYJ, EYJEYJ )
! Outer product of the orientation of particle j along z-direction
CALL OUTER_PRODUCT( EZJ, EZJ, EZJEZJ )

! Semiaxes of particle i
SEMIAI(1) = 0.5D0 * DIAMETER(CI)
SEMIAI(2) = 0.5D0 * DIAMETER(CI)
SEMIAI(3) = 0.5D0 * LENGTH(CI)

! Semiaxes of particle j
SEMIBJ(1) = 0.5D0 * DIAMETER(CJ)
SEMIBJ(2) = 0.5D0 * DIAMETER(CJ)
SEMIBJ(3) = 0.5D0 * LENGTH(CJ)

! Matrix of the semiaxes of particle i
EXIEXI(:,:) = EXIEXI(:,:) / ( SEMIAI(1) * SEMIAI(1) )
EYIEYI(:,:) = EYIEYI(:,:) / ( SEMIAI(2) * SEMIAI(2) )
EZIEZI(:,:) = EZIEZI(:,:) / ( SEMIAI(3) * SEMIAI(3) )
AI(:,:)     = EXIEXI(:,:) + EYIEYI(:,:) + EZIEZI(:,:)
! Degenerate case (sphere)
IF( DABS( SEMIAI(1) - SEMIAI(2) ) < EPSILON( 1.D0 ) .AND. DABS( SEMIAI(1) - SEMIAI(3) ) < EPSILON( 1.D0 ) ) THEN
  AI(:,:) = IDMATRIX(:,:) / ( SEMIAI(1) * SEMIAI(1) )
END IF

! Matrix of the semiaxes of particle j
EXJEXJ(:,:) = EXJEXJ(:,:) / ( SEMIBJ(1) * SEMIBJ(1) )
EYJEYJ(:,:) = EYJEYJ(:,:) / ( SEMIBJ(2) * SEMIBJ(2) )
EZJEZJ(:,:) = EZJEZJ(:,:) / ( SEMIBJ(3) * SEMIBJ(3) )
BJ(:,:)     = EXJEXJ(:,:) + EYJEYJ(:,:) + EZJEZJ(:,:)
! Degenerate case (sphere)
IF( DABS( SEMIBJ(1) - SEMIBJ(2) ) < EPSILON( 1.D0 ) .AND. DABS( SEMIBJ(1) - SEMIBJ(3) ) < EPSILON( 1.D0 ) ) THEN
  BJ(:,:) = IDMATRIX(:,:) / ( SEMIBJ(1) * SEMIBJ(1) )
END IF

! Inverse matrix of the semiaxes of particle i
CALL INVERSE( AI, AIINV )

! Inverse matrix of the semiaxes of particle j
CALL INVERSE( BJ, BJINV )

! *********************************************************************************************** !
! Interpolating function at extremum of the interval [0,1] with λ = 0                             !
! *********************************************************************************************** !

! Extremum (λ = 0)
LAMBDAA = 0.D0

! Auxiliary matrix G
G(:,:) = LAMBDAA * BJINV(:,:) + (1.D0 - LAMBDAA) * AIINV(:,:)

! Inverse matrix G
CALL INVERSE( G, GINV )

! Auxiliary scalar (transpose of vector distance times inverse matrix G times vector distance)
AUX = 0.D0
AUX = AUX + ( RIJ(1) * GINV(1,1) + RIJ(2) * GINV(2,1) + RIJ(3) * GINV(3,1) ) * RIJ(1)
AUX = AUX + ( RIJ(1) * GINV(1,2) + RIJ(2) * GINV(2,2) + RIJ(3) * GINV(3,2) ) * RIJ(2)
AUX = AUX + ( RIJ(1) * GINV(1,3) + RIJ(2) * GINV(2,3) + RIJ(3) * GINV(3,3) ) * RIJ(3)

! First derivative of matrix G with respect to λ
DG_DLAMB(:,:) = ( (1.D0 - LAMBDAA) * (1.D0 - LAMBDAA) * AIINV(:,:) ) - ( LAMBDAA * LAMBDAA * BJINV(:,:) )

! Auxiliary vector (derivative)
DAUX1_DLAMB(1) = ( GINV(1,1) * RIJ(1) + GINV(1,2) * RIJ(2) + GINV(1,3) * RIJ(3) )
DAUX1_DLAMB(2) = ( GINV(2,1) * RIJ(1) + GINV(2,2) * RIJ(2) + GINV(2,3) * RIJ(3) )
DAUX1_DLAMB(3) = ( GINV(3,1) * RIJ(1) + GINV(3,2) * RIJ(2) + GINV(3,3) * RIJ(3) )

! Auxiliary vector (derivative)
DAUX2_DLAMB(1) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,1) + DAUX1_DLAMB(2) * DG_DLAMB(2,1) + DAUX1_DLAMB(3) * DG_DLAMB(3,1) )
DAUX2_DLAMB(2) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,2) + DAUX1_DLAMB(2) * DG_DLAMB(2,2) + DAUX1_DLAMB(3) * DG_DLAMB(3,2) )
DAUX2_DLAMB(3) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,3) + DAUX1_DLAMB(2) * DG_DLAMB(2,3) + DAUX1_DLAMB(3) * DG_DLAMB(3,3) )

! First derivative of the interpolating function with respect to λ (λ = 0)
DFUNC_DLAMB = 0.D0
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,1) + DAUX2_DLAMB(2) * GINV(2,1) + DAUX2_DLAMB(3) * GINV(3,1) ) * RIJ(1)
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,2) + DAUX2_DLAMB(2) * GINV(2,2) + DAUX2_DLAMB(3) * GINV(3,2) ) * RIJ(2)
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,3) + DAUX2_DLAMB(2) * GINV(2,3) + DAUX2_DLAMB(3) * GINV(3,3) ) * RIJ(3)
FA          = DFUNC_DLAMB

! *********************************************************************************************** !
! Interpolating function at extremum of the interval [0,1] with λ = 1                             !
! *********************************************************************************************** !

! Extremum (λ = 1)
LAMBDAB = 1.0D0

! Matrix G
G(:,:) = LAMBDAB * BJINV(:,:) + (1.D0 - LAMBDAB) * AIINV(:,:)

! Inverse Matrix G
CALL INVERSE( G, GINV )

! Auxiliary scalar (transpose of vector distance times inverse of matrix G times vector distance)
AUX = 0.D0
AUX = AUX + ( RIJ(1) * GINV(1,1) + RIJ(2) * GINV(2,1) + RIJ(3) * GINV(3,1) ) * RIJ(1)
AUX = AUX + ( RIJ(1) * GINV(1,2) + RIJ(2) * GINV(2,2) + RIJ(3) * GINV(3,2) ) * RIJ(2)
AUX = AUX + ( RIJ(1) * GINV(1,3) + RIJ(2) * GINV(2,3) + RIJ(3) * GINV(3,3) ) * RIJ(3)

! First derivative of matrix G with respect to λ
DG_DLAMB(:,:) = ( (1.D0 - LAMBDAB) * (1.D0 - LAMBDAB) * AIINV(:,:) ) - ( LAMBDAB * LAMBDAB * BJINV(:,:) )

! Auxiliary vector (derivative)
DAUX1_DLAMB(1) = ( GINV(1,1) * RIJ(1) + GINV(1,2) * RIJ(2) + GINV(1,3) * RIJ(3) )
DAUX1_DLAMB(2) = ( GINV(2,1) * RIJ(1) + GINV(2,2) * RIJ(2) + GINV(2,3) * RIJ(3) )
DAUX1_DLAMB(3) = ( GINV(3,1) * RIJ(1) + GINV(3,2) * RIJ(2) + GINV(3,3) * RIJ(3) )

! Auxiliary vector (derivative)
DAUX2_DLAMB(1) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,1) + DAUX1_DLAMB(2) * DG_DLAMB(2,1) + DAUX1_DLAMB(3) * DG_DLAMB(3,1) )
DAUX2_DLAMB(2) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,2) + DAUX1_DLAMB(2) * DG_DLAMB(2,2) + DAUX1_DLAMB(3) * DG_DLAMB(3,2) )
DAUX2_DLAMB(3) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,3) + DAUX1_DLAMB(2) * DG_DLAMB(2,3) + DAUX1_DLAMB(3) * DG_DLAMB(3,3) )

! First derivative of the interpolating function with respect to λ (λ = 1)
DFUNC_DLAMB = 0.D0
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,1) + DAUX2_DLAMB(2) * GINV(2,1) + DAUX2_DLAMB(3) * GINV(3,1) ) * RIJ(1)
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,2) + DAUX2_DLAMB(2) * GINV(2,2) + DAUX2_DLAMB(3) * GINV(3,2) ) * RIJ(2)
DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,3) + DAUX2_DLAMB(2) * GINV(2,3) + DAUX2_DLAMB(3) * GINV(3,3) ) * RIJ(3)
FB          = DFUNC_DLAMB

! *********************************************************************************************** !
! BISSECTION METHOD                                                                               !
! *********************************************************************************************** !
!  The bissection condition will hardly fail since the interpolating function, S(λ), is concave   !
!  down and has a single maximum in the interval λ ∈ [0,1], which means the derivative of S(λ)    !
!  has opposite signs at the extrema of the interval.                                             !
! *********************************************************************************************** !
!  See Perram and Wertheim, J. Comput. Phys 15, 409-416 (1985), for more information.             !
! *********************************************************************************************** !

! Bissection condition
IF( (FA * FB) < 0.D0 ) THEN

  ! Midpoint λ
  LAMBDAC = 0.5D0 * ( LAMBDAA + LAMBDAB )

  ! Auxiliary matrix G
  G(:,:) = LAMBDAC * BJINV(:,:) + (1.D0 - LAMBDAC) * AIINV(:,:)

  ! Inverse matrix G
  CALL INVERSE( G, GINV )

  ! Auxiliary scalar (transpose of vector distance times inverse matrix G times vector distance)
  AUX = 0.D0
  AUX = AUX + ( RIJ(1) * GINV(1,1) + RIJ(2) * GINV(2,1) + RIJ(3) * GINV(3,1) ) * RIJ(1)
  AUX = AUX + ( RIJ(1) * GINV(1,2) + RIJ(2) * GINV(2,2) + RIJ(3) * GINV(3,2) ) * RIJ(2)
  AUX = AUX + ( RIJ(1) * GINV(1,3) + RIJ(2) * GINV(2,3) + RIJ(3) * GINV(3,3) ) * RIJ(3)

  ! First derivative of matrix G with respect to λ
  DG_DLAMB(:,:) = ( (1.D0 - LAMBDAC) * (1.D0 - LAMBDAC) * AIINV(:,:) ) - ( LAMBDAC * LAMBDAC * BJINV(:,:) )

  ! Auxiliary vector (derivative)
  DAUX1_DLAMB(1) = ( GINV(1,1) * RIJ(1) + GINV(1,2) * RIJ(2) + GINV(1,3) * RIJ(3) )
  DAUX1_DLAMB(2) = ( GINV(2,1) * RIJ(1) + GINV(2,2) * RIJ(2) + GINV(2,3) * RIJ(3) )
  DAUX1_DLAMB(3) = ( GINV(3,1) * RIJ(1) + GINV(3,2) * RIJ(2) + GINV(3,3) * RIJ(3) )

  ! Auxiliary vector (derivative)
  DAUX2_DLAMB(1) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,1) + DAUX1_DLAMB(2) * DG_DLAMB(2,1) + DAUX1_DLAMB(3) * DG_DLAMB(3,1) )
  DAUX2_DLAMB(2) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,2) + DAUX1_DLAMB(2) * DG_DLAMB(2,2) + DAUX1_DLAMB(3) * DG_DLAMB(3,2) )
  DAUX2_DLAMB(3) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,3) + DAUX1_DLAMB(2) * DG_DLAMB(2,3) + DAUX1_DLAMB(3) * DG_DLAMB(3,3) )

  ! First derivative of the interpolating function with respect to λ (λ = λmid)
  DFUNC_DLAMB = 0.D0
  DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,1) + DAUX2_DLAMB(2) * GINV(2,1) + DAUX2_DLAMB(3) * GINV(3,1) ) * RIJ(1)
  DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,2) + DAUX2_DLAMB(2) * GINV(2,2) + DAUX2_DLAMB(3) * GINV(3,2) ) * RIJ(2)
  DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,3) + DAUX2_DLAMB(2) * GINV(2,3) + DAUX2_DLAMB(3) * GINV(3,3) ) * RIJ(3)
  FC          = DFUNC_DLAMB

  ! Stop criterion
  DO WHILE( DABS(FC) >= TOLERANCE )

    ! Bissection criterion
    IF( ( FA * FC ) > 0.D0 )THEN
      LAMBDAA = LAMBDAC
      FA      = FC
    ELSE
      LAMBDAB = LAMBDAC
      FB      = FC
    ENDIF

    ! New midpoint λ
    LAMBDAC = 0.5D0 * ( LAMBDAA + LAMBDAB )

    ! Matrix G
    G(:,:) = LAMBDAC * BJINV(:,:) + (1.D0 - LAMBDAC) * AIINV(:,:)

    ! Inverse matrix G
    CALL INVERSE( G, GINV )

    ! Auxiliary scalar (transpose of vector distance times inverse matrix G times vector distance)
    AUX = 0.D0
    AUX = AUX + ( RIJ(1) * GINV(1,1) + RIJ(2) * GINV(2,1) + RIJ(3) * GINV(3,1) ) * RIJ(1)
    AUX = AUX + ( RIJ(1) * GINV(1,2) + RIJ(2) * GINV(2,2) + RIJ(3) * GINV(3,2) ) * RIJ(2)
    AUX = AUX + ( RIJ(1) * GINV(1,3) + RIJ(2) * GINV(2,3) + RIJ(3) * GINV(3,3) ) * RIJ(3)

    ! First derivative of matrix G with respect to λ
    DG_DLAMB(:,:) = ( (1.D0 - LAMBDAC) * (1.D0 - LAMBDAC) * AIINV(:,:) ) - ( LAMBDAC * LAMBDAC * BJINV(:,:) )

    ! Auxiliary vector (derivative)
    DAUX1_DLAMB(1) = ( GINV(1,1) * RIJ(1) + GINV(1,2) * RIJ(2) + GINV(1,3) * RIJ(3) )
    DAUX1_DLAMB(2) = ( GINV(2,1) * RIJ(1) + GINV(2,2) * RIJ(2) + GINV(2,3) * RIJ(3) )
    DAUX1_DLAMB(3) = ( GINV(3,1) * RIJ(1) + GINV(3,2) * RIJ(2) + GINV(3,3) * RIJ(3) )

    ! Auxiliary vector (derivative)
    DAUX2_DLAMB(1) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,1) + DAUX1_DLAMB(2) * DG_DLAMB(2,1) + DAUX1_DLAMB(3) * DG_DLAMB(3,1) )
    DAUX2_DLAMB(2) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,2) + DAUX1_DLAMB(2) * DG_DLAMB(2,2) + DAUX1_DLAMB(3) * DG_DLAMB(3,2) )
    DAUX2_DLAMB(3) = ( DAUX1_DLAMB(1) * DG_DLAMB(1,3) + DAUX1_DLAMB(2) * DG_DLAMB(2,3) + DAUX1_DLAMB(3) * DG_DLAMB(3,3) )

    ! First derivative of the interpolating function with respect to λ (λ = λmid)
    DFUNC_DLAMB = 0.D0
    DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,1) + DAUX2_DLAMB(2) * GINV(2,1) + DAUX2_DLAMB(3) * GINV(3,1) ) * RIJ(1)
    DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,2) + DAUX2_DLAMB(2) * GINV(2,2) + DAUX2_DLAMB(3) * GINV(3,2) ) * RIJ(2)
    DFUNC_DLAMB = DFUNC_DLAMB + ( DAUX2_DLAMB(1) * GINV(1,3) + DAUX2_DLAMB(2) * GINV(2,3) + DAUX2_DLAMB(3) * GINV(3,3) ) * RIJ(3)
    FC          = DFUNC_DLAMB

  END DO

  ! Interpolating function, S(λ)
  INTFUNC = LAMBDAC * (1.D0 - LAMBDAC) * AUX

  ! Contact distance (Perram-Wertheim method)
  CONTACT_D = INTFUNC * RIJSQ

  ! Non-overlapping criterion
  IF( INTFUNC > 1.D0 ) THEN
    OVERLAP_HER = .FALSE.
  ELSE
    OVERLAP_HER = .TRUE.
  END IF

ELSE

  ! Error
  WRITE( *, "(G0)" ) "The Perram-Wertheim failed! Exiting..."
  CALL EXIT(  )

END IF

RETURN

END SUBROUTINE ELLIPSOID_OVERLAP

! *********************************************************************************************** !
!    This subroutine takes the relative orientations of two molecular spherocylinders i and j     !
!    and the unit vector joining their centers of mass and calculates their contact distance.     !
!           See Vega and Lago, Computers Chem. 18, 55-59 (1993), for more information.            !
! *********************************************************************************************** !
SUBROUTINE SPHEROCYLINDER_OVERLAP( EI, EJ, RIJ, RIJSQ, CI, CJ, DLAMBDAEI, DMUEJ, RVL, PARALLEL, OVERLAP_SPC )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: CI, CJ ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: RIJSQ          ! Squared vector distance between the centers of mass of particles i and j
REAL( KIND= REAL64 )                 :: RIJEI          ! Dot product of vector distance and orientation of particle i
REAL( KIND= REAL64 )                 :: RIJEJ          ! Dot product of vector distance and orientation of particle j
REAL( KIND= REAL64 )                 :: EIEJ           ! Dot product of both orientations (particles i and j)
REAL( KIND= REAL64 )                 :: DLAMBDA, DMU   ! Values that minimize r² (∂r²/∂μ = 0 and ∂r²/∂λ = 0)
REAL( KIND= REAL64 )                 :: CC, AUXI, AUXJ ! Auxiliary variables
REAL( KIND= REAL64 )                 :: RVL            ! Vega-Lago contact distance (variable)
REAL( KIND= REAL64 )                 :: SHORTEST_D     ! Shortest distance
REAL( KIND= REAL64 )                 :: SHORTEST_DSQ   ! Shortest distance (squared)
REAL( KIND= REAL64 ), DIMENSION( 2 ) :: HALFLENGTH     ! Half-length of spherocylinders i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: EI             ! Orientation of particle i
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: EJ             ! Orientation of particle j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: RIJ            ! Vector distance between the centers of mass of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DLAMBDAEI      ! Auxiliar vector (cylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DMUEJ          ! Auxiliar vector (cylinder overlap algorithm)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP_SPC ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PARALLEL    ! Checks whether spherocylinders are parallel or not

! Initialization
OVERLAP_SPC = .FALSE.
PARALLEL    = .FALSE.

! Shortest distance between two spherocylinders
SHORTEST_D   = 0.5D0 * ( DIAMETER(CI) + DIAMETER(CJ) )
SHORTEST_DSQ = SHORTEST_D * SHORTEST_D

! Half length of spherocylinder i
HALFLENGTH(1) = ( 0.5D0 * LENGTH(CI) )
! Half length of spherocylinder j
HALFLENGTH(2) = ( 0.5D0 * LENGTH(CJ) )

! Initial calculation
RIJEI = ( RIJ(1) * EI(1) ) + ( RIJ(2) * EI(2) ) + ( RIJ(3) * EI(3) )
RIJEJ = ( RIJ(1) * EJ(1) ) + ( RIJ(2) * EJ(2) ) + ( RIJ(3) * EJ(3) )
EIEJ  = ( EI(1) * EJ(1) ) + ( EI(2) * EJ(2) ) + ( EI(3) * EJ(3) )
CC    = 1.D0 - ( EIEJ * EIEJ )

! Checking whether the spherocylinders are parallel
IF( DABS( CC ) < 1.D-10 ) THEN
  PARALLEL = .TRUE.
  ! Checking whether the parallel spherocylinders are perpendicular to the intermolecular axis (avoid the indeterminate form 0/0)
  IF( DABS( RIJEI ) >= 1.D-10 .AND. DABS( RIJEJ ) >= 1.D-10 ) THEN
    ! Take the extreme side of particle i
    DLAMBDA = DSIGN( HALFLENGTH(1), RIJEI )
    ! Closest point between particle i and particle j
    DMU = ( DLAMBDA * EIEJ ) - RIJEJ
    ! Take the extreme side of particle j if λ' > L/2
    IF( DABS( DMU ) > HALFLENGTH(2) ) THEN
      DMU = DSIGN( HALFLENGTH(2), DMU )
    END IF
    ! Shortest distance (squared)
    RVL = RIJSQ + ( DLAMBDA * DLAMBDA ) + ( DMU * DMU ) - ( 2.D0 * DLAMBDA * DMU * EIEJ ) + ( 2.D0 * DMU * RIJEJ ) - &
    &     ( 2.D0 * DLAMBDA * RIJEI )
    ! Overlap criterion
    IF( RVL <= SHORTEST_DSQ ) THEN
      OVERLAP_SPC = .TRUE.
    END IF
    ! Vectors to be used in the overlap algorithm of cylinders
    DLAMBDAEI(:) = DLAMBDA * EI(:)
    DMUEJ(:)     = DMU * EJ(:)
    ! Return immediately
    RETURN
  ! Parallel spherocylinders almost orthogonal to the intermolecular axis (avoid the indeterminate form 0/0)
  ELSE
    DLAMBDA = 0.D0
    DMU     = 0.D0
    ! Shortest distance (squared)
    RVL = RIJSQ
    ! Overlap criterion
    IF( RVL <= SHORTEST_DSQ ) THEN
      OVERLAP_SPC = .TRUE.
    END IF
    ! Vectors to be used in the overlap algorithm of cylinders
    DLAMBDAEI(:) = 0.D0
    DMUEJ(:)     = 0.D0
    ! Return immediately
    RETURN
  END IF
END IF

! *********************************************************************************************** !
! STEP 1: Evaluation of (λ’, μ’) according to Equations (3) and (4)                               !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
DLAMBDA = ( RIJEI - ( EIEJ * RIJEJ ) ) / CC
DMU     = ( -RIJEJ + ( EIEJ * RIJEI ) ) / CC

! *********************************************************************************************** !
! STEP 2: Check whether the point (λ’, μ’) is in the rectangle (λ, μ), with λ = μ = [-L/2, L/2].  !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
! Point (λ’, μ’) in the rectangle (λ, μ)
IF( ( DABS( DLAMBDA ) <= HALFLENGTH(1) ) .AND. ( DABS( DMU ) <= HALFLENGTH(2) ) ) THEN
  RVL = RIJSQ + ( DLAMBDA * DLAMBDA ) + ( DMU * DMU ) - ( 2.D0 * DLAMBDA * DMU * EIEJ ) + ( 2.D0 * DMU * RIJEJ ) - &
  &     ( 2.D0 * DLAMBDA * RIJEI )
  ! Overlap criterion
  IF( RVL <= SHORTEST_DSQ ) THEN
    OVERLAP_SPC = .TRUE.
  END IF
  ! Vectors to be used in the overlap algorithm of cylinders
  DLAMBDAEI(:) = DLAMBDA * EI(:)
  DMUEJ(:)     = DMU * EJ(:)
  ! Return immediately
  RETURN
! Point (λ’, μ’) not in the rectangle (λ, μ)
ELSE
  AUXI = DABS( DLAMBDA ) - HALFLENGTH(1) ! Used to determine the region (extreme)
  AUXJ = DABS( DMU ) - HALFLENGTH(2)     ! Used to determine the region (extreme)
END IF

! *********************************************************************************************** !
! STEP 3-7: Shortest distance of the considered extreme (regions 1, 2, 3, or 4) to the line where !
! the other rod is contained                                                                      !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
! Region 1 or 3
IF( AUXI > AUXJ ) THEN
  DLAMBDA = DSIGN( HALFLENGTH(1), DLAMBDA )
  DMU = ( DLAMBDA * EIEJ ) - RIJEJ
  IF( DABS( DMU ) > HALFLENGTH(2) ) THEN
    DMU = DSIGN( HALFLENGTH(2), DMU )
  END IF
! Region 2 or 4
ELSE
  DMU = DSIGN( HALFLENGTH(2), DMU )
  DLAMBDA = ( DMU * EIEJ ) + RIJEI
  IF( DABS( DLAMBDA ) > HALFLENGTH(1) ) THEN
    DLAMBDA = DSIGN( HALFLENGTH(1), DLAMBDA )
  END IF
END IF

! *********************************************************************************************** !
! STEP 8: Evaluate shortest distance (squared)                                                    !
! *********************************************************************************************** !
RVL = RIJSQ + ( DLAMBDA * DLAMBDA ) + ( DMU * DMU ) - ( 2.D0 * DLAMBDA * DMU * EIEJ ) + ( 2.D0 * DMU * RIJEJ ) - &
&     ( 2.D0 * DLAMBDA * RIJEI )
! Overlap criterion
IF( RVL <= SHORTEST_DSQ ) THEN
  OVERLAP_SPC = .TRUE.
END IF

! Vectors to be used in the overlap algorithm of cylinders
DLAMBDAEI(:) = DLAMBDA * EI(:)
DMUEJ(:)     = DMU * EJ(:)

! Return immediately
RETURN

END SUBROUTINE SPHEROCYLINDER_OVERLAP

! *********************************************************************************************** !
! This subroutine takes the relative quaternions/orientations of two molecular cylinders i and j  !
!      as well as the position of their centers of mass and the unit vector joining them and      !
!                        calculates whether the cylinders overlap or not.                         !
!             See Lopes et al., Chem. Phys. 154, 104902 (2021), for more information.             !
! *********************************************************************************************** !
SUBROUTINE CYLINDER_OVERLAP( QI, QJ, EI, EJ, RIJ, RI, RJ, CI, CJ, DLAMBDAEI, DMUEJ, PARALLEL, OVERLAP_CYL )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I, J, K ! Counters
INTEGER( KIND= INT64 ) :: CI, CJ  ! Counters (component)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                    :: LSQ            ! Cylinder geometrical length (squared)
REAL( KIND= REAL64 )                    :: DSQ            ! Cylinder geometrical diameter (squared)
REAL( KIND= REAL64 )                    :: RIJEI          ! Dot product of vector distance and orientation of particle i
REAL( KIND= REAL64 )                    :: RIJEJ          ! Dot product of vector distance and orientation of particle j
REAL( KIND= REAL64 )                    :: EIEJ           ! Dot product of both orientations (particles i and j)
REAL( KIND= REAL64 )                    :: RIJSQ_PARALLEL ! Squared vector distance between particles i and j (parallel)
REAL( KIND= REAL64 )                    :: RIJSQ_ORTHO    ! Squared vector distance between particles i and j (orthogonal)
REAL( KIND= REAL64 )                    :: CC             ! Auxiliary variable
REAL( KIND= REAL64 ), DIMENSION( 2 )    :: DSQ_DISKRIM    ! Diameter of the disk for the disk-rim configuration
REAL( KIND= REAL64 ), DIMENSION( 2 )    :: HALFLENGTH     ! Half length of cylinders i and j
REAL( KIND= REAL64 ), DIMENSION( 2 )    :: HALFDIAMETER   ! Half diameter of cylinders i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DLAMBDAEI      ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DMUEJ          ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EI, EJ         ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RI, RJ         ! Position of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RIJ            ! Vector distance between particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RIJ_PARALLEL   ! Vector distance between particles i and j (parallel)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RIJ_ORTHO      ! Vector distance between particles i and j (orthogonal)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RRIM           ! Position of cylinder rim
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: ERIM           ! Orientation of cylinder rim
REAL( KIND= REAL64 ), DIMENSION( 0: 3 ) :: QI, QJ         ! Quaternion of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 0: 3 ) :: QDISK          ! Quaternion of cylinder disks
REAL( KIND= REAL64 ), DIMENSION( 2, 3 ) :: DI, DJ         ! Position of cylinder disks of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAP_CYL ! Detects overlap between two cylindrical particles : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAPRIM  ! Detects overlap between two particles (rim-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAPDISK ! Detects overlap between two particles (disk-disk configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OVERLAPDRIM ! Detects overlap between two particles (disk-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PARALLEL    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Half length of cylinder i
HALFLENGTH(1)   = 0.5D0 * LENGTH(CI)
! Half length of cylinder j
HALFLENGTH(2)   = 0.5D0 * LENGTH(CJ)
! Half diameter of cylinder i
HALFDIAMETER(1) = 0.5D0 * DIAMETER(CI)
! Half diameter of cylinder j
HALFDIAMETER(2) = 0.5D0 * DIAMETER(CJ)

! Cylindrical length (squared)
LSQ = 0.5D0 * ( LENGTH(CI) + LENGTH(CJ) )
LSQ = LSQ * LSQ
! Cylindrical diameter (squared)
DSQ = 0.5D0 * ( DIAMETER(CI) + DIAMETER(CJ) )
DSQ = DSQ * DSQ

! Diameter of the disks (squared)
DSQ_DISKRIM(1) = DIAMETER(CI) * DIAMETER(CI)
DSQ_DISKRIM(2) = DIAMETER(CJ) * DIAMETER(CJ)

! Initialize logical variables
OVERLAP_CYL  = .FALSE.

! Initial calculation
RIJEI = ( RIJ(1) * EI(1) ) + ( RIJ(2) * EI(2) ) + ( RIJ(3) * EI(3) )
RIJEJ = ( RIJ(1) * EJ(1) ) + ( RIJ(2) * EJ(2) ) + ( RIJ(3) * EJ(3) )
EIEJ  = ( EI(1) * EJ(1) ) + ( EI(2) * EJ(2) ) + ( EI(3) * EJ(3) )
CC    = 1.D0 - ( EIEJ * EIEJ )

! *********************************************************************************************** !
! CASE 1: Parallel cylinders                                                                      !
! *********************************************************************************************** !
IF( PARALLEL ) THEN
  ! Vector distance between cylinders i and j (parallel)
  RIJ_PARALLEL(1) = EI(1) * RIJEI
  RIJ_PARALLEL(2) = EI(2) * RIJEI
  RIJ_PARALLEL(3) = EI(3) * RIJEI
  ! Vector distance between cylinders i and j (orthogonal)
  RIJ_ORTHO(1)    = RIJ(1) - RIJ_PARALLEL(1)
  RIJ_ORTHO(2)    = RIJ(2) - RIJ_PARALLEL(2)
  RIJ_ORTHO(3)    = RIJ(3) - RIJ_PARALLEL(3)
  ! Magnitude of vector distance between cylinders i and j (parallel)
  RIJSQ_PARALLEL  = ( RIJ_PARALLEL(1) * RIJ_PARALLEL(1) ) + ( RIJ_PARALLEL(2) * RIJ_PARALLEL(2) ) + ( RIJ_PARALLEL(3) * &
  &                 RIJ_PARALLEL(3) )
  ! Magnitude of vector distance between cylinders i and j (orthogonal)
  RIJSQ_ORTHO     = ( RIJ_ORTHO(1) * RIJ_ORTHO(1) ) + ( RIJ_ORTHO(2) * RIJ_ORTHO(2) ) + ( RIJ_ORTHO(3) * RIJ_ORTHO(3) )
  ! Overlap criterion
  IF( ( RIJSQ_PARALLEL <= LSQ ) .AND. ( RIJSQ_ORTHO <= DSQ ) ) THEN
    OVERLAP_CYL = .TRUE.
    RETURN
  ELSE
    OVERLAP_CYL = .FALSE.
    RETURN
  END IF
END IF

! *********************************************************************************************** !
! CASE 2: Rim-rim configuration                                                                   !
! *********************************************************************************************** !
CALL RIM_RIM( EI, EJ, RIJ, DLAMBDAEI, DMUEJ, HALFLENGTH, OVERLAPRIM )
IF( OVERLAPRIM ) THEN
  OVERLAP_CYL = .TRUE.
  RETURN
END IF

! *********************************************************************************************** !
! CASE 3: Disk-disk configuration                                                                 !
! *********************************************************************************************** !
! Disk 1 of cylinder i
DI(1,1) = RI(1) + ( EI(1) * HALFLENGTH(1) )
DI(1,2) = RI(2) + ( EI(2) * HALFLENGTH(1) )
DI(1,3) = RI(3) + ( EI(3) * HALFLENGTH(1) )
! Disk 2 of cylinder i
DI(2,1) = RI(1) - ( EI(1) * HALFLENGTH(1) )
DI(2,2) = RI(2) - ( EI(2) * HALFLENGTH(1) )
DI(2,3) = RI(3) - ( EI(3) * HALFLENGTH(1) )
! Disk 1 of cylinder j
DJ(1,1) = RJ(1) + ( EJ(1) * HALFLENGTH(2) )
DJ(1,2) = RJ(2) + ( EJ(2) * HALFLENGTH(2) )
DJ(1,3) = RJ(3) + ( EJ(3) * HALFLENGTH(2) )
! Disk 2 of cylinder j
DJ(2,1) = RJ(1) - ( EJ(1) * HALFLENGTH(2) )
DJ(2,2) = RJ(2) - ( EJ(2) * HALFLENGTH(2) )
DJ(2,3) = RJ(3) - ( EJ(3) * HALFLENGTH(2) )
! Search for overlaps between all pairs of disks from both cylinders
DO I = 1, 2
  DO J = 1, 2
    CALL DISK_DISK( EI, EJ, DI(I,:), DJ(J,:), CC, HALFDIAMETER, OVERLAPDISK )
    IF( OVERLAPDISK ) THEN
      OVERLAP_CYL = .TRUE.
      RETURN
    END IF
  END DO
END DO

! *********************************************************************************************** !
! CASE 4: Disk-rim configuration (modified)                                                       !
! *********************************************************************************************** !
! Rim of cylinder i
RRIM(1)  = RI(1)
RRIM(2)  = RI(2)
RRIM(3)  = RI(3)
! Orientation of cylinder i
ERIM(1)  = EI(1)
ERIM(2)  = EI(2)
ERIM(3)  = EI(3)
! Quaternion of the disks of cylinder j
QDISK(0) = QJ(0)
QDISK(1) = QJ(1)
QDISK(2) = QJ(2)
QDISK(3) = QJ(3)
! Search for overlaps between both disks of cylinder j and the rim of cylinder i
DO K = 1, 2
  CALL DISK_RIM( DJ(K,:), QDISK, HALFDIAMETER(2), RRIM, ERIM, HALFDIAMETER(1), HALFLENGTH(1), OVERLAPDRIM )
  IF( OVERLAPDRIM ) THEN
    OVERLAP_CYL = .TRUE.
    RETURN
  END IF
END DO
! Rim of cylinder j
RRIM(1)  = RJ(1)
RRIM(2)  = RJ(2)
RRIM(3)  = RJ(3)
! Orientation of cylinder j
ERIM(1)  = EJ(1)
ERIM(2)  = EJ(2)
ERIM(3)  = EJ(3)
! Quaternion of the disks of cylinder i
QDISK(0) = QI(0)
QDISK(1) = QI(1)
QDISK(2) = QI(2)
QDISK(3) = QI(3)
! Search for overlaps between both disks of cylinder i and the rim of cylinder j
DO K = 1, 2
  CALL DISK_RIM( DI(K,:), QDISK, HALFDIAMETER(1), RRIM, ERIM, HALFDIAMETER(2), HALFLENGTH(2), OVERLAPDRIM )
  IF( OVERLAPDRIM ) THEN
    OVERLAP_CYL = .TRUE.
    RETURN
  END IF
END DO

RETURN

END SUBROUTINE CYLINDER_OVERLAP

! *********************************************************************************************** !
!                                   Rim-Rim Overlap Algorithm                                     !
! *********************************************************************************************** !
SUBROUTINE RIM_RIM( EI, EJ, RIJ, DLAMBDAEI, DMUEJ, HALFL, OVERLAPRIM )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: PROJI, PROJJ ! Cylindrical projections
REAL( KIND= REAL64 ), DIMENSION( 2 ) :: HALFL        ! Half length of cylinders i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DLAMBDAEI    ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DMUEJ        ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: VECI, VECJ   ! Vector distance between the point of closest approach in one cylinder and the center of the other cylinder
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: EI, EJ       ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: RIJ          ! Vector distance between particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAPRIM ! Detects overlap between two particles (rim-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected

! Initialize logical variable
OVERLAPRIM = .FALSE.

! Vector distance between the point of closest approach on cylinder i and the center of cylinder j
VECI(1) = -RIJ(1) + DLAMBDAEI(1)
VECI(2) = -RIJ(2) + DLAMBDAEI(2)
VECI(3) = -RIJ(3) + DLAMBDAEI(3)

! Vector distance between the point of closest approach on cylinder j and the center of cylinder i
VECJ(1) = RIJ(1) + DMUEJ(1)
VECJ(2) = RIJ(2) + DMUEJ(2)
VECJ(3) = RIJ(3) + DMUEJ(3)

! Projection of the vector distance j along orientation of cylinder i
PROJI = ( VECJ(1) * EI(1) ) + ( VECJ(2) * EI(2) ) + ( VECJ(3) * EI(3) )

! Projection of the vector distance i along orientation of cylinder j
PROJJ = ( VECI(1) * EJ(1) ) + ( VECI(2) * EJ(2) ) + ( VECI(3)* EJ(3) )

! Overlap criterion
IF( ( DABS( PROJI ) <= HALFL(1) ) .AND. ( DABS( PROJJ ) <= HALFL(2) ) ) THEN
  OVERLAPRIM = .TRUE.
  RETURN
END IF

! *********************************************************************************************** !
! -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   OBSERVATION   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* !
! *********************************************************************************************** !
! The overlap criterion must also assess whether the distance between the points of closest       !
! approach falls within the boundaries defined by the combined diameter of both cylinders.        !
! However, given that we perform an initial test involving the circumscribing spherocylinders,    !
! this condition, while essential, is irrelevant in this context.                                 !
! *********************************************************************************************** !

RETURN

END SUBROUTINE RIM_RIM

! *********************************************************************************************** !
!                                  Disk-Disk Overlap Algorithm                                    !
! *********************************************************************************************** !
SUBROUTINE DISK_DISK( EI, EJ, DI, DJ, CC, HALFD, OVERLAPDISK )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                 :: DIJSQ  ! Vector distance between disks i and j (squared)
REAL( KIND= REAL64 )                 :: IDIJSQ ! Magnitude of the shortest distance between disk i and intersection line between planes of the disks
REAL( KIND= REAL64 )                 :: JDIJSQ ! Magnitude of the shortest distance between disk j and intersection line between planes of the disks
REAL( KIND= REAL64 )                 :: SPH_SQ ! Magnitude of the vector distance between the circumscribing spheres (squared)
REAL( KIND= REAL64 )                 :: SEG    ! Segment (sum)
REAL( KIND= REAL64 )                 :: SEGI   ! Segment of cylinder i
REAL( KIND= REAL64 )                 :: SEGJ   ! Segment of cylinder j
REAL( KIND= REAL64 )                 :: CC     ! Auxiliary variable
REAL( KIND= REAL64 )                 :: MODEIJ ! Magnitude of vector eij
REAL( KIND= REAL64 )                 :: PIPJ   ! Projection of dij on eij
REAL( KIND= REAL64 ), DIMENSION( 2 ) :: RSQ    ! Radius of cylinders i and j (squared)
REAL( KIND= REAL64 ), DIMENSION( 2 ) :: HALFD  ! Half diameter of cylinders i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DIJ    ! Vector distance between disks i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: EI, EJ ! Orientation of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: EIJ    ! Cross product of orientations of particles i and j
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: UEIJ   ! Versor of vector eij
REAL( KIND= REAL64 ), DIMENSION( 3 ) :: DI, DJ ! Position of both disks of cylinders (particles) i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAPDISK ! Detects overlap between two particles (disk-disk configuration) : TRUE = overlap detected; FALSE = overlap not detected

! Initialize logical variable
OVERLAPDISK = .FALSE.

! Magnitude of the vector distance between the circumscribing spheres (squared)
SPH_SQ = ( HALFD(1) + HALFD(2) ) * ( HALFD(1) + HALFD(2) )

! Vector distance between both disks
DIJ(1) = DJ(1) - DI(1)
DIJ(2) = DJ(2) - DI(2)
DIJ(3) = DJ(3) - DI(3)

! Magnitude of the vector distance between both disks
DIJSQ  = ( DIJ(1) * DIJ(1) ) + ( DIJ(2) * DIJ(2) ) + ( DIJ(3) * DIJ(3) )

! Non-overlapping condition (spheres circumscribing the disks)
IF( DIJSQ > SPH_SQ ) THEN
  RETURN
END IF

! Magnitude of the shortest distance between disk i and intersection line between planes of the disks
IDIJSQ = ( DIJ(1) * EJ(1) ) + ( DIJ(2) * EJ(2) ) + ( DIJ(3) * EJ(3) )
IDIJSQ = ( IDIJSQ * IDIJSQ )
IDIJSQ = ( IDIJSQ / CC )

! Magnitude of the shortest distance between disk j and intersection line between planes of the disks
JDIJSQ = ( DIJ(1) * EI(1) ) + ( DIJ(2) * EI(2) ) + ( DIJ(3) * EI(3) )
JDIJSQ = ( JDIJSQ * JDIJSQ )
JDIJSQ = ( JDIJSQ / CC )

! Radius (squared)
RSQ(1) = ( HALFD(1) * HALFD(1) )
RSQ(2) = ( HALFD(2) * HALFD(2) )

! Initial overlap condition
IF( ( IDIJSQ < RSQ(1) ) .AND. ( JDIJSQ < RSQ(2) ) ) THEN
  ! Lengths of segments symmetrical to the intersection line between planes of the disks
  SEGI = DSQRT( RSQ(1) - IDIJSQ )
  SEGJ = DSQRT( RSQ(2) - JDIJSQ )
  ! Sum of segments
  SEG = SEGI + SEGJ
  ! Cross product of orientations of disks of cylinders i and j (orientation of intersection line)
  EIJ(1) = ( EI(2) * EJ(3) ) - ( EI(3) * EJ(2) )
  EIJ(2) = ( EI(3) * EJ(1) ) - ( EI(1) * EJ(3) )
  EIJ(3) = ( EI(1) * EJ(2) ) - ( EI(2) * EJ(1) )
  ! Magnitude of the orientation of intersection line
  MODEIJ = ( EIJ(1) * EIJ(1) ) + ( EIJ(2) * EIJ(2) ) + ( EIJ(3) * EIJ(3) )
  MODEIJ = DSQRT( MODEIJ )
  ! Versor of the orientation of intersection line
  UEIJ(1) = ( EIJ(1) / MODEIJ )
  UEIJ(2) = ( EIJ(2) / MODEIJ )
  UEIJ(3) = ( EIJ(3) / MODEIJ )
  ! Projection of the vector distance of both disks along the intersection line between planes of the disks
  PIPJ = ( DIJ(1) * UEIJ(1) ) + ( DIJ(2) * UEIJ(2) ) + ( DIJ(3) * UEIJ(3) )
  ! Overlap criterion
  IF( DABS( PIPJ ) <= SEG ) THEN
    OVERLAPDISK = .TRUE.
    RETURN
  END IF
END IF

RETURN

END SUBROUTINE DISK_DISK

! *********************************************************************************************** !
!                                   Disk-Rim Overlap Algorithm                                    !
! *********************************************************************************************** !
SUBROUTINE DISK_RIM( DK, QDISK, HALFDDISK, RRIM, ERIM, HALFDRIM, HALFLRIM, OVERLAPDRIM )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( KIND= INT64 ) :: I        ! Counter
INTEGER( KIND= INT64 ) :: COUNTER  ! Iteration counter (numerical methods)
INTEGER( KIND= INT64 ) :: COUNTERB ! Cycles of the bissection method
INTEGER( KIND= INT64 ) :: INTERVAL ! Loop over λ intervals
INTEGER( KIND= INT64 ) :: N_POINTS ! Maximum number of intervals

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                    :: DKRRIM_ERIM               ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylinder rim (i or j)
REAL( KIND= REAL64 )                    :: DKRRIM_EXDISK             ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylinder disk along x-direction (j or i)
REAL( KIND= REAL64 )                    :: DKRRIM_EYDISK             ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylinder disk along y-direction (j or i)
REAL( KIND= REAL64 )                    :: EXDISK_ERIM               ! Dot product of the orientation of the cylinder disk (j or i) along x-direction and the orientation of the cylinder rim (i or j)
REAL( KIND= REAL64 )                    :: EYDISK_ERIM               ! Dot product of the orientation of the cylinder disk (j or i) along y-direction and the orientation of the cylinder rim (i or j)
REAL( KIND= REAL64 )                    :: DKURIMSQ                  ! Magnitude of the vector distance between the closest point on the cylinder rim (i or j) and the center of the cylinder disk (j or i)
REAL( KIND= REAL64 )                    :: DAMPING                   ! Damping factor of the Newton-Raphson method
REAL( KIND= REAL64 )                    :: TOLERANCE_NR, TOLERANCE_B ! Numerical tolerance of the Newton-Raphson and bissection methods
REAL( KIND= REAL64 )                    :: FACTOR                    ! Newton-Raphson's factor (function over its derivative)
REAL( KIND= REAL64 )                    :: COSPHI, SINPHI            ! Cossine and sine of the angle φ that defines a point in the circumference of the cylinder disk (j or i)
REAL( KIND= REAL64 )                    :: OPPOSITE_CAT              ! Opposite cathetus
REAL( KIND= REAL64 )                    :: ADJACENT_CAT              ! Adjacent cathetus
REAL( KIND= REAL64 )                    :: HYP                       ! Hypothenuse
REAL( KIND= REAL64 )                    :: TPARALSQ                  ! Magnitude of vector T parallel to the cylinder rim (i or j) (squared)
REAL( KIND= REAL64 )                    :: TORTHOSQ                  ! Magnitude of vector T orthogonal to the cylinder rim (i or j) (squared)
REAL( KIND= REAL64 )                    :: TSQ                       ! Magnitude of the vector T (squared)
REAL( KIND= REAL64 )                    :: FI                        ! Objective function (midpoint)
REAL( KIND= REAL64 )                    :: LAMBDAI                   ! Minimization variable related to a point on the cylinder rim (i or j) (midpoint)
REAL( KIND= REAL64 )                    :: RSQ                       ! Radius of the cylinder (squared)
REAL( KIND= REAL64 )                    :: DSQ                       ! Diameter of the combined cylinder disks (squared)
REAL( KIND= REAL64 )                    :: HALFDDISK                 ! Half diameter of the cylinder disk
REAL( KIND= REAL64 )                    :: HALFDRIM                  ! Half diameter of the cylinder rim
REAL( KIND= REAL64 )                    :: HALFLRIM                  ! Half length of the cylinder rim
REAL( KIND= REAL64 )                    :: PMINDISKSQ                ! Distance between the center of the cylinder disk and the intersection point (squared)
REAL( KIND= REAL64 )                    :: PMINCYLSQ                 ! Distance between the center of the cylinder rim and the intersection point (squared)
REAL( KIND= REAL64 )                    :: ALPHA, BETA, GAMMA        ! Coefficients
REAL( KIND= REAL64 )                    :: CARDANO_Q, CARDANO_R      ! Cardano's coefficients
REAL( KIND= REAL64 )                    :: CARDANO_S, CARDANO_T      ! Cardano's coefficients
REAL( KIND= REAL64 )                    :: THETA, ARGUMENT           ! Cardano's coefficients
REAL( KIND= REAL64 )                    :: DISCRIMINANT              ! Cardano's discriminant
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: ROOTS                     ! Roots of the cubic equation
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: F                         ! Objective function
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DF                        ! Derivative of the objective function with respect to λ
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: LAMBDA                    ! Minimization variable related to a point on the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: RRIM                      ! Position of the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: ERIM                      ! Orientation of the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EXDISK, EYDISK, EZDISK    ! Orientation of the cylinder disks along x-, y-, and z-directions (j or i)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DK                        ! Position of cylinder disks (j or i)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DKRRIM                    ! Vector distance between a cylinder disk (j or i) and a cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: DKURIM                    ! Vector distance between a cylinder disk (j or i) and the closest point on the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: URIM                      ! Closest point on the cylinder rim (i or j) with respect to the center of the cylinder disk (j or i)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: PC                        ! Closest point on the axial axis of the cylinder rim (i or j) with respect to the circumference of the cylinder disk (j or i)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: PD                        ! Closest point on the circumference of the cylinder disk (j or i) with respect to the axial axis of the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: T                         ! Vector distance between the closest point on the circumference of the cylinder disk (j or i) and the center of the cylinder rim (i or j)
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: PK                        ! Point in a plane defined by the cylinder disk
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: LA, LB                    ! Points in a line defined by the cylinder rim
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: INTERSECTION              ! Intersection point between cylinder axis and cylinder disk
REAL( KIND= REAL64 ), DIMENSION( 5 )    :: COEFF_QUARTIC             ! Coefficients of the quartic function
REAL( KIND= REAL64 ), DIMENSION( 4 )    :: COEFF_CARDANO             ! Coefficients of the cubic function
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: CPOINTS                   ! Points where the derivative of the objective function with respect to λ are 0
REAL( KIND= REAL64 ), DIMENSION( 0: 3 ) :: QDISK                     ! Quaternion of the cylinder disks (j or i)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: OVERLAPDRIM   ! Detects overlap between two particles (disk-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: BISECT_METHOD ! Detects whether the bissection method will be used after n iterations of the Newton-Raphson method

! Initialize logical variable
OVERLAPDRIM   = .FALSE.
BISECT_METHOD = .FALSE.

! Initialize counters
COUNTER  = 0
COUNTERB = 0

! Radius of the cylinder rim (squared)
RSQ = ( HALFDRIM * HALFDRIM )

! Diameter of the combined cylinder disks (squared)
DSQ = ( HALFDDISK + HALFDRIM )
DSQ = ( DSQ * DSQ )

! Vector distance between a cylinder disk and a cylinder rim
DKRRIM(1)   = DK(1) - RRIM(1)
DKRRIM(2)   = DK(2) - RRIM(2)
DKRRIM(3)   = DK(3) - RRIM(3)
! Projection of the vector distance between a cylinder disk and a cylinder rim along the orientation of the cylinder rim
DKRRIM_ERIM = ( DKRRIM(1) * ERIM(1) ) + ( DKRRIM(2) * ERIM(2) ) + ( DKRRIM(3) * ERIM(3) )
! Closest point on a cylinder rim with respect to the center of a cylinder disk
URIM(1)     = RRIM(1) + ( DKRRIM_ERIM * ERIM(1) )
URIM(2)     = RRIM(2) + ( DKRRIM_ERIM * ERIM(2) )
URIM(3)     = RRIM(3) + ( DKRRIM_ERIM * ERIM(3) )
! Vector distance between a cylinder disk and the closest point on a cylinder rim
DKURIM(1)   = DK(1) - URIM(1)
DKURIM(2)   = DK(2) - URIM(2)
DKURIM(3)   = DK(3) - URIM(3)
! Magnitude of the vector distance between the closest point on the cylinder rim and the center of the cylinder disk
DKURIMSQ    = ( DKURIM(1) * DKURIM(1) ) + ( DKURIM(2) * DKURIM(2) ) + ( DKURIM(3) * DKURIM(3) )

! *********************************************************************************************** !
! Initial conditions                                                                              !
! *********************************************************************************************** !
! Sphere circumscribing the cylinder disks
IF( DKURIMSQ > DSQ ) THEN
  ! No overlap
  RETURN
! Center of the cylinder disk in the space zone that corresponds to the extension of the cylinder rim
ELSE IF( ( DKURIMSQ <= RSQ ) .AND. ( DABS( DKRRIM_ERIM ) > HALFLRIM ) ) THEN
  ! Disk-disk configuration (already tested)
  RETURN
! Center of the cylinder disk in the space zone that corresponds to the cylinder rim
ELSE IF( ( DKURIMSQ <= RSQ ) .AND. ( DABS( DKRRIM_ERIM ) <= HALFLRIM ) ) THEN
  OVERLAPDRIM = .TRUE.
  RETURN
END IF

! *********************************************************************************************** !
! Other configurations                                                                            !
! *********************************************************************************************** !

! Orientation of the cylinder disk along the x-direction
CALL ACTIVE_TRANSFORMATION( AXISX, QDISK, EXDISK )
! Orientation of the cylinder disk along the y-direction
CALL ACTIVE_TRANSFORMATION( AXISY, QDISK, EYDISK )
! Dot product of the orientation of the cylinder disk along x-direction and the orientation of the cylinder rim
EXDISK_ERIM = ( EXDISK(1) * ERIM(1) ) + ( EXDISK(2) * ERIM(2) ) + ( EXDISK(3) * ERIM(3) )
! Dot product of the orientation of the cylinder disk along y-direction and the orientation of the cylinder rim
EYDISK_ERIM = ( EYDISK(1) * ERIM(1) ) + ( EYDISK(2) * ERIM(2) ) + ( EYDISK(3) * ERIM(3) )
! Projection of the vector distance between a disk and a rim on the orientation of the cylinder disk along x-direction
DKRRIM_EXDISK = ( DKRRIM(1) * EXDISK(1) ) + ( DKRRIM(2) * EXDISK(2) ) + ( DKRRIM(3) * EXDISK(3) )
! Projection of the vector distance between a disk and a rim on the orientation of the cylinder disk along y-direction
DKRRIM_EYDISK = ( DKRRIM(1) * EYDISK(1) ) + ( DKRRIM(2) * EYDISK(2) ) + ( DKRRIM(3) * EYDISK(3) )

! *********************************************************************************************** !
! Preliminary configuration                                                                       !
! *********************************************************************************************** !
!  Calculate a point on the circumference of the cylinder disk that minimizes the objective       !
!  function for λ = L/2 and λ = -L/2 (extremum points on the cylinder rim).                       !
! *********************************************************************************************** !

! *********************************************************************************************** !
! Case 1 (λ = L/2)                                                                                !
! *********************************************************************************************** !
LAMBDA(1)    = HALFLRIM
OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
! Trigonometric relations in the rectangle triangle
COSPHI = ADJACENT_CAT / HYP
SINPHI = OPPOSITE_CAT / HYP
! Closest point on the circumference of the cylinder disk
PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
! Closest point on the axial axis of the cylinder rim
PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
T(1)  = PD(1) - RRIM(1)
T(2)  = PD(2) - RRIM(2)
T(3)  = PD(3) - RRIM(3)
! Magnitude of the vector T (squared)
TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
! Magnitude of vector T parallel to the cylinder rim (squared)
TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
TPARALSQ = TPARALSQ * TPARALSQ
! Magnitude of vector T orthogonal to the cylinder rim (squared)
TORTHOSQ = TSQ - TPARALSQ
! Overlap condition
IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
  OVERLAPDRIM = .TRUE.
  RETURN
END IF

! *********************************************************************************************** !
! Case 2 (λ = -L/2)                                                                               !
! *********************************************************************************************** !
LAMBDA(1)    = -HALFLRIM
OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
! Trigonometric relations in the rectangle triangle
COSPHI = ADJACENT_CAT / HYP
SINPHI = OPPOSITE_CAT / HYP
! Closest point on the circumference of the cylinder disk
PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
! Closest point on the axial axis of the cylinder rim
PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
T(1)  = PD(1) - RRIM(1)
T(2)  = PD(2) - RRIM(2)
T(3)  = PD(3) - RRIM(3)
! Magnitude of the vector T (squared)
TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
! Magnitude of vector T parallel to the cylinder rim (squared)
TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
TPARALSQ = TPARALSQ * TPARALSQ
! Magnitude of vector T orthogonal to the cylinder rim (squared)
TORTHOSQ = TSQ - TPARALSQ
! Overlap condition
IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
  OVERLAPDRIM = .TRUE.
  RETURN
END IF

! *********************************************************************************************** !
! SPECIAL CASE (Mainly for anisomorphic cylinders)                                                !
! *********************************************************************************************** !
! Check if the cylinder axis of particle i intersect any of the cylinder disks of particle j.     !
! *********************************************************************************************** !

! Point in the plane defined by the cylindrical disk
PK = DK + HALFDDISK * EXDISK

! Normal vector of the plane defined by the cylindrical disk
CALL ACTIVE_TRANSFORMATION( AXISZ, QDISK, EZDISK )

! Points in the line defined by the cylindrical axis
LA = RRIM + HALFLRIM * ERIM
LB = RRIM - HALFLRIM * ERIM

! Checking whether the cylindrical axis is orthogonal to the normal vector of the plane (avoid the indeterminate form 0/0)
IF( DABS( DOT_PRODUCT( EZDISK, ERIM ) ) >= 1.D-10 ) THEN

  ! Scale parameter of the line equation that demarks the intersection point between the line and the disk
  LAMBDA(1) = ( DOT_PRODUCT( EZDISK, PK ) - DOT_PRODUCT( EZDISK, LA ) ) / DOT_PRODUCT( EZDISK, (LB - LA) )

  ! Coordinates of the intersection point
  INTERSECTION = LA + LAMBDA(1) * (LB - LA)

  ! Distance between the intersection point and the center of the cylindrical disk (squared)
  PMINDISKSQ = DOT_PRODUCT( (DK - INTERSECTION), (DK - INTERSECTION) )

  ! Distance between the intersection point and the center of the cylindrical axis (squared)
  PMINCYLSQ = DOT_PRODUCT( (RRIM - INTERSECTION), (RRIM - INTERSECTION) )

  ! Checks whether the intersection point is within the limits of both the cylindrical disk and cylindrical axis
  IF( PMINDISKSQ <= (HALFDDISK * HALFDDISK) .AND. PMINCYLSQ <= (HALFLRIM * HALFLRIM) ) THEN
    OVERLAPDRIM = .TRUE.
    RETURN
  END IF

END IF

! *********************************************************************************************** !
! More complex configurations                                                                     !
! *********************************************************************************************** !
!  Calculate a point on the circumference of the cylinder disk and a point on the axial axis of   !
!  the cylinder rim that minimizes the objective function.                                        !
! *********************************************************************************************** !

! Constants
ALPHA = (EXDISK_ERIM * EXDISK_ERIM) + (EYDISK_ERIM * EYDISK_ERIM)
BETA  = (EXDISK_ERIM * DKRRIM_EXDISK) + (EYDISK_ERIM * DKRRIM_EYDISK)
GAMMA = (DKRRIM_EXDISK * DKRRIM_EXDISK) + (DKRRIM_EYDISK * DKRRIM_EYDISK)

! Particles perfectly aligned through their molecular axis (end-to-end configuration => α = 0)
IF( ALPHA > 0.D0 ) THEN

  ! Coefficients of the quartic function (objective function)
  COEFF_QUARTIC(1) = 1.D0
  COEFF_QUARTIC(2) = - (2.D0 / ALPHA) * (BETA + ALPHA * DKRRIM_ERIM)
  COEFF_QUARTIC(3) = (1.D0 / ALPHA) * ( (4.D0 * BETA * DKRRIM_ERIM) + (ALPHA * DKRRIM_ERIM * DKRRIM_ERIM) - &
  &                  (ALPHA * ALPHA * HALFDDISK * HALFDDISK) + GAMMA )
  COEFF_QUARTIC(4) = (2.D0 / ALPHA) * ( (ALPHA * BETA * HALFDDISK * HALFDDISK) - (GAMMA * DKRRIM_ERIM) - &
  &                  (BETA * DKRRIM_ERIM * DKRRIM_ERIM) )
  COEFF_QUARTIC(5) = (1.D0 / ALPHA) * ( (DKRRIM_ERIM * DKRRIM_ERIM * GAMMA) - (BETA * BETA * HALFDDISK * HALFDDISK) )

  ! Coefficients of the cubic function (derivative of the objective function with respect to λ)
  COEFF_CARDANO(1) = 4.D0 * COEFF_QUARTIC(1)
  COEFF_CARDANO(2) = 3.D0 * COEFF_QUARTIC(2)
  COEFF_CARDANO(3) = 2.D0 * COEFF_QUARTIC(3)
  COEFF_CARDANO(4) = COEFF_QUARTIC(4)

  ! *********************************************************************************************** !
  ! Cardano's solution for the cubic function (finding the points where the derivative is zero)     !
  ! *********************************************************************************************** !

  ! Cardano's coefficients
  CARDANO_Q = ( ( 3.D0 * COEFF_CARDANO(1) * COEFF_CARDANO(3) ) - ( COEFF_CARDANO(2) * COEFF_CARDANO(2) ) ) / &
  &           ( 9.D0 * COEFF_CARDANO(1) * COEFF_CARDANO(1) )
  CARDANO_R = ( ( 9.D0 * COEFF_CARDANO(1) * COEFF_CARDANO(2) * COEFF_CARDANO(3) ) - ( 27.D0 *COEFF_CARDANO(1) * COEFF_CARDANO(1) * &
  &           COEFF_CARDANO(4) ) - ( 2.D0 * COEFF_CARDANO(2) * COEFF_CARDANO(2) * COEFF_CARDANO(2) ) ) / &
  &           ( 54.D0 * COEFF_CARDANO(1) * COEFF_CARDANO(1) * COEFF_CARDANO(1) )

  ! Cardano's discriminant
  DISCRIMINANT = CARDANO_Q * CARDANO_Q * CARDANO_Q + CARDANO_R * CARDANO_R

  ! If Δ < 0, all roots are real and unequal (Trigonometric Solution)
  IF( DISCRIMINANT < 0.D0 ) THEN

    ! Trigonometric function
    ARGUMENT = CARDANO_R / DSQRT( - CARDANO_Q * CARDANO_Q * CARDANO_Q )
    IF( ARGUMENT <= -1.D0 ) THEN
      ARGUMENT = -1.D0
    ELSE IF( ARGUMENT >= 1.D0 ) THEN
      ARGUMENT = 1.D0
    END IF
    THETA = DACOS( ARGUMENT )

    ! Roots of the cubic equation
    ROOTS(1) = 2.D0 * DSQRT( - CARDANO_Q ) * DCOS( THETA / 3.D0 ) - COEFF_CARDANO(2) / ( 3.D0 * COEFF_CARDANO(1) )
    ROOTS(2) = 2.D0 * DSQRT( - CARDANO_Q ) * DCOS( (THETA / 3.D0) + (2.D0 * PI / 3.D0) ) - COEFF_CARDANO(2) / &
    &          ( 3.D0 * COEFF_CARDANO(1) )
    ROOTS(3) = 2.D0 * DSQRT( - CARDANO_Q ) * DCOS( (THETA / 3.D0) + (4.D0 * PI / 3.D0) ) - COEFF_CARDANO(2) / &
    &          ( 3.D0 * COEFF_CARDANO(1) )

    ! Sort roots in ascending order
    CPOINTS(1) = MINVAL( ROOTS )
    CPOINTS(3) = MAXVAL( ROOTS )
    DO I = 1, 3
      IF( ROOTS(I) > CPOINTS(1) .AND. ROOTS(I) < CPOINTS(3) ) THEN
        CPOINTS(2) = ROOTS(I)
      END IF
    END DO

  ! If Δ = 0, all roots are real and at least two are equal
  ELSE IF( DABS( DISCRIMINANT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN

    ! Cardano's coefficients
    CARDANO_S = CARDANO_R ** (1.D0 / 3.D0)
    CARDANO_T = CARDANO_S

    ! Roots of the cubic equation
    ROOTS(1) = CARDANO_S + CARDANO_T - COEFF_CARDANO(2) / ( 3.D0 * COEFF_CARDANO(1) )
    ROOTS(2) = - 0.5D0 * (CARDANO_S + CARDANO_T) - COEFF_CARDANO(2) / ( 3.D0 * COEFF_CARDANO(1) )
    ROOTS(3) = ROOTS(2)

    ! Sort roots in ascending order
    CPOINTS(1) = MINVAL( ROOTS )
    CPOINTS(2) = MAXVAL( ROOTS )
    CPOINTS(3) = MAXVAL( ROOTS )

  ! If Δ > 0, one roots is real and two are complex conjugates
  ELSE IF( DISCRIMINANT > 0.D0 ) THEN

    ! Cardano's coefficients
    CARDANO_S = ( CARDANO_R + DSQRT( (CARDANO_Q * CARDANO_Q * CARDANO_Q) + (CARDANO_R * CARDANO_R) ) ) ** (1.D0 / 3.D0)
    CARDANO_T = ( CARDANO_R - DSQRT( (CARDANO_Q * CARDANO_Q * CARDANO_Q) + (CARDANO_R * CARDANO_R) ) ) ** (1.D0 / 3.D0)

    ! Roots of the cubic equation
    ROOTS(1) = CARDANO_S + CARDANO_T - COEFF_CARDANO(2) / ( 3.D0 * COEFF_CARDANO(1) )
    ROOTS(2) = ROOTS(1) ! We are only interested in the real solution
    ROOTS(3) = ROOTS(1) ! We are only interested in the real solution

    ! Sort roots in ascending order
    CPOINTS(1) = ROOTS(1)
    CPOINTS(2) = ROOTS(2)
    CPOINTS(3) = ROOTS(3)

  END IF

  ! Number of intervals
  IF( DISCRIMINANT < 0.D0 ) THEN
    N_POINTS = 4
  ELSE IF( DABS( DISCRIMINANT - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
    N_POINTS = 3
  ELSE IF( DISCRIMINANT > 0.D0 ) THEN
    N_POINTS = 2
  END IF
  ! Single interval
  IF( MINVAL( CPOINTS ) >= HALFLRIM .OR. MAXVAL( CPOINTS ) <= - HALFLRIM ) THEN
    N_POINTS = 1
  END IF

  ! *********************************************************************************************** !
  ! Loops (Initial Guesses/Intervals)                                                               !
  ! *********************************************************************************************** !
  !  First loop  :     -L/2      <= λ0 <= CUBIC_ROOT(1)                                             !
  !  Second loop : CUBIC_ROOT(1) <= λ0 <= CUBIC_ROOT(2)                                             !
  !  Third loop  : CUBIC_ROOT(2) <= λ0 <= CUBIC_ROOT(3)                                             !
  !  Fourth loop : CUBIC_ROOT(3) <= λ0 <=     L/2                                                   !
  ! *********************************************************************************************** !
  INTERVAL_LOOP: DO INTERVAL = 1, N_POINTS

    ! Initial parameters of the Newton-Raphson method
    DAMPING       = 1.D0    ! Damping factor
    COUNTER       = 0       ! Iteration counter
    TOLERANCE_NR  = 1D-10   ! Numerical tolerance

    ! Initial guess, λ0 (Newton-Raphson)
    IF( N_POINTS == 4 ) THEN
      IF( INTERVAL == 1 ) THEN
        LAMBDA(1) = 0.5D0 * ( - HALFLRIM + CPOINTS(1) )
      ELSE IF( INTERVAL == 2 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(1) + CPOINTS(2) )
      ELSE IF( INTERVAL == 3 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(2) + CPOINTS(3) )
      ELSE IF( INTERVAL == 4 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(3) + HALFLRIM )
      END IF
    ELSE IF( N_POINTS == 3 ) THEN
      IF( INTERVAL == 1 ) THEN
        LAMBDA(1) = 0.5D0 * ( - HALFLRIM + CPOINTS(1) )
      ELSE IF( INTERVAL == 2 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(1) + CPOINTS(2) )
      ELSE IF( INTERVAL == 3 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(2) + HALFLRIM )
      END IF
    ELSE IF( N_POINTS == 2 ) THEN
      IF( INTERVAL == 1 ) THEN
        LAMBDA(1) = 0.5D0 * ( - HALFLRIM + CPOINTS(1) )
      ELSE IF( INTERVAL == 2 ) THEN
        LAMBDA(1) = 0.5D0 * ( CPOINTS(1) + HALFLRIM )
      END IF
    ELSE IF( N_POINTS == 1 .AND. MINVAL( CPOINTS ) >= HALFLRIM ) THEN
      LAMBDA(1) = 0.5D0 * ( - HALFLRIM + MINVAL( CPOINTS ) )
    ELSE IF( N_POINTS == 1 .AND. MAXVAL( CPOINTS ) <= - HALFLRIM ) THEN
      LAMBDA(1) = 0.5D0 * ( HALFLRIM + MAXVAL( CPOINTS ) )
    END IF

    ! Objective function
    F(1) = ( COEFF_QUARTIC(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) ) + ( COEFF_QUARTIC(2) * LAMBDA(1) * LAMBDA(1) * &
    &      LAMBDA(1) ) + ( COEFF_QUARTIC(3) * LAMBDA(1) * LAMBDA(1) ) + ( COEFF_QUARTIC(4) * LAMBDA(1) ) + COEFF_QUARTIC(5)

    ! Derivative of the objective function with respect to λ
    DF(1) = ( 4.D0 * COEFF_QUARTIC(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) ) + ( 3.D0 * COEFF_QUARTIC(2) * LAMBDA(1) * LAMBDA(1) ) &
    &       + ( 2.D0 * COEFF_QUARTIC(3) * LAMBDA(1) ) + COEFF_QUARTIC(4)

    ! Newton-Raphson method
    DO WHILE ( DABS( F(1) ) >= TOLERANCE_NR .AND. .NOT. BISECT_METHOD)
      ! Increment factor
      FACTOR = F(1) / DF(1)
      ! New λ point
      LAMBDA(1) = LAMBDA(1) - DAMPING * FACTOR
      ! Recalculate the function with a new λ
      F(1) = ( COEFF_QUARTIC(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) ) + ( COEFF_QUARTIC(2) * LAMBDA(1) * LAMBDA(1) * &
      &         LAMBDA(1) ) + ( COEFF_QUARTIC(3) * LAMBDA(1) * LAMBDA(1) ) + ( COEFF_QUARTIC(4) * LAMBDA(1) ) + COEFF_QUARTIC(5)
      ! Recalculate the derivative of the objective function with a new λ
      DF(1) = ( 4.D0 * COEFF_QUARTIC(1) * LAMBDA(1) * LAMBDA(1) * LAMBDA(1) ) + ( 3.D0 * COEFF_QUARTIC(2) * LAMBDA(1) * & 
      &       LAMBDA(1) ) + ( 2.D0 * COEFF_QUARTIC(3) * LAMBDA(1) ) + COEFF_QUARTIC(4)
      ! Iteration
      COUNTER = COUNTER + 1
      ! Stop criterion (bissection method)
      IF( COUNTER > 25 ) THEN
        BISECT_METHOD = .TRUE.
      END IF
    END DO

    ! Initialize bissection method if the Newton-Raphson method doesn't converge
    IF( BISECT_METHOD ) THEN
      ! Initial parameters of the bissection method
      COUNTERB    = 0      ! Iteration counter
      TOLERANCE_B = 1.D-10 ! Numerical tolerance
      ! Initial guess, λ0 (Newton-Raphson)
      IF( N_POINTS == 4 ) THEN
        IF( INTERVAL == 1 ) THEN
          LAMBDA(2) = - HALFLRIM
          LAMBDA(3) = CPOINTS(1)
        ELSE IF( INTERVAL == 2 ) THEN
          LAMBDA(2) = CPOINTS(1)
          LAMBDA(3) = CPOINTS(2)
        ELSE IF( INTERVAL == 3 ) THEN
          LAMBDA(2) = CPOINTS(2)
          LAMBDA(3) = CPOINTS(3)
        ELSE IF( INTERVAL == 4 ) THEN
          LAMBDA(2) = CPOINTS(3)
          LAMBDA(3) = HALFLRIM
        END IF
      ELSE IF( N_POINTS == 3 ) THEN
        IF( INTERVAL == 1 ) THEN
          LAMBDA(2) = - HALFLRIM
          LAMBDA(3) = CPOINTS(1)
        ELSE IF( INTERVAL == 2 ) THEN
          LAMBDA(2) = CPOINTS(1)
          LAMBDA(3) = CPOINTS(2)
        ELSE IF( INTERVAL == 3 ) THEN
          LAMBDA(2) = CPOINTS(2)
          LAMBDA(3) = HALFLRIM
        END IF
      ELSE IF( N_POINTS == 2 ) THEN
        IF( INTERVAL == 1 ) THEN
          LAMBDA(2) = - HALFLRIM
          LAMBDA(3) = CPOINTS(1)
        ELSE IF( INTERVAL == 2 ) THEN
          LAMBDA(2) = CPOINTS(1)
          LAMBDA(3) = HALFLRIM
        END IF
      ELSE IF( N_POINTS == 1 .AND. MINVAL( CPOINTS ) >= HALFLRIM ) THEN
        LAMBDA(2) = - HALFLRIM
        LAMBDA(3) = MINVAL( CPOINTS )
      ELSE IF( N_POINTS == 1 .AND. MAXVAL( CPOINTS ) <= - HALFLRIM ) THEN
        LAMBDA(2) = MAXVAL( CPOINTS )
        LAMBDA(3) = HALFLRIM
      END IF
      ! Calculate function at λ0
      F(2)  = ( COEFF_QUARTIC(1) * LAMBDA(2) * LAMBDA(2) * LAMBDA(2) * LAMBDA(2) ) + ( COEFF_QUARTIC(2) * LAMBDA(2) * LAMBDA(2) * &
      &       LAMBDA(2) ) + ( COEFF_QUARTIC(3) * LAMBDA(2) * LAMBDA(2) ) + ( COEFF_QUARTIC(4) * LAMBDA(2) ) + COEFF_QUARTIC(5)
      ! Calculate function at λ0
      F(3)  = ( COEFF_QUARTIC(1) * LAMBDA(3) * LAMBDA(3) * LAMBDA(3) * LAMBDA(3) ) + ( COEFF_QUARTIC(2) * LAMBDA(3) * LAMBDA(3) * &
      &       LAMBDA(3) ) + ( COEFF_QUARTIC(3) * LAMBDA(3) * LAMBDA(3) ) + ( COEFF_QUARTIC(4) * LAMBDA(3) ) + COEFF_QUARTIC(5)
      ! Avoid bissection method if initial guess gives already the minimum value
      IF( DABS( F(2) ) < TOLERANCE_B .AND. DABS( F(3) ) >= TOLERANCE_B ) THEN
        LAMBDA(1)    = LAMBDA(2)
        OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
        ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
        HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
        ! Trigonometric relations in the rectangle triangle
        COSPHI = ADJACENT_CAT / HYP
        SINPHI = OPPOSITE_CAT / HYP
        ! Closest point on the circumference of the cylinder disk
        PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
        PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
        PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
        ! Closest point on the axial axis of the cylinder rim
        PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
        PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
        PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
        ! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
        T(1)  = PD(1) - RRIM(1)
        T(2)  = PD(2) - RRIM(2)
        T(3)  = PD(3) - RRIM(3)
        ! Magnitude of the vector T (squared)
        TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
        ! Magnitude of vector T parallel to the cylinder rim (squared)
        TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
        TPARALSQ = TPARALSQ * TPARALSQ
        ! Magnitude of vector T orthogonal to the cylinder rim (squared)
        TORTHOSQ = TSQ - TPARALSQ
        ! Overlap condition
        IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
          OVERLAPDRIM = .TRUE.
          RETURN
        ELSE
          OVERLAPDRIM = .FALSE.
        END IF
        CYCLE INTERVAL_LOOP
      ELSE IF( DABS( F(2) ) >= TOLERANCE_B .AND. DABS( F(3) ) < TOLERANCE_B ) THEN
        LAMBDA(1)    = LAMBDA(3)
        OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
        ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
        HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
        ! Trigonometric relations in the rectangle triangle
        COSPHI = ADJACENT_CAT / HYP
        SINPHI = OPPOSITE_CAT / HYP
        ! Closest point on the circumference of the cylinder disk
        PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
        PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
        PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
        ! Closest point on the axial axis of the cylinder rim
        PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
        PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
        PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
        ! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
        T(1)  = PD(1) - RRIM(1)
        T(2)  = PD(2) - RRIM(2)
        T(3)  = PD(3) - RRIM(3)
        ! Magnitude of the vector T (squared)
        TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
        ! Magnitude of vector T parallel to the cylinder rim (squared)
        TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
        TPARALSQ = TPARALSQ * TPARALSQ
        ! Magnitude of vector T orthogonal to the cylinder rim (squared)
        TORTHOSQ = TSQ - TPARALSQ
        ! Overlap condition
        IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
          OVERLAPDRIM = .TRUE.
          RETURN
        ELSE
          OVERLAPDRIM = .FALSE.
        END IF
        CYCLE INTERVAL_LOOP
      ELSE IF( DABS( F(2) ) < TOLERANCE_B .AND. DABS( F(3) ) < TOLERANCE_B ) THEN
        ! First guess
        LAMBDA(1)    = LAMBDA(2)
        OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
        ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
        HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
        ! Trigonometric relations in the rectangle triangle
        COSPHI = ADJACENT_CAT / HYP
        SINPHI = OPPOSITE_CAT / HYP
        ! Closest point on the circumference of the cylinder disk
        PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
        PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
        PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
        ! Closest point on the axial axis of the cylinder rim
        PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
        PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
        PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
        ! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
        T(1)  = PD(1) - RRIM(1)
        T(2)  = PD(2) - RRIM(2)
        T(3)  = PD(3) - RRIM(3)
        ! Magnitude of the vector T (squared)
        TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
        ! Magnitude of vector T parallel to the cylinder rim (squared)
        TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
        TPARALSQ = TPARALSQ * TPARALSQ
        ! Magnitude of vector T orthogonal to the cylinder rim (squared)
        TORTHOSQ = TSQ - TPARALSQ
        ! Overlap condition
        IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
          OVERLAPDRIM = .TRUE.
          RETURN
        ELSE
          OVERLAPDRIM = .FALSE.
        END IF
        ! Second guess
        LAMBDA(1)    = LAMBDA(3)
        OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
        ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
        HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
        ! Trigonometric relations in the rectangle triangle
        COSPHI = ADJACENT_CAT / HYP
        SINPHI = OPPOSITE_CAT / HYP
        ! Closest point on the circumference of the cylinder disk
        PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
        PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
        PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
        ! Closest point on the axial axis of the cylinder rim
        PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
        PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
        PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
        ! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
        T(1)  = PD(1) - RRIM(1)
        T(2)  = PD(2) - RRIM(2)
        T(3)  = PD(3) - RRIM(3)
        ! Magnitude of the vector T (squared)
        TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
        ! Magnitude of vector T parallel to the cylinder rim (squared)
        TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
        TPARALSQ = TPARALSQ * TPARALSQ
        ! Magnitude of vector T orthogonal to the cylinder rim (squared)
        TORTHOSQ = TSQ - TPARALSQ
        ! Overlap condition
        IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
          OVERLAPDRIM = .TRUE.
          RETURN
        ELSE
          OVERLAPDRIM = .FALSE.
        END IF
        CYCLE INTERVAL_LOOP
      END IF
      ! Bissection condition
      IF( ( F(2) * F(3) ) < 0.D0 ) THEN
        ! Midpoint λ
        LAMBDAI = 0.5D0 * ( LAMBDA(2) + LAMBDA(3) )
        ! Calculate function at the midpoint
        FI  = ( COEFF_QUARTIC(1) * LAMBDAI * LAMBDAI * LAMBDAI * LAMBDAI ) + ( COEFF_QUARTIC(2) * LAMBDAI * LAMBDAI * LAMBDAI ) + &
        &     ( COEFF_QUARTIC(3) * LAMBDAI * LAMBDAI ) + ( COEFF_QUARTIC(4) * LAMBDAI ) + COEFF_QUARTIC(5)
        ! Bissection method
        DO WHILE( DABS(FI) >= TOLERANCE_B )
          ! Bissection criterion
          IF( ( F(2) * FI ) > 0.D0 )THEN
            LAMBDA(2) = LAMBDAI
            F(2)      = FI
          ELSE
            LAMBDA(3) = LAMBDAI
            F(3)      = FI
          ENDIF
          ! New midpoint λ
          LAMBDAI = 0.5D0 * ( LAMBDA(2) + LAMBDA(3) )
          ! Calculate function at the new midpoint
          FI  = ( COEFF_QUARTIC(1) * LAMBDAI * LAMBDAI * LAMBDAI * LAMBDAI ) + ( COEFF_QUARTIC(2) * LAMBDAI * LAMBDAI * LAMBDAI ) &
          &     + ( COEFF_QUARTIC(3) * LAMBDAI * LAMBDAI ) + ( COEFF_QUARTIC(4) * LAMBDAI ) + COEFF_QUARTIC(5)
          ! Iteration
          COUNTERB = COUNTERB + 1
          ! Stop condition
          IF( COUNTERB > 50 ) THEN
            ! Bissection method will not converge (test other interval)
            CYCLE INTERVAL_LOOP
          END IF
        END DO
        ! Root that minimizes objective function
        LAMBDA(1) = LAMBDAI
      ! Bissection method cannot be applied
      ELSE
        ! Bissection method will not converge
        CYCLE INTERVAL_LOOP
      END IF
    END IF

    ! Point that minimizes the objective function is outside the cylinder rim
    IF( LAMBDA(1) < -HALFLRIM .OR. LAMBDA(1) > HALFLRIM) THEN
      CYCLE INTERVAL_LOOP
    END IF

    ! ********************************************************************************************* !
    ! Case 3 (λ = parameter that minimizes the objective function)                                  !
    ! ********************************************************************************************* !
    !   λ = λroot [NEWTON-RAPHSON or BISSECTION]                                                    !
    ! ********************************************************************************************* !
    OPPOSITE_CAT = LAMBDA(1) * EYDISK_ERIM - DKRRIM_EYDISK 
    ADJACENT_CAT = LAMBDA(1) * EXDISK_ERIM - DKRRIM_EXDISK 
    HYP          = DSQRT( (OPPOSITE_CAT * OPPOSITE_CAT) + (ADJACENT_CAT * ADJACENT_CAT) )
    ! Trigonometric relations in the rectangle triangle
    COSPHI = ADJACENT_CAT / HYP
    SINPHI = OPPOSITE_CAT / HYP
    ! Closest point on the circumference of the cylinder disk
    PD(1) = DK(1) + ( HALFDDISK * COSPHI * EXDISK(1) ) + ( HALFDDISK * SINPHI * EYDISK(1) )
    PD(2) = DK(2) + ( HALFDDISK * COSPHI * EXDISK(2) ) + ( HALFDDISK * SINPHI * EYDISK(2) )
    PD(3) = DK(3) + ( HALFDDISK * COSPHI * EXDISK(3) ) + ( HALFDDISK * SINPHI * EYDISK(3) )
    ! Closest point on the axial axis of the cylinder rim
    PC(1) = RRIM(1) + ( LAMBDA(1) * ERIM(1) )
    PC(2) = RRIM(2) + ( LAMBDA(1) * ERIM(2) )
    PC(3) = RRIM(3) + ( LAMBDA(1) * ERIM(3) )
    ! Vector distance between the closest point on the circumference of a cylinder disk and the center of a cylinder rim
    T(1)  = PD(1) - RRIM(1)
    T(2)  = PD(2) - RRIM(2)
    T(3)  = PD(3) - RRIM(3)
    ! Magnitude of the vector T (squared)
    TSQ      = ( T(1) * T(1) ) + ( T(2) * T(2) ) + ( T(3) * T(3) )
    ! Magnitude of vector T parallel to the cylinder rim (squared)
    TPARALSQ = ( ERIM(1) * T(1) ) + ( ERIM(2) * T(2) ) + ( ERIM(3) * T(3) )
    TPARALSQ = TPARALSQ * TPARALSQ
    ! Magnitude of vector T orthogonal to the cylinder rim (squared)
    TORTHOSQ = TSQ - TPARALSQ
    ! Overlap condition
    IF( TPARALSQ <= (HALFLRIM * HALFLRIM) .AND. TORTHOSQ <= (HALFDRIM * HALFDRIM) ) THEN
      OVERLAPDRIM = .TRUE.
      RETURN
    ELSE
      OVERLAPDRIM = .FALSE.
    END IF

  END DO INTERVAL_LOOP

END IF

! *********************************************************************************************** !
! -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  IMPORTANT NOTES  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* !
! *********************************************************************************************** !
!                                                                                                 !
! IF THE ALGORITHM REACHES THIS POINT,                                                            !
!                                                                                                 !
! (1) it is because the root obtained from the numerical methods returns a point in the           !
! circumference of the cylinder disk (j or i) that is outside the cylinder rim (i or j) but       !
! within the limits of its length. RESULT = NO OVERLAP                                            !
!                                                                                                 !
! (2) it is because no root was found in the considered intervals (-L/2 <= λ0 <= L/2), which      !
! means the root that minimizes the objective function is beyond the limits of the length of      !
! the cylinder rim (i or j). This situation will never lead to molecular overlaps except when     !
! the disks of both cylinder intersect each other. However, this situation has already been       !
! tested in the disk-disk configuration. RESULT = NO OVERLAP                                      !
!                                                                                                 !
! *********************************************************************************************** !

! No overlaps
OVERLAPDRIM = .FALSE.

RETURN

END SUBROUTINE DISK_RIM

! *********************************************************************************************** !
!                                    Outer Product of Vectors                                     !
! *********************************************************************************************** !
SUBROUTINE OUTER_PRODUCT( EI, EJ, EIJ )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 ), DIMENSION( 3 )    :: EI, EJ ! Vectors
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: EIJ    ! Cross product of vectors

! First row
EIJ(1,1)  = EI(1) * EJ(1)
EIJ(1,2)  = EI(1) * EJ(2)
EIJ(1,3)  = EI(1) * EJ(3)

! Second row
EIJ(2,1)  = EI(2) * EJ(1)
EIJ(2,2)  = EI(2) * EJ(2)
EIJ(2,3)  = EI(2) * EJ(3)

! Third row
EIJ(3,1)  = EI(3) * EJ(1)
EIJ(3,2)  = EI(3) * EJ(2)
EIJ(3,3)  = EI(3) * EJ(3)

RETURN

END SUBROUTINE OUTER_PRODUCT

! *********************************************************************************************** !
!                                 Inverse Matrix using Cofactors                                  !
! *********************************************************************************************** !
SUBROUTINE INVERSE( A, B )

! Uses one module: global variables
USE GLOBALVAR

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( KIND= REAL64 )                    :: DET, RDET ! Determinant
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: A         ! Matrix
REAL( KIND= REAL64 ), DIMENSION( 3, 3 ) :: B         ! Matrix (inverse)

! *********************************************************************************************** !
! Tranpose matrix of the matrix of cofactors of matrix A                                          !
! *********************************************************************************************** !
B(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
B(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
B(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
B(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
B(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
B(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)
B(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
B(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
B(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)

! *********************************************************************************************** !
! Determinant of matrix A                                                                         !
! *********************************************************************************************** !
DET  = A(1,1) * B(1,1) + A(2,1) * B(1,2) + A(3,1) * B(1,3)
RDET = 0.D0

IF( DABS( DET ) > 0.D0 ) THEN
  RDET = 1.D0 / DET
END IF

! *********************************************************************************************** !
! Inverse of matrix A                                                                             !
! *********************************************************************************************** !
B(:,:) = RDET * B(:,:)

RETURN

END SUBROUTINE INVERSE
