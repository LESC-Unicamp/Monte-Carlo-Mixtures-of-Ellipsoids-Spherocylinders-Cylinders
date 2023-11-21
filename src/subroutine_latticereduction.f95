! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!    to apply a lattice reduction technique to orthogonalize the shape of the simulation box.     ! 
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
! Main References:               D. Gottwald, G. Kahl, C. N. Likos                                !
!                                J. Chem. Phys. 122, 204503 (2005)                                !
!                                     DOI: 10.1063/1.1901585                                      !
!                             --------------------------------------                              !
!                           A. K. Lenstra, H. W. Lenstra Jr., L. Lovász                           !
!                                 Math. Ann. 261, 515-534 (1982)                                  !
!                                     DOI: 10.1007/BF01457454                                     !
!                             --------------------------------------                              !
!                  J. de Graaf, L. Filion, M. Marechal, R. von Roij, M. Dijkstra                  !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!           This subroutine calculates the lattice reduction and reshapes the box size            !
!            See Graaf et al., J. Chem. Phys. 137, 214101 (2012) for more information.            !
! *********************************************************************************************** !
SUBROUTINE LatticeReduction( CurrentBoxLength, CurrentBoxDistortion, LatticeReductionLogical )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER :: nDimension = 3 ! Lattice dimension

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: CurrentBoxDistortion                           ! Distortion parameter (surface-to-volume ratio)
REAL( Kind= Real64 )                 :: AreaXY                                         ! Magnitude of the cross product between vectors vx and vy of the simulation box (area of plane xy)
REAL( Kind= Real64 )                 :: AreaXZ                                         ! Magnitude of the cross product between vectors vx and vz of the simulation box (area of plane xz)
REAL( Kind= Real64 )                 :: AreaYZ                                         ! Magnitude of the cross product between vectors vy and vz of the simulation box (area of plane yz)
REAL( Kind= Real64 )                 :: ParallelepipedVolume                           ! Dot product of vector vx and the cross product of vectors vy and vz (volume)
REAL( Kind= Real64 )                 :: SurfaceArea                                    ! Surface area
REAL( Kind= Real64 )                 :: VectorLengthX, VectorLengthY, VectorLengthZ    ! Magnitude of the box vectors
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BoxVectorX, BoxVectorY, BoxVectorZ             ! Vectors vx, vy, and vz of the simulation box
REAL( Kind= Real64 ), DIMENSION( 3 ) :: CrossProductXY, CrossProductXZ, CrossProductYZ ! Cross product of box vectors
REAL( Kind= Real64 ), DIMENSION( 9 ) :: CurrentBoxLength                               ! Box vector matrix

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: LatticeReductionLogical ! Checks whether a lattice reduction is necessary based on the distortion parameter value

! Box vectors
BoxVectorX = CurrentBoxLength(1:3) ! Vector x
BoxVectorY = CurrentBoxLength(4:6) ! Vector y
BoxVectorZ = CurrentBoxLength(7:9) ! Vector z

! Vector moduli
VectorLengthX = ( BoxVectorX(1) * BoxVectorX(1) ) + ( BoxVectorX(2) * BoxVectorX(2) ) + ( BoxVectorX(3) * BoxVectorX(3) )
VectorLengthX = DSQRT( VectorLengthX )
VectorLengthY = ( BoxVectorY(1) * BoxVectorY(1) ) + ( BoxVectorY(2) * BoxVectorY(2) ) + ( BoxVectorY(3) * BoxVectorY(3) )
VectorLengthY = DSQRT( VectorLengthY )
VectorLengthZ = ( BoxVectorZ(1) * BoxVectorZ(1) ) + ( BoxVectorZ(2) * BoxVectorZ(2) ) + ( BoxVectorZ(3) * BoxVectorZ(3) )
VectorLengthZ = DSQRT( VectorLengthZ )

! Cross product of vx and vy
CrossProductXY(1) = ( BoxVectorX(2) * BoxVectorY(3) ) - ( BoxVectorX(3) * BoxVectorY(2) )
CrossProductXY(2) = ( BoxVectorX(3) * BoxVectorY(1) ) - ( BoxVectorX(1) * BoxVectorY(3) )
CrossProductXY(3) = ( BoxVectorX(1) * BoxVectorY(2) ) - ( BoxVectorX(2) * BoxVectorY(1) )

! Cross product of vx and vz
CrossProductXZ(1) = ( BoxVectorX(2) * BoxVectorZ(3) ) - ( BoxVectorX(3) * BoxVectorZ(2) )
CrossProductXZ(2) = ( BoxVectorX(3) * BoxVectorZ(1) ) - ( BoxVectorX(1) * BoxVectorZ(3) )
CrossProductXZ(3) = ( BoxVectorX(1) * BoxVectorZ(2) ) - ( BoxVectorX(2) * BoxVectorZ(1) )

! Cross product of vy and vz
CrossProductYZ(1) = ( BoxVectorY(2) * BoxVectorZ(3) ) - ( BoxVectorY(3) * BoxVectorZ(2) )
CrossProductYZ(2) = ( BoxVectorY(3) * BoxVectorZ(1) ) - ( BoxVectorY(1) * BoxVectorZ(3) )
CrossProductYZ(3) = ( BoxVectorY(1) * BoxVectorZ(2) ) - ( BoxVectorY(2) * BoxVectorZ(1) )

! Dot product of vx and the cross product of vy and vz
ParallelepipedVolume = DOT_PRODUCT( BoxVectorX, CrossProductYZ )

! Vector moduli
AreaXY = DOT_PRODUCT( CrossProductXY, CrossProductXY )
AreaXY = DSQRT( AreaXY )
AreaXZ = DOT_PRODUCT( CrossProductXZ, CrossProductXZ )
AreaXZ = DSQRT( AreaXZ )
AreaYZ = DOT_PRODUCT( CrossProductYZ, CrossProductYZ )
AreaYZ = DSQRT( AreaYZ )

! Lattice distortion (relation among perimeter × surface area ÷ volume)
CurrentBoxDistortion = (VectorLengthX + VectorLengthY + VectorLengthZ) * (AreaXY + AreaXZ + AreaYZ) / ParallelepipedVolume
CurrentBoxDistortion = CurrentBoxDistortion / 9.D0

! Surface area
SurfaceArea = 2.D0 * AreaXY + 2.D0 * AreaXZ + 2.D0 * AreaYZ

! Verification
IF( CurrentBoxDistortion > MaxBoxDistortion ) THEN
  LatticeReductionLogical = .TRUE.
  ! Lattice reduction method
  IF( LatticeReductionTypeLogical(1) ) THEN
    ! Gottwald method
    CALL GottwaldMethod( BoxVectorX, BoxVectorY, BoxVectorZ, SurfaceArea )
    ! Update box vectors
    CurrentBoxLength(1:3) = BoxVectorX
    CurrentBoxLength(4:6) = BoxVectorY
    CurrentBoxLength(7:9) = BoxVectorZ
  ELSE IF( LatticeReductionTypeLogical(2) ) THEN
    ! Lenstra-Lenstra-Lovász method
    CALL LenstraLenstraLovaszMethod( CurrentBoxLength, nDimension )
  END IF
ELSE
  LatticeReductionLogical = .FALSE.
  RETURN
END IF

RETURN

END SUBROUTINE LatticeReduction

! *********************************************************************************************** !
!          This subroutine applies the Gottwald method to orthogonalize the box vectors           !
!           See Gottwald et al., J. Chem. Phys. 122, 204503 (2005) for more information           !
! *********************************************************************************************** !
SUBROUTINE GottwaldMethod( BoxVectorX, BoxVectorY, BoxVectorZ, MinSurfaceArea )

USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCombination      ! Lattice combination
INTEGER( Kind= Int64 ) :: MinSurfaceAreaLoc ! Array location of minimum surface area

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                     :: AreaXY                                         ! Magnitude of the cross product between vectors vx and vy of the simulation box (area of plane xy)
REAL( Kind= Real64 )                     :: AreaXZ                                         ! Magnitude of the cross product between vectors vx and vz of the simulation box (area of plane xz)
REAL( Kind= Real64 )                     :: AreaYZ                                         ! Magnitude of the cross product between vectors vy and vz of the simulation box (area of plane yz)
REAL( Kind= Real64 )                     :: MinSurfaceArea                                 ! Minimum surface area
REAL( Kind= Real64 ), DIMENSION( 3 )     :: BoxVectorX, BoxVectorY, BoxVectorZ             ! Vectors vx, vy, and vz of the simulation box
REAL( Kind= Real64 ), DIMENSION( 3 )     :: NewBoxVectorX, NewBoxVectorY, NewBoxVectorZ    ! New vectors vx, vy, and vz of the simulation box
REAL( Kind= Real64 ), DIMENSION( 3 )     :: CrossProductXY, CrossProductXZ, CrossProductYZ ! Cross product of box vectors
REAL( Kind= Real64 ), DIMENSION( 12 )    :: SurfaceAreaCombination                         ! Surface area of every lattice combination
REAL( Kind= Real64 ), DIMENSION( 12, 9 ) :: VectorCombination                              ! Lattice combination

! Gottwald orthogonalization process
OrthogonalizationProcess: DO
  ! Lattice combination
  DO iCombination = 1, 12
    CALL LatticeCombination( iCombination, BoxVectorX, BoxVectorY, BoxVectorZ, VectorCombination(iCombination,:) )
    ! New lattice vectors
    NewBoxVectorX = VectorCombination(iCombination,1:3)
    NewBoxVectorY = VectorCombination(iCombination,4:6)
    NewBoxVectorZ = VectorCombination(iCombination,7:9)
    ! Cross product of vx and vy
    CrossProductXY(1) = ( NewBoxVectorX(2) * NewBoxVectorY(3) ) - ( NewBoxVectorX(3) * NewBoxVectorY(2) )
    CrossProductXY(2) = ( NewBoxVectorX(3) * NewBoxVectorY(1) ) - ( NewBoxVectorX(1) * NewBoxVectorY(3) )
    CrossProductXY(3) = ( NewBoxVectorX(1) * NewBoxVectorY(2) ) - ( NewBoxVectorX(2) * NewBoxVectorY(1) )
    ! Cross product of vx and vz
    CrossProductXZ(1) = ( NewBoxVectorX(2) * NewBoxVectorZ(3) ) - ( NewBoxVectorX(3) * NewBoxVectorZ(2) )
    CrossProductXZ(2) = ( NewBoxVectorX(3) * NewBoxVectorZ(1) ) - ( NewBoxVectorX(1) * NewBoxVectorZ(3) )
    CrossProductXZ(3) = ( NewBoxVectorX(1) * NewBoxVectorZ(2) ) - ( NewBoxVectorX(2) * NewBoxVectorZ(1) )
    ! Cross product of vy and vz
    CrossProductYZ(1) = ( NewBoxVectorY(2) * NewBoxVectorZ(3) ) - ( NewBoxVectorY(3) * NewBoxVectorZ(2) )
    CrossProductYZ(2) = ( NewBoxVectorY(3) * NewBoxVectorZ(1) ) - ( NewBoxVectorY(1) * NewBoxVectorZ(3) )
    CrossProductYZ(3) = ( NewBoxVectorY(1) * NewBoxVectorZ(2) ) - ( NewBoxVectorY(2) * NewBoxVectorZ(1) )
    ! Vector moduli
    AreaXY = DOT_PRODUCT( CrossProductXY, CrossProductXY )
    AreaXY = DSQRT( AreaXY )
    AreaXZ = DOT_PRODUCT( CrossProductXZ, CrossProductXZ )
    AreaXZ = DSQRT( AreaXZ )
    AreaYZ = DOT_PRODUCT( CrossProductYZ, CrossProductYZ )
    AreaYZ = DSQRT( AreaYZ )
    ! Surface area
    SurfaceAreaCombination(iCombination) = 2.D0 * AreaXY + 2.D0 * AreaXZ + 2.D0 * AreaYZ
  END DO
  ! Minimum surface area
  IF( MINVAL( SurfaceAreaCombination ) <= MinSurfaceArea ) THEN
    MinSurfaceAreaLoc = MINLOC( SurfaceAreaCombination, Dim= 1 )
    MinSurfaceArea    = MINVAL( SurfaceAreaCombination )
    ! Update box vectors
    BoxVectorX = VectorCombination(MinSurfaceAreaLoc,1:3)
    BoxVectorY = VectorCombination(MinSurfaceAreaLoc,4:6)
    BoxVectorZ = VectorCombination(MinSurfaceAreaLoc,7:9)
  ELSE
    EXIT OrthogonalizationProcess
  END IF
END DO OrthogonalizationProcess

RETURN

END SUBROUTINE GottwaldMethod

! *********************************************************************************************** !
!                       This subroutine calculates the lattice combinations                       !
! *********************************************************************************************** !
SUBROUTINE LatticeCombination( iCombination, BoxVectorX, BoxVectorY, BoxVectorZ, VectorCombination )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCombination ! Lattice combination

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 3 ) :: BoxVectorX, BoxVectorY, BoxVectorZ ! Vectors vx, vy, and vz of the simulation box
REAL( Kind= Real64 ), DIMENSION( 9 ) :: VectorCombination                  ! Lattice combination

! Lattice combination (unchanged determinant)
IF( iCombination == 1 ) THEN
  VectorCombination(1:3) = BoxVectorX + BoxVectorY
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 2 ) THEN
  VectorCombination(1:3) = BoxVectorX - BoxVectorY
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 3 ) THEN
  VectorCombination(1:3) = BoxVectorX + BoxVectorZ
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 4 ) THEN
  VectorCombination(1:3) = BoxVectorX - BoxVectorZ
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 5 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY + BoxVectorX
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 6 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY - BoxVectorX
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 7 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY + BoxVectorZ
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 8 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY - BoxVectorZ
  VectorCombination(7:9) = BoxVectorZ
ELSE IF( iCombination == 9 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ + BoxVectorX
ELSE IF( iCombination == 10 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ - BoxVectorX
ELSE IF( iCombination == 11 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ + BoxVectorY
ELSE IF( iCombination == 12 ) THEN
  VectorCombination(1:3) = BoxVectorX
  VectorCombination(4:6) = BoxVectorY
  VectorCombination(7:9) = BoxVectorZ - BoxVectorY
END IF

RETURN

END SUBROUTINE LatticeCombination

! *********************************************************************************************** !
!  This subroutine applies the Lenstra–Lenstra–Lovász algorithm to orthogonalize the box vectors  !
!       See Lenstra, Lenstra & Lovász, Math. Ann. 261, 515-534 (1982) for more information        !
! *********************************************************************************************** !
SUBROUTINE LenstraLenstraLovaszMethod( BoxMatrix, nDimension )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )                 :: iDimension, jDimension, kDimension ! Counters
INTEGER( Kind= Int64 )                 :: nDimension                         ! Dimension
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: XYZ                                ! XYZ positions
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: TempXYZ                            ! XYZ positions (temporary)

! *********************************************************************************************** !
! REAL VARIABLES (PARAMETER)                                                                      !
! *********************************************************************************************** !
REAL( Kind= Real64 ), PARAMETER :: Delta = 0.75D0 ! Lovász delta

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                                       :: LovaszFactor             ! Lovász factor
REAL( Kind= Real64 )                                       :: LovaszCondition          ! Lovász condition
REAL( Kind= Real64 ), DIMENSION( nDimension * nDimension ) :: BoxMatrix                ! Box vectors
REAL( Kind= Real64 ), DIMENSION( nDimension, nDimension )  :: LatticeBasisVectors      ! Lattice basis vectors
REAL( Kind= Real64 ), DIMENSION( nDimension, nDimension )  :: LatticeBasisVectorsTemp  ! Lattice basis vectors (temporary)
REAL( Kind= Real64 ), DIMENSION( nDimension, nDimension )  :: OldLatticeBasisVectors   ! Lattice basis vectors (old)
REAL( Kind= Real64 ), DIMENSION( nDimension, nDimension )  :: OrthonormalBasisVectorGS ! Gram-Schmidt orthonormal basis vectors
REAL( Kind= Real64 ), DIMENSION( nDimension, nDimension )  :: CoefficientsGS           ! Gram-Schmidt coefficients

! Defining the lattice vectors
DO iDimension = 1, nDimension
  LatticeBasisVectors(:,iDimension) = BoxMatrix((1+nDimension*(iDimension-1)):(nDimension*iDimension))
END DO

! X, Y, and Z positions
XYZ(1) = 1
XYZ(2) = 2
XYZ(3) = 3

! Run LLL algorithm until convergence
LovaszCriterion: DO
  ! Store old lattice vectors
  OldLatticeBasisVectors = LatticeBasisVectors
  ! Initialization (working vector)
  kDimension = 2
  ! Lenstra–Lenstra–Lovász algorithm
  DO WHILE( kDimension <= nDimension )
    ! Initialize/reset Gram-Schmidt basis vectors and coefficients
    DO iDimension = 1, kDimension
      ! Initialization
      OrthonormalBasisVectorGS(:,iDimension) = LatticeBasisVectors(:,iDimension)
      ! Gram-Schimidt operators
      DO jDimension = 1, iDimension - 1
        ! Calculate the necessary Gram-Schmidt coefficients
        CoefficientsGS(iDimension,jDimension) = DOT_PRODUCT( LatticeBasisVectors(:,iDimension), &
        &                                       OrthonormalBasisVectorGS(:,jDimension) ) / &
        &                                       DOT_PRODUCT( OrthonormalBasisVectorGS(:,jDimension), &
        &                                       OrthonormalBasisVectorGS(:,jDimension) )
        ! Calculate the necessary Gram-Schmidt orthonormal basis vectors
        OrthonormalBasisVectorGS(:,iDimension) = OrthonormalBasisVectorGS(:,iDimension) - CoefficientsGS(iDimension,jDimension) * &
        &                                        OrthonormalBasisVectorGS(:,jDimension)
      END DO
    END DO
    ! Calculate the Gram-Schmidt reduction (K > 1)
    DO jDimension = 1, kDimension - 1
      ! Recalculate the necessary Gram-Schmidt coefficients
      CoefficientsGS(kDimension,jDimension) = DOT_PRODUCT( LatticeBasisVectors(:,kDimension), &
      &                                       OrthonormalBasisVectorGS(:,jDimension) ) / &
      &                                       DOT_PRODUCT( OrthonormalBasisVectorGS(:,jDimension), &
      &                                       OrthonormalBasisVectorGS(:,jDimension) )
      ! Compute the reduced lattice vector
      LatticeBasisVectors(:,kDimension) = LatticeBasisVectors(:,kDimension) - NINT( CoefficientsGS(kDimension,jDimension) ) * &
      &                                   LatticeBasisVectors(:,jDimension)
    END DO
    ! Calculate Lovász operators
    LovaszFactor    = Delta - ( CoefficientsGS(kDimension,kDimension-1) * CoefficientsGS(kDimension,kDimension-1) )
    LovaszCondition = LovaszFactor * DOT_PRODUCT( OrthonormalBasisVectorGS(:,kDimension-1), &
    &                 OrthonormalBasisVectorGS(:,kDimension-1) )
    ! Check the Lovász condition
    IF( DOT_PRODUCT( OrthonormalBasisVectorGS(:,kDimension), OrthonormalBasisVectorGS(:,kDimension) ) >= LovaszCondition ) THEN
      ! Iteration
      kDimension = kDimension + 1
    ELSE
      ! Swap working vector and preceding lattice vector
      LatticeBasisVectorsTemp(:,kDimension) = LatticeBasisVectors(:,kDimension)
      LatticeBasisVectors(:,kDimension)     = LatticeBasisVectors(:,kDimension-1)
      LatticeBasisVectors(:,kDimension-1)   = LatticeBasisVectorsTemp(:,kDimension)
      ! Swap XYZ positions
      TempXYZ(kDimension) = XYZ(kDimension)
      XYZ(kDimension)     = XYZ(kDimension-1)
      XYZ(kDimension-1)   = TempXYZ(kDimension)
      ! Return one iteration of K until K = 2
      kDimension = MAX( INT( kDimension ) - 1, 2 )
    END IF
  END DO
  ! Compare new lattice vectors with old lattice vectors
  DO iDimension = 1, nDimension
    DO jDimension = 1, nDimension
      ! Cycle if they differ
      IF( DABS( OldLatticeBasisVectors(jDimension,iDimension) - LatticeBasisVectors(jDimension,iDimension) ) >= &
      &   EPSILON( 1.D0 ) ) CYCLE LovaszCriterion
    END DO
  END DO
  EXIT LovaszCriterion
END DO LovaszCriterion

! Redefining the box vectors
IF( nDimension == 3 ) THEN
  ! Undo LLL swaps (may cause the inversion of the determinant sign)
  BoxMatrix(1:3) = LatticeBasisVectors(:,MINLOC( XYZ, Dim= 1 )) ! X vector
  BoxMatrix(7:9) = LatticeBasisVectors(:,MAXLOC( XYZ, Dim= 1 )) ! Z vector
  DO iDimension = 1, nDimension
    IF( iDimension /= MINLOC( XYZ, Dim= 1 ) .AND. iDimension /= MAXLOC( XYZ, Dim= 1 ) ) THEN
      BoxMatrix(4:6) = LatticeBasisVectors(:,iDimension) ! Y vector
    END IF
  END DO
ELSE IF( nDimension /= 3 ) THEN
  DO iDimension = 1, nDimension
    BoxMatrix((1+nDimension*(iDimension-1)):(nDimension*iDimension)) = LatticeBasisVectors(:,iDimension)
  END DO
END IF

RETURN

END SUBROUTINE LenstraLenstraLovaszMethod