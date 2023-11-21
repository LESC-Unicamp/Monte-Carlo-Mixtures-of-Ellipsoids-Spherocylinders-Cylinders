! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                         to make operations using vectors and matrices.                          !
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
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!               This subroutine calculates the inverse of a matrix using cofactors                !
! *********************************************************************************************** !
SUBROUTINE InverseMatrixCofactorVec( InputMatrix, InverseMatrix, Determinant )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: Determinant        ! Determinant
REAL( Kind= Real64 )                 :: DeterminantInverse ! Inverse of determinant
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InverseMatrix      ! Matrix (inverse)
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InputMatrix        ! Matrix (input)
REAL( Kind= Real64 ), DIMENSION( 9 ) :: tCofactorMatrix    ! Inverse matrix

! Tranpose matrix of the matrix of cofactors
tCofactorMatrix(1) = InputMatrix(5) * InputMatrix(9) - InputMatrix(6) * InputMatrix(8)
tCofactorMatrix(2) = InputMatrix(3) * InputMatrix(8) - InputMatrix(2) * InputMatrix(9)
tCofactorMatrix(3) = InputMatrix(2) * InputMatrix(6) - InputMatrix(3) * InputMatrix(5)
tCofactorMatrix(4) = InputMatrix(6) * InputMatrix(7) - InputMatrix(4) * InputMatrix(9)
tCofactorMatrix(5) = InputMatrix(1) * InputMatrix(9) - InputMatrix(3) * InputMatrix(7)
tCofactorMatrix(6) = InputMatrix(3) * InputMatrix(4) - InputMatrix(1) * InputMatrix(6)
tCofactorMatrix(7) = InputMatrix(4) * InputMatrix(8) - InputMatrix(5) * InputMatrix(7)
tCofactorMatrix(8) = InputMatrix(2) * InputMatrix(7) - InputMatrix(1) * InputMatrix(8)
tCofactorMatrix(9) = InputMatrix(1) * InputMatrix(5) - InputMatrix(2) * InputMatrix(4)

! Determinant of matrix
Determinant = InputMatrix(1) * tCofactorMatrix(1) + InputMatrix(4) * tCofactorMatrix(2) + InputMatrix(7) * tCofactorMatrix(3)

! Inverse of the determinant of matrix
DeterminantInverse = 0.D0
IF( DABS( Determinant ) > 0.D0 ) THEN
  DeterminantInverse = 1.0D0 / Determinant
END IF

! Inverse matrix
InverseMatrix = DeterminantInverse * tCofactorMatrix

RETURN

END SUBROUTINE InverseMatrixCofactorVec

! *********************************************************************************************** !
!                        This subroutine multiplies a matrix and a vector                         !
! *********************************************************************************************** !
SUBROUTINE MatrixVectorMultiplication( InputMatrix, InputVector, OutputVector )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InputMatrix  ! Input matrix
REAL( Kind= Real64 ), DIMENSION( 3 ) :: InputVector  ! Vector (input)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: OutputVector ! Vector (output)

! Multiplication of a matrix and a vector
OutputVector(1) = InputMatrix(1) * InputVector(1) + InputMatrix(4) * InputVector(2) + InputMatrix(7) * InputVector(3)
OutputVector(2) = InputMatrix(2) * InputVector(1) + InputMatrix(5) * InputVector(2) + InputMatrix(8) * InputVector(3)
OutputVector(3) = InputMatrix(3) * InputVector(1) + InputMatrix(6) * InputVector(2) + InputMatrix(9) * InputVector(3)

RETURN

END SUBROUTINE MatrixVectorMultiplication