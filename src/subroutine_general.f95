! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!              This code contains all general subroutines used in the main program.               !
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
!                                           O. K. Smith                                           !
!                           Communications of the ACM, 4(4), 168 (1961)                           !
!                                    DOI: 10.1145/355578.366316                                   !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
!                             --------------------------------------                              !
!                  J. de Graaf, L. Filion, M. Marechal, R. von Roij, M. Dijkstra                  !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
!                             --------------------------------------                              !
!                                D. Gottwald, G. Kahl, C. N. Likos                                !
!                                J. Chem. Phys. 122, 204503 (2005)                                !
!                                     DOI: 10.1063/1.1901585                                      !
!                             --------------------------------------                              !
!                           A. K. Lenstra, H. W. Lenstra Jr., L. Lovász                           !
!                                 Math. Ann. 261, 515-534 (1982)                                  !
!                                     DOI: 10.1007/BF01457454                                     !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!            This subroutine generates a progress bar for the Monte Carlo simulation.             !
! *********************************************************************************************** !
SUBROUTINE ProgressBarMC( iCycle, nCycle, Ensemble )

! Use kind Real64 and Int64, and Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iCycle       ! Counter
INTEGER( Kind= Int64 ) :: nCycle       ! Total/Maximum number of cycles
INTEGER( Kind= Int64 ) :: AuxiliarInt1 ! Auxiliar

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 23 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 03 ) :: Ensemble    ! Ensemble type
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0

! Progress bar (FORMAT)
IF( ( DBLE( iCycle ) / DBLE( nCycle ) ) * 100.D0 < 10.D0 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F4.2)" ) ( DBLE( iCycle ) / DBLE( nCycle ) ) * 100.D0
ELSE IF( ( DBLE( iCycle ) / DBLE( nCycle ) ) * 100.D0 < 100.D0 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ?????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F5.2)" ) ( DBLE( iCycle ) / DBLE( nCycle ) ) * 100.D0
ELSE
  AuxiliarInt1 = 2
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ??????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F6.2)" ) ( DBLE( iCycle ) / DBLE( nCycle ) ) * 100.D0
END IF
ProgressBar((21+AuxiliarInt1):23) = REPEAT( " ", ( (2 - AuxiliarInt1) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 20 + AuxiliarInt1 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(20+AuxiliarInt1+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarMC

! *********************************************************************************************** !
!                 This function converts any string into uppercase (from A to Z)                  !
! *********************************************************************************************** !
SUBROUTINE ToUpper( StringInput, StringLength, StringOutput )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER :: StringLength ! String length
INTEGER :: iCharacter   ! ASCII character code
INTEGER :: iString      ! Counter

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= StringLength ) :: StringInput  ! Length of input string
CHARACTER( LEN= StringLength ) :: StringOutput ! Length of output string

! Character positions
DO iString = 1, StringLength
  ! ASCII character code
  iCharacter = IACHAR( StringInput(iString:iString) )
  ! Convert to uppercase (letters only)
  IF( iCharacter >= IACHAR( "a" ) .AND. iCharacter <= IACHAR( "z" ) ) THEN
    StringOutput(iString:iString) = ACHAR(IACHAR(StringInput(iString:iString))-32)
  ! Do not convert symbols or numbers (special characters included)
  ELSE
    StringOutput(iString:iString) = StringInput(iString:iString)
  END IF
END DO

RETURN

END SUBROUTINE ToUpper