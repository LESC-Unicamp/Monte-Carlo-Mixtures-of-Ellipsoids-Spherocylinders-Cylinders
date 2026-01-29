! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                                to generate pseudorandom numbers.                                !
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
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                       Random number generator based on bitwise operations                       !
!          Please visit <https://github.com/MaginnGroup/Cassandra> for more information.          !
! *********************************************************************************************** !
SUBROUTINE RandomNumberGenBitwise(  )

! Uses one module: global variables
USE GlobalVar
    
IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLE                                                                                !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: bSave

! Bitwise operations
bSave = ISHFT( IEOR( ISHFT( SeedValue(1), 1 ), SeedValue(1) ), -53_INT64 )
SeedValue(1) = IEOR( ISHFT( IAND( SeedValue(1), -2_INT64 ), 10 ), bSave )
bSave = ISHFT( IEOR( ISHFT( SeedValue(2), 24 ), SeedValue(2) ), -50 )
SeedValue(2) = IEOR( ISHFT( IAND( SeedValue(2), -512_INT64 ), 5 ), bSave )
bSave = ISHFT( IEOR( ISHFT( SeedValue(3), 3 ), SeedValue(3) ), -23 )
SeedValue(3) = IEOR( ISHFT( IAND( SeedValue(3), -4096_INT64 ), 29 ), bSave )
bSave = ISHFT( IEOR( ISHFT( SeedValue(4), 5 ), SeedValue(4) ), -24 )
SeedValue(4) = IEOR( ISHFT( IAND( SeedValue(4), -131072_INT64 ), 23 ), bSave )
bSave = ISHFT( IEOR( ISHFT( SeedValue(5), 3 ), SeedValue(5) ) , -33 )
SeedValue(5) = IEOR( ISHFT( IAND( SeedValue(5), -8388608_INT64 ), 8 ), bSave )

! Random number in range [0,1[
RandomNumber = IEOR( IEOR( IEOR( IEOR( SeedValue(1), SeedValue(2) ), SeedValue(3) ), SeedValue(4) ), SeedValue(5) ) * &
&              5.4210108624275221E-20_REAL64 + 0.5_REAL64

! Stop condition
IF( RandomNumber >= 1.0_REAL64 ) WRITE( *, "(G0)" ) "Random number >= 1.0"
IF( RandomNumber >= 1.0_REAL64 ) CALL Exit(  )

RETURN

END SUBROUTINE RandomNumberGenBitwise