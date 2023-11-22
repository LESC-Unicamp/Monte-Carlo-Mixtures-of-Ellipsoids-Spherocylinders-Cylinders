! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!                                to generate pseudorandom numbers.                                !
!                                                                                                 !
! Version number: 1.2.1                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       November 22nd, 2023                                       !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

! *********************************************************************************************** !
!                               Linear congruential generator (LCG)                               !
!    This subroutine generates a random number from the uniform distribution over the range       !
!   0 ≤ x < 1. It does not take any arguments. The number generator seed has an in/out intent,    !
!     i. e., its value is changed every time the RandomNumberGenLCG(  ) subroutine is called.     !
!                  This seed is the heart of the pseudorandom number generator.                   !
!         See Allen and Tildesley, 2nd Edition (2017), Appendix E, for more information.          !
! *********************************************************************************************** !
SUBROUTINE RandomNumberGenLCG(  )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER PARAMETERS                                                                              !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER :: L = 1029
INTEGER( Kind= Int64 ), PARAMETER :: C = 221591
INTEGER( Kind= Int64 ), PARAMETER :: G = 1048576

! Finite modulus arithmetic
SeedValue    = MOD( ( (SeedValue * L) + C ), G ) ! Pseudorandom number generator seed
RandomNumber = DBLE( SeedValue ) / DBLE( G )     ! Pseudorandom number

RETURN

END SUBROUTINE RandomNumberGenLCG