! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!         This code contains the algorithms of the force fields used in the main program          !
!                         to compute the potential energy of the system.                          !
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
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE ForceFields

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!            This subroutine computes the pair potential between particles i and j by             !
!             applying a discrete square-well potential to compute the pair potential             !
! *********************************************************************************************** !
SUBROUTINE SquareWellPotential( SquaredDistance, iComponent, jComponent, PairPotentialEnergy )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Counters (component)
INTEGER( Kind= Int64 ) :: rRange                 ! Counter (range)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                      :: SquaredDistance       ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                      :: SquaredPotentialRange ! Effective range of attraction (squared)
REAL( Kind= Real64 ), DIMENSION( nRange ) :: PairPotentialEnergy   ! Pair potential energy

! Compute pair potential for every attractive range
DO rRange = 1, nRange
  ! Effective range of attraction (squared)
  SquaredPotentialRange = 0.5D0 * ( cPotentialRange(iComponent,rRange) + cPotentialRange(jComponent,rRange) )
  SquaredPotentialRange = SquaredPotentialRange * SquaredPotentialRange
  ! Overlap in the region of attraction
  IF( SquaredDistance <= SquaredPotentialRange ) THEN
    ! Pair potential (reduced units)
    PairPotentialEnergy(rRange) = - 1.D0
  ELSE
    PairPotentialEnergy(rRange) = 0.D0
  END IF
END DO

RETURN

END SUBROUTINE SquareWellPotential

END MODULE ForceFields