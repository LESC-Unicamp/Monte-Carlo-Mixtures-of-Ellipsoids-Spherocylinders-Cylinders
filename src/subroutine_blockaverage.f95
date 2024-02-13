! ############################################################################################### !
!               NVT/NPT-Monte Carlo algorithm for cylindrically-symmetric molecules               !
!                   This code contains all subroutines used in the main program                   !
!  to compute the TPT coefficients and their standard deviations using the block-average method.  !
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

! *********************************************************************************************** !
!     This subroutine uses a block-averaging method to calculate the first- and second-order      !
!         TPT coefficients, the perturbed Helmholtz free energy, and their uncertainties          !
! *********************************************************************************************** !
!                           Programmed by: Luis Fernando Mercier Franco                           !
!                               Modified by: Nathan Barros de Souza                               !
!                     University of Campinas, School of Chemical Engineering                      !
! *********************************************************************************************** !
!         See Allen and Tildesley, 2nd Edition (2017), pages 281-287 for more information         !
! *********************************************************************************************** !
SUBROUTINE BlockAverage( Flag, CoefficientTPT1, CoefficientTPT2, rFreeEnergy, CoefficientTPT1Deviation, &
&                        CoefficientTPT2Deviation, rFreeEnergyDeviation, ExecutionTime )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 )            :: nEquilibrationPoints   ! Equilibration data points in the potential file
INTEGER( Kind= Int64 )            :: nProductionPoints      ! Production data points in the potential file
INTEGER( Kind= Int64 )            :: nBlocks                ! Counter for the number of blocks
INTEGER( Kind= Int64 )            :: nBlockPoints           ! Number of data points in a block
INTEGER( Kind= Int64 )            :: iBlock, jBlock, kBlock ! Counters (block)
INTEGER( Kind= Int64 )            :: cBlock                 ! Counters (block)
INTEGER( Kind= Int64 )            :: iPoint                 ! Counter (data points)
INTEGER( Kind= Int64 ), PARAMETER :: MaxData = 10000        ! Maximum number of block data for linear regression

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: rFreeEnergyMoment1                     ! First moment (average), <exp(-U*/T*)>
REAL( Kind= Real64 ) :: rFreeEnergyMoment2                     ! Second moment, <exp(-2U*/T*)>
REAL( Kind= Real64 ) :: PotentialEnergyMoment1                 ! First moment (average), <U*>
REAL( Kind= Real64 ) :: PotentialEnergyMoment2                 ! Second moment, <U*²>
REAL( Kind= Real64 ) :: PotentialEnergyMoment3                 ! Third moment, <U*³>
REAL( Kind= Real64 ) :: PotentialEnergyMoment4                 ! Fourth moment, <U*⁴>
REAL( Kind= Real64 ) :: MaxExponentialArgument                 ! Largest exponential argument (see Supplementary Material of Abreu and Escobedo, J. Chem. Phys. 124)
REAL( Kind= Real64 ) :: BoltzmannFactorArgument                ! Boltzmann factor argument, -βU
REAL( Kind= Real64 ) :: rFreeEnergyVariance                    ! Variance of the 1st moment, <exp(-2U*/T*)> - <exp(-U*/T*)>²
REAL( Kind= Real64 ) :: PotentialEnergyVariance1               ! Variance of the 1st moment, <U*²> - <U*>²
REAL( Kind= Real64 ) :: PotentialEnergyVariance2               ! Variance of the 2nd moment, <U*⁴> - <U*²>²
REAL( Kind= Real64 ) :: rFreeEnergyBlockAverage                ! Block average (perturbed free energy)
REAL( Kind= Real64 ) :: PotentialEnergyBlockAverage            ! Block average (potential energy)
REAL( Kind= Real64 ) :: rFreeEnergyBlockVariance               ! Block variance (perturbed free energy)
REAL( Kind= Real64 ) :: PotentialEnergyBlockVariance           ! Block variance (potential energy)
REAL( Kind= Real64 ) :: rFreeEnergyVarianceSum                 ! Auxiliary summation variable for the block variance (perturbed free energy)
REAL( Kind= Real64 ) :: PotentialEnergyVarianceSum             ! Auxiliary summation variable for the block variance (potential energy)
REAL( Kind= Real64 ) :: bIntercept, aSlope                     ! Parameters of linear regression
REAL( Kind= Real64 ) :: DeterminationCoefficient               ! Coefficient of determination
REAL( Kind= Real64 ) :: rFreeEnergyStatisticalInefficiency     ! Statistical inefficiency
REAL( Kind= Real64 ) :: PotentialEnergyStatisticalInefficiency ! Statistical inefficiency
REAL( Kind= Real64 ) :: bFreeEnergyVariance                    ! Variance of <exp(-U*/T*)>
REAL( Kind= Real64 ) :: bPotentialEnergyVariance1              ! Variance of <U*>
REAL( Kind= Real64 ) :: bPotentialEnergyVariance2              ! Variance of <U*²>
REAL( Kind= Real64 ) :: Covariance                             ! Covariance between U* and U*²
REAL( Kind= Real64 ) :: Correlation                            ! Correlation between U* and U*²
REAL( Kind= Real64 ) :: ExponentialArgument                    ! Boltzmann factor argument
REAL( Kind= Real64 ) :: rFreeEnergy                            ! Perturbed Helmholtz free energy
REAL( Kind= Real64 ) :: CoefficientTPT1, CoefficientTPT2       ! TPT coefficients
REAL( Kind= Real64 ) :: rFreeEnergyDeviation                   ! Perturbed Helmholtz free energy (standard deviation)
REAL( Kind= Real64 ) :: CoefficientTPT1Deviation               ! TPT coefficient (standard deviation)
REAL( Kind= Real64 ) :: CoefficientTPT2Deviation               ! TPT coefficient (standard deviation)
REAL( Kind= Real64 ) :: InitialTime                            ! Initial time
REAL( Kind= Real64 ) :: FinalTime                              ! Final time
REAL( Kind= Real64 ) :: ExecutionTime                          ! Execution time


! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: bBlockSizeInverse                       ! Inverse of the number of data points in each block
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: bPotentialEnergyStatisticalInefficiency ! Statistical inefficiency per number of blocks (potential energy)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: bFreeEnergyStatisticalInefficiency      ! Statistical inefficiency per number of blocks (free energy)

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: PotentialEnergy ! Potential energy (production only)

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Flag ! Generic true/false flag

! Number of data points
nProductionPoints    = INT( ( DBLE( MaxSimulationCycles ) - DBLE( nEquilibrationCycles ) ) / DBLE( nSavingFrequency ) ) ! Production points
nEquilibrationPoints = INT( DBLE( nEquilibrationCycles ) / DBLE( nSavingFrequency ) ) ! Equilibration points

! Return to main program if number of block data exceeds 'max_data' parameter
IF ( MaxBlocks - MinBlocks >= MaxData ) THEN
  Flag = .TRUE.
  RETURN
END IF

! Initialization of variables
rFreeEnergyMoment1     = 0.D0
rFreeEnergyMoment2     = 0.D0
PotentialEnergyMoment1 = 0.D0
PotentialEnergyMoment2 = 0.D0
PotentialEnergyMoment3 = 0.D0
PotentialEnergyMoment4 = 0.D0

! Allocation
ALLOCATE( PotentialEnergy(nProductionPoints) )
ALLOCATE( bBlockSizeInverse(MaxData) )
ALLOCATE( bFreeEnergyStatisticalInefficiency(MaxData) )
ALLOCATE( bPotentialEnergyStatisticalInefficiency(MaxData) )

! Start subroutine timer
CALL CPU_Time( InitialTime )

! *********************************************************************************************** !
! Initialization of the Largest Exponential Argument                                              !
! *********************************************************************************************** !
!  This corresponds to a method to calculate the natural logarithm of large exponential           !
!  iterations, avoiding mathematical inaccuracies. This is applied to determine the Helmholtz     !
!  free energy of the perturbed system by taking the negative natural logarithm of the canonical  !
!  ensemble average over the configurations visited by the reference fluid.                       !
! *********************************************************************************************** !
!  See Supplementary Material of Abreu and Escobedo (2006) for more information.                  !
! *********************************************************************************************** !
MaxExponentialArgument = 0.D0

! Largest Exponential Argument and Potential Energy
OPEN( Unit= 150, File= "Potential_Energy/"//TRIM( DescriptorDate )//"/Range_"//TRIM( DescriptorRange )//"/"// &
&                      TRIM( DescriptorHour )//"_thermo_η"//TRIM( DescriptorFileThermoVariable )//"_C"// &
&                      TRIM( DescriptorFileComponents )//"_"//TRIM( DescriptorFileGeometry )//".dat", Action= "READ" )
! Skip header
READ ( 150, * )
! Skip equilibration data points (if necessary)
IF( .NOT. PotentialEnergyLogical ) THEN
  DO iPoint = 1, nEquilibrationPoints
    READ ( 150, * )
  END DO
END IF
! Read production data points
DO iPoint = 1, nProductionPoints
  READ ( 150, * ) Dummy, PotentialEnergy(iPoint)
  BoltzmannFactorArgument = - PotentialEnergy(iPoint) / ReducedTemperature
  IF( BoltzmannFactorArgument > MaxExponentialArgument ) THEN
    MaxExponentialArgument = BoltzmannFactorArgument
  END IF
END DO
CLOSE( 150 )

! Iterative process
DO iPoint = 1, nProductionPoints
  ExponentialArgument    = - ( PotentialEnergy(iPoint) / ReducedTemperature ) - MaxExponentialArgument
  rFreeEnergyMoment1     = rFreeEnergyMoment1 + DEXP( ExponentialArgument )
  rFreeEnergyMoment2     = rFreeEnergyMoment2 + DEXP( 2.D0 * ExponentialArgument )
  PotentialEnergyMoment1 = PotentialEnergyMoment1 + ( PotentialEnergy(iPoint) )
  PotentialEnergyMoment2 = PotentialEnergyMoment2 + ( PotentialEnergy(iPoint) * PotentialEnergy(iPoint) )
  PotentialEnergyMoment3 = PotentialEnergyMoment3 + ( PotentialEnergy(iPoint) * PotentialEnergy(iPoint) * PotentialEnergy(iPoint) )
  PotentialEnergyMoment4 = PotentialEnergyMoment4 + ( PotentialEnergy(iPoint) * PotentialEnergy(iPoint) * PotentialEnergy(iPoint) &
  &                        * PotentialEnergy(iPoint) )
END DO

! First moments (average)
rFreeEnergyMoment1     = rFreeEnergyMoment1 / DBLE( nProductionPoints )
PotentialEnergyMoment1 = PotentialEnergyMoment1 / DBLE( nProductionPoints )

! Second moments
rFreeEnergyMoment2     = rFreeEnergyMoment2 / DBLE( nProductionPoints )
PotentialEnergyMoment2 = PotentialEnergyMoment2 / DBLE( nProductionPoints )

! Third moment
PotentialEnergyMoment3 = PotentialEnergyMoment3 / DBLE( nProductionPoints )

! Fourth moment
PotentialEnergyMoment4 = PotentialEnergyMoment4 / DBLE( nProductionPoints )

! Variances of <X>
rFreeEnergyVariance      = rFreeEnergyMoment2 - ( rFreeEnergyMoment1 * rFreeEnergyMoment1 )
PotentialEnergyVariance1 = PotentialEnergyMoment2 - ( PotentialEnergyMoment1 * PotentialEnergyMoment1 )

! Variance of <U*²>
PotentialEnergyVariance2 = PotentialEnergyMoment4 - ( PotentialEnergyMoment2 * PotentialEnergyMoment2 )

! Covariance
Covariance = PotentialEnergyMoment3 - ( PotentialEnergyMoment1 * PotentialEnergyMoment2 )

! Number of blocks
DO nBlocks = MinBlocks, MaxBlocks
  ! Steps in each block for a certain number of blocks
  nBlockPoints = INT( DBLE( nProductionPoints ) / DBLE( nBlocks ) )
  ! Initialization
  rFreeEnergyVarianceSum     = 0.D0
  PotentialEnergyVarianceSum = 0.D0
  jBlock                     = 1
  ! Block average
  DO kBlock = 1, nBlocks
    rFreeEnergyBlockAverage     = 0.D0 ! Reset
    PotentialEnergyBlockAverage = 0.D0 ! Reset
    DO iBlock = 1, nBlockPoints
      ! Block average
      ExponentialArgument = - ( PotentialEnergy(jBlock) / ReducedTemperature ) - MaxExponentialArgument
      rFreeEnergyBlockAverage = rFreeEnergyBlockAverage + DEXP( ExponentialArgument )
      PotentialEnergyBlockAverage = PotentialEnergyBlockAverage + PotentialEnergy(jBlock)
      ! Increment counter
      jBlock = jBlock + 1
    END DO
    rFreeEnergyBlockAverage = rFreeEnergyBlockAverage / DBLE( nBlockPoints )
    rFreeEnergyVarianceSum  = rFreeEnergyVarianceSum + ( rFreeEnergyBlockAverage - rFreeEnergyMoment1 ) ** 2.D0
    PotentialEnergyBlockAverage = PotentialEnergyBlockAverage / DBLE( nBlockPoints )
    PotentialEnergyVarianceSum  = PotentialEnergyVarianceSum + ( PotentialEnergyBlockAverage - PotentialEnergyMoment1 ) ** 2.D0
  END DO
  ! Variance of the block averages
  rFreeEnergyBlockVariance     = rFreeEnergyVarianceSum / DBLE( nBlocks - 1 )
  PotentialEnergyBlockVariance = PotentialEnergyVarianceSum / DBLE( nBlocks - 1 )
  ! Auxiliary counter
  cBlock = nBlocks - MinBlocks + 1
  ! Inverse of block size
  bBlockSizeInverse(cBlock) = 1.D0 / DBLE( nBlockPoints )
  ! Statistical inefficiency
  bFreeEnergyStatisticalInefficiency(cBlock)      = DBLE( nBlockPoints ) * ( rFreeEnergyBlockVariance / rFreeEnergyVariance )
  bPotentialEnergyStatisticalInefficiency(cBlock) = DBLE( nBlockPoints ) * ( PotentialEnergyBlockVariance / &
  &                                                 PotentialEnergyVariance1 )
END DO

! Deallocation
DEALLOCATE( PotentialEnergy )

! Linear regression (perturbed Helmholtz free energy)
CALL LinearFitting( bBlockSizeInverse, bFreeEnergyStatisticalInefficiency, bIntercept, aSlope, DeterminationCoefficient, cBlock )

! Statistical inefficiency (perturbed Helmholtz free energy)
rFreeEnergyStatisticalInefficiency = bIntercept

! Linear regression (TPT coefficients)
CALL LinearFitting( bBlockSizeInverse, bPotentialEnergyStatisticalInefficiency, bIntercept, aSlope, DeterminationCoefficient, &
&                   cBlock )

! Statistical inefficiency (TPT coefficients)
PotentialEnergyStatisticalInefficiency = bIntercept

! Deallocation
DEALLOCATE( bBlockSizeInverse, bFreeEnergyStatisticalInefficiency, bPotentialEnergyStatisticalInefficiency )

! Variances of <X>
bFreeEnergyVariance       = rFreeEnergyVariance * ( rFreeEnergyStatisticalInefficiency / DBLE( nProductionPoints ) )
bPotentialEnergyVariance1 = PotentialEnergyVariance1 * ( PotentialEnergyStatisticalInefficiency / DBLE( nProductionPoints ) )

! Variance of <U*²>
bPotentialEnergyVariance2 = PotentialEnergyVariance2 * ( PotentialEnergyStatisticalInefficiency / DBLE( nProductionPoints ) )

! Covariance between U* and U*²
Covariance = Covariance * ( PotentialEnergyStatisticalInefficiency / DBLE( nProductionPoints ) )

! Correlation between U* and U*²
Correlation = Covariance / ( DSQRT( bPotentialEnergyVariance1 ) ) / ( DSQRT( bPotentialEnergyVariance2 ) )

! TPT Coefficients
CoefficientTPT1 = PotentialEnergyMoment1 / nParticles
CoefficientTPT2 = - 0.5D0 * ( PotentialEnergyMoment2 - ( PotentialEnergyMoment1 * PotentialEnergyMoment1 ) ) / nParticles

! TPT Coefficients (Uncertainty propagation)
CoefficientTPT1Deviation = DSQRT( bPotentialEnergyVariance1 ) / nParticles
CoefficientTPT2Deviation = 0.5D0 * ( DSQRT( ( 4.D0 * bPotentialEnergyVariance1 * PotentialEnergyMoment1 * &
&                          PotentialEnergyMoment1 ) + bPotentialEnergyVariance2 - ( 4.D0 * PotentialEnergyMoment1 * &
&                          Covariance ) ) ) / nParticles

! Perturbed Helmholtz free energy
rFreeEnergy = - ( 1.D0 / nParticles ) * ( DLOG( rFreeEnergyMoment1 ) + MaxExponentialArgument )

! Perturbed Helmholtz free energy (Uncertainty propagation)
rFreeEnergyDeviation = ( 1.D0 / nParticles ) * ( DSQRT( ( bFreeEnergyVariance ) / ( rFreeEnergyMoment1 * rFreeEnergyMoment1 ) ) )

! Finish subroutine timer
CALL CPU_Time( FinalTime )

! Execution time
ExecutionTime = FinalTime - InitialTime

RETURN

END SUBROUTINE BlockAverage

! *********************************************************************************************** !
!                    This subroutine makes a linear regression of x and y data                    !
! *********************************************************************************************** !
!                           Programmed by: Luis Fernando Mercier Franco                           !
!                     University of Campinas, School of Chemical Engineering                      !
! *********************************************************************************************** !
SUBROUTINE LinearFitting( IndependentVar, DependentVar, bIntercept, aSlope, DeterminationCoefficient, nPoints )

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iPoint, nPoints

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                       :: SumX                     ! Sum of X
REAL( Kind= Real64 )                       :: SumY                     ! Sum of Y
REAL( Kind= Real64 )                       :: SumXSquared              ! Sum of X²
REAL( Kind= Real64 )                       :: SumXY                    ! Sum of XY
REAL( Kind= Real64 )                       :: SumYSquared              ! Sum of Y²
REAL( Kind= Real64 )                       :: bIntercept               ! Y-intercept
REAL( Kind= Real64 )                       :: aSlope                   ! Slope
REAL( Kind= Real64 )                       :: DeterminationCoefficient ! Coefficient of determination
REAL( Kind= Real64 ), DIMENSION( nPoints ) :: IndependentVar           ! Independent variable
REAL( Kind= Real64 ), DIMENSION( nPoints ) :: DependentVar             ! Dependent variable

! Initialization
SumX        = 0.D0
SumY        = 0.D0
SumXY       = 0.D0
SumXSquared = 0.D0
SumYSquared = 0.D0

! Iteration
DO iPoint = 1, nPoints
  SumX  = SumX + IndependentVar(iPoint)
  SumY  = SumY + DependentVar(iPoint)
  SumXY = SumXY + ( IndependentVar(iPoint) * DependentVar(iPoint) )
  SumXSquared = SumXSquared + ( IndependentVar(iPoint) * IndependentVar(iPoint) )
  SumYSquared = SumYSquared + ( DependentVar(iPoint) * DependentVar(iPoint) )
END DO

! Slope
aSlope = ( ( DBLE( nPoints ) * SumXY ) - ( SumX * SumY ) ) / ( ( DBLE( nPoints ) * SumXSquared ) - ( SumX * SumX ) )
! Y-intercept
bIntercept = ( SumY - ( aSlope * SumX ) ) / DBLE( nPoints )
! Coefficient of correlation
DeterminationCoefficient = ( ( DBLE( nPoints ) * SumXY ) - ( SumX * SumY ) ) / DSQRT( ( ( DBLE( nPoints ) * SumXSquared ) - &
&                          ( SumX * SumX ) ) * ( ( DBLE( nPoints ) * SumYSquared ) - ( SumY * SumY ) ) )
! Coefficient of determination
DeterminationCoefficient = ( DeterminationCoefficient * DeterminationCoefficient )

RETURN

END SUBROUTINE LinearFitting