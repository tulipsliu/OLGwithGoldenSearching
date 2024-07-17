MODULE PARAMETERS
IMPLICIT NONE
!--------------------------------------------------------------
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(P=13, R=307)
!--------------------------------------------------------------
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COUNTE VARIABLE                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER, PARAMETER  ::               &
    TT = 40,                         &  ! number of transition periods
    JJ = 12,                         &  ! number of years the household lives
    JR = 10,                         &  ! number of years the household retires
    NP = 2,                          &  ! number of persistent shock process values
    NS = 5,                          &  ! number of transitory shock process values
    NA = 100                            ! number of points on the asset grid
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HOUSEHOLD PREFERENCES PARAMETERS            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(DP), PARAMETER ::               &
    gamma = 0.5_dp,                  &  ! intertemporal elasticity of substitution
    beta  = 0.998**5,                &  ! time discout factor (5-years)
    nu    = 0.335,                   &  ! labor force-participation ratio
    egam  = 1.0-gamma  
    
!------------------------------------------------
! household risk process
!------------------------------------------------        

REAL(DP), PARAMETER  ::             &
    sigma_theta = 0.23,             &
    sigma_eps   = 0.0500,           &  ! variance of epsilon
    rho         = 0.9800,           &  ! autocorrelation coefficient
    multiple    = 1.0                  !  
!------------------------------------------------
! production parameters/ TECHNOLOGY
!------------------------------------------------

REAL(DP), PARAMETER  ::                             &
    alpha       = 0.36,                             &  
    delta       = 1.00_dp-(1.00_dp-0.0823_dp)**5,   &   ! Depreciation rate
    Omega       = 1.60
    
!------------------------------------------------
! size of the asset grid
!------------------------------------------------

REAL(DP), PARAMETER  ::              &
    a_l  = 0.0,                      &
    a_u  = 35.00,                    &
    a_grow = 0.01
    
!-------------------------------------------------
! demographic parameters
!-------------------------------------------------

REAL(DP), PARAMETER :: n_p   = (1.0+0.01)**5-1.00

!-------------------------------------------------
! simulation parameters
!-------------------------------------------------

REAL(DP), PARAMETER  :: damp = 0.30, sig = 0.0001

INTEGER,  PARAMETER  :: itermax = 50

INTEGER :: iter

         
!------------------------------------------------------
! Calculating time elapsed
!------------------------------------------------------

REAL(DP)  :: cpu_start_time,  cpu_stop_time


END MODULE PARAMETERS
