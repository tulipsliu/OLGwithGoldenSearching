MODULE GLOBAL
USE PARAMETERS
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! macroeconomic variables                  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(DP) ::  r(0:TT), rn(0:TT), w(0:TT), wn(0:TT), p(0:TT),       &
               KK(0:TT), AA(0:TT), BB(0:TT), LL(0:TT), HH(0:TT),  &
                YY(0:TT), CC(0:TT), II(0:TT), GG(0:TT), INC(0:TT)
                
!----------------------------------------
! government variables
!----------------------------------------

REAL(DP) :: tauc(0:TT), tauw(0:TT), taur(0:TT), taup(0:TT), kappa(0:TT),  &
            gy, by, pen(JJ,0:TT), PP(0:TT), taxrev(4,0:TT)
INTEGER  :: tax(0:TT)            

!----------------------------------------------------
! LSRA variables
!----------------------------------------------------

REAL(DP) :: BA(0:TT) = 0d0, SV(0:TT) = 0d0, lsra_comp, lsra_all, Vstar
LOGICAL  :: lsra_on


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COHORT AGGREGATE VARIABLES  !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(DP) ::  c_coh(JJ, 0:TT), l_coh(JJ, 0:TT), y_coh(JJ, 0:TT), a_coh(JJ, 0:TT)
REAL(DP) ::  v_coh(JJ, 0:TT)= 0d0, VV_coh(JJ, 0:TT) = 0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STOCHASTIC PROCESS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(DP) :: dist_theta(NP), theta(NP),     &
            pi(NS, NS), eta(NS)
 
INTEGER  :: is_initial = 3    ! determine \phi(1,a_v,\theta_p, \eta_s)

!----------------------------------------------
! demographic and other model parameters
!----------------------------------------------

REAL(DP) :: m(JJ, 0:TT), pop(JJ, 0:TT), eff(JJ)

!---------------------------------------------------------------------------------------------------
! HOUSEHOLD   ||  individual variables
!---------------------------------------------------------------------------------------------------

REAL(DP) :: a(0:NA), aplus(JJ, 0:NA, NP, NS, 0:TT),                & ! asset grids with respect to a, and the policy function asset accumulation in next period $t+1$  
     c(JJ, 0:NA, NP, NS, 0:TT), l(JJ, 0:NA, NP, NS, 0:TT)            ! policy functions: consumption c, and labor supply l.
REAL(DP) ::    phi(JJ, 0:NA, NP, NS, 0:TT), VV(JJ, 0:NA, NP, NS, 0:TT) = 0d0  ! Invariant distribution with transition dynamic (policy feform) time index $t$.
     
REAL(DP) :: v(JJ, 0:NA, NP, NS, 0:TT) = 0d0, FLC(JJ,0:TT)

!----------------------------------------------------
! numerical variables
!----------------------------------------------------

REAL(DP) :: RHS(JJ, 0:NA, NP, NS, 0:TT),  & ! Right hand side
            EV(JJ, 0:NA, NP, NS, 0:TT)      ! Expected value function
            
INTEGER  :: ij_com, ia_com, ip_com, is_com, it_com ! communative variables
REAL(DP) :: cons_com, lab_com, DIFF(0:TT)            


END MODULE GLOBAL
