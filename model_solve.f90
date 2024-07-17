MODULE MODEL_SOLVE
USE PARAMETERS
USE GLOBAL
USE auxiliary
!---------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------


CONTAINS

  !---------------------------------------------------------------------------
  ! calculates marginal utility of consumption
  !---------------------------------------------------------------------------
  FUNCTION margu(cons, lab, it)
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  ::    cons,   lab
      INTEGER,  INTENT(IN)  ::    it
      ! outpus
      REAL(DP) :: margu
      
      margu = nu*(cons**nu*(1.00-lab)**(1.00-nu))**egam/(p(it)*cons)
        
  END FUNCTION margu

  
  
  !-------------------------------------------------------------------------------
  ! calculates the value function
  !-------------------------------------------------------------------------------
  FUNCTION valuefunc(a_plus, cons, lab, ij, ip, is, it)

        IMPLICIT NONE
        INTEGER,  INTENT(IN) :: ij, ip, is, it
        REAL(DP), INTENT(IN) :: a_plus, cons, lab
        !LOCAL VARIABLE
        REAL(DP) :: valuefunc, varphi, c_help, l_help
        INTEGER :: ial, iar, itp

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)
        l_help = min(max(lab, 0d0),1d0-1d-10)

        ! get tomorrows year
        itp = year(it, ij, ij+1)

        ! get tomorrows utility
        !CALL linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
        CALL gridlookup1(NA, a, a_plus, ial, varphi)
        iar = ial + 1
        
        
        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        IF(ij < JJ) THEN
            valuefunc = DMAX1(varphi*EV(ij+1, ial, ip, is, itp) + &
                (1d0-varphi)*EV(ij+1, iar, ip, is, itp), 1d-10)**egam/egam
        ENDIF

        ! add todays part and discount
        valuefunc = (c_help**nu*(1.0-l_help)**(1.0-nu))**egam/egam + beta*valuefunc

  END FUNCTION valuefunc
  
  
  
  !---------------------------------------------------------
  ! calculates year at which age ij agent is ij_p
  !--------------------------------------------------------
  FUNCTION year(it, ij, ij_p)

      implicit none
      integer, intent(in) :: it, ij, ij_p
      integer :: year

      year = it + ij_p - ij

      if(it == 0 .or. year <= 0)year = 0
      if(it == TT .or. year >= TT)year = TT

  END FUNCTION year
  
  
  
  !-------------------------------------------------------------------------------------------------
  ! the value function
  !-------------------------------------------------------------------------------------------------
    FUNCTION valuefunc2(x_in)

        IMPLICIT NONE
        ! Inputs
        REAL(DP), INTENT(IN) :: x_in
        ! Outputs
        REAL(DP) :: valuefunc2
        ! Local variables
        REAL(DP) :: a_plus, wage, v_ind, available


        ! calculate tomorrows assets
        a_plus  = x_in

        ! get lsra transfer payment
        v_ind = v(ij_com, ia_com, ip_com, is_com, it_com)

        ! calculate the wage rate
        wage = wn(it_com)*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! calculate available resources
        available = (1.00+rn(it_com))*a(ia_com) + pen(ij_com, it_com) + v_ind

        ! determine labor
        if(ij_com < JR)then
            lab_com = DMIN1( DMAX1( nu + (1d0-nu)*(a_plus-available)/wage, 0.0_dp) , 1d0-1d-10)
        else
            lab_com = 0.00_dp
        endif

        ! calculate consumption
        cons_com = max( (available + wage*lab_com - a_plus)/p(it_com) , 1d-10)

        ! calculates the value function
        valuefunc2 = -(valuefunc(a_plus, cons_com, lab_com, ij_com, ip_com, is_com, it_com))

    END FUNCTION valuefunc2
    
    
    !-----------------------------------------------------------------------------------------------
    ! Golden Search Method
    !-----------------------------------------------------------------------------------------------
    SUBROUTINE goldensolver(xroot, fret, lowerbound, upperbound)
        !Inputs
        REAL(DP), INTENT(IN)   ::  lowerbound, upperbound
        !Outputs
        REAL(DP), INTENT(OUT)   :: fret
        REAL(DP), INTENT(INOUT) :: xroot 
        !Local variables
        REAL(DP) :: eps,  goldenpi   
        REAL(DP) :: lb, ub, xlb, xub, flow, fhigh
        INTEGER  :: giter
        
        EPS   = 0.000001
        giter = 1
        
        ! determine lower bound and upper bound
        lb = lowerbound
        ub = upperbound
        ! golden search pi
        goldenpi = (3.00_dp - SQRT(5.00_dp))/2.00_dp
        
        !----------------------------------------------
        ! Begin Golden Section Searching
        !----------------------------------------------
        DO WHILE( (ABS(ub-lb) > eps) .AND. (giter < 300))
            
            ! calculate x1 and x2 and function values
            xlb = lb + goldenpi*(ub-lb)
            xub = lb + (1.00_dp - goldenpi)*(ub-lb)
            
            flow  = valuefunc2(xlb)
            fhigh = valuefunc2(xub)
            
            !---------------------------
            !Get new values / Updating
            !---------------------------
            IF(flow < fhigh) THEN
              ub    = xub
              xroot = xlb
              fret  = flow
            ELSE  
              lb    = xlb
              xroot = xub
              fret  = fhigh
            ENDIF
            
             giter = giter + 1
        
        END DO
        
        
        
        
    
    END SUBROUTINE goldensolver
  



END MODULE MODEL_SOLVE

