MODULE auxiliary
! ----------------------------------------------------------------------
IMPLICIT NONE
! -----------
INTEGER, PARAMETER, PRIVATE  ::  rk          = SELECTED_REAL_KIND(P=15, R=307)
INTEGER, PARAMETER, PRIVATE  ::  ik9          =  SELECTED_INT_KIND(9)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
CONTAINS
!BL
subroutine markovss(pij,sizex,tol,xss)
    ! M Yum: obtain ergodic distribution

    integer, intent(in) :: sizex
    real(8), intent(in), dimension(sizex,sizex) :: pij
    real(8), intent(in) :: tol
    real(8), intent(out), dimension(sizex) :: xss

    integer :: ix, jx
    real(8) :: dist
    real(8), dimension(sizex) :: xold, xnew

    xold = (1.0d0/dble(sizex))
    xnew = xold
    dist = 1000.0d0

    do while (dist > tol*2.0d0)
        do jx = 1,sizex
            xnew(jx) = 0.0d0
            do ix = 1,sizex
                xnew(jx) = xnew(jx) + pij(ix,jx)*xold(ix)
            end do
        end do
        dist = maxval(abs(log(xnew)-log(xold)))
        xold = xnew
    end do
    xss = xnew

end subroutine markovss
!BL    
! ==========================================================================
!      SUBROUTINE:  linspace [create a linearly-spaced vector] ::Julia K. Thomas
! ==========================================================================
subroutine linspace(lb, ub, gridnum, x)
    integer(ik9):: gridnum, j1
    real(rk):: lb, ub, x(gridnum), y(gridnum-2_ik9)

    intent(in):: lb, ub, gridnum
    intent(out):: x

    ! this subroutine produces a lin-spaced row vector.
    do j1 = 1_ik9, gridnum - 2_ik9
        y(j1) = dble(j1)
    end do

    y = ((ub - lb)/(gridnum-1_ik9))*y + lb
    x(1_ik9) = lb
    x(2_ik9:gridnum-1_ik9) = dble(y)
    x(gridnum) = ub

end subroutine linspace

! ==========================================================================
!      SUBROUTINE:  logspace [create a log-spaced vector]  :: Aubhik Khan
! ==========================================================================
subroutine logspace(lowerbound, upperbound, n, x)

    integer:: n, j1
    real(rk):: lowerbound, upperbound, term0, lb, ub
    real(rk):: x(n), y(n-1)
    intent(in):: lowerbound, upperbound, n
    intent(out):: x

    ! this code departs from matlab in that lower and upper bounds
    ! are the actual bounds on the logarithmically spaced vector.
    term0 = dlog(10.0_rk)
    lb = dlog(lowerbound)/term0
    ub = dlog(upperbound)/term0

    ! lifted directly from the linspace subroutine.
    do j1 = 1, n - 2
        y(j1) = dble(j1)
    end do
    y = ((ub - lb)/(n-1))*y + lb
    x(1) = lb
    x(2:n-1) = y
    x(n) = ub

    ! this is extremely cool, raising a scalar to a vector yields a vector.
    x = 10.0_rk**x
end subroutine logspace
!BL

! ====================================================================
!         FUNCTION: gridlookup: find nearest (at or below) gridloc
! ====================================================================
subroutine gridlookup1(gridnum, xdist, xval, ixlow, iweight)

    integer(ik9):: gridnum, ixhigh, ixlow, ixplace
    real(rk):: xdist(gridnum), xval, iweight
    intent(in):: gridnum, xdist, xval
    intent(out):: ixlow, iweight

    ixhigh = gridnum; ixlow = 1_ik9

    do

        if ((ixhigh	- ixlow).le.1_ik9) then
            exit
        end	if

        ixplace	= (ixhigh +	ixlow)/2_ik9

        if (xdist(ixplace).ge.xval)	then
            ixhigh = ixplace
        else
            ixlow =	ixplace
        end	if

    end	do

    iweight = (xdist(ixlow + 1_ik9) - xval)/(xdist(ixlow + 1_ik9) - xdist(ixlow))

end subroutine gridlookup1
!BL

FUNCTION interp1(X,Y,siz,xi)
	!! Interpolates to find the values of the underlying function Y at the points in the vector or array xi.
	!! X and Y must be one-dimensional arrays with the same dimension siz
	!! xi must be a scalar
	!! if xi exceeds the bounds of X, i.e. requiring extrapolation, then Yi is set to the closest point
	
	! INPUT & OUTPUT
    INTEGER, INTENT(IN) 					:: siz
    REAL(rk), DIMENSION(siz), INTENT(IN)	:: X, Y
    REAL(rk), INTENT(IN)	                :: xi
    REAL(rk)                                :: interp1
    INTEGER								    :: ind
    
	! PROGRAM
	ind = MINLOC( ABS( xi - X ) , 1 )
	outer: IF ( xi >= X(ind) ) THEN
		inner1: IF  ( ind == siz ) THEN
			interp1 = Y(ind)
		ELSE
			interp1 =	&
				( Y(ind+1)*(xi-X(ind)) + Y(ind)*(X(ind+1)-xi) )	&
				/( X(ind+1) - X(ind) )
		END IF inner1
    ELSE ! IF ( xi < X(ind) )
		inner2: IF  ( ind == 1 ) THEN
			interp1 = Y(ind)
		ELSE
			interp1 =	&
				( Y(ind)*(xi-X(ind-1)) + Y(ind-1)*(X(ind)-xi) )	&
				/( X(ind) - X(ind-1) )
		END IF inner2
    END IF outer
	
END FUNCTION interp1   
!BL


!##############################################################################
! SUBROUTINE error
!
! Throws error message and stops program.
!##############################################################################
subroutine error(routine, message)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
        
    ! routine in which error occured
    character(len=*), intent(in) :: routine
        
    ! error message
    character(len=*), intent(in) :: message
        
        
    !##### ROUTINE CODE #######################################################
        
    ! write error message
    write(*,'(/a,a,a,a/)')'ERROR ',routine,': ',message
        
    ! stop program
    stop
    
end subroutine error
!BL

!##############################################################################
! SUBROUTINE grid_Cons_Grow
!
! Constructs a growing grid on [left, right].
!##############################################################################
subroutine grid_Cons_Grow(a, left, right,Num, growth)
       
    implicit none
    !##### INPUT/OUTPUT VARIABLES #############################################
    INTEGER, INTENT(IN) :: Num
    ! left and right interval point
    real*8, intent(in) :: left, right
    ! the array to fill
    real*8, intent(out) :: a(1:Num)       
    
    ! the growth rate of the grid
    real*8 :: growth     
     
    !##### OTHER VARIABLES ####################################################     
    real*8 :: h
    integer :: j, n     
     
    !##### ROUTINE CODE #######################################################
        
    ! get number of grid points
    !n = size(a, 1)-1
    n = Num - 1
    
    ! check for left <= right
    if(left >= right)call error('grid_Cons_Grow', &
        'left interval point greater than right point')

    ! check for growth
    if(growth <= 0d0)call error('grid_Cons_Grow', &
        'growth rate must be greater than zero')

    ! calculate factor
    h = (right-left)/((1d0+growth)**n-1d0)
     
    ! calculate grid value
    a = h*((1d0+growth)**(/(dble(j), j=0,n)/)-1d0)+left

end subroutine grid_Cons_Grow
!BL
END MODULE auxiliary 
