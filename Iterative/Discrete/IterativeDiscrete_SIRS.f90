program IterativeDiscreteSIRS
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of the SIRS model, a stochastic model with two
!          stochastic variables, using the iterative algorithm. 
!	       m and n represent the number of susceptible and infected individuals, respectively. 
!		   m and n can take the values 0, 1,..., NN and 0, 1,..., NN-m, respectively, such that m + n <= NN.
!		   Set of absorbing states: 	S = {(m, 0) : 0 <= m <= NN} (with respect to the infected population)
!   	   Set of non-absorbing states: R = {(m, n) : m = 0, 1,..., NN-1, n = 1, 2,..., NN-m}
!
! Inputs : NN - population size
!		   T  - reproduction ratio
!          K  - number of iteration steps
!		   a  - relaxation factor
! 
! 		   Notes: The transition rates and the initial distribution for the algorithm are computed internally 
!		   as a function of T.  
!     	
! Outputs: q - array containing the QSD Q(m, n) for m = 0, 1,..., NN-1 and n = 1, 2,..., NN (some values are 0)
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
	integer, parameter :: NN = 20, K = 10**3
	double precision, parameter :: T = 0.6d0, a = 0.1d0
	double precision, dimension(-1:NN, 0:NN+1) :: q, q0, fq0
  	call cpu_time(t1)
  	
  	!----------------------------------------------
	! Initial distribution
	!----------------------------------------------
 	if (T.gt.1) then      ! Process biased towards the metastable state 
 		m0 = int(NN/T) 
 		n0 = int(NN * (T - 1) / (2 * T))
 	else if (T.le.1) then ! Process biased towards the absorbing state (m, n) = (NN, 0)
 		m0 = NN - 10
 		n0 = 10
 	end if
 	q(m0, n0) = 1 ! Initial guess for the QSD: Kronecker delta distribution centred at (m, n) = (m0, n0)
	
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K
		q0 = q
  		call iterSIRS(NN, q0, fq0)
  		q = a * q0 + (1 - a) * fq0
		q = q / sum(q)
		if (mod(i, 100).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
	end do 
	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
	open(unit = 1, file = "QSD_SIRS.txt", status = "replace", action = "write")
   	write(1,*) "# m, n, Q(m, n)"
  	do n = 1, NN
  		do m = 0, NN-1
			write(1, *) m, n, q(m, n)
  		end do
  	enddo  	

  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1

contains

	double precision function WI(m, n, NN, T) ! Rate W[(m, n) -> (m-1, n+1)]
        integer, intent(in) :: m, n, NN
        double precision, intent(in):: T
        WI = T * m * n / NN
	end function WI

	double precision function WR(n)			  ! Rate W[(m, n) -> (m, n-1)]
        integer, intent(in) :: n
        WR = n
    end function WR
	
	double precision function WS(m, n, NN)    ! Rate W[(m, n) -> (m+1, n)]
        integer, intent(in) :: m, n, NN
        WS = NN - m - n
	end function WS
	
	subroutine iterSIRS(NN, q, fq)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : NN - population size
	!             q  - current distribution
	! Outputs   : fq - updated distribution after one iteration
	!----------------------------------------------------------------------------------
		integer, intent(in) :: NN
		double precision, dimension(-1:NN, 0:NN+1), intent(in) :: q
		double precision, dimension(-1:NN, 0:NN+1), intent(out) :: fq
		double precision :: JS
		
		JS = sum(q(:,1)) ! Flux into the set of absorbing states S = {(m, 0) : 0 <= m <= NN}
  		do m = 0, NN - 1
  			do n = 1, NN - m
  				fq(m, n) = q(m + 1, n - 1) * WI(m + 1, n - 1, NN, T) + q(m, n + 1) * WR(n + 1) + q0(m - 1, n) * WS(m - 1, n, NN)
  				fq(m, n) = fq(m, n) / (WI(m, n, NN, T) + WR(n) + WS(m, n, NN) - JS)	
  			end do
  		end do	
  		
  	end subroutine iterSIRS
  		
end program IterativeDiscreteSIRS


