program MCDiscreteSIRS
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of the SIRS model, a stochastic model with two
!          stochastic variables, using a Monte Carlo method with resetting at the absorbing boundaries.
!	       m and n represent the number of susceptible and infected individuals, respectively. 
!		   m and n can take the values 0, 1,..., NN and 0, 1,..., NN-m, respectively, such that m + n <= NN.
!		   Set of absorbing states: 	S = {(m, 0) : 0 <= m <= NN} (with respect to the infected population)
!   	   Set of non-absorbing states: R = {(m, n) : m = 0, 1,..., NN-1, n = 1, 2,..., NN-m}
!
! Inputs : NN - population size
!		   TT - reproduction ratio
!          K  - number of simulation steps
! 
! 		   Notes: The transition rates and the initial state of the simulation are computed internally 
!		   as a function of TT. 
!     	
! Outputs: q - array containing the QSD Q(m, n) for m = 0, 1,..., NN-1 and n = 1, 2,..., NN (some values are 0) 
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
	integer, parameter :: NN = 20
	double precision, parameter :: TT = 0.6d0
	integer*8, parameter :: K = int8(10)**int8(9)
	integer*8 :: j, cabs
	double precision, dimension(0:NN-1, 1:NN) :: hist, qnorm
	double precision :: MTE, MTE2, MTE_error
    call cpu_time(t1)    
	call dran_ini(1234) ! Initialisation of the random number generator (program dranxor.f90, available in this repository)

	!----------------------------------------------
	! Initial condition
	!----------------------------------------------
	if (TT.gt.1) then      ! Process biased towards the metastable state 
 		m = int(NN/TT) 
 		n = int(NN * (TT - 1) / (2 * TT))
 	else if (TT.le.1) then ! Process biased towards the absorbing state (m, n) = (NN, 0)
 		m = NN - 10
 		n = 10
 	end if
 	t = 0.0d0		
 	hist = 0.0d0
 	
 	!----------------------------------------------
	! Auxiliary quantities to calculate the mean time to extinction (MTE) along with its error
	!----------------------------------------------
    tabs = 0.0d0	  ! Elapsed time since the last absorption
    tabs_sum = 0.0d0  ! Sum of absorption times
    tabs2_sum = 0.0d0 ! Sum of squared absorption times
    cabs = 0 		  ! Counter of the number of absorptions
 	
 	!----------------------------------------------
	! Simulation of the process conditioned on non-extinction using the Gillespie algorithm + a resetting mechanism
	!----------------------------------------------
    do j = 1, K

		!!! Gillespie step
		m0 = m
		n0 = n
		t0 = t
		call Gillespie_step_SIRS(m, n, dt)
		hist(m0, n0) = hist(m0, n0) + dt ! Update of the histogram of residence times
    	tabs = tabs + dt
    	
     	!!! Resetting after absorption
     	if (n.eq.0) then 
     		call search_list_SIRS(hist / sum(hist), dran_u(), m, n)
			cabs = cabs + 1
    		tabs_sum = tabs_sum + tabs
    		tabs2_sum = tabs2_sum + tabs**2
    		tabs = 0 
     	end if
     	
     	if (mod(j, K/100).eq.0) print '(I12, E20.13)', j, sum((hist/sum(hist))**2) ! Monitoring of the convergence measure

    end do
    
    !----------------------------------------------
	! Export the QSD and MTE to a text file
	!----------------------------------------------    
    open(unit = 1, file = "QSD_SIRS.txt", status = "replace", action = "write")

	!!! Calculation of the QSD
 	hist = hist / sum(hist)
  	write(1,*) "# m, n, Q(m, n)"
  	do n = 1, NN
		do m = 0, NN-1
			write(1, *) m, n, hist(m, n)
		end do
  	end do
  	
  	!!! Calculation of the MTE
	MTE = tabs_sum / cabs
	MTE2 = tabs2_sum / cabs
    MTE_error = sqrt((MTE2 - MTE**2) / cabs) 
  	write(1,*) "# Number of absorptions, MTE, MTE error"
  	write(1,*) cabs, MTE, MTE_error
  	
  	call cpu_time(t2)
	print *, "tCPU (s) = ", t2 - t1
  	
contains

	subroutine Gillespie_step_SIRS(m, n, dt)
	!----------------------------------------------------------------------------------
	! Purpose   : Advances the stochastic process by one Gillespie step. Given the current state
	!             (m, n), the routine selects the next state (m, n) and computes the time dt until
	!             the transition occurs.  
	! Inputs    : m  - current number of S
	!             n  - current number of I
	! Outputs   : m  - updated number of S
	!			  n  - updated number of I
	!			  dt - time required to perform the jump 
	!----------------------------------------------------------------------------------
		integer, intent(inout) :: m, n
		double precision, intent(out) :: dt
  		
  		WI = TT * m * n / NN ! Rate W[(m, n) -> (m-1, n+1)]
  		WR = n				 ! Rate W[(m, n) -> (m, n-1)]
  		WS = NN - m - n	     ! Rate W[(m, n) -> (m+1, n)]
  		W = WI + WR + WS     ! Total escape rate from the current state
  		
  		dt = -log(dran_u()) / W ! Time required to perform the transition 
    	
  		! State after the jump 
  		r = dran_u() * W
  		if (r.lt.WI) then
  			m = m - 1
  			n = n + 1
  		else if (r.lt.WI + WR) then
  			n = n - 1
  		else
  			m = m + 1
  		end if
  		
	end subroutine Gillespie_step_SIRS
	
	subroutine search_list_SIRS(hist_norm, p, m, n)
	!----------------------------------------------------------------------------------
	! Purpose   : Given a 2D array of probabilities "hist_norm" such that sum(hist_norm) = 1
	!             and a probability value "p" in [0, 1], looks for the pair of indices (m, n)
	!             where the cumulative sum of the array first exceeds p.
	! Inputs    : hist_norm, p
	! Outputs   : m, n 
	!----------------------------------------------------------------------------------
		double precision, dimension(0:NN-1, 1:NN), intent (in) :: hist_norm
    	double precision, intent (in) :: p
    	integer, intent(out) :: m, n
    	integer :: S
        
		C = 0.0d0
    	S = 0
    	I = 1
	
		do while (C.lt.p .and. S.le.NN-1)
			if (I.le.NN - S) then
				C = C + hist_norm(S, I)
				if (C.ge.p) exit
				I = I + 1
			else
				S = S + 1
				I = 1
			end if
		end do
		m = S
		n = I

	end subroutine search_list_SIRS
	
end program MCDiscreteSIRS









