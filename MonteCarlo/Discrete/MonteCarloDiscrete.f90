program MCDiscrete
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of a stochastic process with a discrete state
!          space {0, 1,..., NN} using a Monte Carlo method with resetting at the absorbing boundaries. 
!		   The process may have multiple absorbing states.  
!
! Inputs : NN       - upper bound of the state space
!          K        - number of simulation steps
!		   nabslist - array of absorbing states
!
! 		   Notes: 
!	       - In this example, the transition rates and the initial state of the simulation are computed  
!		   internally for the biased voter model, a one-step process with 2 absorbing states (n = 0 and n = NN).
!          - For multi-step processes, all transition rates must be specified, and the subroutine "Gillespie_step"
!     	   must be adapted accordingly.
!
! Outputs: q   - array containing the QSD Q(n) for all non-absorbing states
!          MTE - mean time to extinction (along with its statistical error)	
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
	integer, parameter :: NN = 100
	integer*8, parameter :: K = int8(10)**int8(11)
	integer*8 :: i, cabs
	integer, allocatable :: nabslist(:)
	double precision, dimension(NN) :: hist, Wp, Wm
	double precision :: MTE, MTE2, MTE_error
    call cpu_time(t1)
	call dran_ini(1234) ! Initialisation of the random number generator (program dranxor.f90, available in this repository)

	!----------------------------------------------
	! Transition rates and initial condition. Example: biased voter model
	!----------------------------------------------
    R = 0.9d0									! Bias towards the preferred opinion
    Wp = [(R * i * (NN - i) / NN, i = 1, NN)]	! Array of rates W(n -> n + 1) 
	Wm = [(i * (NN - i) / real(NN), i = 1, NN)] ! Array of rates W(n -> n - 1) 
	nabslist = [0, NN]							! Array of absorbing states
	if (R.lt.1) then	  ! Process biased towards the absorbing state n = 0
		n = 10
	else if (R.gt.1) then ! Process biased towards the absorbing state n = N	
		n = 90
	else if (R.eq.1) then ! There is no preferred opinion     
		n = 50
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
	do i = 1, K
        
    	!!! Gillespie step
    	n0 = n
    	t0 = t
    	call Gillespie_step(n, dt, Wp, Wm)
    	t = t0 + dt
    	hist(n0) = hist(n0) + dt ! Update of the histogram of residence times
    	tabs = tabs + dt
    	
    	!!! Resetting after absorption
		if (any(nabslist == n)) then
    		call search_list_linear_algorithm(hist / sum(hist), dran_u(), n)    		
    		cabs = cabs + 1
    		tabs_sum = tabs_sum + tabs
    		tabs2_sum = tabs2_sum + tabs**2
    		tabs = 0 
    	end if
    	
		if (mod(i, K/100).eq.0) print '(I12, E20.13)', i, sum((hist/sum(hist))**2) ! Monitoring of the convergence measure
     	
	end do

	!----------------------------------------------
	! Export the QSD and MTE to a text file
	!----------------------------------------------
	open(unit = 1, file = "QSD_VM.txt", status = "replace", action = "write")

    !!! Calculation of the QSD
	hist = hist / sum(hist)
  	write(1,*) "# n, Q(n)"
  	do i = 0, NN
  		if (.not. any(nabslist == i)) then ! The state i is non-absorbing
  			write(1,*) i, hist(i)
  		end if
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

	subroutine Gillespie_step(n, dt, Wp, Wm)
	!----------------------------------------------------------------------------------
	! Purpose   : Advances a one-step stochastic process by one Gillespie step. Given the current state n
	!             and the transition rates Wp(n) = W(n -> n + 1) and Wm(n) = W(n -> n−1), the
	!             routine selects the next state n and computes the time dt until the transition occurs.  
	! Inputs    : n  - current state of the process
	!             Wp - array of rates W(n -> n + 1) 
	!			  Wm - array of rates W(n -> n - 1) 
	!             Notes: For multi-step processes, the subroutine must be adapted to include all allowed
	!	          transitions.
	! Outputs   : n  - updated state after the jump
	!			  dt - time required to perform the jump 
	!----------------------------------------------------------------------------------
		integer, intent(inout) :: n
		double precision, dimension (:), intent(in) :: Wp, Wm
		double precision, intent(out) :: dt
  		
  		Wpn = Wp(n)
  		Wmn = Wm(n)
  		W = Wpn + Wmn ! Total escape rate from the current state
  		
  		dt = -log(dran_u()) / W ! Time required to perform the transition 
  		
  		! State after the jump 
  		r = dran_u() * W
  		if (r.lt.Wpn) then
  			n = n + 1
  		else
  			n = n - 1
  		end if
  		
	end subroutine Gillespie_step

	subroutine search_list_linear_algorithm(list, p, index)
	!----------------------------------------------------------------------------------
	! Purpose   : Given a list of probabilities "list" such that sum(list) = 1 and a 
	!             probability value "p" in [0, 1], looks for "index" such that:       
	!             	C(index) >= p and C(j) < p for all j in [1, index[
	!             where C(j) = list(1) + list(2) + ... + list(j)
	! Inputs    : list, p    
	! Outputs   : index 
	!----------------------------------------------------------------------------------
    	double precision, dimension (:), intent (in) :: list
    	double precision, intent (in) :: p
    	integer, intent(out) :: index
        
        j = 1
        C = list(j)
        do while (p.gt.C)
        	j = j + 1
            C  = C + list(j)
        end do
    	index = j
    	
	end subroutine search_list_linear_algorithm

end program MCDiscrete









