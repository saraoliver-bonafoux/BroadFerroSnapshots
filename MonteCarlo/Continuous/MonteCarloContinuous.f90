program MCContinuous
!----------------------------------------------------------------------------------------------------------!
! Purpose: Computes the quasi-stationary distribution (QSD) of a stochastic process with a continuous state
!          space [0, L] using a Monte Carlo method with resetting at the absorbing boundaries. 
!		   We consider the following two types of systems:
!		   - flag = 0: absorbing boundary at x = 0 and reflecting boundary at x = L
!		   - flag = 1: absorbing boundaries at both x = 0 and x = L
!		   The state space is discretised into subintervals of length dx, which we write as Ii = [(i-1)*dx, i*dx]
!	       for i = 1,..., N, with N = L / dx. 	
!
! Inputs : L    - upper boundary of the state space
!	       flag - system type (0 = absorbing/reflecting, 1 = absorbing/absorbing)
!          K    - number of simulation steps
!          dx   - spatial step size
!          dt   - time step size
!
! 		   Notes: In this example, the drift and diffusion functions and the initial state of the simulation 
!		   are computed internally for the biased random walk (corresponding to flag = 0).
!
! Outputs: q   - array containing the QSD Q(xi) evaluated at xi = dx / 2 + (i - 1) * dx, with i = 1,..., N
!          MTE - mean time to extinction (along with its statistical error)	
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
	integer*8, parameter :: K = int8(10)**int8(11) 
	integer, parameter :: flag = 0
 	integer*8 :: i, cabs 
	double precision, parameter :: L = 1.0d0, dx = 1 * 10.0d0**(-3), dt = 1 * 10.0d0**(-4)
	double precision, allocatable :: hist(:) 
 	double precision :: MTE, MTE2, MTE_error
    call cpu_time(t1)
	call dran_ini(1234) ! Initialisation of the random number generator (program dranxor.f90, available in this repository)

	N = int(L / dx)
    allocate(hist(N))
    sqrtdt = sqrt(dt) 
	
	!----------------------------------------------
	! Parameters and initial state. Example: Biased random walk
	!----------------------------------------------
    D = 0.4d0
    f = -1.0d0
    x = 0.1d0
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
	! Simulation of the process conditioned on non-extinction using the Euler-Maruyama scheme + a resetting mechanism
	!----------------------------------------------     
    do i = 1, K

		!!! Euler-Maruyama step     
     	x0 = x
     	t0 = t
     	x = x0 + A(x0) * dt + sqrt(B(x0)) * sqrtdt * dran_gbmw() 
		t = t0 + dt
		i0 = int(x0/dx) + 1 ! Index of the bin of the previous position     	
     	hist(i0) = hist(i0) + 1.0d0
     	tabs = tabs + dt 

		!!! Absorption / Reflection
     	if ((flag.eq.0 .and. x.le.0) .or. (flag.eq.1 .and. (x.le.0 .or. x.ge.L))) then 	! Absorption 
    		call search_list_linear_algorithm(hist / sum(hist), dran_u(), index) ! We randomly select an interval according to the histogram of residence times
    		x = (index - 1) * dx + dran_u() * dx ! We place the process in a position uniformly sampled in the selected interval	
    		cabs = cabs + 1
    		tabs_sum = tabs_sum + tabs
    		tabs2_sum = tabs2_sum + tabs**2
    		tabs = 0 
    	else if (flag.eq.0 .and. x.ge.L) then ! Reflection
    		x = 2.0d0 * L - x
        end if
        
        if (mod(i, K/100).eq.0) print '(I12, E20.13)', i, sum((hist/sum(hist))**2) ! Monitoring of the convergence measure
    
    end do	
 	 	
 	!----------------------------------------------
	! Export the QSD and MTE to a text file
	!----------------------------------------------
	open(unit = 1, file = "QSD_BRW.txt", status = "replace", action = "write")
	
	!!! Calculation of the QSD
 	hist = hist / (K * dx)	
 	write(1,*) "# x, Q(x)"  	  	
    do i = 1, N
		write(1,*) dx / 2 + (i - 1) * dx, hist(i)
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

	double precision function A(x) ! Drift function
        double precision, intent(in) :: x
        A = f 
	end function A
	
	double precision function B(x) ! Diffusion function
        double precision, intent(in) :: x
        B = 2 * D 
	end function B
	
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

end program MCContinuous









