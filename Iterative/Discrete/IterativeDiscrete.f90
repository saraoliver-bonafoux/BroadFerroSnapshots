program IterativeDiscrete 
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of a one-step stochastic process with a discrete
!          state space {0, 1,..., N} and an absorbing state at n = 0 using the iterative algorithm.
!
! Inputs : N  - upper bound of the state space
!          K  - number of iteration steps
!		   s  - relaxation factor		
!		   Wp - array of rates W(n -> n + 1) 
!		   Wm - array of rates W(n -> n - 1) 
!
!		   Notes: In this example, the transition rates and the initial distribution for the algorithm are  
!		   computed internally for a branching process.
!
! Outputs: q - array containing the QSD Q(n) for n = 1,..., N
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
  	integer, parameter :: N = 150, K = 10**4
  	double precision, parameter :: s = 0.1d0 
  	double precision, dimension(N) :: q, q0, fq0, Wp, Wm
  	call cpu_time(t1)

	!----------------------------------------------
	! Transition rates and initial distribution. Example: branching process
	!----------------------------------------------
	R = 0.8d0 					! Reproduction ratio 
    Wp = [(R * i, i = 1, N)] 	! Birth rates
    Wp(N) = 0					! No transitions out of the state space
	Wm = [(i, i = 1, N)] 		! Death rates
	n0 = 10 	  				! State near the absorbing state n = 0
	q(n0) = 1.0d0 				! Initial guess for the QSD: Kronecker delta distribution centred at n = n0
				
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K
		q0 = q
     	call iter(N, q0, fq0, Wp, Wm)
     	q = s * q0 + (1 - s) * fq0 
     	q = q / sum(q) ! Normalisation of the probability distribution 
     	if (mod(i, 1000).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
  	end do
  	  	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
  	open(unit = 1, file = "QSD_Branching.txt", status = "replace", action = "write")
  	do i = 1, N
  		write(1,*) i, q(i)
  	end do
  	
  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1
   	
contains

	subroutine iter(N, q, fq, Wp, Wm)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : N  - upper bound of the state space
	!             q  - current distribution
	!             Wp - birth rates
	!             Wm - death rates
	! Outputs   : fq - updated distribution after one iteration
	!----------------------------------------------------------------------------------
		integer, intent(in) :: N
		double precision, dimension (N), intent(in) :: q, Wp, Wm
		double precision, dimension (N), intent(out) :: fq
		double precision :: JS
  
		JS = q(1) * Wm(1) ! Flux into the absorbing state n = 0 
		
		fq(1) = q(2) * Wm(2) / (Wp(1) + Wm(1) - JS)
		do j = 2, N-1
			fq(j) = (q(j+1) * Wm(j+1) + q(j-1) * Wp(j-1)) / (Wp(j) + Wm(j) - JS)
		end do
		fq(N) = q(N-1) * Wp(N-1) / (Wp(N) + Wm(N) - JS)

	end subroutine iter
	
end program IterativeDiscrete
