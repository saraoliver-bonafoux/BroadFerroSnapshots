program IterativeDiscreteVM
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of the biased voter model, which has two absorbing 
!          states, using the iterative algorithm. The process has a discrete state space {0, 1,..., N} and two
!	       absorbing states: n = 0 and n = N.
!
! Inputs : N - population size
! 		   R - bias towards the preferred opinion	
!          K - number of iteration steps
!		   s - relaxation factor	
!
! 		   Notes: The transition rates and the initial distribution for the algorithm are computed internally 
!		   as a function of R.
!
! Outputs: q - array containing the QSD Q(n) for n = 1,..., N-1
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
  	integer, parameter :: N = 100, K = 5 * 10**4
  	double precision, parameter :: R = 0.9d0, s = 0.1d0
  	double precision, dimension(N-1) :: q, q0, fq0, Wp, Wm
  	call cpu_time(t1)

	!----------------------------------------------
	! Transition rates and initial distribution
	!----------------------------------------------
    Wp = [(R * i * (N - i) / N, i = 1, N-1)]   ! Array of rates W(n -> n + 1) 
	Wm = [(i * (N - i) / real(N), i = 1, N-1)] ! Array of rates W(n -> n - 1) 
	if (R.lt.1) then      ! Process biased towards the absorbing state n = 0	
		n0 = 10
		q(n0) = 1.0d0        ! Initial guess for the QSD: Kronecker delta distribution centred at n = n0
	else if (R.gt.1) then ! Process biased towards the absorbing state n = N	
		n0 = 90
		q(n0) = 1.0d0        ! Initial guess for the QSD: Kronecker delta distribution centred at n = n0       
	else if (R.eq.1) then ! There is no preferred opinion     
		q = 1.0d0 / (N - 1)  ! Initial guess for the QSD: uniform distribution
	end if
				
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K
		q0 = q
     	call iterVM(N, q0, fq0, Wp, Wm)
     	q = s * q0 + (1 - s) * fq0 
     	q = q / sum(q) ! Normalisation of the probability distribution 
     	if (mod(i, 1000).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
  	end do
  	  	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
  	open(unit = 1, file = "QSD_VM.txt", status = "replace", action = "write")
  	do i = 1, N-1
  		write(1,*) i, q(i)
  	end do
  	
  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1
   	
contains

	subroutine iterVM(N, q, fq, Wp, Wm)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : N  - population size
	!             q  - current distribution
	!             Wp - array of rates W(n -> n+1)
	!             Wm - array of rates W(n -> n-1)
	! Outputs   : fq - updated distribution after one iteration
	!----------------------------------------------------------------------------------
		integer, intent(in) :: N
		double precision, dimension (N-1), intent(in) :: q, Wp, Wm
		double precision, dimension (N-1), intent(out) :: fq
		double precision :: JS
  
		JS = q(1) * Wm(1) + q(N-1) * Wp(N-1) ! Flux into the absorbing states n = 0 and n = N

		fq(1) = q(2) * Wm(2) / (Wp(1) + Wm(1) - JS)
		do j = 2, N-2
			fq(j) = (q(j+1) * Wm(j+1) + q(j-1) * Wp(j-1)) / (Wp(j) + Wm(j) - JS)
		end do
		fq(N-1) = q(N-2) * Wp(N-2) / (Wp(N-1) + Wm(N-1) - JS)

	end subroutine iterVM

end program IterativeDiscreteVM
