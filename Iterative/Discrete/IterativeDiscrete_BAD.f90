program IterativeDiscreteBAD 
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of the branching-annihilation-decay process, 
!          which involves two-step transitions, using the iterative algorithm. The process has a discrete 
!          state space {0, 1,..., N} and an absorbing state at n = 0.
!
! Inputs : N - upper bound of the state space
!          V - volume of the system
!		   R - reproduction ratio
!          K - number of iteration steps
!		   s - relaxation factor
!
! 		   Notes: The transition rates and the initial distribution for the algorithm are computed internally
!		   as a function of R.
!
! Outputs: q - array containing the QSD Q(n) for n = 1,..., N
!----------------------------------------------------------------------------------------------------------
	implicit double precision(a-h,o-z)
  	integer, parameter :: N = 150, V = 250, K = 10**4
  	double precision, parameter :: R = 0.9d0 , s = 0.1d0
  	double precision, dimension(N) :: q, q0, fq0, Wp, Wm, Wm2
  	call cpu_time(t1)

	!----------------------------------------------
	! Transition rates and initial distribution
	!----------------------------------------------
    Wp = [(R * i, i = 1, N)]  					  ! Array of rates W(n -> n + 1)
    Wp(N) = 0									  ! No transitions out of the state space
	Wm = [(i, i = 1, N)]						  ! Array of rates W(n -> n - 1)
	Wm2 = [(R * i * (i - 1) / (2 * V), i = 1, N)] ! Array of rates W(n -> n - 2)
	if (R.gt.1) then ! Process biased towards the metastable state
		n0 = nint(1 + N * (1 - 1/R))				 
	else             ! Process biased towards the absorbing state n = 0
		n0 = 10										 
	end if
	q(n0) = 1.0d0 	 ! Initial guess for the QSD: Kronecker delta distribution centred at n = n0
	
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K
		q0 = q
     	call iterBAD(N, q0, fq0, Wp, Wm, Wm2)
     	q = s * q0 + (1 - s) * fq0 
     	q = q / sum(q) ! Normalisation of the probability distribution 
     	if (mod(i, 1000).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
  	end do
  	  	  	  	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
  	open(unit = 1, file = "QSD_BAD.txt", status = "replace", action = "write")
  	do i = 1, N
  		write(1,*) i, q(i)
  	end do
  	
  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1
   	
contains

	subroutine iterBAD(N, q, fq, Wp, Wm, Wm2)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : N   - upper bound of the state space
	!             q   - current distribution
	!             Wp  - birth rates
	!             Wm  - death rates
	!			  Wm2 - annihilation rates
	! Outputs   : fq  - updated distribution after one iteration
	!----------------------------------------------------------------------------------
		integer, intent(in) :: N
		double precision, dimension (N), intent(in) :: q, Wp, Wm, Wm2
		double precision, dimension (N), intent(out) :: fq
		double precision :: JS
  
		JS = q(1) * Wm(1) + q(2) * Wm2(2) ! Flux into the absorbing state n = 0 

		fq(1) = (q(2) * Wm(2) + q(3) * Wm2(3)) / (Wp(1) + Wm(1) - JS)
		do j = 2, N-1
			fq(j) = (q(j+1) * Wm(j+1) + q(j-1) * Wp(j-1) + q(j+2) * Wm2(j+2)) / (Wp(j) + Wm(j) + Wm2(j) - JS)
		end do
		fq(N-1) = (q(N) * Wm(N) + q(N-2) * Wp(N-2)) / (Wp(N-1) + Wm(N-1) + Wm2(N-1) - JS)
		fq(N) = q(N-1) * Wp(N-1) / (Wp(N) + Wm(N) - JS)

	end subroutine iterBAD

end program IterativeDiscreteBAD
