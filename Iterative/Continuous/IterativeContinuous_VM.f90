program IterativeContinuousVM
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of the continuous biased voter model using the 
!		   iterative algorithm. The process has a continuous state space [0, L], with both boundaries x = 0
!	       and x = L absorbing.       
!          The state space is discretised into subintervals of length dx, so that xi = i * dx for
!	       i = 0,..., N, with N = L / dx. 
!	
! Inputs : NN - population size
!		   R  - bias towards the preferred opinion	
!	       L  - upper boundary of the state space (always equal to 1 in this system)
!		   K  - number of iteration steps
!		   s  - relaxation factor
!		   dx - spatial step size	
!		
! 		   Notes: The drift and diffusion functions and the initial distribution for the algorithm are computed
!		   internally as a function of R.
!
! Outputs: q - array containing the QSD Q(xi) for i = 1,..., N-1
!----------------------------------------------------------------------------------------------------------	
	implicit double precision(a-h,o-z)
  	integer*8, parameter :: NN = 100, K = 1.3 * 10**7, L = 1
  	double precision, parameter :: R = 0.9d0, s = 0.1d0, dx = 5 * 10.0d0**(-4) 
  	double precision, allocatable :: q(:), q0(:), fq0(:)
  	call cpu_time(t1)
  	
  	N = int(L / dx) ! Number of mesh points
    allocate(q(N-1))
	allocate(q0(N-1))
	allocate(fq0(N-1))
	
	!----------------------------------------------
	! Initial distribution
	!----------------------------------------------
	if (R.lt.1) then      ! Process biased towards the absorbing boundary x = 0
		x0 = 0.1d0
		q(x0/dx) = 1.0d0      ! Initial guess for the QSD: Kronecker delta distribution centred at x = x0
	else if (R.gt.1) then ! Process biased towards the absorbing boundary x = L
		n0 = 0.9d0
		q(x0/dx) = 1.0d0      ! Initial guess for the QSD: Kronecker delta distribution centred at x = x0       
	else if (R.eq.1) then ! There is no preferred opinion     
		q = 1.0d0             ! Initial guess for the QSD: uniform distribution
	end if
	q = q / (sum(q) * dx) ! Normalisation of the probability density function
	
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K
		q0 = q 		
  		call iterVM(N, L, q0, fq0)
		q = s * q0 + (1 - s) * fq0
		q = q / (sum(q) * dx) ! Normalisation of the probability density function
     	if (mod(i, 10**5).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
	end do    
	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
	open(unit = 1, file = "QSD_VMC.txt", status = "replace", action = "write")
	do i = 1, N-1
		write(1,*) i * dx, q(i)
   	end do
   	
  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1
	
contains
 
	double precision function A(x) ! Drift function
        double precision, intent(in) :: x
		A = (R - 1) * x * (1 - x)
	end function A
	
	double precision function B(x) ! Diffusion function
        double precision, intent(in) :: x
        B = (R + 1) * x * (1 - x) / NN
	end function B
 
	subroutine iterVM(N, L, q, fq)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : N  - number of mesh points
	!		      L  - upper boundary of the state space
	!             q  - current distribution
	! Outputs   : fq - updated distribution after one iteration
	!----------------------------------------------------------------------------------
  		integer, intent(in) :: N
  		double precision, dimension (N-1), intent(in) :: q
  		double precision, dimension (N-1), intent(out) :: fq
  		
		do j = 1, N-1
		    x = j * dx
		    xpdx = x + dx
		    xmdx = x - dx
			fq(j) = (q(j+1) * (B(xpdx) - dx * A(xpdx)) + q(j-1) * (B(xmdx) + dx * A(xmdx))) / (2 * B(x) - dx * (B(L-dx) * q(N-1) + B(dx) * q(1)))
  		end do

	end subroutine iterVM
	
end program IterativeContinuousVM


