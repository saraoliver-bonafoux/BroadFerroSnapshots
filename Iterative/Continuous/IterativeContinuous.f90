program IterativeContinuous
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of a stochastic process with a continuous state
!          space [0, L] using the iterative algorithm. The boundary x = 0 is absorbing, while the boundary
!	       x = L is reflecting.
!          The state space is discretised into subintervals of length dx, so that xi = i * dx for
!	       i = 0,..., N, with N = L / dx. 	
!
! Inputs : L    - upper boundary of the state space
!          K    - number of iteration steps
!		   s    - relaxation factor
!		   dx   - spatial step size			
!		   A(x) - drift function
!		   B(x) - diffusion function
! 
! 		   Notes: In this example, the drift and diffusion functions correspond to the Ornstein-Uhlenbeck process.
!
! Outputs: q - array containing the QSD Q(xi) for i = 1,..., N
!----------------------------------------------------------------------------------------------------------	
    implicit double precision(a-h,o-z)
	integer*8, parameter :: K = 1.8 * 10**7
  	double precision, parameter :: L = 3.0d0, s = 0.1d0, dx = 5 * 10.0d0**(-4)
  	double precision, allocatable :: q(:), q0(:), fq0(:)
	double precision :: mu
  	call cpu_time(t1)
  	
  	N = int(L / dx) ! Number of mesh points
    allocate(q(N))
	allocate(q0(N))
	allocate(fq0(N))

	!----------------------------------------------
	! Parameters and initial distribution. Example: Ornstein-Uhlenbeck process
	!----------------------------------------------
	D = 0.3d0
    f = 2.0d0
    mu = 0.0d0
	x0 = 0.1d0	          ! Point near the absorbing boundary x = 0
	q(int(x0/dx)) = 1.0d0 ! Initial guess for the QSD: Kronecker delta distribution centred at x = x0
	q = q / (sum(q) * dx) ! Normalisation of the probability density function
	
	!----------------------------------------------
	! Iterative scheme
	!----------------------------------------------
	do i = 1, K	
		q0 = q 		
  		call iter(N, L, q0, fq0)
		q = s * q0 + (1 - s) * fq0
		q = q / (sum(q) * dx) ! Normalisation of the probability density function
     	if (mod(i, 10**5).eq.0) print '(I8, E20.13)', i, sqrt(sum((q - q0)**2.0d0)) ! Monitoring of the convergence measure
	end do    
	
	!----------------------------------------------
	! Export the QSD to a text file
	!----------------------------------------------
	open(unit = 1, file = "QSD_Ornstein.txt", status = "replace", action = "write")
	do i = 1, N
		write(1,*) i * dx, q(i)
   	end do
   	
  	call cpu_time(t2)
   	print *, "tCPU (s) = ", t2 - t1
    
contains
 
	double precision function A(x) ! Drift function
        double precision, intent(in) :: x
        A = f * (mu - x) 
	end function A
	
	double precision function B(x) ! Diffusion function
        double precision, intent(in) :: x
        B = 2 * D 
	end function B
 
	subroutine iter(N, L, q, fq)
	!----------------------------------------------------------------------------------
	! Purpose   : Performs one iteration step of the algorithm.
	! Inputs    : N  - number of mesh points
	!		      L  - upper boundary of the state space
	!             q  - current distribution
	! Outputs   : fq - updated distribution after one iteration
	!----------------------------------------------------------------------------------
  		integer, intent(in) :: N
  		double precision, intent(in) :: L
  		double precision, dimension (N), intent(in) :: q
  		double precision, dimension (N), intent(out) :: fq
  		double precision :: JS
  		
  		JS = B(dx) * q(1) / (2 * dx) ! Flux into the absorbing boundary x = 0
  		dxB1q1 = 2 * dx**2 * JS 	 ! Auxiliary quantity
  		
		fq(1) = q(2) * (B(2*dx) - dx * A(2*dx)) / (2 * B(dx) - dxB1q1)
		do j = 2, N-1
		    x = j * dx
		    xpdx = x + dx
		    xmdx = x - dx
			fq(j) = (q(j+1) * (B(xpdx) - dx * A(xpdx)) + q(j-1) * (B(xmdx) + dx * A(xmdx))) / (2 * B(x) - dxB1q1)
  		end do
  		fq(N) = B(L - dx) / (B(L) - 2 * dx * A(L)) * q(N - 1)

	end subroutine iter	
	
end program IterativeContinuous


