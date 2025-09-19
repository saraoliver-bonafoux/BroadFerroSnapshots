program MCContinuous2DBM
!----------------------------------------------------------------------------------------------------------
! Purpose: Computes the quasi-stationary distribution (QSD) of a 2D Brownian motion in the (x, y) plane using  
!          a Monte Carlo method with resetting at the absorbing boundaries. 
! 		   The radial variable is confined to Rin <= r <= Rout. Rin is absorbing, and Rout is reflecting. 
!
!          2D Histogram (x, y): 
!          - Simulation domain: [-Rext, Rext] × [-Rext, Rext], with Rext >= Rout.
!          - The domain is discretised into subintervals of length h along both axes, giving a uniform n × n grid with 
!            n = nint(2 * Rext / h). 
!          - For each grid cell (i, j), the coordinates of its bottom-left corner are given by:
!            x = -Rext + (j - 1) * h, i = 1,..., n
!		     y = Rext - i * h,        j = 1,..., n
! 		   - Example: Rext = 1.5, h = 0.5 -> n = 6, hist_xy = 6 × 6.
!                     Rows (i) decrease in y from top to bottom.
!                     Columns (j) increase in x from left to right.
!
!          1D Histogram (r):
!          - The radial domain r \in [Rin, Rout] is discretised into subintervals of length dr:
!            Ik = [Rin + (k-1)*dr, Rin + k*dr] for k = 1,..., nr, with nr = nint((Rout - Rin)/dr)
!                
! Inputs : Rin  - inner radius (absorbing barrier) 
!          Rout - outer radius (reflecting barrier)
!          Rext - size of the square domain for the 2D QSD (Rext >= Rout)
!		   h    - spatial step size for the 2D histogram (x, y)
!          dr   - spatial step size for the 1D radial histogram
!          dt   - time step size
!          D    - diffusion coefficient
!
! Outputs: hist_xy - 2D QSD histogram Q(x, y)
!          hist_r  - 1D QSD histogram Q(r)
!          MTE     - mean time to extinction (along with its statistical error)	 
!----------------------------------------------------------------------------------------------------------
  	implicit double precision(a-h,o-z)
	integer*8, parameter :: K = int8(10)**int8(11)
	double precision, parameter :: Rin = 0.3d0, Rout = 1.2d0, Rext = 1.3d0, D = 0.5d0, h = 5.0d0 * 10.0d0**(-3), dr = 1.0d0 * 10.0d0**(-3), dt = 1 * 10.0d0**(-4)
	double precision, allocatable :: hist_xy(:,:), hist_r(:)
	integer*8 :: cabs, i 
	double precision :: MTE, MTE2, MTE_error
	call cpu_time(t1)
	call dran_ini(1234) ! Initialisation of the random number generator (program dranxor.f90, available in this repository)
	
	n = nint(2 * Rext / h)       ! Dimension of the 2D histogram (x, y)
	allocate(hist_xy(n, n))
	
	nr = nint((Rout - Rin) / dr) ! Dimension of the 1D radial histogram
    allocate(hist_r(nr))
    
    !----------------------------------------------
	! Auxiliary quantities and initial state
	!----------------------------------------------
	Ddt = D * dt
	sqrt2Ddt = sqrt(2 * Ddt)
	call sample_xy_donut(Rin, Rout, x, y) ! Initial condition randomly sampled uniformly within the annular domain 
	t = 0.0d0
	hist_xy = 0.0d0
	hist_r = 0.0d0
	
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
     	y0 = y
     	r0 = sqrt(x0**2 + y0**2)
     	t0 = t
     	x = x0 - x0 / r0**2 * Ddt + sqrt2Ddt * dran_gbmw()
     	y = y0 - y0 / r0**2 * Ddt + sqrt2Ddt * dran_gbmw()
     	r = sqrt(x**2 + y**2)
     	t = t0 + dt
		tabs = tabs + dt 
     	! Update of the 2D and 1D histograms
     	i0 = int((Rext - y0) / h) + 1  ! y-index of the bin of the previous position 
		j0 = int((x0 + Rext) / h) + 1  ! x-index of the bin of the previous position 
     	hist_xy(i0, j0) = hist_xy(i0, j0) + 1
     	k0 = int((r0 - Rin) / dr) + 1  ! Index of the bin of the previous position in 1D (Radial coordinate)
     	hist_r(k0) = hist_r(k0) + 1
     	
     	!!! Absorption / Reflection
     	if (r.ge.Rout) then     ! Reflection
     		call reflection(x0, y0, x, y, Rout)
     	else if (r.le.Rin) then ! Absorption
  			call search_list_linear_algorithm(reshape(hist_xy / sum(hist_xy), [n**2]), dran_u(), index) ! We randomly select an interval (in 2D) according to the histogram of residence times
  			ii = mod(index - 1, n) + 1 ! y-index of the selected bin
			jj = (index - 1) / n + 1   ! x-index of the selected bin 
     		call sample_xy_bin(ii, jj, Rin, Rout, x, y) ! We place the process in a position uniformly sampled within the selected interval
     		cabs = cabs + 1
    		tabs_sum = tabs_sum + tabs
    		tabs2_sum = tabs2_sum + tabs**2
    		tabs = 0 
     	end if
     	
     	if (mod(i, K/100).eq.0) print '(I12, E20.13)', i, sum((hist_xy/sum(hist_xy))**2) ! Monitoring of the convergence measure
    	
     end do
    
 	!----------------------------------------------
	! Export the QSD to a text file + Print the value obtained for the MTE 
	!----------------------------------------------
 	open(unit = 1, file = "QSD_RBM_xy.txt", status = "replace", action = "write")
	open(unit = 2, file = "QSD_RBM_r.txt", status = "replace", action = "write")
 	 	 	
 	!!! Calculation of the QSD in 2D 
 	hist_xy = hist_xy / (sum(hist_xy) * h**2)
 	
 	write(1,*) "# x, y, Q(x, y)"
   	do kk = 1, n**2
   		ii = mod(kk - 1, n) + 1
		jj = (kk - 1) / n + 1
		xc = -Rext + (jj - 1) * h + h / 2
     	yc = Rext - ii * h + h / 2
     	write(1,*) xc, yc, hist_xy(ii, jj)
	end do
	
	!!! Calculation of the QSD in 1D (Radial component)
	hist_r = hist_r / sum(hist_r)
	hist_r = hist_r / dr
	
	do kk = 1, nr
		write(2, *), Rin + (kk - 1) * dr + dr / 2, hist_r(kk)
	end do

	!!! Calculation of the MTE
	MTE = tabs_sum / cabs
	MTE2 = tabs2_sum / cabs
    MTE_error = sqrt((MTE2 - MTE**2) / cabs) 
	print *, "# Number of absorptions, MTE, MTE error"
  	print *, cabs, MTE, MTE_error
  	
  	call cpu_time(t2)
 	print *, "tCPU (s) = ", t2 - t1


contains

	recursive subroutine sample_xy_donut(Rin, Rout, x, y)
	!----------------------------------------------------------------------------------
	! Purpose   : Randomly samples a point (x, y) uniformly within an annular domain  
	!			  bounded by Rin < r < Rout. 
	! Inputs    : Rin, Rout  
	! Outputs   : x, y 
	! Notes     : Recursive call ensures the sampled point is inside the annular domain.
	!----------------------------------------------------------------------------------
  		double precision, intent (in) :: Rin, Rout
  		double precision, intent (out) :: x, y  		
  		
     	x = Rout * (2 * dran_u() - 1)
     	y = Rout * (2 * dran_u() - 1)
     	r = sqrt(x**2 + y**2)
     	
  		if (r.gt.Rout .or. r.lt.Rin) then 
     		call sample_xy_donut(Rin, Rout, x, y)
  		endif
  		
	end subroutine sample_xy_donut

	recursive subroutine sample_xy_bin(ii, jj, Rin, Rout, x, y)
	!----------------------------------------------------------------------------------
	! Purpose   : Randomly samples a position (x, y) uniformly within a specific 2D grid
	!             cell (ii, jj), constrained to lie within the annular domain Rin <= r <= Rout.
	! Inputs    : ii, jj, Rin, Rout  
	! Outputs   : x, y 
	! Notes     : Recursive call ensures the sampled point is inside the annular domain.
	!----------------------------------------------------------------------------------
		double precision, intent (in) :: Rin, Rout
  		integer, intent(in) :: ii, jj
  		double precision, intent (out) :: x, y
  		  		
  		x_lc = -Rext + (jj - 1) * h ! x-coordinate of the bottom-left corner of the bin
     	y_lc = Rext - ii * h        ! y-coordinate of the bottom-left corner of the bin
     	x = x_lc + h * dran_u()
     	y = y_lc + h * dran_u()
     	r = sqrt(x**2 + y**2)
     	
  		if (r.gt.Rout .or. r.lt.Rin) then 
     		call sample_xy_bin(ii, jj, Rin, Rout, x, y)
  		endif
  		
	end subroutine sample_xy_bin	
	
	recursive subroutine reflection(x0, y0, x, y, R)
	!----------------------------------------------------------------------------------	
	! Purpose   : Given an initial point (x0, y0) that attempts to move to (x, y) outside a circle
	!             of radius R, reflects (x, y) inside the circle, allowing for multiple reflections. 
	! Inputs    : (x0, y0) - initial point (inside the circle)
	!             (x, y)   - final point   (outside the circle)
	!             R        - circle radius
	! Outputs   : (x, y)   - reflected point (inside the circle)
	! Notes     : - Recursive call ensures the reflected point is inside the circle.
	! 			  - See the document "ReflectionCircle.pdf", available in this repository, for an
	!             explanation of this implementation of reflecting boundary conditions.
	!----------------------------------------------------------------------------------
  		double precision, intent(in) :: x0, y0, R
  		double precision, intent(in out) :: x, y
  		double precision :: lambda
  		
  		dx = x - x0
  		dy = y - y0
  		
  		!!! Coordinates of point a (Point on the boundary of the circle on the line connecting the initial and final points)
  		dd = dx**2 + dy**2 						    ! Auxiliary quantity
  		det = dd * R**2 - (x0 * dy - y0 * dx)**2    ! Auxiliary quantity
  		lambda = (-(x0 * dx + y0 * dy) + sqrt(det)) / dd
  		ax = x0 + lambda * dx
  		ay = y0 + lambda * dy
  		
  		!!! Coordinates of the reflected point
  		prod = (x * ax + y * ay) / R**2		        ! Auxiliary quantity
  		x = x + 2.0d0 * (1.0d0 - prod) * ax
  		y = y + 2.0d0 * (1.0d0 - prod) * ay
  		
  		!!! Multiple reflections
  		! If the reflected point lies outside the circle, we apply another reflection with (ax, ay) as initial point, 
  		! and (x, y) (reflected point) as final point.
  		if (x**2 + y**2.gt.R**2) then 
     		call reflection(ax, ay, x, y, R)
  		end if
  		
	end subroutine reflection
	
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
	
end program MCContinuous2DBM


