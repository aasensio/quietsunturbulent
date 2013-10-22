module geodesiclmModule
use globalModule
use mathsModule, only : fvoigt
	real(kind=8), allocatable :: xvar(:), fvec(:), fjac(:,:), dtd(:,:)
	integer :: m, n, mode, niters, nfev, njev, naev, maxiters, maxfev, maxjev, converged, print_level, info
	integer :: maxaev, print_unit, imethod, iaccel, ibroyden, ibold, whichLine
	logical :: analytic_jac, analytic_avv, center_diff
	real(kind=8) :: eps, h1, h2, maxlam, artol, cgoal, gtol, xtol, xrtol, ftol, frtol, initial_factor
	real(kind=8) :: factoraccept, factorreject, avmax
		
contains

!-----------------------------------------------------------------------
! Solve the optimization problem using a geodesic Levenberg-Marquardt algorithm
!-----------------------------------------------------------------------
	subroutine maximumLikelihoodGeodesicLM(nParameters, pars, line, stddev)
	integer :: nParameters, line
	real(kind=8) :: pars(nParameters), stddev

		n = nParameters
		m = lineList%transition(line)%nLambda
				
		whichLine = line
		
		if (allocated(xvar)) then			
			deallocate(xvar)
			allocate(xvar(n))
			
			deallocate(fvec)
			allocate(fvec(m))
			
			deallocate(fjac)
			allocate(fjac(m,n))
			
			deallocate(dtd)
			allocate(dtd(n,n))
		else
			allocate(xvar(n))
			allocate(fvec(m))
			allocate(fjac(m,n))
			allocate(dtd(n,n))
		endif		
		
		dtd = 1.d-3
		
		analytic_jac = .FALSE.
		analytic_Avv = .FALSE.
		center_diff = .FALSE.
		
		xvar = pars
		
		eps = 1.d-5
		h1 = 0.1d0
		h2 = 0.1d0
		mode = 0
					
		maxiters = 200
		maxfev = 0
		maxjev = 0
		maxlam = -1
		artol = 1.d-5
		cgoal = 1.d-8
		gtol = 1.d-8
		xtol = 1.d-8
		xrtol = 1.d-8
		frtol = 1.d-8
		print_level = 0
		print_unit = 6
		imethod = 0
		initial_factor = 0.1d0
		factoraccept = 10
		factorreject = 10
		avmax = 2
		info = 0
		ibroyden = 0
		
		iaccel = 1
	! 	ibold = 1     ! Type of acceptance
				
		call geolevmar(func, jacobian, Avv, xvar, fvec, fjac, n, m, callback, info,&
					analytic_jac, analytic_Avv, center_diff, eps, h1, h2,&
					dtd, mode, niters, nfev, njev, naev,&
					maxiters, maxfev, maxjev, maxaev, maxlam,&
					artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,&
					converged, print_level, print_unit,&
					imethod, iaccel, ibold, ibroyden,&
					initial_factor, factoraccept, factorreject, avmax)
					
		pars = xvar
		
		stddev = sqrt(sum(fvec**2) / m)
								
	end subroutine maximumLikelihoodGeodesicLM

!-----------------------------------------------------------------------
! Forward problem
!-----------------------------------------------------------------------
	subroutine func(m, n, x, fvec)
	integer :: m, n
	real(kind=8) :: x(n), fvec(m), f, g(n)
	real(kind=8) :: beta0, etal, deltaL, damping, prof(m)
	
		beta0 = x(1)
		etal = x(2)
		deltaL = x(3)
		damping = x(4)
		
		if (beta0 < 0 .or. etal < 0 .or. deltaL < 0 .or. damping < 0) then
			fvec = 1.d10			
		else
		
			prof = fvoigt(damping, (lineList%transition(whichLine)%lambda - lineList%transition(whichLine)%lambdaNew) / deltaL) / sqrt(PI)
		
			fvec = (1.d0 + beta0 / (1.d0+etal*prof)) / (1.d0 + beta0)	
		
			fvec = (fvec - lineList%transition(whichLine)%observed)
		endif
															
	end subroutine func

!-----------------------------------------------------------------------
! Jacobian
!-----------------------------------------------------------------------
	subroutine jacobian(m, n, x, fjac)
	integer :: m, n
	real(kind=8) :: x(n), fjac(m,n), f, g(n)
				
							
	end subroutine jacobian

!-----------------------------------------------------------------------
! Directional derivative
!-----------------------------------------------------------------------
	subroutine Avv(m, n, x, v, acc)
	integer :: m, n
	real(kind=8) :: x(n), v(n), acc(m)
		
		
	end subroutine Avv

!-----------------------------------------------------------------------
! Callback
!-----------------------------------------------------------------------
	subroutine callback(m,n,x,fvec,info)
	integer :: m, n, info
	real(kind=8) :: x(n), fvec(m)
		return
	end subroutine callback

end module geodesiclmModule