module adMaths
implicit none
contains

!------------------------------------------------------------------
! Calculates log(Gamma(x))
! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
! http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
!------------------------------------------------------------------
	function alngam ( xvalue, ifault )

	  real ( kind = 8 ) alngam
	  real ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
	  integer ( kind = 4 ) ifault
	  real ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
	    -2.66685511495D+00, &
	    -24.4387534237D+00, &
	    -21.9698958928D+00, &
	     11.1667541262D+00, &
	     3.13060547623D+00, &
	     0.607771387771D+00, &
	     11.9400905721D+00, &
	     31.4690115749D+00, &
	     15.2346874070D+00 /)
	  real ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
	    -78.3359299449D+00, &
	    -142.046296688D+00, &
	     137.519416416D+00, &
	     78.6994924154D+00, &
	     4.16438922228D+00, &
	     47.0668766060D+00, &
	     313.399215894D+00, &
	     263.505074721D+00, &
	     43.3400022514D+00 /)
	  real ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
	    -2.12159572323D+05, &
	     2.30661510616D+05, &
	     2.74647644705D+04, &
	    -4.02621119975D+04, &
	    -2.29660729780D+03, &
	    -1.16328495004D+05, &
	    -1.46025937511D+05, &
	    -2.42357409629D+04, &
	    -5.70691009324D+02 /)
	  real ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
	     0.279195317918525D+00, &
	     0.4917317610505968D+00, &
	     0.0692910599291889D+00, &
	     3.350343815022304D+00, &
	     6.012459259764103D+00 /)
	  real ( kind = 8 ) x
	  real ( kind = 8 ) x1
	  real ( kind = 8 ) x2
	  real ( kind = 8 ), parameter :: xlge = 5.10D+05
	  real ( kind = 8 ), parameter :: xlgst = 1.0D+30
	  real ( kind = 8 ) xvalue
	  real ( kind = 8 ) y

	  x = xvalue
	  alngam = 0.0D+00
	!
	!  Check the input.
	!
	  if ( xlgst <= x ) then
	    ifault = 2
	    return
	  end if

	  if ( x <= 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  ifault = 0
	!
	!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
	!
	  if ( x < 1.5D+00 ) then

	    if ( x < 0.5D+00 ) then

	      alngam = - log ( x )
	      y = x + 1.0D+00
	!
	!  Test whether X < machine epsilon.
	!
	      if ( y == 1.0D+00 ) then
	        return
	      end if

	    else

	      alngam = 0.0D+00
	      y = x
	      x = ( x - 0.5D+00 ) - 0.5D+00

	    end if

	    alngam = alngam + x * (((( &
	        r1(5)   * y &
	      + r1(4) ) * y &
	      + r1(3) ) * y &
	      + r1(2) ) * y &
	      + r1(1) ) / (((( &
	                  y &
	      + r1(9) ) * y &
	      + r1(8) ) * y &
	      + r1(7) ) * y &
	      + r1(6) )

	    return

	  end if
	!
	!  Calculation for 1.5 <= X < 4.0.
	!
	  if ( x < 4.0D+00 ) then

	    y = ( x - 1.0D+00 ) - 1.0D+00

	    alngam = y * (((( &
	        r2(5)   * x &
	      + r2(4) ) * x &
	      + r2(3) ) * x &
	      + r2(2) ) * x &
	      + r2(1) ) / (((( &
	                  x &
	      + r2(9) ) * x &
	      + r2(8) ) * x &
	      + r2(7) ) * x &
	      + r2(6) )
	!
	!  Calculation for 4.0 <= X < 12.0.
	!
	  else if ( x < 12.0D+00 ) then

	    alngam = (((( &
	        r3(5)   * x &
	      + r3(4) ) * x &
	      + r3(3) ) * x &
	      + r3(2) ) * x &
	      + r3(1) ) / (((( &
	                  x &
	      + r3(9) ) * x &
	      + r3(8) ) * x &
	      + r3(7) ) * x &
	      + r3(6) )
	!
	!  Calculation for 12.0 <= X.
	!
	  else

	    y = log ( x )
	    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

	    if ( x <= xlge ) then

	      x1 = 1.0D+00 / x
	      x2 = x1 * x1

	      alngam = alngam + x1 * ( ( &
	             r4(3)   * &
	        x2 + r4(2) ) * &
	        x2 + r4(1) ) / ( ( &
	        x2 + r4(5) ) * &
	        x2 + r4(4) )

	    end if

	  end if

	  return
	end function alngam

!*****************************************************************************80
!
!! DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!  Modified:
!
!    03 June 2013
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, 1976, pages 315-317.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the digamma function.
!    0 < X.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) DIGAMA, the value of the digamma function at X.
!
! http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
	function digama ( x, ifault )


	implicit none

	real ( kind = 8 ), parameter :: euler_mascheroni = 0.57721566490153286060D+00
	real ( kind = 8 ) digama
	integer ( kind = 4 ) ifault
	real ( kind = 8 ) r
	real ( kind = 8 ) x
	real ( kind = 8 ) x2
	!
	!  Check the input.
	!
		if ( x <= 0.0D+00 ) then
			digama = 0.0D+00
			ifault = 1
			return
		end if
		!
		!  Initialize.
		!
		ifault = 0
		x2 = x
		digama = 0.0D+00
		!
		!  Approximation for small argument.
		!
		if ( x2 <= 0.00001D+00 ) then
			digama = - euler_mascheroni - 1.0D+00 / x2
			return
		end if
		!
		!  Reduce to DIGAMA(X + N).
		!
		do while ( x2 < 8.5D+00 )
			digama = digama - 1.0D+00 / x2
			x2 = x2 + 1.0D+00
		end do
		!
		!  Use Stirling's (actually de Moivre's) expansion.
		!
		r = 1.0D+00 / x2
		digama = digama + log ( x2 ) - 0.5D+00 * r
		r = r * r
		digama = digama &
			- r * ( 1.0D+00 / 12.0D+00 &
			- r * ( 1.0D+00 / 120.0D+00 &
			- r *   1.0D+00 / 252.0D+00 ) )

		return
	end function digama
	
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function voigtAutomaticDiff(daInput, dv)
	real*8 :: da, daInput
	real*8 :: dv, voigtAutomaticDiff(2)
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
		
! There is a problem for very small dampings. We force it to be 0 using a threshold
		da = daInput
		if (daInput < 1.d-3) da = 0.d0

		z = cmplx(dv, da)
		t = cmplx(da, -dv)
		s = dabs(dv) + da
		u = t*t


		if (s >= 15.d0) then
			w4 = t * 0.5641896d0 / (0.5d0+u)
		elseif (s >= 5.5) then
			w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		elseif (da >= 0.195d0*dabs(dv)-0.176d0) then
			w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
				(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
		else
			w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
				u*(1.320522d0-u*0.56419d0))))))
			v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
				u*(61.57037d0-u*(1.841439d0-u)))))))
			w4 = cexp(u) - w4/v4
		endif
		voigtAutomaticDiff(1) = dble(w4)
		voigtAutomaticDiff(2) = aimag(w4)

	end function voigtAutomaticDiff

end module adMaths