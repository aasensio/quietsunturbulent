module mathsModule
use globalModule, only : PHC, PK, PHK, PI, PME, PH, PC, PH2C3, PHC2, EV_ERG, BOHR_R, UMA
use dualClass
implicit none
	
contains

		
!-------------------------------------------------------------------
! Calculate the natural logarithm of the Gamma function
!-------------------------------------------------------------------
	function gammln(xx)
	real(kind=8) :: gammln, xx
	integer :: j 
	real(kind=8) :: ser,tmp,x,y
	real(kind=8) :: cof(6) = (/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,&
		.1208650973866179d-2,-.5395239384953d-5/)
	real(kind=8) :: stp = 2.5066282746310005d0
		
		x = xx 
		y = x 
		tmp = x+5.5d0 
		tmp = (x+0.5d0)*log(tmp)-tmp 
		ser = 1.000000000190015d0 
		do j=1,6 
			y = y+1.d0 
			ser = ser + cof(j)/y 
		enddo 
		gammln = tmp + log(stp*ser/x)
	end function gammln

!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function close_to(wave_total, wave)
	real(kind=8) :: wave_total(:), wave
	real(kind=8), allocatable :: diff(:)
	integer :: n, i, location(1)
	integer :: close_to

		do i = 1, size(wave_total)
			if (wave_total(i) > wave) then
				close_to = i
				return
			endif
		enddo

	end function close_to

!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function extremes(lambda, leftLambda, rightLambda)
	real(kind=8) :: lambda(:), leftLambda, rightLambda
	integer :: left, right, extremes(2)


		left = close_to(lambda, leftLambda)
		right = close_to(lambda, rightLambda)

		extremes(1) = min(left,right)
		extremes(2) = max(left,right)

	end function extremes

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt(size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n

		n = size(dv)
		do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold
			if (da(i) < 1.d-3) da(i) = 0.d0
			z = cmplx(dv(i), da(i))
			t = cmplx(da(i), -dv(i))
			s = dabs(dv(i)) + da(i)
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da(i) >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			vecfvoigt(i) = dble(w4)
		enddo

	end function vecfvoigt
	
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function fvoigt(daInput, dv)
	real*8 :: da, daInput
	real*8 :: dv(:), fvoigt(size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n

		n = size(dv)
		
! There is a problem for very small dampings. We force it to be 0 using a threshold
		da = daInput
		if (daInput < 1.d-3) da = 0.d0
		do i = 1, n
			z = cmplx(dv(i), da)
			t = cmplx(da, -dv(i))
			s = dabs(dv(i)) + da
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			fvoigt(i) = dble(w4)
		enddo

	end function fvoigt
		
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function fvoigtScalar(daInput, dv)
	real*8 :: da, daInput
	real*8 :: dv, fvoigtScalar
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
		fvoigtScalar = dble(w4)

	end function fvoigtScalar


	!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	subroutine voigtZeeman(da, dv, H, L)
	real*8 :: da
	real*8 :: dv(:), H(size(dv)), L(size(dv))
	complex :: w4(size(dv)), z(size(dv)), t(size(dv)), u(size(dv)), v4(size(dv))
	real*8 :: s(size(dv))
	integer :: i, n
                
		n = size(dv)
				
		if (da < 1.d-3) then
			da = 0.d0
		endif
		
		z = cmplx(dv,da)
		t = cmplx(da, -dv)
		s = dabs(dv) + da
		u = t*t
		
		where(s >= 15.d0)
			w4 = t * 0.5641896d0 / (0.5d0+u)
		endwhere
		where(s < 15.d0 .and. s >= 5.5d0)
			w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da >= 0.195d0*dabs(dv)-0.176d0)
			w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da < 0.195d0*dabs(dv)-0.176d0)
			w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
			v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
			w4 = cexp(u) - w4/v4
		endwhere
		
		H = dble(w4)
		L = aimag(w4)
				
	end subroutine voigtZeeman
		
!-----------------------------------------------------------------------
! Return a sigmoid function
!-----------------------------------------------------------------------
	elemental function sigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = lower + (upper-lower) / (1.d0 + exp(-x))
	end function sigmoid
	
!-----------------------------------------------------------------------
! Return the inverse of the sigmoid function
!-----------------------------------------------------------------------
	elemental function invSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = log( (lower - x) / (x - upper) )
	end function invSigmoid
	
!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	elemental function diffSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = (upper-lower) * exp(-x) / (1.d0+exp(-x))**2
	end function diffSigmoid	

end module mathsModule