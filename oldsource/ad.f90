module dualClass
implicit none

	type, public :: dual	
		real(kind=8) :: rp
		real(kind=8), dimension(dualSize) :: ip   ! Real and infinitesimal part
	contains
		procedure :: initIndependent
		procedure :: initDependent
	end type dual
	
! Operator overload
	interface operator(+)
		module procedure dual_dual_add, dual_real_add, real_dual_add
	end interface
	
	interface operator(-)
		module procedure dual_dual_subs, dual_real_subs, real_dual_subs, dual_sign_change
	end interface
	
	interface operator(*)
		module procedure dual_dual_mult, dual_real_mult, real_dual_mult
	end interface
	
	interface operator(/)
		module procedure dual_dual_div, dual_real_div, real_dual_div
	end interface
	
	interface operator(**)
		module procedure dual_pow, dual_pow_int
	end interface
	
	interface exp
		module procedure dual_exp
	end interface
	
	interface sin
		module procedure dual_sin
	end interface
	
	interface cos
		module procedure dual_cos
	end interface
	
	interface log
		module procedure dual_log
	end interface
	
	interface tanh
		module procedure dual_tanh
	end interface
	
	contains

! Creation of the dual number
	subroutine initIndependent(this, nDerivatives, whichDerivative)
	class(dual), intent(inout) :: this
	integer :: nDerivatives, whichDerivative
			
		this%ip = 0.d0
		this%ip(whichDerivative) = 1.d0
		
	end subroutine initIndependent
	
	subroutine initDependent(this, nDerivatives)
	class(dual), intent(inout) :: this
	integer :: nDerivatives
		
		this%ip = 0.d0		
		
	end subroutine initDependent
	
! Functions for defining the basic operations

!---------------------------
! Multiplication of two dual numbers
	function dual_dual_mult(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp * b%rp
		c%ip = a%rp * b%ip + b%rp * a%ip
	end function dual_dual_mult
	
! Multiplication of a dual and a real
	function dual_real_mult(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp * b
		c%ip = a%ip * b
	end function dual_real_mult
	
! Multiplication of a real and a dual
	function real_dual_mult(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp * a
		c%ip = b%ip * a
	end function real_dual_mult
	
!---------------------------	
! Division of two dual numbers
	function dual_dual_div(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp / b%rp
		c%ip = (b%rp * a%ip - a%rp * b%ip) / b%rp**2
	end function dual_dual_div
	
! Division of a dual and a real number
	function dual_real_div(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp / b
		c%ip = a%ip / b
	end function dual_real_div
	
! Division of a real number and a dual
	function real_dual_div(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = a / b%rp
		c%ip = -a * b%ip / b%rp**2
	end function real_dual_div	

!---------------------------
! Addition of two dual numbers
	function dual_dual_add(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp + b%rp
		c%ip = a%ip + b%ip
	end function dual_dual_add
	
! Addition of a dual and a real number
	function dual_real_add(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp + b
		c%ip = a%ip
	end function dual_real_add
	
! Addition of a real number and a dual
	function real_dual_add(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp + a
		c%ip = b%ip
	end function real_dual_add	

!---------------------------
! Substraction of two dual numbers
	function dual_dual_subs(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp - b%rp
		c%ip = a%ip - b%ip
	end function dual_dual_subs
	
! Substraction of a dual and a real number
	function dual_real_subs(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp - b
		c%ip = a%ip
	end function dual_real_subs
	
! Substraction of a real number and a dual
	function real_dual_subs(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp - a
		c%ip = b%ip
	end function real_dual_subs
	
! Change of sign of a dual number
	function dual_sign_change(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = -a%rp
		c%ip = -a%ip
	end function dual_sign_change

!---------------------------
! Power	
	function dual_pow(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function dual_pow
	
! Power	
	function dual_pow_int(a,b) result(c)
		type (dual), intent(in) :: a
		integer, intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function dual_pow_int

!---------------------------
! Exponential	
	function dual_exp(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		real(kind=8) :: t
		t = exp(a%rp)
		c%rp = t
		c%ip = a%ip * t
	end function dual_exp	

!---------------------------
! Sin	
	function dual_sin(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = sin(a%rp)
		c%ip = a%ip * cos(a%rp)
	end function dual_sin

!---------------------------	
! Cos
	function dual_cos(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = cos(a%rp)
		c%ip = -a%ip * sin(a%rp)
	end function dual_cos

!---------------------------
! ln
	function dual_log(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = log(a%rp)
		c%ip = a%ip / a%rp
	end function dual_log
	
!---------------------------
! Hyperbolic tangent
	function dual_tanh(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		real(kind=8) :: t
		c%rp = tanh(a%rp)
		c%ip = 1.d0/cosh(a%rp)**2 * a%ip
	end function dual_tanh
	
	
end module dualClass